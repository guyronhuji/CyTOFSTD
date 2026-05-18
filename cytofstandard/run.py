"""Run class for cytofstandard."""

import yaml
import json
import importlib
import importlib.util
import re
import inspect
import time
import sys
import pandas as pd
import anndata
import numpy as np
from pathlib import Path
from datetime import datetime
from typing import Any

from cytofstandard.exceptions import (
    IngestionError,
    CytofStandardError,
    RunValidationError,
    RunNotIngestedError,
    MetadataValidationError,
    MarkerValidationError,
)
from cytofstandard.storage import (
    ensure_directory,
    write_yaml,
    write_parquet,
    read_parquet,
    compute_sha256,
    copy_file,
    file_exists,
    directory_exists,
    set_zarr_parts_writable,
)
from cytofstandard.uuid_utils import generate_cell_uuid
from cytofstandard.validation import (
    validate_file_list,
    validate_sample_metadata,
    validate_marker_consistency,
)
from cytofstandard.io import read_fcs, read_csv
from cytofstandard.provenance import (
    ProvenanceLogger,
    log_ingestion_started,
    log_run_ingested,
    log_ingestion_failed,
)


class Run:
    """CytOF run for managing ingestion and data access."""

    def __init__(
        self,
        project,
        run_id: str,
        run_config: dict,
    ):
        """Initialize run (internal use - use Project.add_run() or get_run()).

        Args:
            project: Parent Project instance
            run_id: Run identifier
            run_config: Run configuration dictionary
        """
        self.project = project
        self.run_id = run_id
        self.run_config = run_config
        self.path = Path(project.path) / "runs" / run_id
        self._adata = None
        self.markers_all: list[str] = []
        self.markers_core_histones: list[str] = []
        self.markers_identity: list[str] = []
        self.markers_functional: list[str] = []
        self.markers_qc: list[str] = []
        self.markers_intracellular: list[str] = []
        self.markers_extracellular: list[str] = []
        self.markers_by_class: dict[str, list[str]] = {}
        self._run_provenance_logger = ProvenanceLogger(
            self.path / "logs" / "provenance.jsonl"
        )
        self._initialize_marker_variables_from_metadata()

    def ingest(
        self,
        files: list[str],
        sample_metadata: str,
        copy_raw: bool = True,
        strict_markers: bool = True,
        allow_extra_markers: bool = False,
        drop_columns: list[str] | None = None,
    ) -> None:
        """Ingest files into this run.

        Args:
            files: List of file paths to ingest
            sample_metadata: Path to sample metadata CSV/Parquet file
            copy_raw: Whether to copy raw files to project
            strict_markers: Whether to fail on unknown markers
            allow_extra_markers: Whether to allow extra markers
            drop_columns: Column names to remove before marker processing

        Raises:
            MetadataValidationError if metadata validation fails
            MarkerValidationError if marker validation fails
            IngestionError if ingestion fails
        """
        # Log ingestion start
        log_ingestion_started(
            self.project._provenance_logger,
            self.project.project_id,
            self.run_id,
            len(files),
        )
        self._log_run_event(
            "run_ingestion_started",
            {
                "n_files": int(len(files)),
                "strict_markers": bool(strict_markers),
                "allow_extra_markers": bool(allow_extra_markers),
                "drop_columns": list(drop_columns or []),
            },
        )

        drop_columns = drop_columns or []

        # Validate files exist
        for f in files:
            if not Path(f).exists():
                raise ValueError(f"File not found: {f}")

        # Load sample metadata
        sample_metadata_path = Path(sample_metadata)
        if sample_metadata_path.suffix == ".csv":
            sample_df = pd.read_csv(sample_metadata_path)
        elif sample_metadata_path.suffix == ".parquet":
            sample_df = pd.read_parquet(sample_metadata_path)
        else:
            sample_df = pd.read_csv(sample_metadata_path)

        # Validate sample metadata
        validate_sample_metadata(sample_df)
        validate_file_list(files, sample_df)

        # Validate sample IDs are unique (per run)
        from cytofstandard.validation import validate_sample_id_unique

        validate_sample_id_unique(sample_df)

        # Copy raw files if requested
        raw_dir = self.path / "raw"
        if copy_raw:
            ensure_directory(raw_dir)
            file_mappings = {}
            for f in files:
                src = Path(f)
                dst = raw_dir / src.name
                copy_file(str(src), str(dst))
                file_mappings[str(src)] = str(dst)
        else:
            file_mappings = {f: f for f in files}

        # Get marker registry
        marker_registry = self.project.get_marker_registry()

        # Read and process each file
        all_data = []
        file_manifest = []
        all_marker_mappings = []
        all_csv_meta_cols = set()

        for f in files:
            src_path = Path(f)
            stored_path = file_mappings[str(src_path)]
            file_hash = compute_sha256(stored_path)

            # Read file
            if src_path.suffix.lower() == ".fcs":
                data_df, marker_metadata = read_fcs(stored_path)
            elif src_path.suffix.lower() == ".csv":
                data_df, marker_metadata = read_csv(stored_path)
            else:
                raise ValueError(f"Unsupported file format: {src_path.suffix}")

            if drop_columns:
                normalized_drop = {str(col).strip() for col in drop_columns}
                columns_to_drop = [
                    col
                    for col in data_df.columns
                    if col in normalized_drop
                    or self._sanitize_observed_marker_name(col) in normalized_drop
                ]
                if columns_to_drop:
                    data_df = data_df.drop(columns=columns_to_drop)

            numeric_columns = data_df.select_dtypes(include=["number"]).columns.tolist()
            if not numeric_columns:
                raise IngestionError(
                    f"No numeric marker columns remain in {src_path.name} after applying drop_columns."
                )

            # Store CSV metadata columns for later
            csv_meta_cols = [c for c in data_df.columns if c not in numeric_columns]
            all_csv_meta_cols.update(csv_meta_cols)

            # Standardize marker names using sanitized channel labels
            observed_markers = [c for c in data_df.columns if c in numeric_columns]
            sanitized_markers = [
                self._sanitize_observed_marker_name(c) for c in observed_markers
            ]
            marker_mapping = marker_registry.standardize_marker_names(
                sanitized_markers,
                strict=strict_markers,
            )
            marker_mapping["sanitized_marker_name"] = marker_mapping[
                "original_marker_name"
            ]
            marker_mapping["original_marker_name"] = observed_markers

            # Add file info to marker mapping
            marker_mapping["run_id"] = self.run_id
            marker_mapping["file_name"] = src_path.name
            all_marker_mappings.append(marker_mapping)

            # Subset data to standard markers that were matched
            mapped = marker_mapping[
                marker_mapping["mapping_status"].isin(
                    ["matched_standard", "matched_alias"]
                )
            ]
            standard_markers = mapped["standard_marker_name"].tolist()
            rename_map = dict(
                zip(
                    mapped["original_marker_name"].tolist(),
                    mapped["standard_marker_name"].tolist(),
                )
            )

            # Rename observed numeric columns to standard names before subsetting.
            data_numeric = data_df[numeric_columns].copy()
            data_standardized = data_numeric.rename(columns=rename_map)
            data_subset = data_standardized[standard_markers].copy()

            # Add cell metadata
            # Get sample info for this file
            file_name_in_metadata = src_path.name
            sample_info = sample_df[
                sample_df["file_name"] == file_name_in_metadata
            ].iloc[0]

            n_events = len(data_subset)
            cell_uuids = [
                generate_cell_uuid(
                    self.project.project_id,
                    self.run_id,
                    file_hash,
                    i,
                )
                for i in range(n_events)
            ]

            cell_metadata = pd.DataFrame(
                {
                    "cell_uuid": cell_uuids,
                    "project_id": self.project.project_id,
                    "run_id": self.run_id,
                    "sample_id": sample_info["sample_id"],
                    "line_id": sample_info["line_id"],
                    "source_file": src_path.name,
                    "source_file_hash": file_hash,
                    "event_index": range(n_events),
                }
            )

            # Add optional metadata fields
            optional_fields = [
                "condition",
                "replicate_id",
                "batch_id",
                "barcoding_id",
                "notes",
            ]
            for field in optional_fields:
                if field in sample_info:
                    cell_metadata[field] = sample_info[field]

            # Add CSV metadata if present (with prefix)
            for col in csv_meta_cols:
                cell_metadata[f"csv_meta_{col}"] = data_df[col].values

            # Combine data and metadata
            data_subset.index = cell_metadata.index
            data_with_meta = pd.concat([data_subset, cell_metadata], axis=1)

            all_data.append(data_with_meta)

            # Add to file manifest
            file_manifest.append(
                {
                    "run_id": self.run_id,
                    "source_path": str(src_path),
                    "stored_path": stored_path,
                    "file_name": src_path.name,
                    "file_size_bytes": src_path.stat().st_size,
                    "file_hash_sha256": file_hash,
                    "file_type": src_path.suffix.lower().lstrip("."),
                    "n_events": n_events,
                    "n_markers": len(standard_markers),
                    "ingested_at": datetime.utcnow().isoformat(),
                }
            )

        # Validate marker consistency across files
        common_markers = validate_marker_consistency(
            all_marker_mappings,
            strict=strict_markers,
            allow_extra_markers=allow_extra_markers,
        )

        # Build AnnData object
        if not all_data:
            raise IngestionError("No data was ingested")

        combined = pd.concat(all_data, ignore_index=True)

        # Separate matrix from metadata
        obs_cols = [
            "cell_uuid",
            "project_id",
            "run_id",
            "sample_id",
            "line_id",
            "source_file",
            "source_file_hash",
            "event_index",
        ]
        optional_field_cols = set(optional_fields)
        csv_meta_cols_list = [f"csv_meta_{c}" for c in all_csv_meta_cols]
        non_marker_cols = set(obs_cols) | optional_field_cols | set(csv_meta_cols_list)

        marker_cols = [c for c in combined.columns if c not in non_marker_cols]
        X = combined[marker_cols].values.astype(np.float32)

        # Create obs DataFrame
        obs = combined[obs_cols].copy()
        for field in optional_fields:
            if field in combined.columns:
                obs[field] = combined[field].values
        obs.index = obs["cell_uuid"].astype(str)
        obs.index.name = None

        # Build var DataFrame
        var = []
        for marker in marker_cols:
            info = marker_registry.get_marker_info(marker)
            entry = {
                "standard_marker_name": marker,
                "marker_class": info.get("marker_class", "unknown")
                if info
                else "unknown",
            }
            if info:
                for field in [
                    "intra_extra",
                    "is_core_histone",
                    "is_identity_marker",
                    "is_functional_marker",
                    "is_qc_marker",
                    "is_intracellular_marker",
                    "is_extracellular_marker",
                    "canonical_metal",
                    "canonical_mass",
                    "description",
                ]:
                    if field in info:
                        entry[field] = info[field]

            # Derive intra/extra boolean flags from intra_extra when available.
            intra_extra_value = entry.get("intra_extra", None)
            if isinstance(intra_extra_value, str):
                token = intra_extra_value.strip().lower()
                if token.startswith("intra"):
                    entry["is_intracellular_marker"] = True
                    entry["is_extracellular_marker"] = False
                elif token.startswith("extra"):
                    entry["is_intracellular_marker"] = False
                    entry["is_extracellular_marker"] = True
            var.append(entry)

        var_df = pd.DataFrame(var)
        var_df.index = var_df["standard_marker_name"].astype(str)

        # Normalize bool-like object columns so AnnData/Zarr does not try to
        # write mixed object arrays as variable-length UTF-8 and fail on bools.
        for col in var_df.columns:
            series = var_df[col]
            if series.dtype != object:
                continue

            non_null = series.dropna()
            if non_null.empty:
                continue

            def _is_bool_like(value):
                if isinstance(value, (bool, np.bool_)):
                    return True
                if isinstance(value, str) and value.strip().lower() in {
                    "true",
                    "false",
                }:
                    return True
                return False

            if non_null.map(_is_bool_like).all():
                var_df[col] = series.map(
                    lambda value: pd.NA
                    if pd.isna(value)
                    else (
                        value
                        if isinstance(value, (bool, np.bool_))
                        else value.strip().lower() == "true"
                    )
                ).astype("boolean")

        # Create AnnData
        adata = anndata.AnnData(
            X=X,
            obs=obs,
            var=var_df,
        )
        adata.layers["raw"] = X

        # Store in uns
        adata.uns["project"] = {
            "project_id": self.project.project_id,
            "project_name": self.project.project_name,
        }
        adata.uns["run"] = {
            "run_id": self.run_id,
            "run_name": self.run_config.get("run_name", ""),
            "panel_id": self.run_config.get("panel_id", ""),
            "acquisition_date": self.run_config.get("acquisition_date", ""),
            "instrument": self.run_config.get("instrument", ""),
            "operator": self.run_config.get("operator", ""),
        }
        adata.uns["ingestion"] = {
            "created_at": datetime.utcnow().isoformat(),
            "package_version": "0.1.0",
            "strict_markers": strict_markers,
            "drop_columns": drop_columns,
            "n_files": len(files),
            "n_cells": len(adata),
            "n_markers": len(marker_cols),
        }

        # Store metadata hashes
        sample_hash = compute_sha256(sample_metadata)
        standard_hash = compute_sha256(
            self.project.metadata_dir / "standard_markers.parquet"
        )
        alias_path = self.project.path / self.project.marker_alias_file
        alias_hash = compute_sha256(alias_path)

        adata.uns["metadata_hashes"] = {
            "sample_metadata_sha256": sample_hash,
            "standard_markers_sha256": standard_hash,
            "marker_aliases_sha256": alias_hash,
        }

        # Write Zarr
        zarr_path = self.path / "processed" / f"{self.run_id}.zarr"
        adata.write_zarr(zarr_path)

        # Save file manifest
        manifest_path = self.path / "metadata" / "file_manifest.parquet"
        manifest_df = pd.DataFrame(file_manifest)
        write_parquet(str(manifest_path), manifest_df)

        # Save marker mappings
        all_marker_df = pd.concat(all_marker_mappings, ignore_index=True)
        marker_mapping_path = self.path / "metadata" / "marker_mapping.parquet"
        write_parquet(str(marker_mapping_path), all_marker_df)

        # Save sample table
        samples_path = self.path / "metadata" / "samples.parquet"
        write_parquet(str(samples_path), sample_df)

        # Generate and save lines table
        lines_df = (
            sample_df.groupby("sample_id")
            .agg(
                {
                    "line_id": "first",
                }
            )
            .reset_index()
        )
        lines_df["n_samples"] = 1
        lines_df["sample_ids"] = lines_df["sample_id"]
        if "condition" in sample_df.columns:
            lines_df["conditions"] = (
                sample_df.groupby("sample_id")["condition"]
                .apply(lambda x: ";".join(x.dropna().astype(str)))
                .values
            )
        else:
            lines_df["conditions"] = ""
        lines_path = self.path / "metadata" / "lines.parquet"
        write_parquet(str(lines_path), lines_df)

        # Update run status
        self.run_config["status"] = "ingested"
        write_yaml(self.path / "run.yaml", self.run_config)

        # Update run registry
        self.project._run_registry.loc[
            self.project._run_registry["run_id"] == self.run_id, "status"
        ] = "ingested"
        self.project._save_run_registry()

        # Log successful ingestion
        log_run_ingested(
            self.project._provenance_logger,
            self.project.project_id,
            self.run_id,
            len(adata),
            len(marker_cols),
            len(files),
        )
        self._log_run_event(
            "run_ingested",
            {
                "n_cells": int(len(adata)),
                "n_markers": int(len(marker_cols)),
                "n_files": int(len(files)),
            },
        )

        self._adata = adata
        self._update_marker_variables(adata)

    def read_adata(self, backed: bool = False):
        """Read the ingested AnnData object.

        Args:
            backed: Whether to use backed mode (currently ignored)

        Returns:
            AnnData object

        Raises:
            RunNotIngestedError if run has not been ingested
        """
        zarr_path = self.zarr_path()

        if not zarr_path.exists():
            raise RunNotIngestedError(
                f"Run {self.run_id} is registered but has no ingested AnnData-Zarr object at: "
                f"{zarr_path}\n\nCall run.ingest(...) first."
            )

        if self._adata is None:
            import anndata

            self._adata = anndata.read_zarr(zarr_path)
            self._update_marker_variables(self._adata)

        return self._adata

    def to_dataframe(
        self,
        fields: list[str],
        layer: str = "X",
    ) -> pd.DataFrame:
        """Return a DataFrame from selected marker and obs fields.

        Args:
            fields: Ordered list of field names to include. Each name must be
                either a marker in `adata.var_names` or an `adata.obs` column.
            layer: Expression layer used for marker values (`X` or layer key).

        Returns:
            DataFrame indexed by `adata.obs_names` with columns in `fields` order.
        """
        if not fields:
            raise ValueError("fields cannot be empty")

        adata = self.read_adata()
        matrix = self._matrix_from_layer(adata, layer)
        var_names = adata.var_names.astype(str).tolist()
        marker_to_idx = {name: i for i, name in enumerate(var_names)}

        missing = [
            field
            for field in fields
            if field not in marker_to_idx and field not in adata.obs.columns
        ]
        if missing:
            raise ValueError(
                f"Requested fields not found as markers or obs columns: {missing}"
            )

        out = pd.DataFrame(index=adata.obs_names)
        for field in fields:
            if field in marker_to_idx:
                out[field] = matrix[:, marker_to_idx[field]]
            else:
                out[field] = adata.obs[field].values

        return out

    def _log_run_event(
        self, event_type: str, data: dict[str, Any] | None = None
    ) -> None:
        """Write provenance event to the per-run JSONL log."""
        payload = {
            "project_id": self.project.project_id,
            "run_id": self.run_id,
        }
        if data:
            payload.update(data)
        self._run_provenance_logger.log(event_type=event_type, data=payload)

    def save(self, adata: anndata.AnnData | None = None) -> None:
        """Persist current run AnnData to disk.

        Use this after external modifications to `run.read_adata()` output.

        Args:
            adata: Optional AnnData object to save. If omitted, saves the current
                in-memory `self._adata`.

        Raises:
            RunNotIngestedError: If run has no zarr and no adata was provided.
        """
        if adata is not None:
            self._adata = adata

        if self._adata is None:
            if not self.zarr_path().exists():
                raise RunNotIngestedError(
                    f"Run {self.run_id} has no AnnData to save. "
                    "Call run.ingest(...) first or pass adata=..."
                )
            self._adata = self.read_adata()

        self._adata.write_zarr(self.zarr_path())
        self._update_marker_variables(self._adata)
        self._log_run_event(
            "run_saved",
            {
                "n_cells": int(self._adata.n_obs),
                "n_markers": int(self._adata.n_vars),
            },
        )

        # Ensure run/project status reflects persisted AnnData.
        if self.run_config.get("status") != "ingested":
            self.run_config["status"] = "ingested"
            write_yaml(self.path / "run.yaml", self.run_config)

            if self.run_id in self.project._run_registry["run_id"].values:
                self.project._run_registry.loc[
                    self.project._run_registry["run_id"] == self.run_id,
                    "status",
                ] = "ingested"
                self.project._save_run_registry()

    def save_adata(self, adata: anndata.AnnData | None = None) -> None:
        """Alias for `save()` for explicit naming in notebooks."""
        self.save(adata=adata)

    def rename(self, new_run_name: str) -> None:
        """Rename this run via project metadata (run_name only)."""
        self.project.rename_run(self.run_id, new_run_name)
        self.run_config["run_name"] = str(new_run_name).strip()

    def _initialize_marker_variables_from_metadata(self) -> None:
        """Populate marker variables for loaded runs without reading full AnnData."""
        marker_mapping_path = self.path / "metadata" / "marker_mapping.parquet"
        if not marker_mapping_path.exists():
            return

        try:
            mapping_df = read_parquet(str(marker_mapping_path))
        except Exception:
            return

        if "standard_marker_name" not in mapping_df.columns:
            return

        matched_df = mapping_df
        if "mapping_status" in mapping_df.columns:
            matched_df = mapping_df[
                mapping_df["mapping_status"].isin(["matched_standard", "matched_alias"])
            ]

        markers = (
            matched_df["standard_marker_name"]
            .dropna()
            .astype(str)
            .drop_duplicates()
            .tolist()
        )
        if not markers:
            return

        registry_df = self.project.get_marker_registry().standard_markers.copy()
        if "standard_marker_name" not in registry_df.columns:
            return

        registry_df["standard_marker_name"] = registry_df[
            "standard_marker_name"
        ].astype(str)
        marker_set = set(markers)
        var_df = registry_df[
            registry_df["standard_marker_name"].isin(marker_set)
        ].copy()
        if var_df.empty:
            return

        marker_order = {marker: i for i, marker in enumerate(markers)}
        var_df["__order"] = var_df["standard_marker_name"].map(marker_order)
        var_df = var_df.sort_values("__order").drop(columns=["__order"])
        var_df.index = var_df["standard_marker_name"].astype(str)
        self._update_marker_variables_from_var_df(var_df)

    @staticmethod
    def _bool_col_to_marker_list(var_df: pd.DataFrame, col: str) -> list[str]:
        """Return marker names where a boolean-like var column is true."""
        if col not in var_df.columns:
            return []

        series = var_df[col]
        truthy = series.map(
            lambda value: False
            if pd.isna(value)
            else (
                bool(value)
                if isinstance(value, (bool, np.bool_))
                else str(value).strip().lower() in {"true", "1", "yes", "y"}
            )
        )
        return var_df.index[truthy].astype(str).tolist()

    def _update_marker_variables(self, adata: anndata.AnnData | None = None) -> None:
        """Populate run marker variables by category for the current AnnData."""
        if adata is None:
            if self._adata is None:
                self.markers_all = []
                self.markers_core_histones = []
                self.markers_identity = []
                self.markers_functional = []
                self.markers_qc = []
                self.markers_intracellular = []
                self.markers_extracellular = []
                self.markers_by_class = {}
                return
            adata = self._adata

        var_df = adata.var.copy()
        var_df.index = var_df.index.astype(str)
        self._update_marker_variables_from_var_df(var_df)

    def _update_marker_variables_from_var_df(self, var_df: pd.DataFrame) -> None:
        """Populate run marker variables from an AnnData var-like DataFrame."""
        var_df = var_df.copy()
        var_df.index = var_df.index.astype(str)

        self.markers_all = var_df.index.tolist()
        self.markers_core_histones = self._bool_col_to_marker_list(
            var_df, "is_core_histone"
        )
        self.markers_identity = self._bool_col_to_marker_list(
            var_df, "is_identity_marker"
        )
        self.markers_functional = self._bool_col_to_marker_list(
            var_df, "is_functional_marker"
        )
        self.markers_qc = self._bool_col_to_marker_list(var_df, "is_qc_marker")

        self.markers_intracellular = self._bool_col_to_marker_list(
            var_df, "is_intracellular_marker"
        )
        self.markers_extracellular = self._bool_col_to_marker_list(
            var_df, "is_extracellular_marker"
        )

        # Prefer explicit intra/extra annotation from marker registry when present.
        if "intra_extra" in var_df.columns:
            intra_extra_series = (
                var_df["intra_extra"]
                .astype("string")
                .fillna("")
                .str.strip()
                .str.lower()
            )
            intra_from_registry = (
                var_df.index[intra_extra_series.str.startswith("intra")]
                .astype(str)
                .tolist()
            )
            extra_from_registry = (
                var_df.index[intra_extra_series.str.startswith("extra")]
                .astype(str)
                .tolist()
            )

            if intra_from_registry:
                self.markers_intracellular = intra_from_registry
            if extra_from_registry:
                self.markers_extracellular = extra_from_registry

        marker_classes: dict[str, list[str]] = {}
        if "marker_class" in var_df.columns:
            class_series = (
                var_df["marker_class"].astype("string").fillna("unknown").astype(str)
            )
            for class_name in sorted(class_series.unique().tolist()):
                marker_classes[class_name] = var_df.index[
                    class_series == class_name
                ].tolist()

            # Heuristic fallbacks when explicit intra/extra flags are unavailable.
            if not self.markers_intracellular:
                intra_classes = {
                    "histone",
                    "functional",
                    "ptm",
                    "intracellular",
                    "nuclear",
                }
                self.markers_intracellular = [
                    marker
                    for marker, class_name in class_series.items()
                    if class_name.strip().lower() in intra_classes
                ]
            if not self.markers_extracellular:
                extra_classes = {"identity", "extracellular", "surface"}
                self.markers_extracellular = [
                    marker
                    for marker, class_name in class_series.items()
                    if class_name.strip().lower() in extra_classes
                ]

        self.markers_by_class = marker_classes

    @staticmethod
    def _matrix_from_layer(adata: anndata.AnnData, layer: str) -> np.ndarray:
        """Return dense expression matrix for selected layer."""
        # Support common alias/typo used in notebooks.
        if layer == "normlized":
            if "normlized" in adata.layers:
                layer = "normlized"
            elif "normalized" in adata.layers:
                layer = "normalized"

        if layer in {"X", "x"}:
            matrix = adata.X
        else:
            if layer not in adata.layers:
                raise ValueError(
                    f"Layer '{layer}' not found. Available layers: {list(adata.layers.keys())}"
                )
            matrix = adata.layers[layer]

        if hasattr(matrix, "toarray"):
            matrix = matrix.toarray()
        return np.asarray(matrix, dtype=np.float32)

    def set_x_from_layer(self, layer: str, inplace: bool = True) -> None:
        """Overwrite `adata.X` with values from a selected layer.

        This operation is intentionally direct: no backup/history is created.

        Args:
            layer: Source layer name (or "X").
            inplace: If True, persist the updated AnnData to run zarr.
        """
        adata = self.read_adata()
        adata.X = self._matrix_from_layer(adata, layer).copy()
        self._adata = adata
        self._log_run_event(
            "x_set_from_layer",
            {
                "source_layer": layer,
                "n_cells": int(adata.n_obs),
                "n_markers": int(adata.n_vars),
                "inplace": bool(inplace),
            },
        )

        if inplace:
            self.save()

    @staticmethod
    def _load_mlx_umap_modules(module_name: str = "mlx_umap"):
        """Import mlx_umap main module and its knn helper submodule."""
        try:
            umap_module = importlib.import_module(module_name)
            knn_module = importlib.import_module(f"{module_name}._knn")
        except ImportError as exc:
            raise ImportError(
                "compute_umap requires mlx-umap (module 'mlx_umap'). "
                "Install/make it importable and retry."
            ) from exc

        if not hasattr(umap_module, "UMAP"):
            raise ImportError(f"{module_name} does not expose UMAP class")
        if not hasattr(knn_module, "compute_knn"):
            raise ImportError(f"{module_name}._knn does not expose compute_knn")

        return umap_module, knn_module

    @staticmethod
    def _knn_to_graph_arrays(knn_idx: np.ndarray, knn_dist: np.ndarray, n_nodes: int):
        """Convert knn arrays to undirected weighted graph arrays.

        Implements the same conversion logic used in the MetalUMAP notebook:
        - remove self-edges
        - distance -> similarity via exp(-d/sigma)
        - undirected edge dedup with max weight
        """
        rows = np.repeat(np.arange(n_nodes, dtype=np.int32), knn_idx.shape[1])
        cols = knn_idx.reshape(-1).astype(np.int32)
        d = knn_dist.reshape(-1).astype(np.float32)

        mask = rows != cols
        rows = rows[mask]
        cols = cols[mask]
        d = d[mask]

        sigma = float(np.median(d[d > 0])) if np.any(d > 0) else 1.0
        w = np.exp(-d / (sigma + 1e-8)).astype(np.float32)

        u = np.minimum(rows, cols)
        v = np.maximum(rows, cols)
        keys = u.astype(np.int64) * n_nodes + v.astype(np.int64)
        uniq, inv = np.unique(keys, return_inverse=True)

        w_max = np.zeros(uniq.shape[0], dtype=np.float32)
        np.maximum.at(w_max, inv, w)

        d_min = np.full(uniq.shape[0], np.inf, dtype=np.float32)
        np.minimum.at(d_min, inv, d)

        uu = (uniq // n_nodes).astype(np.int32)
        vv = (uniq % n_nodes).astype(np.int32)

        return uu, vv, w_max, d_min, sigma

    @staticmethod
    def _upsert_nested_uns_dict(adata: anndata.AnnData, key: str) -> dict:
        """Return mutable nested dict under adata.uns[key], creating if needed."""
        if key not in adata.uns or not isinstance(adata.uns[key], dict):
            adata.uns[key] = {}
        return adata.uns[key]

    def _clear_embedding_artifacts(
        self, adata: anndata.AnnData, embedding_name: str
    ) -> None:
        """Remove existing embedding artifacts for overwrite semantics."""
        embeddings = self._upsert_nested_uns_dict(adata, "embeddings")
        prev = embeddings.get(embedding_name, {})
        if not isinstance(prev, dict):
            prev = {}

        for container_name, key_field in [
            ("obsm", "embedding_key"),
            ("obsm", "knn_indices_key"),
            ("obsm", "knn_distances_key"),
            ("obsp", "connectivities_key"),
            ("obsp", "distances_key"),
        ]:
            artifact_key = prev.get(key_field)
            if artifact_key and artifact_key in getattr(adata, container_name):
                del getattr(adata, container_name)[artifact_key]

        if embedding_name in embeddings:
            del embeddings[embedding_name]

    def compute_umap(
        self,
        markers: list[str],
        source_layer: str = "X",
        embedding_name: str = "X_umap",
        n_neighbors: int = 15,
        n_components: int = 2,
        min_dist: float = 0.1,
        metric: str = "euclidean",
        random_state: int = 42,
        knn_method: str = "brute",
        verbose: bool = False,
        module_name: str = "mlx_umap",
        inplace: bool = True,
    ) -> dict[str, Any]:
        """Compute UMAP embedding from selected markers using mlx-umap.

        Also computes/stores KNN and graph artifacts required for downstream Leiden.
        Existing artifacts with the same `embedding_name` are overwritten.
        """
        import inspect as _inspect
        import time as _time

        if not markers:
            raise ValueError("markers cannot be empty")
        if n_neighbors <= 1:
            raise ValueError("n_neighbors must be > 1")
        if n_components <= 0:
            raise ValueError("n_components must be > 0")

        adata = self.read_adata()
        umap_module, knn_module = self._load_mlx_umap_modules(module_name)

        matrix = self._matrix_from_layer(adata, source_layer)
        var_names = adata.var_names.astype(str).tolist()
        marker_to_idx = {name: i for i, name in enumerate(var_names)}

        missing = [marker for marker in markers if marker not in marker_to_idx]
        if missing:
            raise ValueError(f"Markers not found for embedding: {missing}")

        marker_idx = [marker_to_idx[marker] for marker in markers]
        x_embed = np.asarray(matrix[:, marker_idx], dtype=np.float32)

        self._clear_embedding_artifacts(adata, embedding_name)

        if verbose:
            print(
                f"[compute_umap] embedding='{embedding_name}' layer='{source_layer}' "
                f"cells={adata.n_obs} markers={len(markers)}"
            )

        umap_kwargs = {
            "n_components": n_components,
            "n_neighbors": n_neighbors,
            "min_dist": min_dist,
            "metric": metric,
            "random_state": random_state,
            "knn_method": knn_method,
        }
        try:
            params = _inspect.signature(umap_module.UMAP).parameters
            if "verbose" in params:
                umap_kwargs["verbose"] = verbose
        except (TypeError, ValueError):
            pass

        t0 = _time.perf_counter()
        umap_model = umap_module.UMAP(**umap_kwargs)
        embedding = np.asarray(umap_model.fit_transform(x_embed), dtype=np.float32)
        umap_sec = float(_time.perf_counter() - t0)

        t1 = _time.perf_counter()
        knn_idx, knn_dist = knn_module.compute_knn(
            x_embed,
            n_neighbors,
            method=knn_method,
            random_state=random_state,
            verbose=verbose,
        )
        knn_idx = np.asarray(knn_idx, dtype=np.int32)
        knn_dist = np.asarray(knn_dist, dtype=np.float32)
        knn_sec = float(_time.perf_counter() - t1)

        uu, vv, weights, dists, sigma = self._knn_to_graph_arrays(
            knn_idx,
            knn_dist,
            adata.n_obs,
        )

        try:
            from scipy import sparse as sp
        except ImportError as exc:
            raise ImportError(
                "compute_umap requires scipy for sparse graph storage"
            ) from exc

        edge_u = np.concatenate([uu, vv])
        edge_v = np.concatenate([vv, uu])
        conn_vals = np.concatenate([weights, weights]).astype(np.float32)
        dist_vals = np.concatenate([dists, dists]).astype(np.float32)

        connectivities = sp.csr_matrix(
            (conn_vals, (edge_u, edge_v)),
            shape=(adata.n_obs, adata.n_obs),
            dtype=np.float32,
        )
        distances = sp.csr_matrix(
            (dist_vals, (edge_u, edge_v)),
            shape=(adata.n_obs, adata.n_obs),
            dtype=np.float32,
        )

        embedding_key = embedding_name
        knn_indices_key = f"{embedding_name}_knn_indices"
        knn_distances_key = f"{embedding_name}_knn_distances"
        connectivities_key = f"{embedding_name}_connectivities"
        distances_key = f"{embedding_name}_distances"

        adata.obsm[embedding_key] = embedding
        adata.obsm[knn_indices_key] = knn_idx
        adata.obsm[knn_distances_key] = knn_dist
        adata.obsp[connectivities_key] = connectivities
        adata.obsp[distances_key] = distances

        embeddings_uns = self._upsert_nested_uns_dict(adata, "embeddings")
        metadata = {
            "timestamp": datetime.utcnow().isoformat(),
            "method": "mlx_umap",
            "module": module_name,
            "verbose": bool(verbose),
            "source_layer": source_layer,
            "markers": list(markers),
            "n_neighbors": int(n_neighbors),
            "n_components": int(n_components),
            "min_dist": float(min_dist),
            "metric": str(metric),
            "random_state": int(random_state),
            "knn_method": str(knn_method),
            "sigma": float(sigma),
            "embedding_key": embedding_key,
            "knn_indices_key": knn_indices_key,
            "knn_distances_key": knn_distances_key,
            "connectivities_key": connectivities_key,
            "distances_key": distances_key,
            "n_cells": int(adata.n_obs),
            "n_markers": int(len(markers)),
            "umap_sec": umap_sec,
            "knn_sec": knn_sec,
        }
        embeddings_uns[embedding_name] = metadata

        if verbose:
            print(
                f"[compute_umap] done embedding='{embedding_name}' "
                f"umap_sec={umap_sec:.3f} knn_sec={knn_sec:.3f}"
            )

        self._adata = adata
        self._log_run_event(
            "umap_computed",
            {
                "embedding_name": embedding_name,
                "source_layer": source_layer,
                "n_neighbors": int(n_neighbors),
                "n_components": int(n_components),
                "min_dist": float(min_dist),
                "metric": str(metric),
                "knn_method": str(knn_method),
                "n_markers": int(len(markers)),
                "n_cells": int(adata.n_obs),
            },
        )
        if inplace:
            self.save()

        return metadata

    def cluster_leiden(
        self,
        embedding_name: str,
        cluster_key: str | None = None,
        resolution: float = 1.0,
        n_iterations: int = 2,
        beta: float = 0.01,
        objective_function: str = "modularity",
        seed: int = 42,
        verbose: bool = False,
        inplace: bool = True,
    ) -> dict[str, Any]:
        """Cluster cells with Leiden using graph artifacts from an embedding.

        Follows the Leiden settings used in the MetalUMAP notebook.
        Existing clusterings with the same key are overwritten.
        """
        import time as _time

        adata = self.read_adata()
        embeddings_uns = self._upsert_nested_uns_dict(adata, "embeddings")

        if embedding_name not in embeddings_uns or not isinstance(
            embeddings_uns[embedding_name], dict
        ):
            raise ValueError(
                f"Embedding '{embedding_name}' not found. Run compute_umap first."
            )
        embed_meta = embeddings_uns[embedding_name]

        connectivities_key = embed_meta.get(
            "connectivities_key", f"{embedding_name}_connectivities"
        )
        if connectivities_key not in adata.obsp:
            raise ValueError(
                f"Connectivity graph '{connectivities_key}' not found in adata.obsp"
            )

        cluster_key = cluster_key or f"{embedding_name}_leiden"

        try:
            import igraph as ig
        except ImportError as exc:
            raise ImportError(
                "cluster_leiden requires python-igraph. Install with: pip install igraph"
            ) from exc

        try:
            from scipy import sparse as sp
        except ImportError as exc:
            raise ImportError(
                "cluster_leiden requires scipy for sparse graph handling"
            ) from exc

        conn = adata.obsp[connectivities_key]
        if sp.issparse(conn):
            tri = sp.triu(conn, k=1).tocoo()
            rows = tri.row.astype(np.int32)
            cols = tri.col.astype(np.int32)
            weights = tri.data.astype(np.float32)
        else:
            dense = np.asarray(conn, dtype=np.float32)
            rows, cols = np.where(np.triu(dense, k=1) > 0)
            weights = dense[rows, cols].astype(np.float32)

        edges = list(zip(rows.tolist(), cols.tolist()))
        graph = ig.Graph(n=adata.n_obs, edges=edges, directed=False)
        graph.es["weight"] = weights.tolist()

        if verbose:
            print(
                f"[cluster_leiden] embedding='{embedding_name}' "
                f"cluster_key='{cluster_key}' edges={len(edges)} resolution={resolution}"
            )

        t0 = _time.perf_counter()
        part = graph.community_leiden(
            objective_function=objective_function,
            weights="weight",
            resolution=resolution,
            n_iterations=n_iterations,
            beta=beta,
        )
        leiden_sec = float(_time.perf_counter() - t0)
        labels = np.asarray(part.membership, dtype=np.int32)
        adata.obs[cluster_key] = pd.Categorical(labels.astype(str))

        n_clusters = int(len(set(labels.tolist())))
        clusterings_uns = self._upsert_nested_uns_dict(adata, "clusterings")
        metadata = {
            "timestamp": datetime.utcnow().isoformat(),
            "method": "igraph_leiden",
            "embedding_name": embedding_name,
            "connectivities_key": connectivities_key,
            "cluster_key": cluster_key,
            "resolution": float(resolution),
            "n_iterations": int(n_iterations),
            "beta": float(beta),
            "objective_function": str(objective_function),
            "seed": int(seed),
            "verbose": bool(verbose),
            "n_cells": int(adata.n_obs),
            "n_clusters": n_clusters,
            "leiden_sec": leiden_sec,
        }
        clusterings_uns[cluster_key] = metadata

        if verbose:
            print(
                f"[cluster_leiden] done cluster_key='{cluster_key}' "
                f"n_clusters={n_clusters} leiden_sec={leiden_sec:.3f}"
            )

        self._adata = adata
        self._log_run_event(
            "leiden_clustered",
            {
                "embedding_name": embedding_name,
                "cluster_key": cluster_key,
                "resolution": float(resolution),
                "n_iterations": int(n_iterations),
                "beta": float(beta),
                "objective_function": str(objective_function),
                "n_clusters": int(n_clusters),
                "n_cells": int(adata.n_obs),
            },
        )
        if inplace:
            self.save()

        return metadata

    def _resolve_permcell_module(
        self,
        permcell_module: Any | None = None,
        module_name: str = "PermCell_Smooth",
    ):
        """Resolve PermCell module without importing it.

        Resolution order:
        1) explicit module object (`permcell_module`)
        2) already-imported module in `sys.modules[module_name]`
        """
        if permcell_module is not None:
            module = permcell_module
        else:
            module = sys.modules.get(module_name)
            if module is None:
                raise ImportError(
                    "run_permcell requires a pre-imported PermCell module. "
                    "Import it first (or pass permcell_module=...) and retry."
                )

        required = ["gaussian_smooth_all_torch", "sipsic_like_scores_v3"]
        missing = [name for name in required if not hasattr(module, name)]
        if missing:
            raise ImportError(
                f"PermCell module is missing required functions: {missing}"
            )
        return module

    def run_permcell(
        self,
        signatures: dict[str, object],
        source_layer: str = "X",
        positions_key: str = "X_umap",
        result_prefix: str = "permcell",
        smoothed_key: str | None = None,
        compute_unsmoothed: bool = True,
        # smoothing params
        bandwidth: float = -1,
        k: int | None = 64,
        radius: float | None = None,
        chunk_size: int = 2048,
        device: str | None = None,
        # scoring params
        n_perm: int = 2000,
        seed: int = 0,
        exclude_set: bool = True,
        two_sided: bool = False,
        abs_variant: bool = True,
        exact_max_combinations: int = 50_000,
        progress: bool = True,
        normalize_set_weights: str | None = None,
        use_sparse_W: bool = False,
        prefer_permutation: bool = True,
        perm_batch: int = 1024,
        # module resolution (no import is performed)
        permcell_module: Any | None = None,
        module_name: str = "PermCell_Smooth",
        module_path: str | None = None,
        inplace: bool = True,
    ) -> dict[str, pd.DataFrame]:
        """Run PermCell smoothing + scoring and store results in AnnData.

        The permanent outputs are stored in obsm/uns (no layers are modified).
        Existing outputs with the same `result_prefix` are overwritten.

        Returns:
            Dict of result DataFrames keyed by: z, p, zdir (and zabs/pabs when enabled).
        """
        if not signatures:
            raise ValueError("signatures cannot be empty")

        adata = self.read_adata()
        if positions_key not in adata.obsm:
            raise ValueError(
                f"positions_key '{positions_key}' not found in adata.obsm."
            )

        positions = np.asarray(adata.obsm[positions_key], dtype=np.float32)
        if positions.ndim != 2 or positions.shape[1] < 2:
            raise ValueError(
                f"positions_key '{positions_key}' must have shape (n_cells, >=2)."
            )
        positions_2d = positions[:, :2]

        matrix = self._matrix_from_layer(adata, source_layer)
        # Resolve module inline to avoid dependence on helper availability
        # in long-lived notebook objects.
        if permcell_module is not None:
            module = permcell_module
        else:
            module = sys.modules.get(module_name)
            if module is None:
                raise ImportError(
                    "run_permcell requires a pre-imported PermCell module. "
                    "Import it first (or pass permcell_module=...) and retry."
                )

        required = ["gaussian_smooth_all_torch", "sipsic_like_scores_v3"]
        missing = [name for name in required if not hasattr(module, name)]
        if missing:
            raise ImportError(
                f"PermCell module is missing required functions: {missing}"
            )

        if smoothed_key is None:
            smoothed_key = f"{result_prefix}_smoothed"

        smoothed = module.gaussian_smooth_all_torch(
            data=np.asarray(matrix, dtype=np.float32),
            positions=positions_2d,
            bandwidth=bandwidth,
            k=k,
            radius=radius,
            chunk_size=chunk_size,
            device=device,
        )
        smoothed = np.asarray(smoothed, dtype=np.float32)
        if smoothed.shape != matrix.shape:
            raise ValueError(
                "PermCell smoothing output shape mismatch: "
                f"expected {matrix.shape}, got {smoothed.shape}"
            )

        adata.obsm[smoothed_key] = smoothed.astype(np.float32)

        def _score_matrix(input_matrix: np.ndarray) -> dict[str, pd.DataFrame]:
            tmp = anndata.AnnData(
                X=np.asarray(input_matrix, dtype=np.float32),
                obs=adata.obs.copy(),
                var=adata.var.copy(),
            )
            scores = module.sipsic_like_scores_v3(
                adata=tmp,
                marker_sets=signatures,
                layer="X",
                n_perm=n_perm,
                seed=seed,
                exclude_set=exclude_set,
                two_sided=two_sided,
                abs_variant=abs_variant,
                exact_max_combinations=exact_max_combinations,
                progress=progress,
                normalize_set_weights=normalize_set_weights,
                use_sparse_W=use_sparse_W,
                prefer_permutation=prefer_permutation,
                perm_batch=perm_batch,
            )

            if abs_variant:
                z_df, p_df, zabs_df, pabs_df, zdir_df = scores
                return {
                    "z": z_df,
                    "p": p_df,
                    "zabs": zabs_df,
                    "pabs": pabs_df,
                    "zdir": zdir_df,
                }

            z_df, p_df, zdir_df = scores
            return {
                "z": z_df,
                "p": p_df,
                "zdir": zdir_df,
            }

        smoothed_scores = _score_matrix(smoothed)
        raw_scores = _score_matrix(matrix) if compute_unsmoothed else None

        signature_names = [str(name) for name in smoothed_scores["z"].columns.tolist()]

        def _to_obsm(df: pd.DataFrame, key: str):
            adata.obsm[key] = df.loc[adata.obs_names].to_numpy(
                dtype=np.float32, copy=True
            )

        def _store_score_block(
            prefix: str, score_block: dict[str, pd.DataFrame]
        ) -> dict[str, str | None]:
            z_key = f"{prefix}_z"
            p_key = f"{prefix}_p"
            zdir_key = f"{prefix}_zdir"
            zabs_key = f"{prefix}_zabs"
            pabs_key = f"{prefix}_pabs"

            _to_obsm(score_block["z"], z_key)
            _to_obsm(score_block["p"], p_key)
            _to_obsm(score_block["zdir"], zdir_key)

            if "zabs" in score_block:
                _to_obsm(score_block["zabs"], zabs_key)
            elif zabs_key in adata.obsm:
                del adata.obsm[zabs_key]

            if "pabs" in score_block:
                _to_obsm(score_block["pabs"], pabs_key)
            elif pabs_key in adata.obsm:
                del adata.obsm[pabs_key]

            return {
                "z": z_key,
                "p": p_key,
                "zdir": zdir_key,
                "zabs": zabs_key if "zabs" in score_block else None,
                "pabs": pabs_key if "pabs" in score_block else None,
            }

        result_keys_smoothed = _store_score_block(result_prefix, smoothed_scores)
        result_keys_raw = (
            _store_score_block(f"{result_prefix}_raw", raw_scores)
            if raw_scores is not None
            else None
        )

        permcell_uns = self._upsert_nested_uns_dict(adata, "permcell")
        permcell_uns[result_prefix] = {
            "timestamp": datetime.utcnow().isoformat(),
            "module_name": getattr(module, "__name__", module_name),
            "module_path": getattr(module, "__file__", None),
            "requested_module_path": module_path,
            "source_layer": source_layer,
            "positions_key": positions_key,
            "smoothed_key": smoothed_key,
            "result_prefix": result_prefix,
            "signature_names": signature_names,
            "result_keys_smoothed": result_keys_smoothed,
            "result_keys_raw": result_keys_raw,
            "params": {
                "bandwidth": float(bandwidth),
                "k": None if k is None else int(k),
                "radius": None if radius is None else float(radius),
                "chunk_size": int(chunk_size),
                "device": device,
                "n_perm": int(n_perm),
                "seed": int(seed),
                "exclude_set": bool(exclude_set),
                "two_sided": bool(two_sided),
                "abs_variant": bool(abs_variant),
                "exact_max_combinations": int(exact_max_combinations),
                "progress": bool(progress),
                "normalize_set_weights": normalize_set_weights,
                "use_sparse_W": bool(use_sparse_W),
                "prefer_permutation": bool(prefer_permutation),
                "perm_batch": int(perm_batch),
                "compute_unsmoothed": bool(compute_unsmoothed),
            },
            "n_cells": int(adata.n_obs),
            "n_markers": int(adata.n_vars),
            "n_signatures": int(len(signature_names)),
        }

        self._adata = adata
        self._log_run_event(
            "permcell_scored",
            {
                "result_prefix": result_prefix,
                "source_layer": source_layer,
                "positions_key": positions_key,
                "smoothed_key": smoothed_key,
                "compute_unsmoothed": bool(compute_unsmoothed),
                "n_signatures": int(len(signature_names)),
                "n_cells": int(adata.n_obs),
            },
        )
        if inplace:
            self.save()

        result = {
            "z": smoothed_scores["z"],
            "p": smoothed_scores["p"],
            "zdir": smoothed_scores["zdir"],
        }
        if "zabs" in smoothed_scores:
            result["zabs"] = smoothed_scores["zabs"]
        if "pabs" in smoothed_scores:
            result["pabs"] = smoothed_scores["pabs"]
        if raw_scores is not None:
            result["raw_z"] = raw_scores["z"]
            result["raw_p"] = raw_scores["p"]
            result["raw_zdir"] = raw_scores["zdir"]
            if "zabs" in raw_scores:
                result["raw_zabs"] = raw_scores["zabs"]
            if "pabs" in raw_scores:
                result["raw_pabs"] = raw_scores["pabs"]
        return result

    def permcell_to_adata(
        self,
        result_prefix: str = "permcell",
        score: str = "z",
        smoothed: bool = True,
        embedding_key: str | None = None,
    ) -> anndata.AnnData:
        """Build an AnnData view of PermCell results.

        Args:
            result_prefix: Prefix used in `run_permcell`.
            score: One of: z, p, zdir, zabs, pabs.
            smoothed: If True use smoothed PermCell results, else raw results.
            embedding_key: Optional embedding to copy into returned `.obsm`.
                If None, uses the run's stored positions_key for that PermCell run.

        Returns:
            AnnData with:
              - X: selected PermCell score matrix
              - var_names: signature names
              - obs: copied run obs
              - obsm: selected embedding (if available)
        """
        adata = self.read_adata()
        permcell_uns = adata.uns.get("permcell", {})
        if not isinstance(permcell_uns, dict) or result_prefix not in permcell_uns:
            raise ValueError(
                f"PermCell result_prefix '{result_prefix}' not found in adata.uns['permcell']"
            )

        meta = permcell_uns[result_prefix]
        if not isinstance(meta, dict):
            raise ValueError(f"Invalid PermCell metadata for prefix '{result_prefix}'")

        key_map_name = "result_keys_smoothed" if smoothed else "result_keys_raw"
        key_map = meta.get(key_map_name)
        if key_map is None:
            if smoothed and "result_keys" in meta:
                key_map = meta.get("result_keys")
            else:
                raise ValueError(
                    f"No {'smoothed' if smoothed else 'raw'} PermCell keys for prefix '{result_prefix}'"
                )

        if not isinstance(key_map, dict):
            raise ValueError(
                f"Invalid key map in PermCell metadata for '{result_prefix}'"
            )

        if score not in key_map:
            raise ValueError(
                f"score '{score}' not available. Available: {sorted(key_map.keys())}"
            )

        score_key = key_map.get(score)
        if score_key is None or score_key not in adata.obsm:
            raise ValueError(
                f"PermCell score '{score}' key missing in obsm: {score_key}"
            )

        signature_names = meta.get("signature_names", [])
        score_matrix = np.asarray(adata.obsm[score_key], dtype=np.float32)
        if len(signature_names) != score_matrix.shape[1]:
            signature_names = [f"sig_{i}" for i in range(score_matrix.shape[1])]

        out = anndata.AnnData(
            X=score_matrix.copy(),
            obs=adata.obs.copy(),
            var=pd.DataFrame(index=np.asarray(signature_names, dtype=str)),
        )

        if embedding_key is None:
            embedding_key = meta.get("positions_key")
        if embedding_key and embedding_key in adata.obsm:
            out.obsm[embedding_key] = np.asarray(
                adata.obsm[embedding_key], dtype=np.float32
            )

        out.uns["permcell_source"] = {
            "result_prefix": result_prefix,
            "score": score,
            "smoothed": bool(smoothed),
            "score_key": score_key,
            "embedding_key": embedding_key,
        }
        return out

    @staticmethod
    def _sanitize_observed_marker_name(name: Any) -> str:
        """Strip common CyTOF metal/isotope channel prefixes/suffixes from marker names."""
        marker = str(name).strip()

        # Examples:
        #   116Cd_H3K9me3 -> H3K9me3
        #   116Cd-H3K9me3 -> H3K9me3
        #   116Cd H3K9me3 -> H3K9me3
        marker = re.sub(r"^\d{2,3}[A-Za-z]{1,3}[\s_\-:.]+", "", marker)

        # Also handle trailing isotope tokens:
        #   H3K9me3_116Cd -> H3K9me3
        marker = re.sub(r"[\s_\-:.]+\d{2,3}[A-Za-z]{1,3}$", "", marker)

        return marker.strip()

    @staticmethod
    def _resolve_gate_bound(bound: Any, values: np.ndarray) -> float | None:
        """Resolve a numeric or percentile gate bound."""
        if bound is None:
            return None

        if isinstance(bound, (int, float, np.number)):
            return float(bound)

        if isinstance(bound, str):
            token = bound.strip().lower()
            if token.startswith("p"):
                percentile_str = token[1:]
                if not percentile_str:
                    raise ValueError(
                        f"Invalid percentile bound '{bound}'. Use format like 'p1' or 'p99.5'."
                    )
                percentile = float(percentile_str)
                if percentile < 0 or percentile > 100:
                    raise ValueError(
                        f"Percentile bound '{bound}' must be between 0 and 100."
                    )
                return float(np.percentile(values, percentile))

        raise ValueError(
            f"Unsupported gate bound '{bound}'. Use numeric values or percentile strings like 'p5'."
        )

    @staticmethod
    def _parse_gate_spec(gate_spec: Any) -> tuple[Any, Any]:
        """Parse gate specification into (lower, upper)."""
        if isinstance(gate_spec, dict):
            return gate_spec.get("lower"), gate_spec.get("upper")

        if isinstance(gate_spec, (list, tuple)) and len(gate_spec) == 2:
            return gate_spec[0], gate_spec[1]

        raise ValueError(
            "Each marker gate must be a dict with 'lower'/'upper' keys or a 2-item tuple/list."
        )

    def plot_marker_histograms(
        self,
        markers: list[str],
        layer: str = "X",
        cofactor: float = 5.0,
        fill: bool = False,
        stat: str = "density",
        element: str = "step",
        bins: int | str = "auto",
    ):
        """Plot per-sample histograms for selected markers.

        Args:
            markers: Marker names to plot (must exist in `adata.var_names`).
            layer: Expression layer to plot from (`X` or a layer key).
            cofactor: Cofactor used in arcsinh transform.
            fill: Passed to seaborn.histplot.
            stat: Passed to seaborn.histplot.
            element: Passed to seaborn.histplot.
            bins: Passed to seaborn.histplot.

        Returns:
            Tuple of (figure, axes) from matplotlib.
        """
        if cofactor <= 0:
            raise ValueError("cofactor must be > 0")
        if not markers:
            raise ValueError("markers cannot be empty")

        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
        except ImportError as exc:
            raise ImportError(
                "plot_marker_histograms requires matplotlib and seaborn. "
                "Install with: pip install matplotlib seaborn"
            ) from exc

        adata = self.read_adata()
        matrix = self._matrix_from_layer(adata, layer)
        var_names = adata.var_names.astype(str).tolist()
        marker_to_index = {name: idx for idx, name in enumerate(var_names)}

        for marker in markers:
            if marker not in marker_to_index:
                raise ValueError(
                    f"Marker '{marker}' not found. Available markers include: {var_names[:10]}"
                )

        if "sample_id" in adata.obs.columns:
            sample_ids = adata.obs["sample_id"].astype(str)
        else:
            sample_ids = pd.Series(["all_cells"] * adata.n_obs, index=adata.obs.index)

        samples = sorted(sample_ids.unique().tolist())
        n_rows = len(samples)
        n_cols = len(markers)

        fig, axes = plt.subplots(
            n_rows,
            n_cols,
            figsize=(4.0 * n_cols, 2.8 * n_rows),
            squeeze=False,
        )

        for row_idx, sample in enumerate(samples):
            sample_mask = sample_ids.values == sample
            for col_idx, marker in enumerate(markers):
                ax = axes[row_idx, col_idx]
                marker_values = matrix[:, marker_to_index[marker]][sample_mask]
                transformed = np.arcsinh(marker_values / cofactor)
                sns.histplot(
                    x=transformed,
                    fill=fill,
                    stat=stat,
                    element=element,
                    bins=bins,
                    ax=ax,
                )
                ax.set_title(f"{sample} | {marker}")
                ax.set_xlabel(f"arcsinh({marker}/{cofactor})")
                ax.set_ylabel(stat)

        fig.tight_layout()
        return fig, axes

    def _resolve_field_values(
        self,
        adata: anndata.AnnData,
        field: str,
        layer: str,
    ) -> tuple[np.ndarray, str]:
        """Return per-cell numeric values for a marker or numeric obs column.

        Args:
            adata: AnnData object.
            field: Marker name (in `adata.var_names`) or numeric obs column.
            layer: Expression layer used when `field` is a marker.

        Returns:
            (values, kind) where kind is "marker" or "obs".
        """
        var_names = adata.var_names.astype(str).tolist()
        if field in var_names:
            matrix = self._matrix_from_layer(adata, layer)
            idx = var_names.index(field)
            return np.asarray(matrix[:, idx], dtype=np.float64), "marker"

        if field in adata.obs.columns:
            series = adata.obs[field]
            values = pd.to_numeric(series, errors="coerce").to_numpy(dtype=np.float64)
            if np.all(np.isnan(values)):
                raise ValueError(
                    f"obs column '{field}' is not numeric and cannot be aggregated or plotted."
                )
            return values, "obs"

        raise ValueError(
            f"Field '{field}' not found in adata.var_names or adata.obs columns."
        )

    @staticmethod
    def _resolve_groupby_series(
        adata: anndata.AnnData,
        groupby,
    ) -> pd.Series:
        """Resolve `groupby` (str or list[str]) to a per-cell group label series."""
        if isinstance(groupby, str):
            if groupby not in adata.obs.columns:
                raise ValueError(
                    f"groupby '{groupby}' not in obs columns: "
                    f"{list(adata.obs.columns)}"
                )
            return adata.obs[groupby].astype(str).rename(groupby)

        if isinstance(groupby, (list, tuple)):
            cols = list(groupby)
            for col in cols:
                if col not in adata.obs.columns:
                    raise ValueError(
                        f"groupby column '{col}' not in obs columns: "
                        f"{list(adata.obs.columns)}"
                    )
            joined = adata.obs[cols].astype(str).agg(" | ".join, axis=1)
            return joined.rename(" | ".join(cols))

        raise TypeError("groupby must be a str or a list/tuple of str.")

    def plot_heatmap(
        self,
        fields: list[str],
        groupby,
        layer: str = "X",
        agg: str = "mean",
        standard_scale: str | None = None,
        cmap: str | None = None,
        center: float | None = None,
        annot: bool = False,
        fmt: str = ".2f",
        order: list[str] | None = None,
        figsize: tuple[float, float] | None = None,
        ax=None,
        heatmap_kwargs: dict | None = None,
    ):
        """Plot a heatmap of aggregated field values grouped by an obs column.

        Fields can be markers (in `adata.var_names`) or numeric `adata.obs`
        columns. `groupby` can be a single obs column or a list of obs columns
        for composite grouping (e.g. `["line_id", "condition"]`).

        Args:
            fields: Markers or numeric obs columns to include as heatmap rows.
            groupby: Single obs column name or a list of obs column names.
            layer: Layer used when a field is a marker (default `X`).
            agg: `"mean"` or `"median"`.
            standard_scale: `None`, `"row"`, or `"column"`. Applies z-score over
                rows (fields) or columns (groups) after aggregation.
            cmap: Matplotlib colormap. Defaults to `viridis` (raw) or `RdBu_r`
                (when `standard_scale` is set).
            center: Value at which to center the colormap. Defaults to 0 when
                `standard_scale` is set, otherwise None.
            annot: If True, write the aggregated value in each cell.
            fmt: Format string used when `annot=True`.
            order: Optional explicit group order along the x-axis.
            figsize: Optional figure size.
            ax: Optional matplotlib axes to draw into.
            heatmap_kwargs: Extra keyword arguments forwarded to
                `seaborn.heatmap` (e.g. `linewidths`, `linecolor`,
                `cbar_kws`, `vmin`, `vmax`). Explicit arguments above take
                precedence over keys in this dict.

        Returns:
            Tuple of (figure, axes).
        """
        if not fields:
            raise ValueError("fields cannot be empty")
        if agg not in {"mean", "median"}:
            raise ValueError("agg must be 'mean' or 'median'")
        if standard_scale not in {None, "row", "column"}:
            raise ValueError("standard_scale must be None, 'row', or 'column'")

        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
        except ImportError as exc:
            raise ImportError(
                "plot_heatmap requires matplotlib and seaborn. "
                "Install with: pip install matplotlib seaborn"
            ) from exc

        adata = self.read_adata()
        group_series = self._resolve_groupby_series(adata, groupby).reset_index(drop=True)

        columns = {}
        for f in fields:
            values, _ = self._resolve_field_values(adata, f, layer)
            columns[f] = values

        long_df = pd.DataFrame(columns)
        long_df["__group__"] = group_series.values

        grouped = long_df.groupby("__group__", observed=True)[fields]
        table = grouped.mean() if agg == "mean" else grouped.median()

        if order is not None:
            missing = [g for g in order if g not in table.index]
            if missing:
                raise ValueError(f"order contains unknown groups: {missing}")
            table = table.loc[order]
        else:
            table = table.sort_index()

        # Rows = fields (markers/obs), columns = groups
        matrix = table.T

        scaled = matrix.astype(float).copy()
        if standard_scale == "row":
            row_mean = scaled.mean(axis=1)
            row_std = scaled.std(axis=1).replace(0, np.nan)
            scaled = scaled.sub(row_mean, axis=0).div(row_std, axis=0).fillna(0.0)
        elif standard_scale == "column":
            col_mean = scaled.mean(axis=0)
            col_std = scaled.std(axis=0).replace(0, np.nan)
            scaled = scaled.sub(col_mean, axis=1).div(col_std, axis=1).fillna(0.0)

        if cmap is None:
            cmap = "RdBu_r" if standard_scale is not None else "viridis"
        if center is None and standard_scale is not None:
            center = 0.0

        if figsize is None:
            figsize = (
                max(4.0, 0.6 * scaled.shape[1] + 2.0),
                max(3.0, 0.4 * scaled.shape[0] + 2.0),
            )

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.figure

        sns.heatmap(
            scaled,
            ax=ax,
            cmap=cmap,
            center=center,
            annot=annot,
            fmt=fmt,
            cbar=True,
            **{
                k: v for k, v in (heatmap_kwargs or {}).items()
                if k not in {"data", "ax", "cmap", "center", "annot", "fmt"}
            },
        )

        group_label = (
            groupby if isinstance(groupby, str) else " | ".join(groupby)
        )
        scale_suffix = f" ({standard_scale}-z)" if standard_scale else ""
        ax.set_xlabel(group_label)
        ax.set_ylabel("field")
        ax.set_title(f"{agg}{scale_suffix} by {group_label}")
        fig.tight_layout()
        return fig, ax

    def plot_boxplot(
        self,
        field: str,
        groupby,
        layer: str = "X",
        order: list[str] | None = None,
        comparisons=None,
        test: str = "mannwhitney",
        multitest: str | None = "holm",
        show_points: bool = False,
        max_points: int = 2000,
        point_alpha: float = 0.4,
        point_size: float = 2.0,
        palette=None,
        figsize: tuple[float, float] | None = None,
        ax=None,
        bracket_color: str = "black",
        bracket_linewidth: float = 1.0,
        bracket_fontsize: float = 11.0,
        ns_label: str = "ns",
        random_state: int = 0,
        boxplot_kwargs: dict | None = None,
        stripplot_kwargs: dict | None = None,
    ):
        """Boxplot for a single field grouped by an obs column, with significance brackets.

        The field can be a marker (in `adata.var_names`) or a numeric obs
        column. `groupby` can be a single obs column or a list of obs columns
        for composite grouping.

        Args:
            field: Marker name or numeric obs column to plot.
            groupby: Single obs column or list of obs columns.
            layer: Layer used when `field` is a marker (default `X`).
            order: Optional explicit group order along the x-axis.
            comparisons: Pairs to test. Accepts:
                - `None` or `"all"`: all unordered pairs.
                - `"adjacent"`: only neighbours in `order`.
                - List of `(group_a, group_b)` tuples for explicit pairs.
            test: `"mannwhitney"` (default), `"ttest"`, or `"welch"`.
            multitest: `None`, `"holm"`, or `"bonferroni"`.
            show_points: Overlay a stripplot of per-cell points (subsampled).
            max_points: Maximum number of stripplot points (subsampled).
            point_alpha: Alpha for overlaid points.
            point_size: Size for overlaid points.
            palette: Seaborn palette name or color list for the boxes.
            figsize: Optional figure size.
            ax: Optional matplotlib axes to draw into.
            bracket_color: Color used for significance brackets.
            bracket_linewidth: Line width used for significance brackets.
            bracket_fontsize: Font size used for significance labels.
            ns_label: Label used for non-significant comparisons.
            random_state: Seed used when subsampling points.
            boxplot_kwargs: Extra keyword arguments forwarded to
                `seaborn.boxplot` (e.g. `width`, `linewidth`, `notch`,
                `whis`, `saturation`). Explicit arguments above take
                precedence over keys in this dict.
            stripplot_kwargs: Extra keyword arguments forwarded to
                `seaborn.stripplot` when `show_points=True` (e.g. `jitter`,
                `dodge`).

        Returns:
            Tuple of (figure, axes).
        """
        if test not in {"mannwhitney", "ttest", "welch"}:
            raise ValueError("test must be 'mannwhitney', 'ttest', or 'welch'")
        if multitest not in {None, "holm", "bonferroni"}:
            raise ValueError("multitest must be None, 'holm', or 'bonferroni'")

        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
        except ImportError as exc:
            raise ImportError(
                "plot_boxplot requires matplotlib and seaborn. "
                "Install with: pip install matplotlib seaborn"
            ) from exc

        adata = self.read_adata()
        values, _ = self._resolve_field_values(adata, field, layer)
        group_series = self._resolve_groupby_series(adata, groupby).reset_index(drop=True)

        plot_df = pd.DataFrame(
            {
                field: values,
                "__group__": group_series.values,
            }
        ).dropna(subset=[field, "__group__"])

        if order is None:
            order = sorted(plot_df["__group__"].unique().tolist())
        else:
            missing = [g for g in order if g not in plot_df["__group__"].unique()]
            if missing:
                raise ValueError(f"order contains unknown groups: {missing}")

        plot_df = plot_df[plot_df["__group__"].isin(order)]
        plot_df["__group__"] = pd.Categorical(
            plot_df["__group__"], categories=order, ordered=True
        )

        if figsize is None:
            figsize = (max(4.0, 0.7 * len(order) + 2.0), 4.5)

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.figure

        sns.boxplot(
            data=plot_df,
            x="__group__",
            y=field,
            order=order,
            ax=ax,
            palette=palette,
            showfliers=not show_points,
            **{
                k: v for k, v in (boxplot_kwargs or {}).items()
                if k not in {
                    "data", "x", "y", "order", "ax", "palette", "showfliers"
                }
            },
        )

        if show_points and len(plot_df) > 0:
            if len(plot_df) > max_points:
                sampled = plot_df.sample(n=max_points, random_state=random_state)
            else:
                sampled = plot_df
            sns.stripplot(
                data=sampled,
                x="__group__",
                y=field,
                order=order,
                color="black",
                alpha=point_alpha,
                size=point_size,
                ax=ax,
                **{
                    k: v for k, v in (stripplot_kwargs or {}).items()
                    if k not in {
                        "data", "x", "y", "order", "ax",
                        "color", "alpha", "size",
                    }
                },
            )

        group_label = (
            groupby if isinstance(groupby, str) else " | ".join(groupby)
        )
        ax.set_xlabel(group_label)
        ax.set_ylabel(field)
        ax.set_title(f"{field} by {group_label}")

        pairs = self._build_comparison_pairs(order, comparisons)
        if pairs:
            values_by_group = {
                g: plot_df.loc[plot_df["__group__"] == g, field].to_numpy()
                for g in order
            }
            pvalues = self._compute_pairwise_pvalues(
                values_by_group, pairs, test=test
            )
            if multitest is not None:
                pvalues = self._adjust_pvalues(pvalues, method=multitest)

            self._draw_significance_brackets(
                ax,
                pairs=pairs,
                group_positions={g: i for i, g in enumerate(order)},
                pvalues=pvalues,
                data_values=plot_df[field].to_numpy(),
                color=bracket_color,
                linewidth=bracket_linewidth,
                fontsize=bracket_fontsize,
                ns_label=ns_label,
            )

        fig.tight_layout()
        return fig, ax

    @staticmethod
    def _build_comparison_pairs(order, comparisons):
        """Return list of (group_a, group_b) pairs to test for significance."""
        if comparisons is None or comparisons == "all":
            return [
                (order[i], order[j])
                for i in range(len(order))
                for j in range(i + 1, len(order))
            ]
        if comparisons == "adjacent":
            return [(order[i], order[i + 1]) for i in range(len(order) - 1)]

        pairs = []
        for pair in comparisons:
            if len(pair) != 2:
                raise ValueError(
                    "Each comparison must be a (group_a, group_b) pair."
                )
            a, b = pair
            if a not in order or b not in order:
                raise ValueError(
                    f"comparison ({a}, {b}) contains unknown groups."
                )
            pairs.append((a, b))
        return pairs

    @staticmethod
    def _compute_pairwise_pvalues(values_by_group, pairs, test):
        """Compute p-values for the requested pair tests using scipy.stats."""
        try:
            from scipy import stats
        except ImportError as exc:
            raise ImportError(
                "Significance tests require scipy. Install with: pip install scipy"
            ) from exc

        pvals = []
        for a, b in pairs:
            x = np.asarray(values_by_group.get(a, np.array([])), dtype=float)
            y = np.asarray(values_by_group.get(b, np.array([])), dtype=float)
            x = x[~np.isnan(x)]
            y = y[~np.isnan(y)]
            if len(x) < 2 or len(y) < 2:
                pvals.append(float("nan"))
                continue
            if test == "mannwhitney":
                _, p = stats.mannwhitneyu(x, y, alternative="two-sided")
            elif test == "ttest":
                _, p = stats.ttest_ind(x, y, equal_var=True)
            elif test == "welch":
                _, p = stats.ttest_ind(x, y, equal_var=False)
            else:
                raise ValueError(f"Unknown test: {test}")
            pvals.append(float(p))
        return pvals

    @staticmethod
    def _adjust_pvalues(pvalues, method="holm"):
        """Apply Holm or Bonferroni correction across pairwise p-values."""
        pvals = np.asarray(pvalues, dtype=float)
        valid_mask = ~np.isnan(pvals)
        n_valid = int(valid_mask.sum())
        if n_valid == 0:
            return pvals.tolist()

        if method == "bonferroni":
            adjusted = pvals.copy()
            adjusted[valid_mask] = np.minimum(pvals[valid_mask] * n_valid, 1.0)
            return adjusted.tolist()

        if method == "holm":
            adjusted = pvals.copy()
            valid_idx = np.where(valid_mask)[0]
            order_local = np.argsort(pvals[valid_idx])
            running_max = 0.0
            for rank, local_idx in enumerate(order_local):
                orig = valid_idx[local_idx]
                value = pvals[orig] * (n_valid - rank)
                value = min(value, 1.0)
                running_max = max(running_max, value)
                adjusted[orig] = running_max
            return adjusted.tolist()

        raise ValueError(f"Unknown method: {method}")

    @staticmethod
    def _pvalue_to_stars(p, ns_label="ns"):
        """Convert a p-value to a star annotation."""
        if p is None:
            return ns_label
        try:
            p_f = float(p)
        except (TypeError, ValueError):
            return ns_label
        if np.isnan(p_f):
            return ns_label
        if p_f <= 1e-4:
            return "****"
        if p_f <= 1e-3:
            return "***"
        if p_f <= 1e-2:
            return "**"
        if p_f <= 5e-2:
            return "*"
        return ns_label

    @staticmethod
    def _draw_significance_brackets(
        ax,
        pairs,
        group_positions,
        pvalues,
        data_values,
        color="black",
        linewidth=1.0,
        fontsize=11.0,
        ns_label="ns",
    ):
        """Draw stacked significance brackets with star annotations."""
        if not pairs:
            return

        y_min, y_max = ax.get_ylim()
        finite_values = np.asarray(data_values, dtype=float)
        finite_values = finite_values[np.isfinite(finite_values)]
        if finite_values.size:
            data_max = float(np.max(finite_values))
        else:
            data_max = y_max

        span = data_max - y_min if data_max > y_min else max(abs(data_max), 1.0)
        base = data_max + 0.06 * span
        step = 0.08 * span
        tick = 0.015 * span

        indexed_pairs = sorted(
            enumerate(pairs),
            key=lambda ip: abs(
                group_positions[ip[1][1]] - group_positions[ip[1][0]]
            ),
        )

        occupied: list[tuple[float, float, float]] = []
        for orig_idx, (a, b) in indexed_pairs:
            x1 = group_positions[a]
            x2 = group_positions[b]
            lo, hi = min(x1, x2), max(x1, x2)

            overlap_y = [
                oy for ox_lo, ox_hi, oy in occupied
                if not (hi < ox_lo or lo > ox_hi)
            ]
            y = max(overlap_y) + step if overlap_y else base
            occupied.append((lo, hi, y))

            ax.plot(
                [x1, x1, x2, x2],
                [y, y + tick, y + tick, y],
                lw=linewidth,
                color=color,
                clip_on=False,
            )
            label = Run._pvalue_to_stars(pvalues[orig_idx], ns_label=ns_label)
            ax.text(
                (x1 + x2) / 2.0,
                y + tick,
                label,
                ha="center",
                va="bottom",
                color=color,
                fontsize=fontsize,
                clip_on=False,
            )

        if occupied:
            new_top = max(oy for _, _, oy in occupied) + 3 * step
            ax.set_ylim(y_min, max(y_max, new_top))

    def qc_gate(
        self,
        gates: dict[str, Any],
        layer: str = "X",
        inplace: bool = True,
    ) -> pd.Series:
        """Apply marker QC gates and optionally persist filtered cells.

        Args:
            gates: Dict mapping marker -> gate specification.
                Gate spec can be:
                - {'lower': value, 'upper': value}
                - (lower, upper)
                where each bound is numeric, None, or percentile string (e.g., 'p1').
            layer: Expression layer used for gating.
            inplace: If True, persist filtered AnnData back to run zarr path.

        Returns:
            Boolean pass mask indexed by original `obs_names`.
        """
        adata = self.read_adata()
        matrix = self._matrix_from_layer(adata, layer)

        if not gates:
            pass_mask = pd.Series(True, index=adata.obs_names, name="pass_qc")
            self._log_run_event(
                "qc_gated",
                {
                    "layer": layer,
                    "inplace": bool(inplace),
                    "n_before": int(adata.n_obs),
                    "n_after": int(adata.n_obs),
                    "n_removed": 0,
                    "gates": {},
                },
            )
            return pass_mask

        var_names = adata.var_names.astype(str).tolist()
        marker_to_index = {name: idx for idx, name in enumerate(var_names)}
        keep_mask = np.ones(adata.n_obs, dtype=bool)
        resolved_gates: dict[str, dict[str, float | None]] = {}

        for marker, gate_spec in gates.items():
            if marker not in marker_to_index:
                raise ValueError(
                    f"Marker '{marker}' not found. Available markers include: {var_names[:10]}"
                )

            lower_spec, upper_spec = self._parse_gate_spec(gate_spec)
            marker_values = matrix[:, marker_to_index[marker]].astype(float)

            lower_bound = self._resolve_gate_bound(lower_spec, marker_values)
            upper_bound = self._resolve_gate_bound(upper_spec, marker_values)

            marker_mask = np.ones(adata.n_obs, dtype=bool)
            if lower_bound is not None:
                marker_mask &= marker_values >= lower_bound
            if upper_bound is not None:
                marker_mask &= marker_values <= upper_bound

            keep_mask &= marker_mask
            resolved_gates[marker] = {
                "lower": lower_bound,
                "upper": upper_bound,
            }

        pass_mask = pd.Series(keep_mask, index=adata.obs_names, name="pass_qc")

        if not inplace:
            self._log_run_event(
                "qc_gated",
                {
                    "layer": layer,
                    "inplace": False,
                    "n_before": int(adata.n_obs),
                    "n_after": int(keep_mask.sum()),
                    "n_removed": int((~keep_mask).sum()),
                    "gates": resolved_gates,
                },
            )
            return pass_mask

        # Preserve ungated values for retained cells if raw layer does not exist yet.
        if "raw" not in adata.layers:
            adata.layers["raw"] = self._matrix_from_layer(adata, "X").copy()

        qc_record = {
            "timestamp": datetime.utcnow().isoformat(),
            "layer": layer,
            "gates": resolved_gates,
            "n_before": int(adata.n_obs),
            "n_after": int(keep_mask.sum()),
            "n_removed": int((~keep_mask).sum()),
        }

        if "qc" not in adata.uns or not isinstance(adata.uns["qc"], dict):
            adata.uns["qc"] = {}
        history = adata.uns["qc"].get("history", [])
        if not isinstance(history, list):
            history = []
        history.append(json.dumps(qc_record, sort_keys=True))
        adata.uns["qc"]["history"] = history
        adata.uns["qc"]["latest"] = qc_record

        filtered = adata[keep_mask].copy()
        filtered.write_zarr(self.zarr_path())
        self._adata = filtered
        self._update_marker_variables(filtered)
        self._log_run_event(
            "qc_gated",
            {
                "layer": layer,
                "inplace": True,
                "n_before": int(adata.n_obs),
                "n_after": int(keep_mask.sum()),
                "n_removed": int((~keep_mask).sum()),
                "gates": resolved_gates,
            },
        )

        return pass_mask

    @staticmethod
    def _load_cytof_transform_module(module_name: str = "cytof_transform"):
        """Import and return the external cytof_transform module."""
        try:
            return importlib.import_module(module_name)
        except ImportError as exc:
            raise ImportError(
                "Normalization requires the external 'cytof_transform' package/module. "
                "Install or make it importable, then retry."
            ) from exc

    def normalize_with_cytof_transform(
        self,
        control_markers: list[str],
        markers_to_correct: list[str],
        source_layer: str = "raw",
        corrected_layer: str = "normalized",
        z_layer: str = "normalized_z",
        groupby_col: str = "sample_id",
        input_is_arcsinh: bool = False,
        arcsinh_cofactor: float = 5.0,
        anchor_to_median: bool = True,
        zscore: bool = True,
        module_name: str = "cytof_transform",
        inplace: bool = True,
    ) -> dict[str, Any]:
        """Normalize markers using external cytof_transform, per sample/line group.

        Args:
            control_markers: Core markers (e.g. histones) used to estimate technical factor.
            markers_to_correct: Markers to normalize.
            source_layer: Layer used as normalization input.
            corrected_layer: Output layer for corrected (asinh-space) values.
            z_layer: Output layer for z-scored corrected values.
            groupby_col: Obs column used for per-group normalization (e.g. sample_id/line_id).
            input_is_arcsinh: If True, source layer is already arcsinh-transformed.
            arcsinh_cofactor: Cofactor used when transforming source data with arcsinh.
            anchor_to_median: Passed to cytof_transform config.
            zscore: Passed to cytof_transform config.
            module_name: Module name for importing cytof_transform.
            inplace: If True, persist updates to run zarr.

        Returns:
            Summary dictionary of normalization outputs and settings.
        """
        if not control_markers:
            raise ValueError("control_markers cannot be empty")
        if not markers_to_correct:
            raise ValueError("markers_to_correct cannot be empty")
        if arcsinh_cofactor <= 0:
            raise ValueError("arcsinh_cofactor must be > 0")

        adata = self.read_adata()
        module = self._load_cytof_transform_module(module_name)

        matrix = self._matrix_from_layer(adata, source_layer)
        var_names = adata.var_names.astype(str).tolist()
        asinh_matrix = (
            matrix if input_is_arcsinh else np.arcsinh(matrix / arcsinh_cofactor)
        )
        asinh_df = pd.DataFrame(asinh_matrix, index=adata.obs_names, columns=var_names)

        missing_ctrl = [m for m in control_markers if m not in asinh_df.columns]
        if missing_ctrl:
            raise ValueError(f"Control markers not found in data: {missing_ctrl}")
        missing_norm = [m for m in markers_to_correct if m not in asinh_df.columns]
        if missing_norm:
            raise ValueError(f"markers_to_correct not found in data: {missing_norm}")
        if groupby_col not in adata.obs.columns:
            raise ValueError(
                f"groupby_col '{groupby_col}' is missing from obs. "
                f"Available columns: {adata.obs.columns.tolist()}"
            )

        groups = adata.obs[groupby_col].astype(str)
        unique_groups = sorted(groups.unique().tolist())

        corrected_df = asinh_df.copy()
        z_df = asinh_df.copy()
        tech_factor = pd.Series(np.nan, index=adata.obs_names, dtype=float)
        gamma_by_group: dict[str, dict[str, float]] = {}
        alpha_by_group: dict[str, dict[str, float]] = {}

        for group in unique_groups:
            idx = groups == group
            group_df = asinh_df.loc[idx]
            if group_df.shape[0] < 2:
                raise ValueError(
                    f"Group '{group}' has only {group_df.shape[0]} cells; "
                    "at least 2 cells are required for normalization."
                )

            cfg = module.CytofTransformConfig(
                control_markers=control_markers,
                markers_to_correct=markers_to_correct,
                use_compartments=False,
                n_pcs_for_T=1,
                anchor_to_median=anchor_to_median,
                zscore=zscore,
                line_col=None,
            )
            result = module.cytof_transform_global(group_df, cfg)

            corrected_df.loc[idx, :] = result.corrected.loc[group_df.index, :]
            z_df.loc[idx, :] = result.residuals_z.loc[group_df.index, :]
            tech_factor.loc[group_df.index] = result.tech_factor.loc[
                group_df.index
            ].astype(float)
            gamma_by_group[group] = {k: float(v) for k, v in result.gamma.items()}
            alpha_by_group[group] = {k: float(v) for k, v in result.alpha.items()}

        adata.layers[corrected_layer] = corrected_df.values.astype(np.float32)
        adata.layers[z_layer] = z_df.values.astype(np.float32)
        adata.obs["norm_tech_factor"] = tech_factor.values

        if "normalization" not in adata.uns or not isinstance(
            adata.uns["normalization"], dict
        ):
            adata.uns["normalization"] = {}

        summary = {
            "timestamp": datetime.utcnow().isoformat(),
            "method": "cytof_transform_global",
            "module": module_name,
            "groupby_col": groupby_col,
            "groups": unique_groups,
            "source_layer": source_layer,
            "corrected_layer": corrected_layer,
            "z_layer": z_layer,
            "input_is_arcsinh": bool(input_is_arcsinh),
            "arcsinh_cofactor": float(arcsinh_cofactor),
            "control_markers": list(control_markers),
            "markers_to_correct": list(markers_to_correct),
            "anchor_to_median": bool(anchor_to_median),
            "zscore": bool(zscore),
            "n_before": int(adata.n_obs),
            "n_after": int(adata.n_obs),
            "gamma_by_group": gamma_by_group,
            "alpha_by_group": alpha_by_group,
        }

        history = adata.uns["normalization"].get("history", [])
        if not isinstance(history, list):
            history = []
        history.append(json.dumps(summary, sort_keys=True))
        adata.uns["normalization"]["history"] = history
        adata.uns["normalization"]["latest"] = summary

        if inplace:
            adata.write_zarr(self.zarr_path())
            self._adata = adata

        self._log_run_event(
            "run_normalized",
            {
                "method": "cytof_transform_global",
                "groupby_col": groupby_col,
                "groups": list(unique_groups),
                "source_layer": source_layer,
                "corrected_layer": corrected_layer,
                "z_layer": z_layer,
                "n_cells": int(adata.n_obs),
                "n_markers": int(adata.n_vars),
            },
        )

        return summary

    def plot_normalization_tech_factor_qc(
        self,
        control_markers: list[str],
        layer: str = "raw",
        group_value: str | None = None,
        groupby_col: str = "sample_id",
        input_is_arcsinh: bool = False,
        arcsinh_cofactor: float = 5.0,
        module_name: str = "cytof_transform",
    ):
        """Plot technical-factor QC using cytof_transform.plot_tech_factor_qc."""
        adata = self.read_adata()
        module = self._load_cytof_transform_module(module_name)

        matrix = self._matrix_from_layer(adata, layer)
        var_names = adata.var_names.astype(str).tolist()
        asinh_matrix = (
            matrix if input_is_arcsinh else np.arcsinh(matrix / arcsinh_cofactor)
        )
        asinh_df = pd.DataFrame(asinh_matrix, index=adata.obs_names, columns=var_names)

        if group_value is not None:
            if groupby_col not in adata.obs.columns:
                raise ValueError(
                    f"groupby_col '{groupby_col}' is missing from obs. "
                    f"Available columns: {adata.obs.columns.tolist()}"
                )
            idx = adata.obs[groupby_col].astype(str) == str(group_value)
            asinh_df = asinh_df.loc[idx]

        return module.plot_tech_factor_qc(
            asinh_data=asinh_df,
            control_markers=control_markers,
            tech_factor=None,
            tech_name="tech1",
        )

    def plot_normalization_marker_correlations_qc(
        self,
        pre_layer: str = "raw",
        post_layer: str = "normalized",
        group_value: str | None = None,
        groupby_col: str = "sample_id",
        input_pre_is_arcsinh: bool = False,
        arcsinh_cofactor: float = 5.0,
        top_n: int = 25,
        module_name: str = "cytof_transform",
    ):
        """Plot marker-tech correlations pre/post normalization for a run or group."""
        adata = self.read_adata()
        module = self._load_cytof_transform_module(module_name)

        if "norm_tech_factor" not in adata.obs.columns:
            raise ValueError(
                "norm_tech_factor not found in obs. Run normalize_with_cytof_transform first."
            )

        var_names = adata.var_names.astype(str).tolist()
        pre_matrix = self._matrix_from_layer(adata, pre_layer)
        pre_asinh = (
            pre_matrix
            if input_pre_is_arcsinh
            else np.arcsinh(pre_matrix / arcsinh_cofactor)
        )
        pre_df = pd.DataFrame(pre_asinh, index=adata.obs_names, columns=var_names)

        post_matrix = self._matrix_from_layer(adata, post_layer)
        post_df = pd.DataFrame(post_matrix, index=adata.obs_names, columns=var_names)

        tech_factor = adata.obs["norm_tech_factor"].astype(float)

        if group_value is not None:
            if groupby_col not in adata.obs.columns:
                raise ValueError(
                    f"groupby_col '{groupby_col}' is missing from obs. "
                    f"Available columns: {adata.obs.columns.tolist()}"
                )
            idx = adata.obs[groupby_col].astype(str) == str(group_value)
            pre_df = pre_df.loc[idx]
            post_df = post_df.loc[idx]
            tech_factor = tech_factor.loc[idx]

        return module.plot_marker_correlations_qc(
            asinh_pre=pre_df,
            asinh_post=post_df,
            tech_factor=tech_factor,
            top_n=top_n,
        )

    def plot_normalization_gamma_qc(
        self,
        group_value: str,
        marker_groups: dict[str, list[str]] | None = None,
        module_name: str = "cytof_transform",
    ):
        """Plot gamma values for a normalized group using cytof_transform.plot_gamma_qc."""
        adata = self.read_adata()
        module = self._load_cytof_transform_module(module_name)

        norm_info = adata.uns.get("normalization", {})
        latest = norm_info.get("latest", {}) if isinstance(norm_info, dict) else {}
        gamma_by_group = (
            latest.get("gamma_by_group", {}) if isinstance(latest, dict) else {}
        )

        if str(group_value) not in gamma_by_group:
            raise ValueError(
                f"No gamma values found for group '{group_value}'. "
                "Run normalize_with_cytof_transform first."
            )

        return module.plot_gamma_qc(
            gamma=gamma_by_group[str(group_value)],
            marker_groups=marker_groups,
        )

    @staticmethod
    def _balanced_reference_indices(
        groups: pd.Series,
        random_state: int = 0,
    ) -> tuple[np.ndarray, dict[str, int], int]:
        """Return balanced per-group indices and summary counts."""
        group_series = groups.astype(str)
        counts = group_series.value_counts().sort_index()
        count_dict = {str(group): int(count) for group, count in counts.items()}

        if len(counts) <= 1:
            return (
                np.arange(group_series.shape[0]),
                count_dict,
                int(group_series.shape[0]),
            )

        target_n = int(counts.min())
        rng = np.random.default_rng(random_state)
        selected: list[np.ndarray] = []

        for group in sorted(count_dict.keys()):
            group_idx = np.where(group_series.values == group)[0]
            if len(group_idx) == target_n:
                selected.append(group_idx)
            else:
                selected.append(rng.choice(group_idx, size=target_n, replace=False))

        return np.concatenate(selected), count_dict, target_n

    def zscore_markers_balanced(
        self,
        source_layer: str = "normlized",
        output_layer: str = "zscore",
        groupby_col: str = "sample_id",
        random_state: int = 0,
        inplace: bool = True,
    ) -> dict[str, Any]:
        """Z-score all markers using a balanced subsample across sample IDs.

        If multiple groups exist in `groupby_col`, z-score parameters (mean/std)
        are estimated from an equal-size subsample per group to avoid bias.

        Args:
            source_layer: Layer to z-score (`X` or layer key).
            output_layer: Target layer name for z-scored values.
            groupby_col: Obs column used for balancing (default: sample_id).
            random_state: RNG seed for balanced subsampling.
            inplace: If True, persist updated AnnData to run zarr.

        Returns:
            Summary dictionary with balancing and z-score metadata.
        """
        adata = self.read_adata()
        matrix = self._matrix_from_layer(adata, source_layer)

        if matrix.shape[1] != adata.n_vars:
            raise ValueError("Source matrix width does not match number of markers.")

        if groupby_col in adata.obs.columns:
            groups = adata.obs[groupby_col]
        else:
            groups = pd.Series(["all_cells"] * adata.n_obs, index=adata.obs.index)

        ref_idx, counts, target_n = self._balanced_reference_indices(
            groups,
            random_state=random_state,
        )

        reference = matrix[ref_idx, :]
        mean = reference.mean(axis=0)
        std = reference.std(axis=0)
        std = np.where(std == 0, 1.0, std)

        z = (matrix - mean) / std
        adata.layers[output_layer] = z.astype(np.float32)

        mean_col = f"{output_layer}_mean"
        std_col = f"{output_layer}_std"
        adata.var[mean_col] = mean.astype(np.float32)
        adata.var[std_col] = std.astype(np.float32)

        summary = {
            "timestamp": datetime.utcnow().isoformat(),
            "method": "zscore_markers_balanced",
            "source_layer": source_layer,
            "output_layer": output_layer,
            "groupby_col": groupby_col,
            "group_counts": counts,
            "balanced_target_per_group": int(target_n),
            "n_reference_cells": int(reference.shape[0]),
            "n_cells": int(adata.n_obs),
            "n_markers": int(adata.n_vars),
            "random_state": int(random_state),
        }

        if "zscore" not in adata.uns or not isinstance(adata.uns["zscore"], dict):
            adata.uns["zscore"] = {}
        history = adata.uns["zscore"].get("history", [])
        if not isinstance(history, list):
            history = []
        history.append(json.dumps(summary, sort_keys=True))
        adata.uns["zscore"]["history"] = history
        adata.uns["zscore"]["latest"] = summary

        if inplace:
            adata.write_zarr(self.zarr_path())
            self._adata = adata

        self._log_run_event(
            "run_zscored",
            {
                "source_layer": source_layer,
                "output_layer": output_layer,
                "groupby_col": groupby_col,
                "group_counts": counts,
                "balanced_target_per_group": int(target_n),
                "n_reference_cells": int(reference.shape[0]),
                "n_cells": int(adata.n_obs),
                "n_markers": int(adata.n_vars),
            },
        )

        return summary

    def create_subset_run(
        self,
        new_run_id: str,
        sample_ids: list[str] | None = None,
        line_ids: list[str] | None = None,
        run_name: str | None = None,
        notes: str | None = None,
    ) -> "Run":
        """Create a new run from a subset of samples and/or lines in this run.

        Args:
            new_run_id: New run ID to create.
            sample_ids: Sample IDs to keep.
            line_ids: Line IDs to keep.
            run_name: Optional name for new run.
            notes: Optional notes for new run metadata.

        Returns:
            Newly created and persisted subset Run.
        """
        if not sample_ids and not line_ids:
            raise ValueError("Provide at least one of sample_ids or line_ids.")

        adata = self.read_adata()
        mask = np.ones(adata.n_obs, dtype=bool)

        if sample_ids is not None:
            if "sample_id" not in adata.obs.columns:
                raise ValueError("sample_id column not found in run.obs")
            sample_set = {str(sample) for sample in sample_ids}
            mask &= adata.obs["sample_id"].astype(str).isin(sample_set).values

        if line_ids is not None:
            if "line_id" not in adata.obs.columns:
                raise ValueError("line_id column not found in run.obs")
            line_set = {str(line) for line in line_ids}
            mask &= adata.obs["line_id"].astype(str).isin(line_set).values

        if int(mask.sum()) == 0:
            raise ValueError("Subset selection produced zero cells.")

        subset_adata = adata[mask].copy()
        subset_adata.uns["derived_from"] = {
            "project_id": self.project.project_id,
            "source_run_id": self.run_id,
            "sample_ids": [str(s) for s in sample_ids]
            if sample_ids is not None
            else None,
            "line_ids": [str(l) for l in line_ids] if line_ids is not None else None,
            "created_at": datetime.utcnow().isoformat(),
        }

        new_run = self.project.add_run(
            run_id=new_run_id,
            run_name=run_name
            or f"{self.run_config.get('run_name', self.run_id)} subset",
            panel_id=self.run_config.get("panel_id"),
            acquisition_date=self.run_config.get("acquisition_date"),
            instrument=self.run_config.get("instrument"),
            operator=self.run_config.get("operator"),
            notes=notes,
        )

        new_run._adata = subset_adata
        new_run.save()
        self._log_run_event(
            "subset_run_created",
            {
                "new_run_id": new_run_id,
                "n_cells_source": int(adata.n_obs),
                "n_cells_subset": int(subset_adata.n_obs),
                "sample_ids": [str(s) for s in sample_ids]
                if sample_ids is not None
                else None,
                "line_ids": [str(l) for l in line_ids]
                if line_ids is not None
                else None,
            },
        )
        new_run._log_run_event(
            "subset_run_initialized",
            {
                "source_run_id": self.run_id,
                "n_cells": int(subset_adata.n_obs),
            },
        )
        return new_run

    def zarr_path(self) -> Path:
        """Get the path to the Zarr file.

        Returns:
            Path to Zarr file
        """
        return self.path / "processed" / f"{self.run_id}.zarr"

    def lock_zarr_parts(
        self,
        parts: list[str] | None = None,
        strict: bool = True,
    ) -> list[str]:
        """Make selected Zarr store parts read-only.

        Args:
            parts: Relative paths in the run Zarr store (for example,
                ``["layers/raw", "obs"]``). If None, lock the full store.
            strict: Whether to raise if a requested part does not exist.

        Returns:
            Normalized part paths that were locked. ``"."`` means full store.
        """
        self.require_ingested()
        locked_parts = set_zarr_parts_writable(
            self.zarr_path(),
            parts=parts,
            writable=False,
            strict=strict,
        )
        self._log_run_event(
            "run_zarr_parts_locked",
            {
                "parts": locked_parts,
            },
        )
        return locked_parts

    def unlock_zarr_parts(
        self,
        parts: list[str] | None = None,
        strict: bool = True,
    ) -> list[str]:
        """Make selected Zarr store parts owner-writable again.

        Args:
            parts: Relative paths in the run Zarr store (for example,
                ``["layers/raw", "obs"]``). If None, unlock the full store.
            strict: Whether to raise if a requested part does not exist.

        Returns:
            Normalized part paths that were unlocked. ``"."`` means full store.
        """
        self.require_ingested()
        unlocked_parts = set_zarr_parts_writable(
            self.zarr_path(),
            parts=parts,
            writable=True,
            strict=strict,
        )
        self._log_run_event(
            "run_zarr_parts_unlocked",
            {
                "parts": unlocked_parts,
            },
        )
        return unlocked_parts

    def is_ingested(self) -> bool:
        """Check if the run has been ingested.

        Returns:
            True if run is ingested
        """
        return self.zarr_path().exists()

    def require_ingested(self) -> None:
        """Require that the run has been ingested.

        Raises:
            RunNotIngestedError if run has not been ingested
        """
        if not self.is_ingested():
            raise RunNotIngestedError(
                f"Run {self.run_id} is registered but has no ingested AnnData-Zarr object at: "
                f"{self.zarr_path()}\n\nCall run.ingest(...) first."
            )

    @property
    def status(self) -> str:
        """Get the run status.

        Returns:
            Status string (registered, ingested, failed_ingestion)
        """
        if self.run_id in self.project._run_registry["run_id"].values:
            return self.project._run_registry[
                self.project._run_registry["run_id"] == self.run_id
            ]["status"].iloc[0]
        return "unknown"


__all__ = ["Run"]
