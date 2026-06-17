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
    ZarrLockedError,
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
    get_locked_zarr_parts,
)
from cytofstandard.uuid_utils import generate_cell_uuid
from cytofstandard.validation import (
    validate_file_list,
    validate_sample_metadata,
    validate_marker_consistency,
)
from cytofstandard.io import read_fcs, read_csv, read_parquet as read_parquet_file
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
        common_markers_only: bool = False,
        drop_columns: list[str] | None = None,
        show_marker_coverage: bool = False,
    ) -> None:
        """Ingest files into this run.

        Args:
            files: List of file paths to ingest
            sample_metadata: Path to sample metadata CSV/Parquet file
            copy_raw: Whether to copy raw files to project
            strict_markers: Whether to fail on unknown markers
            allow_extra_markers: Whether to allow extra markers
            common_markers_only: If True, drop markers that are not present in
                all files for this run.
            drop_columns: Column names to remove before marker processing
            show_marker_coverage: If True, display a heatmap after ingestion
                showing which standard markers are present or absent in each
                file.  Only markers that are missing from at least one file are
                shown; if all markers are present in all files the plot is
                skipped and a message is printed.

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
                "common_markers_only": bool(common_markers_only),
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

        # Detect extra columns in the metadata CSV (beyond required + known optional).
        # These will be carried through to adata.obs unchanged.
        _known_sample_cols = {
            "file_name", "sample_id", "line_id",
            "condition", "replicate_id", "batch_id", "barcoding_id", "notes",
        }
        extra_sample_cols = [c for c in sample_df.columns if c not in _known_sample_cols]

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
            elif src_path.suffix.lower() == ".parquet":
                data_df, marker_metadata = read_parquet_file(stored_path)
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

            # Add extra sample metadata columns (any beyond the known set)
            for col in extra_sample_cols:
                cell_metadata[col] = sample_info[col]

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

        if show_marker_coverage:
            self._plot_marker_coverage(all_marker_mappings, sample_df)

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
        non_marker_cols = (
            set(obs_cols) | optional_field_cols | set(csv_meta_cols_list) | set(extra_sample_cols)
        )

        marker_cols = sorted(c for c in combined.columns if c not in non_marker_cols)
        if common_markers_only:
            marker_cols = [c for c in marker_cols if c in common_markers]
            if not marker_cols:
                raise IngestionError(
                    "No common markers found across files after applying common_markers_only=True."
                )
        X = combined[marker_cols].values.astype(np.float32)

        # Create obs DataFrame
        obs = combined[obs_cols].copy()
        for field in optional_fields:
            if field in combined.columns:
                obs[field] = combined[field].values
        for col in extra_sample_cols:
            obs[col] = combined[col].values
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
            "package_version": "0.1.2",
            "strict_markers": strict_markers,
            "common_markers_only": common_markers_only,
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
        self._write_zarr(adata, zarr_path)

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

    def read_adata(self, backed: bool = False, force: bool = False):
        """Read the ingested AnnData object.

        Args:
            backed: Whether to use backed mode (currently ignored)
            force: Re-read from disk even if an in-memory copy is cached.
                Use this after external modifications (e.g. gating via the
                Streamlit app) to pick up changes written by another process.

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

        if self._adata is None or force:
            import anndata

            self._adata = anndata.read_zarr(zarr_path)
            self._update_marker_variables(self._adata)

        return self._adata

    def ingestion_summary(self) -> pd.DataFrame:
        """Return a DataFrame summarising every file and sample ingested into this run.

        Columns always present:

        - ``file_name`` — original filename
        - ``sample_id`` — sample identifier from the metadata sheet
        - ``line_id`` — cell-line identifier
        - ``n_events_file`` — cells read from the file at ingest time
        - ``n_cells_obs`` — cells currently in the zarr (after any QC gating)
        - ``file_type`` — ``"csv"``, ``"fcs"``, ``"parquet"``, …
        - ``file_size_bytes`` — original file size
        - ``ingested_at`` — UTC timestamp of ingestion

        Additional metadata columns from the sample metadata sheet (e.g.
        ``"condition"``, ``"replicate_id"``) are included when available.

        Raises:
            RunNotIngestedError: If the run has not been ingested yet.
        """
        self.require_ingested()

        manifest_path = self.path / "metadata" / "file_manifest.parquet"
        if not manifest_path.exists():
            raise RunNotIngestedError(
                f"Run '{self.run_id}' has no file manifest. "
                "Re-ingest the run to generate one."
            )

        manifest = read_parquet(str(manifest_path)).rename(
            columns={"n_events": "n_events_file"}
        )

        adata = self.read_adata()
        obs = adata.obs

        # Internal / cell-level columns that are never sample-level metadata
        _skip = {
            "cell_uuid", "project_id", "run_id", "sample_id", "line_id",
            "source_file", "source_file_hash", "event_index",
        }

        # Only include columns where every cell within a given file shares the
        # same value — i.e. genuinely sample/file-level metadata (e.g. condition)
        # rather than per-cell computed values (e.g. norm_tech_factor, cluster IDs).
        def _is_file_level(col: str) -> bool:
            return (
                obs.groupby("source_file_hash", observed=True)[col]
                .nunique()
                .le(1)
                .all()
            )

        extra_meta_cols = [
            c for c in obs.columns
            if c not in _skip and _is_file_level(c)
        ]

        # Build one row per file from obs: sample_id, line_id, extra metadata
        # Join key: source_file_hash (obs) ↔ file_hash_sha256 (manifest)
        obs_per_file = (
            obs[["source_file_hash", "sample_id", "line_id"] + extra_meta_cols]
            .drop_duplicates(subset=["source_file_hash"])
            .rename(columns={"source_file_hash": "file_hash_sha256"})
            .reset_index(drop=True)
        )

        # Per-sample cell counts after any QC gating
        cell_counts = (
            obs.groupby("sample_id", observed=True)
            .size()
            .rename("n_cells_obs")
            .reset_index()
        )

        # Merge everything onto the manifest (one row per ingested file)
        summary = manifest.merge(obs_per_file, on="file_hash_sha256", how="left")
        summary = summary.merge(cell_counts, on="sample_id", how="left")

        # Drop internal/redundant columns
        drop_cols = {"run_id", "stored_path", "file_hash_sha256"}
        summary = summary.drop(columns=[c for c in drop_cols if c in summary.columns])

        # Preferred column order
        front = [
            "file_name", "sample_id", "line_id",
            "n_events_file", "n_cells_obs",
            "file_type", "file_size_bytes", "ingested_at",
        ]
        ordered = [c for c in front if c in summary.columns]
        rest = [c for c in summary.columns if c not in ordered]
        return summary[ordered + rest].reset_index(drop=True)

    def subsample_by_group(
        self,
        groupby_col: str,
        n_per_group: int | None = None,
        random_state: int = 0,
        replace: bool = False,
    ) -> anndata.AnnData:
        """Return a balanced subsample AnnData by group.

        Args:
            groupby_col: Obs column to balance on.
            n_per_group: Number of cells per group. If None, uses the
                smallest group size.
            random_state: RNG seed for sampling.
            replace: Sample with replacement if True.

        Returns:
            Subsampled AnnData copy.
        """
        adata = self.read_adata()
        if groupby_col not in adata.obs.columns:
            raise ValueError(
                f"groupby_col '{groupby_col}' not found in obs. "
                f"Available columns: {adata.obs.columns.tolist()}"
            )

        groups = adata.obs[groupby_col]
        idx, _, _ = self._balanced_group_indices(
            groups,
            n_per_group=n_per_group,
            random_state=random_state,
            replace=replace,
        )
        return adata[idx].copy()

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

        self._write_zarr(self._adata)
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
            self._try_save_inplace()

    @staticmethod
    def _load_umap_module():
        """Import standard umap-learn module."""
        try:
            import umap as umap_module
        except ImportError as exc:
            raise ImportError(
                "compute_umap requires umap-learn. "
                "Install it with: pip install umap-learn"
            ) from exc
        if not hasattr(umap_module, "UMAP"):
            raise ImportError("umap module does not expose UMAP class")
        return umap_module

    @staticmethod
    def _compute_knn_sklearn(
        x: np.ndarray, n_neighbors: int, metric: str
    ) -> tuple[np.ndarray, np.ndarray]:
        """Compute KNN using sklearn NearestNeighbors.

        Returns (knn_idx, knn_dist) each of shape (n_cells, effective_k),
        excluding self-neighbors. effective_k is clamped to n_cells - 1.
        """
        try:
            from sklearn.neighbors import NearestNeighbors
        except ImportError as exc:
            raise ImportError(
                "compute_umap requires scikit-learn for KNN. "
                "Install it with: pip install scikit-learn"
            ) from exc
        n_cells = x.shape[0]
        effective_k = min(n_neighbors, n_cells - 1)
        nn = NearestNeighbors(n_neighbors=effective_k + 1, metric=metric, algorithm="auto")
        nn.fit(x)
        dists, indices = nn.kneighbors(x)
        return indices[:, 1:].astype(np.int32), dists[:, 1:].astype(np.float32)

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
        subsample_size: int | None = 50_000,
        chunk_size: int = 50_000,
        verbose: bool = False,
        inplace: bool = True,
    ) -> dict[str, Any]:
        """Compute UMAP embedding from selected markers using umap-learn.

        Fits on a random subsample of up to `subsample_size` cells, then
        transforms all cells in chunks of `chunk_size`. Set `subsample_size=None`
        to fit on the full dataset. Also computes/stores KNN and graph artifacts
        required for downstream Leiden. Existing artifacts with the same
        `embedding_name` are overwritten.
        """
        import time as _time

        if not markers:
            raise ValueError("markers cannot be empty")
        if n_neighbors <= 1:
            raise ValueError("n_neighbors must be > 1")
        if n_components <= 0:
            raise ValueError("n_components must be > 0")

        adata = self.read_adata()
        umap_module = self._load_umap_module()

        matrix = self._matrix_from_layer(adata, source_layer)
        var_names = adata.var_names.astype(str).tolist()
        marker_to_idx = {name: i for i, name in enumerate(var_names)}

        missing = [marker for marker in markers if marker not in marker_to_idx]
        if missing:
            raise ValueError(f"Markers not found for embedding: {missing}")

        marker_idx = [marker_to_idx[marker] for marker in markers]
        x_embed = np.asarray(matrix[:, marker_idx], dtype=np.float32)
        n_cells = x_embed.shape[0]

        self._clear_embedding_artifacts(adata, embedding_name)

        if verbose:
            print(
                f"[compute_umap] embedding='{embedding_name}' layer='{source_layer}' "
                f"cells={n_cells} markers={len(markers)}"
            )

        umap_model = umap_module.UMAP(
            n_components=n_components,
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            metric=metric,
            random_state=random_state,
            verbose=verbose,
        )

        t0 = _time.perf_counter()
        if subsample_size is not None and subsample_size < n_cells:
            rng = np.random.default_rng(random_state)
            fit_idx = rng.choice(n_cells, size=subsample_size, replace=False)
            if verbose:
                print(
                    f"[compute_umap] fit_transform on {subsample_size} cells, "
                    f"transform rest ({n_cells - subsample_size}) in chunks of {chunk_size}"
                )
            embedding = np.empty((n_cells, n_components), dtype=np.float32)
            embedding[fit_idx] = umap_model.fit_transform(x_embed[fit_idx])
            rest_idx = np.setdiff1d(np.arange(n_cells), fit_idx)
            for start in range(0, len(rest_idx), chunk_size):
                chunk = rest_idx[start : start + chunk_size]
                embedding[chunk] = umap_model.transform(x_embed[chunk])
        else:
            embedding = umap_model.fit_transform(x_embed).astype(np.float32)
            subsample_size = None
        umap_sec = float(_time.perf_counter() - t0)

        t1 = _time.perf_counter()
        knn_idx, knn_dist = self._compute_knn_sklearn(x_embed, n_neighbors, metric)
        knn_sec = float(_time.perf_counter() - t1)

        uu, vv, weights, dists, sigma = self._knn_to_graph_arrays(
            knn_idx,
            knn_dist,
            n_cells,
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
            shape=(n_cells, n_cells),
            dtype=np.float32,
        )
        distances = sp.csr_matrix(
            (dist_vals, (edge_u, edge_v)),
            shape=(n_cells, n_cells),
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
            "method": "umap",
            "verbose": bool(verbose),
            "source_layer": source_layer,
            "markers": list(markers),
            "n_neighbors": int(n_neighbors),
            "n_components": int(n_components),
            "min_dist": float(min_dist),
            "metric": str(metric),
            "random_state": int(random_state),
            "subsample_size": int(subsample_size) if subsample_size is not None else None,
            "chunk_size": int(chunk_size),
            "sigma": float(sigma),
            "embedding_key": embedding_key,
            "knn_indices_key": knn_indices_key,
            "knn_distances_key": knn_distances_key,
            "connectivities_key": connectivities_key,
            "distances_key": distances_key,
            "n_cells": int(n_cells),
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
                "subsample_size": int(subsample_size) if subsample_size is not None else None,
                "n_markers": int(len(markers)),
                "n_cells": int(n_cells),
            },
        )
        if inplace:
            self._try_save_inplace()

        return metadata

    def compute_umap_balanced(
        self,
        markers: list[str],
        source_layer: str = "X",
        embedding_name: str = "X_umap",
        groupby_col: str = "sample_id",
        n_per_group: int | None = None,
        replace: bool = False,
        n_neighbors: int = 15,
        n_components: int = 2,
        min_dist: float = 0.1,
        metric: str = "euclidean",
        random_state: int = 42,
        chunk_size: int = 50_000,
        verbose: bool = False,
        inplace: bool = True,
    ) -> dict[str, Any]:
        """Compute UMAP with fit on balanced subsample, then transform all cells.

        The balanced subsample is defined by `groupby_col`, using `n_per_group`
        cells per group (or the smallest group size if None). The UMAP model
        is fit on the subsample, then all cells are transformed in chunks of
        `chunk_size`. KNN and graph artifacts are computed on the full dataset
        to support downstream Leiden clustering.
        """
        import time as _time

        if not markers:
            raise ValueError("markers cannot be empty")
        if n_neighbors <= 1:
            raise ValueError("n_neighbors must be > 1")
        if n_components <= 0:
            raise ValueError("n_components must be > 0")

        adata = self.read_adata()
        if groupby_col not in adata.obs.columns:
            raise ValueError(
                f"groupby_col '{groupby_col}' not found in obs. "
                f"Available columns: {adata.obs.columns.tolist()}"
            )

        umap_module = self._load_umap_module()

        matrix = self._matrix_from_layer(adata, source_layer)
        var_names = adata.var_names.astype(str).tolist()
        marker_to_idx = {name: i for i, name in enumerate(var_names)}

        missing = [marker for marker in markers if marker not in marker_to_idx]
        if missing:
            raise ValueError(f"Markers not found for embedding: {missing}")

        marker_idx = [marker_to_idx[marker] for marker in markers]
        x_embed = np.asarray(matrix[:, marker_idx], dtype=np.float32)
        n_cells = x_embed.shape[0]

        self._clear_embedding_artifacts(adata, embedding_name)

        group_series = adata.obs[groupby_col]
        fit_idx, counts, target_n = self._balanced_group_indices(
            group_series,
            n_per_group=n_per_group,
            random_state=random_state,
            replace=replace,
        )

        if verbose:
            print(
                f"[compute_umap_balanced] embedding='{embedding_name}' layer='{source_layer}' "
                f"cells={n_cells} markers={len(markers)} groups={len(counts)} "
                f"n_per_group={target_n}"
            )

        umap_model = umap_module.UMAP(
            n_components=n_components,
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            metric=metric,
            random_state=random_state,
            verbose=verbose,
        )

        t0 = _time.perf_counter()
        embedding = np.empty((n_cells, n_components), dtype=np.float32)
        embedding[fit_idx] = umap_model.fit_transform(x_embed[fit_idx])
        rest_idx = np.setdiff1d(np.arange(n_cells), fit_idx)
        for start in range(0, len(rest_idx), chunk_size):
            chunk = rest_idx[start : start + chunk_size]
            embedding[chunk] = umap_model.transform(x_embed[chunk])
        umap_sec = float(_time.perf_counter() - t0)

        t1 = _time.perf_counter()
        knn_idx, knn_dist = self._compute_knn_sklearn(x_embed, n_neighbors, metric)
        knn_sec = float(_time.perf_counter() - t1)

        uu, vv, weights, dists, sigma = self._knn_to_graph_arrays(
            knn_idx,
            knn_dist,
            n_cells,
        )

        try:
            from scipy import sparse as sp
        except ImportError as exc:
            raise ImportError(
                "compute_umap_balanced requires scipy for sparse graph storage"
            ) from exc

        edge_u = np.concatenate([uu, vv])
        edge_v = np.concatenate([vv, uu])
        conn_vals = np.concatenate([weights, weights]).astype(np.float32)
        dist_vals = np.concatenate([dists, dists]).astype(np.float32)

        connectivities = sp.csr_matrix(
            (conn_vals, (edge_u, edge_v)),
            shape=(n_cells, n_cells),
            dtype=np.float32,
        )
        distances = sp.csr_matrix(
            (dist_vals, (edge_u, edge_v)),
            shape=(n_cells, n_cells),
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
            "method": "umap_balanced",
            "verbose": bool(verbose),
            "source_layer": source_layer,
            "markers": list(markers),
            "n_neighbors": int(n_neighbors),
            "n_components": int(n_components),
            "min_dist": float(min_dist),
            "metric": str(metric),
            "random_state": int(random_state),
            "chunk_size": int(chunk_size),
            "sigma": float(sigma),
            "embedding_key": embedding_key,
            "knn_indices_key": knn_indices_key,
            "knn_distances_key": knn_distances_key,
            "connectivities_key": connectivities_key,
            "distances_key": distances_key,
            "n_cells": int(n_cells),
            "n_markers": int(len(markers)),
            "umap_sec": umap_sec,
            "knn_sec": knn_sec,
            "balanced_groupby_col": groupby_col,
            "balanced_n_per_group": int(target_n),
            "balanced_group_counts": counts,
            "balanced_n_groups": int(len(counts)),
            "balanced_n_cells": int(len(fit_idx)),
        }
        embeddings_uns[embedding_name] = metadata

        if verbose:
            print(
                f"[compute_umap_balanced] done embedding='{embedding_name}' "
                f"umap_sec={umap_sec:.3f} knn_sec={knn_sec:.3f}"
            )

        self._adata = adata
        self._log_run_event(
            "umap_balanced_computed",
            {
                "embedding_name": embedding_name,
                "source_layer": source_layer,
                "groupby_col": groupby_col,
                "n_per_group": int(target_n),
                "n_groups": int(len(counts)),
                "n_markers": int(len(markers)),
                "n_cells": int(n_cells),
            },
        )
        if inplace:
            self._try_save_inplace()

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
            self._try_save_inplace()

        return metadata

    def cluster_dbscan(
        self,
        embedding_name: str,
        cluster_key: str | None = None,
        eps: float = 0.5,
        min_samples: int = 5,
        metric: str = "euclidean",
        verbose: bool = False,
        inplace: bool = True,
    ) -> dict[str, Any]:
        """Cluster cells with DBSCAN on the coordinates of a stored embedding.

        Operates directly on adata.obsm[embedding_name] (e.g. 2-D UMAP coordinates).
        Noise points (DBSCAN label -1) are retained as the category ``"-1"``.
        Existing clusterings with the same key are overwritten.

        Args:
            embedding_name: Key in adata.obsm to use as input coordinates.
            cluster_key: Column name written to adata.obs. Defaults to
                ``"{embedding_name}_dbscan"``.
            eps: Maximum distance between two samples to be considered neighbours.
            min_samples: Minimum neighbours for a point to be a core point.
            metric: Distance metric passed to sklearn DBSCAN.
            verbose: Print progress and result summary.
            inplace: Persist to zarr immediately after clustering.

        Returns:
            Metadata dict with parameters, timing, n_clusters, and n_noise.
        """
        import time as _time

        try:
            from sklearn.cluster import DBSCAN
        except ImportError as exc:
            raise ImportError(
                "cluster_dbscan requires scikit-learn. Install with: pip install scikit-learn"
            ) from exc

        adata = self.read_adata()

        if embedding_name not in adata.obsm:
            raise ValueError(
                f"Embedding '{embedding_name}' not found in adata.obsm. "
                f"Available: {list(adata.obsm.keys())}"
            )

        cluster_key = cluster_key or f"{embedding_name}_dbscan"
        coords = np.asarray(adata.obsm[embedding_name], dtype=np.float32)

        if verbose:
            print(
                f"[cluster_dbscan] embedding='{embedding_name}' "
                f"cluster_key='{cluster_key}' eps={eps} min_samples={min_samples} "
                f"metric='{metric}' cells={adata.n_obs}"
            )

        t0 = _time.perf_counter()
        labels = DBSCAN(eps=eps, min_samples=min_samples, metric=metric).fit_predict(coords)
        dbscan_sec = float(_time.perf_counter() - t0)

        n_noise = int(np.sum(labels == -1))
        n_clusters = int(len(set(labels.tolist())) - (1 if -1 in labels else 0))

        adata.obs[cluster_key] = pd.Categorical(labels.astype(str))

        clusterings_uns = self._upsert_nested_uns_dict(adata, "clusterings")
        metadata = {
            "timestamp": datetime.utcnow().isoformat(),
            "method": "dbscan",
            "embedding_name": embedding_name,
            "cluster_key": cluster_key,
            "eps": float(eps),
            "min_samples": int(min_samples),
            "metric": str(metric),
            "verbose": bool(verbose),
            "n_cells": int(adata.n_obs),
            "n_clusters": n_clusters,
            "n_noise": n_noise,
            "dbscan_sec": dbscan_sec,
        }
        clusterings_uns[cluster_key] = metadata

        if verbose:
            print(
                f"[cluster_dbscan] done cluster_key='{cluster_key}' "
                f"n_clusters={n_clusters} n_noise={n_noise} dbscan_sec={dbscan_sec:.3f}"
            )

        self._adata = adata
        self._log_run_event(
            "dbscan_clustered",
            {
                "embedding_name": embedding_name,
                "cluster_key": cluster_key,
                "eps": float(eps),
                "min_samples": int(min_samples),
                "metric": str(metric),
                "n_clusters": n_clusters,
                "n_noise": n_noise,
                "n_cells": int(adata.n_obs),
            },
        )
        if inplace:
            self._try_save_inplace()

        return metadata

    @staticmethod
    def _knn_to_jaccard_graph_arrays(
        knn_idx: np.ndarray,
        n_nodes: int,
        min_jaccard: float = 0.0,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Build Jaccard-weighted undirected graph arrays from KNN indices.

        For each directed KNN edge (i -> j), computes the Jaccard similarity
        between the neighbor sets of i and j.  Self-loops are excluded.
        The result is symmetrized by taking the max Jaccard for each pair.

        Args:
            knn_idx: (n_nodes, k) integer array of neighbor indices.
            n_nodes: Total number of nodes.
            min_jaccard: Edges with Jaccard < this value are dropped.

        Returns:
            (rows, cols, jaccard_weights) — upper-triangular undirected edges.
        """
        knn_idx = np.asarray(knn_idx, dtype=np.int32)
        k = knn_idx.shape[1]

        # Build neighbor sets as frozensets (excluding self).
        neighbor_sets = []
        for i in range(n_nodes):
            nbrs = set(int(x) for x in knn_idx[i] if int(x) != i)
            neighbor_sets.append(nbrs)

        # Iterate over all directed edges once; accumulate Jaccard per
        # undirected pair using a dict keyed by (lo, hi).
        jaccard_dict: dict[tuple[int, int], float] = {}
        for i in range(n_nodes):
            for j_raw in knn_idx[i]:
                j = int(j_raw)
                if j == i:
                    continue
                lo, hi = (i, j) if i < j else (j, i)
                if (lo, hi) in jaccard_dict:
                    continue  # already computed
                si = neighbor_sets[i] | {i}
                sj = neighbor_sets[j] | {j}
                inter = len(si & sj)
                union = len(si | sj)
                jac = inter / union if union > 0 else 0.0
                if jac >= min_jaccard:
                    jaccard_dict[(lo, hi)] = jac

        if not jaccard_dict:
            empty = np.empty(0, dtype=np.int32)
            return empty, empty, np.empty(0, dtype=np.float32)

        pairs = list(jaccard_dict.keys())
        rows = np.array([p[0] for p in pairs], dtype=np.int32)
        cols = np.array([p[1] for p in pairs], dtype=np.int32)
        weights = np.array([jaccard_dict[p] for p in pairs], dtype=np.float32)
        return rows, cols, weights

    def cluster_leiden_jaccard(
        self,
        embedding_name: str,
        cluster_key: str | None = None,
        jaccard_connectivities_key: str | None = None,
        min_jaccard: float = 0.0,
        resolution: float = 1.0,
        n_iterations: int = 2,
        beta: float = 0.01,
        objective_function: str = "modularity",
        seed: int = 42,
        verbose: bool = False,
        inplace: bool = True,
    ) -> dict[str, Any]:
        """Cluster cells using a PhenoGraph-style Jaccard graph with Leiden.

        Equivalent to PhenoGraph but uses Leiden instead of Louvain:
            1. KNN graph (pre-computed by `compute_umap`)
            2. Edge weights replaced by Jaccard similarity of neighbor sets
            3. Leiden community detection on the Jaccard-weighted graph

        The KNN indices are read from the `obsm` artifacts stored by
        `compute_umap`.  No recomputation of neighbors is performed.

        Args:
            embedding_name: Name of the embedding whose KNN artifacts to use
                (must have been computed via `compute_umap`).
            cluster_key: Key under which cluster labels are stored in
                `adata.obs`. Defaults to
                `f"{embedding_name}_jaccard_leiden"`.
            jaccard_connectivities_key: Key used to store the Jaccard
                connectivity matrix in `adata.obsp`. Defaults to
                `f"{embedding_name}_jaccard_connectivities"`.
            min_jaccard: Edges with Jaccard similarity below this threshold
                are pruned before clustering.  Useful to remove very weak
                connections (default 0.0 keeps all edges).
            resolution: Leiden resolution parameter.
            n_iterations: Number of Leiden iterations.
            beta: Leiden randomness parameter.
            objective_function: `"modularity"` or `"CPM"`.
            seed: Random seed passed to Leiden.
            verbose: If True, print progress to stdout.
            inplace: If True, persist updated AnnData to Zarr after
                clustering.

        Returns:
            Metadata dict stored under `adata.uns["clusterings"][cluster_key]`.
        """
        import time as _time

        adata = self.read_adata()
        embeddings_uns = self._upsert_nested_uns_dict(adata, "embeddings")

        if embedding_name not in embeddings_uns or not isinstance(
            embeddings_uns[embedding_name], dict
        ):
            raise ValueError(
                f"Embedding '{embedding_name}' not found in uns. "
                "Run compute_umap first."
            )
        embed_meta = embeddings_uns[embedding_name]

        knn_indices_key = embed_meta.get(
            "knn_indices_key", f"{embedding_name}_knn_indices"
        )
        if knn_indices_key not in adata.obsm:
            raise ValueError(
                f"KNN indices '{knn_indices_key}' not found in adata.obsm. "
                "Re-run compute_umap to regenerate KNN artifacts."
            )

        cluster_key = cluster_key or f"{embedding_name}_jaccard_leiden"
        jaccard_connectivities_key = (
            jaccard_connectivities_key
            or f"{embedding_name}_jaccard_connectivities"
        )

        try:
            import igraph as ig
        except ImportError as exc:
            raise ImportError(
                "cluster_leiden_jaccard requires python-igraph. "
                "Install with: pip install igraph"
            ) from exc

        try:
            from scipy import sparse as sp
        except ImportError as exc:
            raise ImportError(
                "cluster_leiden_jaccard requires scipy. "
                "Install with: pip install scipy"
            ) from exc

        knn_idx = np.asarray(adata.obsm[knn_indices_key], dtype=np.int32)
        n_nodes = adata.n_obs

        if verbose:
            print(
                f"[cluster_leiden_jaccard] embedding='{embedding_name}' "
                f"n_cells={n_nodes} k={knn_idx.shape[1]} "
                f"min_jaccard={min_jaccard}"
            )

        t_jac = _time.perf_counter()
        rows, cols, weights = self._knn_to_jaccard_graph_arrays(
            knn_idx, n_nodes, min_jaccard=min_jaccard
        )
        jac_sec = float(_time.perf_counter() - t_jac)

        if verbose:
            print(
                f"[cluster_leiden_jaccard] jaccard edges={len(rows)} "
                f"jac_sec={jac_sec:.3f}"
            )

        # Symmetrize into sparse matrix.
        edge_u = np.concatenate([rows, cols])
        edge_v = np.concatenate([cols, rows])
        sym_w = np.concatenate([weights, weights]).astype(np.float32)
        jaccard_conn = sp.csr_matrix(
            (sym_w, (edge_u, edge_v)),
            shape=(n_nodes, n_nodes),
            dtype=np.float32,
        )
        adata.obsp[jaccard_connectivities_key] = jaccard_conn

        # Build igraph from upper-triangular edges.
        edges = list(zip(rows.tolist(), cols.tolist()))
        graph = ig.Graph(n=n_nodes, edges=edges, directed=False)
        graph.es["weight"] = weights.tolist()

        if verbose:
            print(
                f"[cluster_leiden_jaccard] running Leiden "
                f"resolution={resolution} n_iterations={n_iterations}"
            )

        t_leiden = _time.perf_counter()
        part = graph.community_leiden(
            objective_function=objective_function,
            weights="weight",
            resolution=resolution,
            n_iterations=n_iterations,
            beta=beta,
        )
        leiden_sec = float(_time.perf_counter() - t_leiden)

        labels = np.asarray(part.membership, dtype=np.int32)
        adata.obs[cluster_key] = pd.Categorical(labels.astype(str))

        n_clusters = int(len(set(labels.tolist())))
        clusterings_uns = self._upsert_nested_uns_dict(adata, "clusterings")
        metadata = {
            "timestamp": datetime.utcnow().isoformat(),
            "method": "jaccard_leiden",
            "embedding_name": embedding_name,
            "knn_indices_key": knn_indices_key,
            "jaccard_connectivities_key": jaccard_connectivities_key,
            "cluster_key": cluster_key,
            "min_jaccard": float(min_jaccard),
            "resolution": float(resolution),
            "n_iterations": int(n_iterations),
            "beta": float(beta),
            "objective_function": str(objective_function),
            "seed": int(seed),
            "verbose": bool(verbose),
            "n_cells": int(n_nodes),
            "n_clusters": n_clusters,
            "jac_sec": jac_sec,
            "leiden_sec": leiden_sec,
        }
        clusterings_uns[cluster_key] = metadata

        if verbose:
            print(
                f"[cluster_leiden_jaccard] done cluster_key='{cluster_key}' "
                f"n_clusters={n_clusters} leiden_sec={leiden_sec:.3f}"
            )

        self._adata = adata
        self._log_run_event(
            "leiden_jaccard_clustered",
            {
                "embedding_name": embedding_name,
                "cluster_key": cluster_key,
                "min_jaccard": float(min_jaccard),
                "resolution": float(resolution),
                "n_iterations": int(n_iterations),
                "n_clusters": int(n_clusters),
                "n_cells": int(n_nodes),
            },
        )
        if inplace:
            self._try_save_inplace()

        return metadata

    def compute_flowsom(
        self,
        markers: list[str],
        source_layer: str = "X",
        grid_size: int | tuple[int, int] = 10,
        n_iterations: int = 100,
        neighborhood_radius: int = 1,
        learning_rate: float = 0.5,
        cluster_method: str = "leiden",
        cluster_resolution: float = 1.0,
        n_meta_clusters: int | None = None,
        random_state: int = 42,
        verbose: bool = False,
        inplace: bool = True,
        use_gpu: bool = True,
    ) -> dict[str, Any]:
        """Compute FlowSOM clustering.

        Uses the external ``flowsom`` package (no local reimplementation).
        Results are copied back into this run's AnnData object:

        - ``adata.obs['X_flowsom']``: FlowSOM metaclusters
        - ``adata.obsm['X_flowsom_node']``: SOM node assignment per cell
        - ``adata.varm['flowsom_weights']``: SOM codebook (markers x nodes)

        Args:
            markers: Markers to use for clustering (must exist in `adata.var_names`).
            source_layer: Layer used for expression values (`X` or layer key).
            grid_size: SOM grid size (10x10) or (rows, cols) tuple.
            n_iterations: FlowSOM training passes (`rlen`).
            neighborhood_radius: Kept for compatibility (stored in metadata).
            learning_rate: Initial FlowSOM learning rate (`alpha[0]`).
            cluster_method: Kept for compatibility. FlowSOM metaclustering is
                always used as final labels.
            cluster_resolution: Backward-compatible hint used to derive default
                ``n_meta_clusters`` when not provided.
            n_meta_clusters: Explicit number of FlowSOM metaclusters.
            random_state: RNG seed.
            verbose: If True, print progress to stdout.
            inplace: If True, persist results to run zarr.
            use_gpu: If True, request MPS when available (best-effort; only
                used when the installed FlowSOM backend exposes a ``device``
                parameter).

        Returns:
            Metadata dict with FlowSOM metadata stored under
            `adata.uns["clusterings"]["X_flowsom"]` and `adata.uns["flowsom"]`.
        """
        import time as _time
        import warnings as _warnings

        try:
            import flowsom
        except ImportError as exc:
            raise ImportError(
                "compute_flowsom requires the 'flowsom' package. "
                "Install with: pip install flowsom"
            ) from exc

        try:
            from flowsom.models import BatchFlowSOMEstimator, FlowSOMEstimator
        except ImportError as exc:
            raise ImportError(
                "Unable to import FlowSOM estimator classes from 'flowsom.models'."
            ) from exc

        adata = self.read_adata()
        matrix = self._matrix_from_layer(adata, source_layer)
        var_names = adata.var_names.astype(str).tolist()

        marker_to_idx = {name: i for i, name in enumerate(var_names)}
        missing = [m for m in markers if m not in marker_to_idx]
        if missing:
            raise ValueError(f"Markers not found for FlowSOM: {missing}")
        marker_idx = [marker_to_idx[m] for m in markers]
        x_embed = np.asarray(matrix[:, marker_idx], dtype=np.float32, order="C")

        n_cells, _ = x_embed.shape

        if isinstance(grid_size, int):
            grid_rows = grid_cols = grid_size
        else:
            grid_rows, grid_cols = grid_size

        if grid_rows <= 0 or grid_cols <= 0:
            raise ValueError("grid_size dimensions must be > 0")
        if n_iterations <= 0:
            raise ValueError("n_iterations must be > 0")
        if learning_rate <= 0:
            raise ValueError("learning_rate must be > 0")

        n_nodes = grid_rows * grid_cols

        if cluster_method not in {"leiden", "flowsom", "metaclustering"}:
            raise ValueError(
                "cluster_method must be one of: 'leiden', 'flowsom', 'metaclustering'"
            )

        if n_meta_clusters is None:
            n_meta_clusters = max(2, int(round(cluster_resolution * 10)))
        n_meta_clusters = int(np.clip(n_meta_clusters, 2, n_nodes))

        device_requested = "cpu"
        device_used = "cpu"
        mps_available = False
        mps_enabled = False
        model_kwargs: dict[str, Any] = {}
        model_class = BatchFlowSOMEstimator

        if use_gpu:
            try:
                import torch

                mps_available = bool(
                    hasattr(torch.backends, "mps") and torch.backends.mps.is_available()
                )
            except Exception:
                mps_available = False

            if mps_available:
                device_requested = "mps"
                for candidate in [BatchFlowSOMEstimator, FlowSOMEstimator]:
                    try:
                        params = inspect.signature(candidate.__init__).parameters
                    except (TypeError, ValueError):
                        params = {}
                    if "device" in params:
                        model_class = candidate
                        model_kwargs["device"] = "mps"
                        mps_enabled = True
                        device_used = "mps"
                        break

                if not mps_enabled:
                    _warnings.warn(
                        "MPS is available but the installed flowsom backend does not expose a device parameter; running FlowSOM on CPU.",
                        UserWarning,
                        stacklevel=2,
                    )

        alpha_end = max(0.01, float(learning_rate) * 0.2)
        alpha = (float(learning_rate), alpha_end)

        fs_adata = anndata.AnnData(
            X=x_embed,
            obs=adata.obs.copy(),
            var=pd.DataFrame(index=np.asarray(markers, dtype=str)),
        )
        fs_adata.obs_names = adata.obs_names.copy()

        if verbose:
            print(
                f"[compute_flowsom] cells={n_cells} markers={len(markers)} "
                f"grid={grid_rows}x{grid_cols} rlen={n_iterations} "
                f"n_meta_clusters={n_meta_clusters} model={model_class.__name__} "
                f"device={device_used}"
            )

        t0 = _time.perf_counter()
        fsom = flowsom.FlowSOM(
            fs_adata,
            cols_to_use=list(markers),
            n_clusters=n_meta_clusters,
            model=model_class,
            xdim=grid_rows,
            ydim=grid_cols,
            rlen=n_iterations,
            alpha=alpha,
            seed=random_state,
            **model_kwargs,
        )
        fit_sec = float(_time.perf_counter() - t0)

        cell_data = fsom.mudata["cell_data"]
        cluster_data = fsom.mudata["cluster_data"]

        node_labels = np.asarray(cell_data.obs["clustering"], dtype=np.int32)
        metacluster_labels = cell_data.obs["metaclustering"].astype(str)

        adata.obsm["X_flowsom_node"] = node_labels.reshape(-1, 1).astype(np.int32)

        cluster_key = "X_flowsom"
        adata.obs[cluster_key] = pd.Categorical(metacluster_labels)

        if "codes" in cluster_data.obsm:
            codes = np.asarray(cluster_data.obsm["codes"], dtype=np.float32)
        else:
            codes = np.asarray(cluster_data.X, dtype=np.float32)

        # AnnData varm must have first dimension == adata.n_vars.
        # FlowSOM codes only cover the selected marker subset, so place those
        # rows into a full marker-by-node matrix and leave other markers as NaN.
        full_weights = np.full((adata.n_vars, codes.shape[0]), np.nan, dtype=np.float32)
        full_weights[np.asarray(marker_idx, dtype=np.int32), :] = codes.T.astype(np.float32)
        adata.varm["flowsom_weights"] = full_weights

        if "grid" in cluster_data.obsm:
            adata.uns["flowsom_grid"] = np.asarray(cluster_data.obsm["grid"], dtype=np.int32)

        n_clusters = int(pd.Series(metacluster_labels).nunique())

        clusterings_uns = self._upsert_nested_uns_dict(adata, "clusterings")
        metadata = {
            "timestamp": datetime.utcnow().isoformat(),
            "method": "flowsom_package",
            "model_class": f"{model_class.__module__}.{model_class.__name__}",
            "source_layer": source_layer,
            "markers": list(markers),
            "grid_size": {"rows": grid_rows, "cols": grid_cols},
            "n_iterations": int(n_iterations),
            "alpha": [float(alpha[0]), float(alpha[1])],
            "neighborhood_radius": float(neighborhood_radius),
            "learning_rate": float(learning_rate),
            "cluster_method_requested": str(cluster_method),
            "cluster_method": "flowsom_metaclustering",
            "cluster_resolution": float(cluster_resolution),
            "n_meta_clusters": int(n_meta_clusters),
            "cluster_key": cluster_key,
            "node_key": "X_flowsom_node",
            "n_cells": int(n_cells),
            "n_markers": int(len(markers)),
            "n_nodes": int(codes.shape[0]),
            "n_clusters": n_clusters,
            "fit_sec": fit_sec,
            "use_gpu": bool(use_gpu),
            "mps_available": bool(mps_available),
            "mps_requested": bool(device_requested == "mps"),
            "mps_enabled": bool(mps_enabled),
            "device_used": device_used,
        }
        clusterings_uns[cluster_key] = metadata

        flowsom_uns = self._upsert_nested_uns_dict(adata, "flowsom")
        flowsom_uns["latest"] = {
            "method": "flowsom_package",
            "grid_size": {"rows": grid_rows, "cols": grid_cols},
            "n_iterations": int(n_iterations),
            "neighborhood_radius": float(neighborhood_radius),
            "learning_rate": float(learning_rate),
            "random_state": int(random_state),
            "cluster_resolution": float(cluster_resolution),
            "n_meta_clusters": int(n_meta_clusters),
            "n_nodes": int(codes.shape[0]),
            "n_clusters": n_clusters,
            "n_cells_assigned": int(n_cells),
            "fit_sec": fit_sec,
            "device_used": device_used,
            "mps_available": bool(mps_available),
            "mps_enabled": bool(mps_enabled),
        }

        if verbose:
            print(
                f"[compute_flowsom] done cluster_key='{cluster_key}' "
                f"n_clusters={n_clusters} fit_sec={fit_sec:.3f} device={device_used}"
            )

        self._adata = adata
        self._log_run_event(
            "flowsom_clustered",
            {
                "grid_rows": int(grid_rows),
                "grid_cols": int(grid_cols),
                "n_iterations": int(n_iterations),
                "learning_rate": float(learning_rate),
                "cluster_resolution": float(cluster_resolution),
                "n_meta_clusters": int(n_meta_clusters),
                "n_clusters": n_clusters,
                "n_cells": int(n_cells),
                "n_markers": int(len(markers)),
                "device_used": device_used,
                "mps_enabled": bool(mps_enabled),
            },
        )
        if inplace:
            self._try_save_inplace()

        return metadata

    def annotate_clusters(
        self,
        cluster_key: str,
        annotation_map: dict[str, str],
        output_key: str | None = None,
        inplace: bool = True,
    ) -> pd.Series:
        """Map cluster labels to named cell types.

        Adds a new obs column with the mapped names.  Unmapped labels are
        passed through unchanged so partial mappings are allowed.

        Args:
            cluster_key: Obs column containing cluster labels (e.g.
                ``"my_umap_leiden"``).
            annotation_map: Dict from cluster label (str) to cell-type name.
                Keys are compared against string-cast cluster values, so
                integer labels like ``0`` should be passed as ``"0"``.
            output_key: Name for the new obs column.  Defaults to
                ``"{cluster_key}_annotated"``.
            inplace: If True, persist updated AnnData to Zarr.

        Returns:
            pd.Series of annotated labels indexed by cell names.

        Raises:
            ValueError: If ``cluster_key`` is not in ``adata.obs``.
        """
        import warnings

        adata = self.read_adata()
        if cluster_key not in adata.obs.columns:
            raise ValueError(
                f"cluster_key '{cluster_key}' not found in adata.obs. "
                f"Available columns: {list(adata.obs.columns)}"
            )

        output_key = output_key or f"{cluster_key}_annotated"
        src = adata.obs[cluster_key].astype(str)

        missing_keys = set(annotation_map.keys()) - set(src.unique())
        if missing_keys:
            warnings.warn(
                f"annotation_map keys not found in '{cluster_key}': {missing_keys}",
                UserWarning,
                stacklevel=2,
            )

        annotated = src.map(annotation_map).fillna(src)
        adata.obs[output_key] = annotated.values

        n_mapped = int((annotated != src).sum())
        self._log_run_event(
            "clusters_annotated",
            {
                "cluster_key": cluster_key,
                "output_key": output_key,
                "n_mapped": n_mapped,
                "annotation_map": annotation_map,
            },
        )

        if inplace:
            self._try_save_inplace(adata)

        return annotated

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

        def _store_z_to_obs(obs_prefix: str, score_block: dict[str, pd.DataFrame]) -> list[str]:
            z_df = score_block["z"].loc[adata.obs_names]
            col_prefix = f"{obs_prefix}_"
            existing = [col for col in adata.obs.columns if col.startswith(col_prefix)]
            if existing:
                adata.obs = adata.obs.drop(columns=existing)
            obs_cols = []
            for col in z_df.columns:
                col_name = f"{col_prefix}{col}"
                adata.obs[col_name] = z_df[col].to_numpy(dtype=np.float32, copy=True)
                obs_cols.append(col_name)
            return obs_cols

        result_keys_smoothed = _store_score_block(result_prefix, smoothed_scores)
        result_keys_raw = (
            _store_score_block(f"{result_prefix}_raw", raw_scores)
            if raw_scores is not None
            else None
        )
        _store_z_to_obs("PC", smoothed_scores)
        if raw_scores is not None:
            _store_z_to_obs("PC_R", raw_scores)

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
            self._try_save_inplace()

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

    def _plot_marker_coverage(
        self,
        all_marker_mappings: list[pd.DataFrame],
        sample_df: pd.DataFrame,
    ) -> None:
        """Display a heatmap of standard-marker presence/absence across files.

        Only markers missing from at least one file are shown.
        If every marker is present in every file, prints a message instead.
        """
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
        except ImportError as exc:
            raise ImportError(
                "show_marker_coverage requires matplotlib and seaborn."
            ) from exc

        matched_statuses = {"matched_standard", "matched_alias"}

        fname_to_sample: dict[str, str] = {}
        if "sample_id" in sample_df.columns and "file_name" in sample_df.columns:
            fname_to_sample = dict(
                zip(sample_df["file_name"], sample_df["sample_id"])
            )

        file_order: list[str] = []
        file_markers: dict[str, set[str]] = {}
        for mapping_df in all_marker_mappings:
            fname = mapping_df["file_name"].iloc[0]
            label = fname_to_sample.get(fname, fname)
            matched = mapping_df[mapping_df["mapping_status"].isin(matched_statuses)]
            file_markers[label] = set(matched["standard_marker_name"].dropna().tolist())
            file_order.append(label)

        all_standard = sorted(
            {m for markers in file_markers.values() for m in markers}
        )

        presence = pd.DataFrame(
            {
                label: [m in file_markers[label] for m in all_standard]
                for label in file_order
            },
            index=all_standard,
        )

        coverage_df = presence[~presence.all(axis=1)]

        if coverage_df.empty:
            print(
                f"[{self.run_id}] Marker coverage: all {len(all_standard)} standard markers "
                "are present in every file — nothing to show."
            )
            return

        n_markers, n_files = coverage_df.shape
        fig, ax = plt.subplots(
            figsize=(max(4, n_files * 0.8 + 2), max(3, n_markers * 0.35 + 1.5))
        )
        # red = absent (0), green = present (1)
        sns.heatmap(
            coverage_df.astype(int),
            annot=False,
            cmap=sns.color_palette(["#d62728", "#2ca02c"], as_cmap=False),
            vmin=0,
            vmax=1,
            linewidths=0.4,
            linecolor="white",
            cbar=False,
            ax=ax,
        )
        ax.set_title(
            f"{self.run_id} — marker coverage\n"
            f"({n_markers} marker(s) missing in ≥1 file)  "
            "green = present  |  red = absent",
            fontsize=10,
        )
        ax.set_xlabel("Sample")
        ax.set_ylabel("Standard marker")
        ax.tick_params(axis="x", rotation=45)
        ax.tick_params(axis="y", rotation=0)
        fig.tight_layout()
        plt.show()

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
            cmap = "seismic"
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

    def compare_groups(
        self,
        field: str,
        groupby,
        layer: str = "X",
        comparisons: str | list[tuple[str, str]] | None = "all",
        method: str = "ttest",
        equal_var: bool = True,
        order: list[str] | None = None,
        multitest: str | None = "bh",
    ) -> pd.DataFrame:
        """Compute pairwise group comparisons for a marker or numeric obs field.

        Args:
            field: Marker name (in `adata.var_names`) or numeric obs column.
            groupby: Single obs column or list of obs columns for composite grouping.
            layer: Layer used when `field` is a marker (default `X`).
            comparisons: `"all"` (default), `"adjacent"`, explicit list of
                (group_a, group_b) pairs, or None to skip comparisons.
            method: `"ttest"` or `"wald"`.
            equal_var: When method is `"ttest"`, use pooled variance (True)
                or Welch correction (False).
            order: Optional explicit group order; otherwise sorted unique groups.
            multitest: Optional p-value correction (`"bh"` or `"bonferroni"`).

        Returns:
            DataFrame with columns:
                group_a, group_b, n_a, n_b, mean_a, mean_b, var_a, var_b,
                stat, p_value, p_adj, method
        """
        if method not in {"ttest", "wald"}:
            raise ValueError("method must be 'ttest' or 'wald'")
        if multitest not in {None, "bh", "bonferroni"}:
            raise ValueError("multitest must be None, 'bh', or 'bonferroni'")

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

        pairs = self._build_comparison_pairs(order, comparisons)
        if not pairs:
            return pd.DataFrame(
                columns=[
                    "group_a",
                    "group_b",
                    "n_a",
                    "n_b",
                    "mean_a",
                    "mean_b",
                    "var_a",
                    "var_b",
                    "stat",
                    "p_value",
                    "p_adj",
                    "method",
                ]
            )

        try:
            from scipy import stats
        except ImportError as exc:
            raise ImportError(
                "compare_groups requires scipy. Install with: pip install scipy"
            ) from exc

        results = []
        pvals = []
        for a, b in pairs:
            x = plot_df.loc[plot_df["__group__"] == a, field].to_numpy(dtype=float)
            y = plot_df.loc[plot_df["__group__"] == b, field].to_numpy(dtype=float)
            x = x[~np.isnan(x)]
            y = y[~np.isnan(y)]

            n_a = int(len(x))
            n_b = int(len(y))
            mean_a = float(np.mean(x)) if n_a > 0 else float("nan")
            mean_b = float(np.mean(y)) if n_b > 0 else float("nan")
            var_a = float(np.var(x, ddof=1)) if n_a > 1 else float("nan")
            var_b = float(np.var(y, ddof=1)) if n_b > 1 else float("nan")

            stat = float("nan")
            p_value = float("nan")

            if n_a >= 2 and n_b >= 2:
                if method == "ttest":
                    stat, p_value = stats.ttest_ind(x, y, equal_var=equal_var)
                    stat = float(stat)
                    p_value = float(p_value)
                else:
                    # Wald z-test on difference of means
                    se = np.sqrt(var_a / n_a + var_b / n_b)
                    if se > 0 and np.isfinite(se):
                        stat = float((mean_a - mean_b) / se)
                        p_value = float(2.0 * stats.norm.sf(abs(stat)))

            results.append(
                {
                    "group_a": a,
                    "group_b": b,
                    "n_a": n_a,
                    "n_b": n_b,
                    "mean_a": mean_a,
                    "mean_b": mean_b,
                    "var_a": var_a,
                    "var_b": var_b,
                    "stat": stat,
                    "p_value": p_value,
                    "method": method,
                }
            )
            pvals.append(p_value)

        p_adj = None
        if multitest is not None:
            p_adj = self._adjust_pvalues(pvals, method=multitest)

        out = pd.DataFrame(results)
        if p_adj is not None:
            out["p_adj"] = p_adj
        else:
            out["p_adj"] = np.nan
        return out

    def differential_abundance(
        self,
        cluster_key: str,
        groupby: str,
        comparisons: str | list[tuple[str, str]] | None = "all",
        method: str = "fisher",
        multitest: str | None = "bh",
        order: list[str] | None = None,
        plot: bool = False,
        figsize: tuple[float, float] | None = None,
        ax=None,
    ) -> pd.DataFrame | tuple[pd.DataFrame, tuple]:
        """Test whether cluster proportions differ between groups.

        For each cluster and each pair of groups, the observed cell counts are
        compared using Fisher's exact test or a chi-squared test.  P-values are
        optionally corrected across all clusters × all pairs.

        Args:
            cluster_key: Obs column with cluster labels (e.g.
                ``"my_umap_leiden"``).
            groupby: Obs column defining the groups to compare (e.g.
                ``"line_id"`` or ``"condition"``).
            comparisons: ``"all"`` (every pair), ``"adjacent"`` (neighbours in
                ``order``), explicit list of ``(group_a, group_b)`` tuples, or
                ``None`` (skip — returns proportions without p-values).
            method: ``"fisher"`` (Fisher's exact) or ``"chi2"``
                (chi-squared contingency).
            multitest: P-value correction: ``"bh"`` (Benjamini-Hochberg),
                ``"bonferroni"``, or ``None``.
            order: Explicit group order.  Defaults to sorted unique groups.
            plot: If True, also return a stacked-bar proportion plot.
            figsize: Figure size when ``plot=True``.
            ax: Existing axes to draw on (single axes, ``plot=True`` only).

        Returns:
            DataFrame with columns:
                ``cluster``, ``group_a``, ``group_b``,
                ``n_a``, ``n_b``, ``N_a``, ``N_b``,
                ``freq_a``, ``freq_b``, ``p_value``, ``p_adj``

            When ``plot=True``: ``(DataFrame, (fig, ax))``.

        Raises:
            ValueError: If ``cluster_key`` or ``groupby`` are not in
                ``adata.obs``.
        """
        if method not in {"fisher", "chi2"}:
            raise ValueError("method must be 'fisher' or 'chi2'")
        if multitest not in {None, "bh", "bonferroni"}:
            raise ValueError("multitest must be None, 'bh', or 'bonferroni'")

        adata = self.read_adata()
        if cluster_key not in adata.obs.columns:
            raise ValueError(
                f"cluster_key '{cluster_key}' not found in adata.obs. "
                f"Available columns: {list(adata.obs.columns)}"
            )
        if groupby not in adata.obs.columns:
            raise ValueError(
                f"groupby '{groupby}' not found in adata.obs. "
                f"Available columns: {list(adata.obs.columns)}"
            )

        ct = pd.crosstab(adata.obs[groupby], adata.obs[cluster_key].astype(str))
        groups = order if order is not None else sorted(ct.index.tolist())
        ct = ct.loc[groups]
        clusters = ct.columns.tolist()

        pairs = self._build_comparison_pairs(groups, comparisons)

        try:
            from scipy import stats
        except ImportError as exc:
            raise ImportError(
                "differential_abundance requires scipy. "
                "Install with: pip install scipy"
            ) from exc

        rows = []
        pvals = []
        for group_a, group_b in pairs:
            N_a = int(ct.loc[group_a].sum())
            N_b = int(ct.loc[group_b].sum())
            for c in clusters:
                n_a = int(ct.loc[group_a, c])
                n_b = int(ct.loc[group_b, c])
                freq_a = n_a / N_a if N_a > 0 else float("nan")
                freq_b = n_b / N_b if N_b > 0 else float("nan")
                table = [[n_a, N_a - n_a], [n_b, N_b - n_b]]
                if method == "fisher":
                    _, p_value = stats.fisher_exact(table)
                else:
                    _, p_value, _, _ = stats.chi2_contingency(table)
                p_value = float(p_value)
                rows.append(
                    {
                        "cluster": c,
                        "group_a": group_a,
                        "group_b": group_b,
                        "n_a": n_a,
                        "n_b": n_b,
                        "N_a": N_a,
                        "N_b": N_b,
                        "freq_a": freq_a,
                        "freq_b": freq_b,
                        "p_value": p_value,
                    }
                )
                pvals.append(p_value)

        results_df = pd.DataFrame(
            rows,
            columns=[
                "cluster",
                "group_a",
                "group_b",
                "n_a",
                "n_b",
                "N_a",
                "N_b",
                "freq_a",
                "freq_b",
                "p_value",
            ],
        )

        if pvals and multitest is not None:
            results_df["p_adj"] = self._adjust_pvalues(pvals, method=multitest)
        else:
            results_df["p_adj"] = np.nan

        self._log_run_event(
            "differential_abundance",
            {
                "cluster_key": cluster_key,
                "groupby": groupby,
                "method": method,
                "n_clusters": len(clusters),
                "n_pairs": len(pairs),
            },
        )

        if not plot:
            return results_df

        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
        except ImportError as exc:
            raise ImportError(
                "differential_abundance with plot=True requires matplotlib and seaborn. "
                "Install with: pip install matplotlib seaborn"
            ) from exc

        props = ct.div(ct.sum(axis=1), axis=0)
        long_df = props.reset_index().melt(
            id_vars=groupby, var_name="cluster", value_name="proportion"
        )

        if ax is None:
            fig, ax = plt.subplots(
                figsize=figsize or (max(6, 0.6 * len(clusters) + 2), 4)
            )
        else:
            fig = ax.figure

        sns.barplot(data=long_df, x="cluster", y="proportion", hue=groupby, ax=ax)
        ax.set_xlabel("Cluster")
        ax.set_ylabel("Proportion")
        ax.legend(title=groupby, bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0)

        if pairs and not results_df.empty:
            y_max = props.values.max()
            step = y_max * 0.08
            pair_offsets: dict[str, int] = {}
            for _, row in results_df.iterrows():
                c = row["cluster"]
                stars = self._pvalue_to_stars(row["p_adj"] if not np.isnan(row["p_adj"]) else row["p_value"])
                if stars == "ns":
                    continue
                cluster_x_positions = {cl: i for i, cl in enumerate(clusters)}
                x_pos = cluster_x_positions.get(c, 0)
                offset_idx = pair_offsets.get(c, 0)
                y_pos = y_max + step * (offset_idx + 1)
                ax.text(
                    x_pos,
                    y_pos,
                    stars,
                    ha="center",
                    va="bottom",
                    fontsize=10,
                )
                pair_offsets[c] = offset_idx + 1

        fig.tight_layout()
        return results_df, (fig, ax)

    def match_clusterings(
        self,
        key_a: str,
        key_b: str,
        score_mode: str = "jaccard",
        n_permutations: int = 1000,
        random_state: int = 0,
        plot: str | bool = False,
        matrix: str = "score",
        annot: bool = True,
        figsize: tuple[float, float] | None = None,
    ) -> dict:
        """Compare two obs columns (clusterings or any categorical labels).

        Performs both many-to-one matching (every A label gets its best B
        label and vice versa) and one-to-one optimal matching via the
        Hungarian algorithm.  Hypergeometric enrichment and permutation-based
        significance are computed for each match.

        Args:
            key_a: First obs column (e.g. ``"CL"``).
            key_b: Second obs column (e.g. ``"sample_id"`` or ``"CL_ID"``).
            score_mode: Matching criterion — ``"jaccard"`` (default),
                ``"a"`` (fraction of A in B), ``"b"`` (fraction of B from A),
                or ``"overlap"`` / ``None`` (raw count).
            n_permutations: Permutations for global significance test.
                Set to 0 to skip.
            random_state: Random seed.
            plot: ``False`` (default, no plot), ``"heatmap"`` (seaborn heatmap,
                no dendrogram), or ``"clustermap"`` (seaborn clustermap with
                row/column dendrograms).
            matrix: Which values to display — ``"score"`` (default, the
                scoring matrix used for matching) or ``"overlap"`` (raw cell
                counts from the contingency table).
            annot: Annotate each cell with its value.  Default ``True``.
            figsize: Figure size passed to the plot.  When ``None`` a
                sensible default is chosen based on the matrix dimensions.

        Returns:
            When ``plot=False``: the report ``dict``.
            When ``plot`` is set: ``(report_dict, figure)`` where *figure* is
            the ``matplotlib.figure.Figure`` for a heatmap or the
            ``seaborn.ClusterGrid`` for a clustermap.

            Report dict keys:

            - ``A_to_B`` — every A label matched to its best B label
            - ``B_to_A`` — every B label matched to its best A label
            - ``one_to_one_matches`` — Hungarian optimal assignment
            - ``contingency`` — raw overlap crosstab
            - ``score_matrix`` — scoring matrix used for matching
            - ``B_receives_from_A``, ``A_receives_from_B`` — split/merge structure
            - ``ari``, ``nmi`` — global agreement metrics
            - ``global_many_to_one_score_A_to_B``, ``_B_to_A``, ``_one_to_one``
            - ``permutation_p_A_to_B``, ``_B_to_A``, ``_one_to_one``

        Raises:
            ValueError: If either key is not in ``adata.obs``, or ``plot``/
                ``matrix`` are invalid values.
        """
        try:
            from scipy.optimize import linear_sum_assignment
            from scipy.stats import hypergeom
        except ImportError as exc:
            raise ImportError(
                "match_clusterings requires scipy. Install with: pip install scipy"
            ) from exc
        try:
            from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
        except ImportError as exc:
            raise ImportError(
                "match_clusterings requires scikit-learn. "
                "Install with: pip install scikit-learn"
            ) from exc
        try:
            from statsmodels.stats.multitest import multipletests
        except ImportError as exc:
            raise ImportError(
                "match_clusterings requires statsmodels. "
                "Install with: pip install statsmodels"
            ) from exc

        adata = self.read_adata()
        for key in (key_a, key_b):
            if key not in adata.obs.columns:
                raise ValueError(
                    f"'{key}' not found in adata.obs. "
                    f"Available columns: {list(adata.obs.columns)}"
                )

        labels_a = adata.obs[key_a].astype(str).reset_index(drop=True)
        labels_b = adata.obs[key_b].astype(str).reset_index(drop=True)

        n = len(labels_a)
        rng = np.random.default_rng(random_state)

        def _make_score_matrix(a, b):
            tab = pd.crosstab(a, b)
            overlap = tab.to_numpy(dtype=float)
            size_a = tab.sum(axis=1).to_numpy(dtype=float)[:, None]
            size_b = tab.sum(axis=0).to_numpy(dtype=float)[None, :]
            if score_mode is None or score_mode == "overlap":
                vals = overlap
            elif score_mode == "a":
                vals = np.divide(overlap, size_a, out=np.zeros_like(overlap), where=size_a != 0)
            elif score_mode == "b":
                vals = np.divide(overlap, size_b, out=np.zeros_like(overlap), where=size_b != 0)
            elif score_mode == "jaccard":
                union = size_a + size_b - overlap
                vals = np.divide(overlap, union, out=np.zeros_like(overlap), where=union != 0)
            else:
                raise ValueError(
                    "score_mode must be one of None, 'overlap', 'jaccard', 'a', or 'b'."
                )
            return tab, pd.DataFrame(vals, index=tab.index, columns=tab.columns)

        def _add_pair_stats(df):
            df = df.copy()
            expected, fold, pvals = [], [], []
            for _, row in df.iterrows():
                a_sz, b_sz, k = int(row["size_A"]), int(row["size_B"]), int(row["overlap"])
                exp = a_sz * b_sz / n
                expected.append(exp)
                fold.append(k / exp if exp > 0 else float("nan"))
                pvals.append(float(hypergeom.sf(k - 1, n, b_sz, a_sz)))
            df["expected_overlap"] = expected
            df["fold_enrichment"] = fold
            df["hypergeom_p"] = pvals
            if len(df) > 0:
                df["hypergeom_q"] = multipletests(df["hypergeom_p"], method="fdr_bh")[1]
            else:
                df["hypergeom_q"] = []
            return df

        contingency, score_matrix = _make_score_matrix(labels_a, labels_b)

        def _many_to_one(src_axis, src_col, dst_col, sort_frac_col):
            best = score_matrix.idxmax(axis=src_axis)
            rows = []
            for src, dst in best.items():
                if src_axis == 1:
                    a_cl, b_cl = src, dst
                else:
                    a_cl, b_cl = dst, src
                ov = contingency.loc[a_cl, b_cl]
                sa = contingency.loc[a_cl, :].sum()
                sb = contingency.loc[:, b_cl].sum()
                sc = score_matrix.loc[a_cl, b_cl]
                rows.append({
                    src_col: src,
                    dst_col: dst,
                    "score": sc,
                    "overlap": ov,
                    "size_A": sa,
                    "size_B": sb,
                    "frac_A_in_B": ov / sa if sa > 0 else float("nan"),
                    "frac_B_from_A": ov / sb if sb > 0 else float("nan"),
                })
            df = pd.DataFrame(rows)
            df = _add_pair_stats(df)
            return df.sort_values(
                ["score", sort_frac_col, "fold_enrichment", "overlap"],
                ascending=False,
            ).reset_index(drop=True)

        A_to_B = _many_to_one(1, "cluster_A", "matched_cluster_B", "frac_A_in_B")
        B_to_A = _many_to_one(0, "cluster_B", "matched_cluster_A", "frac_B_from_A")

        row_ind, col_ind = linear_sum_assignment(-score_matrix.to_numpy())
        one_to_one = pd.DataFrame({
            "cluster_A": score_matrix.index[row_ind],
            "cluster_B": score_matrix.columns[col_ind],
            "score": score_matrix.to_numpy()[row_ind, col_ind],
            "overlap": contingency.to_numpy()[row_ind, col_ind],
            "size_A": contingency.sum(axis=1).to_numpy()[row_ind],
            "size_B": contingency.sum(axis=0).to_numpy()[col_ind],
        })
        one_to_one["frac_A_in_B"] = one_to_one["overlap"] / one_to_one["size_A"]
        one_to_one["frac_B_from_A"] = one_to_one["overlap"] / one_to_one["size_B"]
        one_to_one = _add_pair_stats(one_to_one)
        one_to_one = one_to_one.sort_values(
            ["score", "fold_enrichment", "overlap"], ascending=False
        ).reset_index(drop=True)

        ari = float(adjusted_rand_score(labels_a, labels_b))
        nmi = float(normalized_mutual_info_score(labels_a, labels_b))

        g_a2b = float(A_to_B["score"].sum())
        g_b2a = float(B_to_A["score"].sum())
        g_121 = float(one_to_one["score"].sum())

        null_a2b = null_b2a = null_121 = None
        p_a2b = p_b2a = p_121 = None

        if n_permutations and n_permutations > 0:
            null_a2b = np.empty(n_permutations)
            null_b2a = np.empty(n_permutations)
            null_121 = np.empty(n_permutations)
            for i in range(n_permutations):
                perm_b = labels_b.iloc[rng.permutation(n)].reset_index(drop=True)
                _, ps = _make_score_matrix(labels_a, perm_b)
                null_a2b[i] = ps.max(axis=1).sum()
                null_b2a[i] = ps.max(axis=0).sum()
                r, c = linear_sum_assignment(-ps.to_numpy())
                null_121[i] = ps.to_numpy()[r, c].sum()
            p_a2b = float((1 + np.sum(null_a2b >= g_a2b)) / (n_permutations + 1))
            p_b2a = float((1 + np.sum(null_b2a >= g_b2a)) / (n_permutations + 1))
            p_121 = float((1 + np.sum(null_121 >= g_121)) / (n_permutations + 1))

        B_receives_from_A = (
            A_to_B.groupby("matched_cluster_B")["cluster_A"]
            .apply(list)
            .reset_index()
            .rename(columns={"matched_cluster_B": "cluster_B",
                              "cluster_A": "A_clusters_mapping_to_B"})
        )
        B_receives_from_A["n_A_clusters"] = B_receives_from_A[
            "A_clusters_mapping_to_B"
        ].apply(len)

        A_receives_from_B = (
            B_to_A.groupby("matched_cluster_A")["cluster_B"]
            .apply(list)
            .reset_index()
            .rename(columns={"matched_cluster_A": "cluster_A",
                              "cluster_B": "B_clusters_mapping_to_A"})
        )
        A_receives_from_B["n_B_clusters"] = A_receives_from_B[
            "B_clusters_mapping_to_A"
        ].apply(len)

        self._log_run_event(
            "clusterings_matched",
            {
                "key_a": key_a,
                "key_b": key_b,
                "score_mode": score_mode,
                "n_permutations": n_permutations,
                "ari": ari,
                "nmi": nmi,
            },
        )

        report = {
            "A_to_B": A_to_B,
            "B_to_A": B_to_A,
            "contingency": contingency,
            "score_matrix": score_matrix,
            "one_to_one_matches": one_to_one,
            "unmatched_A_one_to_one": sorted(
                set(contingency.index) - set(one_to_one["cluster_A"])
            ),
            "unmatched_B_one_to_one": sorted(
                set(contingency.columns) - set(one_to_one["cluster_B"])
            ),
            "B_receives_from_A": B_receives_from_A,
            "A_receives_from_B": A_receives_from_B,
            "ari": ari,
            "nmi": nmi,
            "global_many_to_one_score_A_to_B": g_a2b,
            "global_many_to_one_score_B_to_A": g_b2a,
            "global_one_to_one_score": g_121,
            "permutation_p_A_to_B": p_a2b,
            "permutation_p_B_to_A": p_b2a,
            "permutation_p_one_to_one": p_121,
            "permutation_null_A_to_B": null_a2b,
            "permutation_null_B_to_A": null_b2a,
            "permutation_null_one_to_one": null_121,
        }

        if not plot:
            return report

        if plot not in ("heatmap", "clustermap"):
            raise ValueError("plot must be False, 'heatmap', or 'clustermap'.")
        if matrix not in ("score", "overlap"):
            raise ValueError("matrix must be 'score' or 'overlap'.")

        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
        except ImportError as exc:
            raise ImportError(
                "Plotting match_clusterings requires matplotlib and seaborn."
            ) from exc

        plot_df = score_matrix if matrix == "score" else contingency.astype(float)
        nr, nc = plot_df.shape
        default_fs = (max(4, nc * 0.6 + 2), max(3, nr * 0.5 + 1.5))
        fmt = ".2f" if matrix == "score" else "g"
        title = (
            f"{key_a} vs {key_b}  —  "
            f"{matrix} matrix  |  ARI={ari:.3f}  NMI={nmi:.3f}"
        )

        if plot == "heatmap":
            fig, ax = plt.subplots(figsize=figsize or default_fs)
            sns.heatmap(
                plot_df,
                annot=annot,
                fmt=fmt,
                cmap="YlOrRd",
                linewidths=0.3,
                ax=ax,
            )
            ax.set_xlabel(key_b)
            ax.set_ylabel(key_a)
            ax.set_title(title)
            fig.tight_layout()
            return report, fig

        # clustermap
        grid = sns.clustermap(
            plot_df,
            annot=annot,
            fmt=fmt,
            cmap="YlOrRd",
            linewidths=0.3,
            figsize=figsize or (max(5, nc * 0.6 + 3), max(4, nr * 0.5 + 2.5)),
        )
        grid.ax_heatmap.set_xlabel(key_b)
        grid.ax_heatmap.set_ylabel(key_a)
        grid.figure.suptitle(title, y=1.02)
        return report, grid

    def plot_cluster_composition(
        self,
        cluster_key: str,
        groupby: str,
        normalize: str = "cluster",
        palette=None,
        figsize: tuple[float, float] | None = None,
        legend_kwargs: dict | None = None,
        ax=None,
    ) -> tuple:
        """Stacked bar chart of label composition.

        Two modes controlled by ``normalize``:

        - ``"cluster"`` *(default)* — one bar per cluster, showing the
          fraction of cells in that cluster that come from each group (e.g.
          sample).  Answers: *"what samples make up each cluster?"*
        - ``"group"`` — one bar per group, showing the fraction of cells in
          that group assigned to each cluster.  Answers: *"how are each
          sample's cells distributed across clusters?"*

        Args:
            cluster_key: Obs column with cluster labels.
            groupby: Obs column with group labels (e.g. ``"sample_id"``).
            normalize: ``"cluster"`` or ``"group"``.
            palette: Colour palette passed to seaborn.  ``None`` uses the
                default categorical palette.
            figsize: Figure size.  Defaults to a width proportional to the
                number of bars.
            legend_kwargs: Extra kwargs forwarded to ``ax.legend()``.
            ax: Existing axes to draw on.

        Returns:
            ``(fig, ax)`` tuple.

        Raises:
            ValueError: If ``cluster_key`` or ``groupby`` are not in obs,
                or ``normalize`` is not ``"cluster"`` or ``"group"``.
        """
        if normalize not in {"cluster", "group"}:
            raise ValueError("normalize must be 'cluster' or 'group'")

        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
        except ImportError as exc:
            raise ImportError(
                "plot_cluster_composition requires matplotlib and seaborn."
            ) from exc

        adata = self.read_adata()
        for key in (cluster_key, groupby):
            if key not in adata.obs.columns:
                raise ValueError(
                    f"'{key}' not found in adata.obs. "
                    f"Available columns: {list(adata.obs.columns)}"
                )

        ct = pd.crosstab(adata.obs[cluster_key].astype(str),
                         adata.obs[groupby].astype(str))

        if normalize == "cluster":
            # rows = clusters, cols = groups → normalise across columns (groups)
            props = ct.div(ct.sum(axis=1), axis=0)
            xlabel = cluster_key
            legend_title = groupby
            hue_labels = ct.columns.tolist()
        else:
            # rows = clusters, cols = groups → normalise down rows (clusters per group)
            props = ct.div(ct.sum(axis=0), axis=1)
            # Transpose so x-axis = groups
            props = props.T
            xlabel = groupby
            legend_title = cluster_key
            hue_labels = props.columns.tolist()

        n_bars = len(props)
        n_hues = len(hue_labels)

        if palette is None:
            palette = sns.color_palette("tab20" if n_hues > 10 else "tab10", n_hues)
        elif isinstance(palette, str):
            palette = sns.color_palette(palette, n_hues)

        if ax is None:
            fig, ax = plt.subplots(
                figsize=figsize or (max(4, 0.55 * n_bars + 1.5), 5)
            )
        else:
            fig = ax.figure

        bottoms = np.zeros(n_bars)
        x = np.arange(n_bars)
        for i, hue in enumerate(hue_labels):
            vals = props[hue].to_numpy()
            ax.bar(x, vals, bottom=bottoms, label=str(hue),
                   color=palette[i % len(palette)], width=0.8)
            bottoms += vals

        ax.set_xticks(x)
        ax.set_xticklabels(props.index.tolist(), rotation=45, ha="right")
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Proportion")
        ax.set_ylim(0, 1)

        lkw = {"title": legend_title, "bbox_to_anchor": (1.02, 1),
                "loc": "upper left", "borderaxespad": 0}
        if legend_kwargs:
            lkw.update(legend_kwargs)
        ax.legend(**lkw)

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
        multitest: str | None = "bh",
        show_points: bool = False,
        show_outliers: bool = True,
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
        significance_thresholds: list[tuple[float, str]] | None = None,
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
                - `None` (default): no significance brackets are drawn.
                - `"all"`: all unordered pairs.
                - `"adjacent"`: only neighbours in `order`.
                - List of `(group_a, group_b)` tuples for explicit pairs.
            test: `"mannwhitney"` (default), `"ttest"`, or `"welch"`.
            multitest: `None`, `"bh"`, or `"bonferroni"`.
            show_points: Overlay a stripplot of per-cell points (subsampled).
            show_outliers: Whether to render boxplot outlier markers
                (maps to seaborn's `showfliers`). Defaults to True. Set to
                False to clean up dense plots (or pass `"showfliers": False`
                via `boxplot_kwargs`, which takes precedence).
            max_points: Maximum number of stripplot points (subsampled).
            point_alpha: Alpha for overlaid points.
            point_size: Size for overlaid points.
            palette: Seaborn palette name, color list, or dict. When None,
                a distinct color per group is generated automatically using
                the `"tab10"` palette.
            figsize: Optional figure size.
            ax: Optional matplotlib axes to draw into.
            bracket_color: Color used for significance brackets.
            bracket_linewidth: Line width used for significance brackets.
            bracket_fontsize: Font size used for significance labels.
            ns_label: Label used for non-significant comparisons.
            significance_thresholds: Ordered list of `(p_threshold, label)`
                tuples used to convert each (adjusted) p-value into a star
                label. The first tuple whose `p_threshold` is `>= p` wins.
                Provide thresholds from strictest to loosest, for example
                `[(1e-4, "****"), (1e-3, "***"), (1e-2, "**"), (5e-2, "*")]`
                (the default). Any p-value larger than every threshold is
                labelled with `ns_label`.
            random_state: Seed used when subsampling points.
            boxplot_kwargs: Extra keyword arguments forwarded to
                `seaborn.boxplot` (e.g. `width`, `linewidth`, `notch`,
                `whis`, `saturation`, `showfliers`). Keys here override the
                explicit arguments above.
            stripplot_kwargs: Extra keyword arguments forwarded to
                `seaborn.stripplot` when `show_points=True` (e.g. `jitter`,
                `dodge`).

        Returns:
            Tuple of (figure, axes).
        """
        if test not in {"mannwhitney", "ttest", "welch"}:
            raise ValueError("test must be 'mannwhitney', 'ttest', or 'welch'")
        if multitest not in {None, "bh", "bonferroni"}:
            raise ValueError("multitest must be None, 'bh', or 'bonferroni'")

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

        # Default palette = one distinct color per group (auto-coloring).
        if palette is None:
            palette = sns.color_palette("tab10", n_colors=len(order))

        boxplot_extra = {
            k: v for k, v in (boxplot_kwargs or {}).items()
            if k not in {
                "data", "x", "y", "hue", "order", "ax", "palette", "legend",
            }
        }
        boxplot_call = dict(
            data=plot_df,
            x="__group__",
            y=field,
            hue="__group__",
            order=order,
            ax=ax,
            palette=palette,
            legend=False,
            showfliers=show_outliers,
        )
        boxplot_call.update(boxplot_extra)

        try:
            sns.boxplot(**boxplot_call)
        except TypeError:
            # Older seaborn versions (<0.13) do not accept hue/legend on boxplot
            # without an explicit hue variable; fall back to the simple form.
            for key in ("hue", "legend"):
                boxplot_call.pop(key, None)
            sns.boxplot(**boxplot_call)

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
                significance_thresholds=significance_thresholds,
            )

        fig.tight_layout()
        return fig, ax

    def vendi_score(
        self,
        k: int = 15,
        markers: list[str] | None = None,
        layer: str | None = None,
        use_rep: str | None = None,
        metric: str = "cosine",
        obs_key: str = "vendi_score",
        groupby: str | list[str] | None = None,
        sigma: float | None = None,
        n_bins: int = 10,
        n_reps: int = 1,
        m: int | None = None,
        random_state: int = 0,
        return_eigenvalues: bool = False,
        inplace: bool = True,
    ) -> np.ndarray | pd.DataFrame | tuple[pd.DataFrame, dict[str, np.ndarray]]:
        """Compute Vendi score over local k-NN neighborhoods or whole groups.

        **Per-cell mode** (``groupby=None``): for each cell i, features of its
        k-nearest-neighbor neighborhood N_k(i) are binned with
        ``KBinsDiscretizer`` (uniform) then scored with
        ``vendi_score.vendi.score_dual``.  With rarefaction (``n_reps > 1``),
        results are bootstrap-averaged and 95 % CI bounds are also stored in
        ``adata.obs`` as ``{obs_key}_ci_low`` / ``{obs_key}_ci_high``.

        **Per-group mode** (``groupby`` set): for each unique value (or
        combination of values when ``groupby`` is a list) in the chosen
        ``adata.obs`` column(s), all cells belonging to that group are
        rarefied, binned, and scored.  Returns a DataFrame indexed by group
        label (MultiIndex when ``groupby`` is a list) with columns
        ``vendi_score``, ``ci_low``, ``ci_high``, stored in
        ``adata.uns["vendi"][obs_key]``.

        Args:
            k: Number of nearest neighbors per cell (per-cell mode only).
            markers: Subset of marker names (``adata.var_names``) to use as
                features. ``None`` uses all markers. Ignored when ``use_rep``
                is set.
            layer: AnnData layer to use as the feature matrix. ``None`` uses
                ``adata.X``.
            use_rep: Key in ``adata.obsm`` to use instead of marker expression
                (e.g. ``"X_umap"``). When provided, ``layer`` is ignored.
            metric: ``"cosine"`` uses ``score_dual`` on binned features;
                ``"rbf"`` builds a Gaussian kernel matrix and calls
                ``score_K``.
            obs_key: Storage key: ``adata.obs`` column (per-cell) or
                ``adata.uns["vendi"]`` key (per-group).
            groupby: ``adata.obs`` column(s) to group by (e.g. ``"leiden"``,
                ``["condition", "cell_cycle"]``). When a list is given, one
                score is produced per unique combination of values.
            sigma: Bandwidth for the RBF kernel. Defaults to the median
                pairwise distance in the data.
            n_bins: Number of uniform bins for ``KBinsDiscretizer``.
            n_reps: Bootstrap repetitions. ``1`` scores once without
                bootstrapping.
            m: Subsample size per repetition. Defaults to pool size
                (``k + 1`` per-cell, group size per-group).
            random_state: Seed for the bootstrap RNG.
            return_eigenvalues: Per-group mode only. When ``True``, also
                compute the eigenvalue spectrum of the dual kernel for each
                group (using the full binned pool, no subsampling) and return
                ``(DataFrame, {group: eigenvalues_array})``. Eigenvalues are
                also stored in ``adata.uns["vendi"][obs_key + "_eigenvalues"]``.
            inplace: If ``True``, persist the updated AnnData to disk.

        Returns:
            Per-cell: ``np.ndarray`` of shape ``(n_cells,)`` stored in
            ``adata.obs[obs_key]``.

            Per-group (``return_eigenvalues=False``): ``pd.DataFrame`` indexed
            by group label, stored in ``adata.uns["vendi"][obs_key]``.

            Per-group (``return_eigenvalues=True``): ``(DataFrame, dict)``
            where the dict maps each group label to a sorted
            ``np.ndarray`` of eigenvalues (ascending).

        Raises:
            ValueError: If ``metric`` is invalid, ``groupby`` column is
                missing, or the requested ``layer``/``use_rep`` does not exist.
            RunNotIngestedError: If the run has not been ingested yet.
        """
        if metric not in {"cosine", "rbf"}:
            raise ValueError("metric must be 'cosine' or 'rbf'")

        try:
            from sklearn.preprocessing import KBinsDiscretizer
            from vendi_score.vendi import score_dual, score_K
        except ImportError as exc:
            raise ImportError(
                "vendi_score requires vendi-score and scikit-learn. "
                "Install with: pip install vendi-score scikit-learn"
            ) from exc

        self.require_ingested()
        adata = self.read_adata()

        # --- Extract feature matrix ---
        if use_rep is not None:
            if use_rep not in adata.obsm:
                raise ValueError(
                    f"use_rep '{use_rep}' not found in adata.obsm. "
                    f"Available: {list(adata.obsm.keys())}"
                )
            X = np.asarray(adata.obsm[use_rep], dtype=np.float64)
        else:
            if layer is not None and layer not in adata.layers:
                raise ValueError(
                    f"layer '{layer}' not found in adata.layers. "
                    f"Available: {list(adata.layers.keys())}"
                )
            raw = adata.layers[layer] if layer is not None else adata.X
            X = np.asarray(
                raw.toarray() if hasattr(raw, "toarray") else raw,
                dtype=np.float64,
            )
            if markers is not None:
                var_names = adata.var_names.tolist()
                missing = [mk for mk in markers if mk not in var_names]
                if missing:
                    raise ValueError(
                        f"Markers not found in adata.var_names: {missing}"
                    )
                col_idx = [var_names.index(mk) for mk in markers]
                X = X[:, col_idx]

        rng = np.random.default_rng(random_state)
        binner = KBinsDiscretizer(
            n_bins=n_bins, strategy="uniform", encode="ordinal"
        )

        def _score_pool(pool: np.ndarray, m_size: int) -> float:
            import warnings as _w
            sub = pool[rng.choice(len(pool), m_size, replace=False)] if m_size < len(pool) else pool
            with _w.catch_warnings():
                _w.filterwarnings("ignore", message=".*constant.*", category=UserWarning)
                MM = binner.fit_transform(sub)
            if metric == "cosine":
                return float(score_dual(MM, normalize=True))
            diff = MM[:, np.newaxis, :] - MM[np.newaxis, :, :]
            sq_dist = np.sum(diff**2, axis=-1)
            K = np.exp(-sq_dist / (2.0 * sigma**2))
            return float(score_K(K))

        def _run_reps(pool: np.ndarray, m_size: int) -> tuple[float, float, float]:
            if n_reps == 1:
                return _score_pool(pool, m_size), float("nan"), float("nan")
            vals = np.array([_score_pool(pool, m_size) for _ in range(n_reps)])
            lo, hi = np.quantile(vals, [0.025, 0.975])
            return float(vals.mean()), float(lo), float(hi)

        def _eigenvalues_pool(pool: np.ndarray) -> np.ndarray:
            """Eigenvalues of the dual kernel on the full binned pool (ascending)."""
            import scipy.linalg as _la
            from sklearn.preprocessing import normalize as _normalize
            MM = binner.fit_transform(pool)
            if metric == "cosine":
                MM = _normalize(MM, axis=1)
                S = MM.T @ MM / len(MM)
            else:
                diff = MM[:, np.newaxis, :] - MM[np.newaxis, :, :]
                sq_dist = np.sum(diff**2, axis=-1)
                S = np.exp(-sq_dist / (2.0 * sigma**2)) / len(MM)
            return _la.eigvalsh(S)

        # ------------------------------------------------------------------ #
        # Per-group mode                                                       #
        # ------------------------------------------------------------------ #
        if groupby is not None:
            groupby_cols = [groupby] if isinstance(groupby, str) else list(groupby)
            missing_cols = [c for c in groupby_cols if c not in adata.obs.columns]
            if missing_cols:
                raise ValueError(
                    f"groupby columns not found in adata.obs: {missing_cols}. "
                    f"Available: {list(adata.obs.columns)}"
                )

            if metric == "rbf" and sigma is None:
                sample_idx = rng.choice(len(X), min(2000, len(X)), replace=False)
                Xs = X[sample_idx]
                flat_dists = np.linalg.norm(
                    Xs[:, np.newaxis, :] - Xs[np.newaxis, :, :], axis=-1
                ).ravel()
                flat_dists = flat_dists[flat_dists > 0]
                sigma = float(np.median(flat_dists)) if len(flat_dists) else 1.0

            if len(groupby_cols) == 1:
                groups = adata.obs[groupby_cols[0]]
                group_labels = (
                    groups.cat.categories.tolist()
                    if hasattr(groups, "cat")
                    else sorted(groups.unique().tolist())
                )
            else:
                groups = adata.obs[groupby_cols].apply(tuple, axis=1)
                group_labels = sorted(groups.unique().tolist())

            try:
                from tqdm.auto import tqdm as _tqdm
            except ImportError:
                def _tqdm(iterable, **_):
                    return iterable

            rows: dict = {}
            eigenvalues: dict[str, np.ndarray] = {}
            for label in _tqdm(group_labels, desc="Vendi score (groups)", unit="group"):
                pool = X[(groups == label).values]
                if len(pool) == 0:
                    continue
                m_size = min(m, len(pool)) if m is not None else len(pool)
                if m_size < 1:
                    continue
                mean_v, lo, hi = _run_reps(pool, m_size)
                rows[label] = {"vendi_score": mean_v, "ci_low": lo, "ci_high": hi}
                if return_eigenvalues:
                    # use string key for eigenvalues storage (tuples aren't JSON-safe)
                    eig_key = "|".join(str(v) for v in label) if isinstance(label, tuple) else str(label)
                    eigenvalues[eig_key] = _eigenvalues_pool(pool)

            result_df = pd.DataFrame.from_dict(rows, orient="index")
            if len(groupby_cols) == 1:
                result_df.index.name = groupby_cols[0]
            else:
                result_df.index = pd.MultiIndex.from_tuples(
                    result_df.index.tolist(), names=groupby_cols
                )

            # Serialize for uns storage: convert tuple keys to pipe-separated strings
            def _uns_key(lbl) -> str:
                return "|".join(str(v) for v in lbl) if isinstance(lbl, tuple) else str(lbl)

            vendi_uns = adata.uns.setdefault("vendi", {})
            vendi_uns[obs_key] = {_uns_key(lbl): v for lbl, v in rows.items()}
            if return_eigenvalues:
                vendi_uns[f"{obs_key}_eigenvalues"] = {
                    k: v.tolist() for k, v in eigenvalues.items()
                }

            if inplace:
                self._try_save_inplace(adata)
            else:
                self._adata = adata

            if return_eigenvalues:
                return result_df, eigenvalues
            return result_df

        # ------------------------------------------------------------------ #
        # Per-cell mode                                                        #
        # ------------------------------------------------------------------ #
        n_cells = X.shape[0]
        if k >= n_cells:
            raise ValueError(
                f"k={k} must be less than the number of cells ({n_cells})"
            )

        try:
            from scipy.spatial import KDTree
        except ImportError as exc:
            raise ImportError("vendi_score requires scipy") from exc

        # For cosine KNN, normalize rows so Euclidean distance ≡ cosine distance.
        if metric == "cosine":
            norms = np.linalg.norm(X, axis=1, keepdims=True)
            norms[norms == 0] = 1.0
            X_feat = X / norms
        else:
            X_feat = X

        tree = KDTree(X_feat)
        query_dists, knn_idx = tree.query(X_feat, k=k + 1)

        if metric == "rbf" and sigma is None:
            knn_dists = query_dists[:, 1:]
            median_dist = float(np.median(knn_dists[knn_dists > 0]))
            sigma = median_dist if median_dist > 0 else 1.0

        neigh_size = k + 1
        m_eff = min(m, neigh_size) if m is not None else neigh_size
        store_ci = n_reps > 1

        try:
            from tqdm.auto import tqdm as _tqdm
        except ImportError:
            def _tqdm(iterable, **_):
                return iterable

        scores = np.empty(n_cells, dtype=np.float64)
        if store_ci:
            ci_low = np.empty(n_cells, dtype=np.float64)
            ci_high = np.empty(n_cells, dtype=np.float64)

        for i in _tqdm(range(n_cells), desc="Vendi score (cells)", unit="cell"):
            mean_v, lo, hi = _run_reps(X[knn_idx[i]], m_eff)
            scores[i] = mean_v
            if store_ci:
                ci_low[i] = lo
                ci_high[i] = hi

        adata.obs[obs_key] = scores
        if store_ci:
            adata.obs[f"{obs_key}_ci_low"] = ci_low
            adata.obs[f"{obs_key}_ci_high"] = ci_high

        if inplace:
            self._try_save_inplace(adata)
        else:
            self._adata = adata

        return scores

    # ------------------------------------------------------------------ #
    # Random Matrix Theory spectrum                                        #
    # ------------------------------------------------------------------ #

    @staticmethod
    def _rmt_one_group(
        X,
        matrix,
        sigma_sq,
        compute_eigvecs=False,
    ):
        """Eigenvalue spectrum and Marchenko-Pastur bounds for one cell group.

        Args:
            X: Expression matrix (n_cells x p_markers), already centred and
                optionally standardised.
            matrix: ``"marker_cov"`` for the p x p covariance matrix;
                ``"cell_gram"`` for the n x n Gram matrix X @ X.T / p.
            sigma_sq: Noise variance for the MP distribution.  ``None``
                estimates it as ``trace(C) / p``.
            compute_eigvecs: When ``True``, include eigenvectors under
                ``"eigvecs"`` (columns in descending eigenvalue order).
                Only meaningful for ``"marker_cov"``.

        Returns:
            Dict with keys ``eigenvalues`` (descending ndarray), ``eigvecs``
            (optional p x p ndarray), ``lambda_max_mp``, ``lambda_min_mp``,
            ``n_signal``, ``q``, ``sigma_sq``, ``n``, ``p``.
        """
        import scipy.linalg as _la

        n, p = X.shape
        if matrix == "marker_cov":
            C = X.T @ X / n
            if compute_eigvecs:
                eigvals_asc, eigvecs_asc = _la.eigh(C)
                eigvals = eigvals_asc[::-1].copy()
                eigvecs = eigvecs_asc[:, ::-1].copy()
            else:
                eigvals = _la.eigvalsh(C)[::-1].copy()
                eigvecs = None
            q = p / n
            s2 = sigma_sq if sigma_sq is not None else float(np.trace(C) / p)
        else:
            G = X @ X.T / p
            eigvals = _la.eigvalsh(G)[::-1].copy()
            eigvecs = None
            q = p / n
            effective_scale = n / p
            s2 = (sigma_sq if sigma_sq is not None else 1.0) * effective_scale

        lam_max = s2 * (1.0 + np.sqrt(q)) ** 2
        lam_min = s2 * max(0.0, 1.0 - np.sqrt(q)) ** 2
        n_signal = int(np.sum(eigvals > lam_max))

        out = dict(
            eigenvalues=eigvals,
            lambda_max_mp=float(lam_max),
            lambda_min_mp=float(lam_min),
            n_signal=n_signal,
            q=float(q),
            sigma_sq=float(s2),
            n=int(n),
            p=int(p),
        )
        if compute_eigvecs and eigvecs is not None:
            out["eigvecs"] = eigvecs
        return out

    def compute_rmt_spectrum(
        self,
        markers=None,
        layer=None,
        groupby=None,
        matrix="marker_cov",
        standardize=True,
        n_cells=10_000,
        sigma_sq=None,
        random_state=0,
        uns_key="rmt_spectrum",
        plot=False,
        inplace=True,
    ):
        """Compute the marker eigenvalue spectrum and compare to the Marchenko-Pastur null.

        Uses Random Matrix Theory (RMT) to identify which principal components
        of the marker covariance (or cell Gram) matrix represent genuine
        biological signal vs. sampling noise.  Eigenvalues above the MP upper
        bound ``lambda_max = sigma2 * (1 + sqrt(q))^2`` are classed as signal.

        The Marchenko-Pastur law states that for a random n x p matrix with
        i.i.d. entries of variance sigma2, the eigenvalues of the sample
        covariance ``C = X.T @ X / n`` lie within
        ``[sigma2*(1-sqrt(q))^2, sigma2*(1+sqrt(q))^2]`` with ``q = p/n``.

        Args:
            markers: Subset of marker names to use.  ``None`` uses all markers.
            layer: AnnData layer to use.  ``None`` uses ``adata.X``.
            groupby: ``adata.obs`` column to group by.  When set, the spectrum
                is computed separately per group.
            matrix: ``"marker_cov"`` (default) builds the p x p marker
                covariance matrix.  ``"cell_gram"`` builds the n x n cell Gram
                matrix ``X @ X.T / p`` on a subsampled cell set.
            standardize: If ``True``, centre and scale each marker to zero mean
                and unit variance before analysis.  For ``"marker_cov"`` this
                yields a correlation matrix (sigma2 = 1 under the null).
            n_cells: Maximum cells per group.  Randomly subsampled when larger.
            sigma_sq: Noise variance for the MP distribution.  ``None``
                auto-estimates as ``trace(C) / p``.
            random_state: Seed for subsampling.
            uns_key: Key under which results are stored in ``adata.uns``.
            plot: If ``True``, return ``(result, fig)`` with a scree plot.
            inplace: If ``True``, persist the updated AnnData to disk.

        Returns:
            Result dict (or ``(dict, Figure)`` when ``plot=True``).  Contains
            ``eigenvalues``, ``lambda_max_mp``, ``lambda_min_mp``,
            ``n_signal``, ``q``, ``sigma_sq``, ``n``, ``p``, ``markers``,
            ``matrix``.  With ``groupby``, nested under ``"groups"``.

        Raises:
            ValueError: Invalid parameter values or missing data.
            RunNotIngestedError: Run has not been ingested.
        """
        if matrix not in {"marker_cov", "cell_gram"}:
            raise ValueError("matrix must be 'marker_cov' or 'cell_gram'")

        self.require_ingested()
        adata = self.read_adata()

        if layer is not None and layer not in {"X", "x"} and layer not in adata.layers:
            raise ValueError(
                f"layer '{layer}' not found. Available: {list(adata.layers.keys())}"
            )
        raw = adata.layers[layer] if layer is not None else adata.X
        X_full = np.asarray(
            raw.toarray() if hasattr(raw, "toarray") else raw,
            dtype=np.float64,
        )
        if markers is not None:
            var_names = adata.var_names.tolist()
            missing = [mk for mk in markers if mk not in var_names]
            if missing:
                raise ValueError(f"Markers not found in adata.var_names: {missing}")
            col_idx = [var_names.index(mk) for mk in markers]
            X_full = X_full[:, col_idx]
            marker_names = list(markers)
        else:
            marker_names = adata.var_names.tolist()

        rng = np.random.default_rng(random_state)

        def _prep(X):
            if n_cells is not None and len(X) > n_cells:
                idx = rng.choice(len(X), n_cells, replace=False)
                X = X[idx]
            X = X - X.mean(axis=0)
            if standardize:
                std = X.std(axis=0)
                std[std == 0] = 1.0
                X = X / std
            return X

        def _one(X):
            res = self._rmt_one_group(_prep(X), matrix=matrix, sigma_sq=sigma_sq)
            return {**res, "eigenvalues": res["eigenvalues"].tolist()}

        if groupby is not None:
            if groupby not in adata.obs.columns:
                raise ValueError(
                    f"groupby '{groupby}' not found in adata.obs. "
                    f"Available: {list(adata.obs.columns)}"
                )
            groups = adata.obs[groupby]
            group_labels = (
                groups.cat.categories.tolist()
                if hasattr(groups, "cat")
                else sorted(groups.unique().tolist())
            )
            group_results = {}
            for label in group_labels:
                pool = X_full[(groups == label).values]
                if len(pool) < 2:
                    continue
                group_results[label] = _one(pool)
            result = {
                "matrix": matrix,
                "groupby": groupby,
                "markers": marker_names,
                "groups": group_results,
            }
        else:
            res = _one(X_full)
            result = {"matrix": matrix, "markers": marker_names, **res}

        adata.uns[uns_key] = result

        if inplace:
            self._try_save_inplace(adata)
        else:
            self._adata = adata

        if not plot:
            return result

        try:
            import matplotlib.pyplot as plt
        except ImportError as exc:
            raise ImportError("plot=True requires matplotlib") from exc

        def _draw_one(ax, gres, title):
            eigvals = np.asarray(gres["eigenvalues"])
            lam_max = gres["lambda_max_mp"]
            lam_min = gres["lambda_min_mp"]
            ranks = np.arange(1, len(eigvals) + 1)
            colors = ["#E06C2B" if v > lam_max else "#AAAAAA" for v in eigvals]
            ax.bar(ranks, eigvals, color=colors, width=0.8, zorder=2)
            ax.axhline(lam_max, color="#333333", linestyle="--", linewidth=1.2,
                       label=f"MP upper bound ({lam_max:.3f})")
            if lam_min > 0:
                ax.axhline(lam_min, color="#888888", linestyle=":", linewidth=1.0,
                           label=f"MP lower bound ({lam_min:.3f})")
            ax.set_xlabel("Component rank")
            ax.set_ylabel("Eigenvalue")
            ax.set_title(f"{title}  (n_signal = {gres['n_signal']})")
            ax.legend(fontsize=8)
            ax.set_xlim(0.5, len(eigvals) + 0.5)

        if groupby is not None:
            items = list(result["groups"].items())
            ncols = min(3, len(items))
            nrows = (len(items) + ncols - 1) // ncols
            fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows),
                                     squeeze=False)
            for idx, (label, gres) in enumerate(items):
                _draw_one(axes[idx // ncols][idx % ncols], gres, str(label))
            for idx in range(len(items), nrows * ncols):
                axes[idx // ncols][idx % ncols].set_visible(False)
            q_repr = list(result["groups"].values())[0]["q"] if result["groups"] else float("nan")
        else:
            fig, ax = plt.subplots(figsize=(max(6, len(result["eigenvalues"]) * 0.4), 4))
            _draw_one(ax, result, "All cells")
            q_repr = result["q"]

        fig.suptitle(f"RMT eigenvalue spectrum — {matrix}  (q = {q_repr:.4f})", fontsize=10)
        fig.tight_layout()
        return result, fig

    def bootstrap_rmt_spectrum(
        self,
        markers=None,
        layer=None,
        groupby=None,
        matrix="marker_cov",
        standardize=True,
        n_cells=10_000,
        frac=0.8,
        n_bootstrap=200,
        random_state=0,
        uns_key="rmt_bootstrap",
        plot=False,
        inplace=True,
    ):
        """Bootstrap stability of the RMT eigenvalue spectrum.

        Repeatedly subsamples a fraction of cells and recomputes the
        eigenvalue spectrum to quantify estimation uncertainty.

        - **Eigenvalue distributions**: spread of each rank's eigenvalue
          across resamples — tight = well-estimated, wide = noisy.
        - **Eigenvector stability** (``matrix="marker_cov"`` only): mean
          absolute cosine similarity between the reference eigenvector and
          each bootstrap replicate.  Near 1 = reproducible; near 0 = unstable.

        These complement the MP boundary from :meth:`compute_rmt_spectrum`:
        RMT says "above the noise floor"; bootstrap says "stably estimated".

        Args:
            markers: Subset of marker names to use.  ``None`` uses all.
            layer: AnnData layer to use.  ``None`` uses ``adata.X``.
            groupby: ``adata.obs`` column to group by.
            matrix: ``"marker_cov"`` or ``"cell_gram"``.
            standardize: Centre and scale each marker before analysis.
            n_cells: Maximum cell pool size per group before bootstrapping.
            frac: Fraction of the pool drawn without replacement per replicate.
            n_bootstrap: Number of bootstrap resamples.
            random_state: Seed for the RNG.
            uns_key: Key under which results are stored in ``adata.uns``.
            plot: If ``True``, return ``(result, fig)`` with a two-panel
                figure: eigenvalue distributions and eigenvector stability.
            inplace: If ``True``, persist the updated AnnData to disk.

        Returns:
            Result dict (or ``(dict, Figure)`` when ``plot=True``).  Contains
            ``eigenvalue_ref``, ``eigenvalue_matrix`` (n_bootstrap x p),
            ``eigenvalue_mean``, ``eigenvalue_std``, ``eigenvalue_ci_low``,
            ``eigenvalue_ci_high``, ``eigenvector_stability`` (marker_cov
            only), ``lambda_max_mp_ref``, ``lambda_max_mp_distribution``,
            ``n_signal_ref``, ``markers``, ``matrix``.
            With ``groupby``, nested under ``"groups"``.

        Raises:
            ValueError: Invalid parameter values or missing data.
            RunNotIngestedError: Run has not been ingested.
        """
        if matrix not in {"marker_cov", "cell_gram"}:
            raise ValueError("matrix must be 'marker_cov' or 'cell_gram'")
        if not (0.0 < frac < 1.0):
            raise ValueError("frac must be in (0, 1)")

        self.require_ingested()
        adata = self.read_adata()

        if layer is not None and layer not in {"X", "x"} and layer not in adata.layers:
            raise ValueError(
                f"layer '{layer}' not found. Available: {list(adata.layers.keys())}"
            )
        raw = adata.layers[layer] if layer is not None else adata.X
        X_full = np.asarray(
            raw.toarray() if hasattr(raw, "toarray") else raw,
            dtype=np.float64,
        )
        if markers is not None:
            var_names = adata.var_names.tolist()
            missing = [mk for mk in markers if mk not in var_names]
            if missing:
                raise ValueError(f"Markers not found in adata.var_names: {missing}")
            col_idx = [var_names.index(mk) for mk in markers]
            X_full = X_full[:, col_idx]
            marker_names = list(markers)
        else:
            marker_names = adata.var_names.tolist()

        rng = np.random.default_rng(random_state)

        def _prep_pool(X):
            if n_cells is not None and len(X) > n_cells:
                idx = rng.choice(len(X), n_cells, replace=False)
                X = X[idx]
            X = X - X.mean(axis=0)
            if standardize:
                std = X.std(axis=0)
                std[std == 0] = 1.0
                X = X / std
            return X

        def _bootstrap_one(X):
            import scipy.linalg as _la

            n, p = X.shape
            n_sub = max(2, int(np.floor(frac * n)))

            ref = self._rmt_one_group(
                X, matrix=matrix, sigma_sq=None,
                compute_eigvecs=(matrix == "marker_cov"),
            )
            ref_eigvals = ref["eigenvalues"]
            ref_eigvecs = ref.get("eigvecs")
            lam_max_ref = ref["lambda_max_mp"]
            sigma_sq_ref = ref["sigma_sq"]

            n_comp = p if matrix == "marker_cov" else min(n_sub, p)

            boot_eigvals = np.empty((n_bootstrap, n_comp), dtype=np.float64)
            boot_stab = (
                np.zeros(n_comp, dtype=np.float64)
                if matrix == "marker_cov" else None
            )
            boot_lam_max = np.empty(n_bootstrap, dtype=np.float64)

            for b in range(n_bootstrap):
                idx_b = rng.choice(n, n_sub, replace=False)
                Xb = X[idx_b]

                if matrix == "marker_cov":
                    Cb = Xb.T @ Xb / n_sub
                    eigvals_b_asc, eigvecs_b_asc = _la.eigh(Cb)
                    eigvals_b = eigvals_b_asc[::-1]
                    eigvecs_b = eigvecs_b_asc[:, ::-1]
                    boot_eigvals[b] = eigvals_b[:n_comp]
                    for r in range(n_comp):
                        dot = float(ref_eigvecs[:, r] @ eigvecs_b[:, r])
                        boot_stab[r] += abs(dot)
                    q_b = p / n_sub
                    boot_lam_max[b] = sigma_sq_ref * (1.0 + np.sqrt(q_b)) ** 2
                else:
                    Gb = Xb @ Xb.T / p
                    eigvals_b = _la.eigvalsh(Gb)[::-1]
                    boot_eigvals[b] = eigvals_b[:n_comp]
                    q_b = p / n_sub
                    s2_sub = 1.0 * (n_sub / p)
                    boot_lam_max[b] = s2_sub * (1.0 + np.sqrt(q_b)) ** 2

            eigvec_stability = boot_stab / n_bootstrap if boot_stab is not None else None

            ci_low = np.percentile(boot_eigvals, 2.5, axis=0)
            ci_high = np.percentile(boot_eigvals, 97.5, axis=0)

            out = dict(
                eigenvalue_ref=ref_eigvals[:n_comp].tolist(),
                eigenvalue_matrix=boot_eigvals.tolist(),
                eigenvalue_mean=boot_eigvals.mean(axis=0).tolist(),
                eigenvalue_std=boot_eigvals.std(axis=0).tolist(),
                eigenvalue_ci_low=ci_low.tolist(),
                eigenvalue_ci_high=ci_high.tolist(),
                lambda_max_mp_ref=float(lam_max_ref),
                lambda_max_mp_distribution=boot_lam_max.tolist(),
                n_signal_ref=ref["n_signal"],
                q_ref=float(ref["q"]),
                n=int(n),
                p=int(p),
                n_sub=int(n_sub),
            )
            if eigvec_stability is not None:
                out["eigenvector_stability"] = eigvec_stability.tolist()
            return out

        if groupby is not None:
            if groupby not in adata.obs.columns:
                raise ValueError(
                    f"groupby '{groupby}' not found in adata.obs. "
                    f"Available: {list(adata.obs.columns)}"
                )
            groups = adata.obs[groupby]
            group_labels = (
                groups.cat.categories.tolist()
                if hasattr(groups, "cat")
                else sorted(groups.unique().tolist())
            )
            group_results = {}
            for label in group_labels:
                pool = X_full[(groups == label).values]
                if len(pool) < 4:
                    continue
                group_results[label] = _bootstrap_one(_prep_pool(pool))
            result = {
                "matrix": matrix,
                "groupby": groupby,
                "markers": marker_names,
                "n_bootstrap": n_bootstrap,
                "frac": frac,
                "groups": group_results,
            }
        else:
            bres = _bootstrap_one(_prep_pool(X_full))
            result = {
                "matrix": matrix,
                "markers": marker_names,
                "n_bootstrap": n_bootstrap,
                "frac": frac,
                **bres,
            }

        adata.uns[uns_key] = result

        if inplace:
            self._try_save_inplace(adata)
        else:
            self._adata = adata

        if not plot:
            return result

        try:
            import matplotlib.pyplot as plt
        except ImportError as exc:
            raise ImportError("plot=True requires matplotlib") from exc

        def _draw_group(axes_pair, gres, title):
            ax_eig, ax_stab = axes_pair
            eigvals_mat = np.asarray(gres["eigenvalue_matrix"])
            ref_vals = np.asarray(gres["eigenvalue_ref"])
            ci_low_arr = np.asarray(gres["eigenvalue_ci_low"])
            lam_max_ref = gres["lambda_max_mp_ref"]
            lam_max_dist = np.asarray(gres["lambda_max_mp_distribution"])
            n_comp = eigvals_mat.shape[1]
            ranks = np.arange(1, n_comp + 1)

            for r in range(n_comp):
                col = "#E06C2B" if ci_low_arr[r] > lam_max_ref else "#AAAAAA"
                vp = ax_eig.violinplot(
                    eigvals_mat[:, r], positions=[r + 1], widths=0.7,
                    showmedians=False, showextrema=False,
                )
                for pc in vp["bodies"]:
                    pc.set_facecolor(col)
                    pc.set_alpha(0.7)
            ax_eig.scatter(ranks, ref_vals, color="#222222", s=18, zorder=5,
                           label="Reference eigenvalue")
            ax_eig.axhline(lam_max_ref, color="#333333", linestyle="--",
                           linewidth=1.2, label=f"MP bound ({lam_max_ref:.3f})")
            mp_med = float(np.median(lam_max_dist))
            mp_std = float(np.std(lam_max_dist))
            ax_eig.axhspan(mp_med - mp_std, mp_med + mp_std, color="#CCCCCC",
                           alpha=0.3, label="MP bound +/- 1 SD")
            ax_eig.set_xlabel("Component rank")
            ax_eig.set_ylabel("Eigenvalue")
            ax_eig.set_title(f"{title} — eigenvalue distributions")
            ax_eig.set_xlim(0.5, n_comp + 0.5)
            ax_eig.legend(fontsize=7)

            if "eigenvector_stability" in gres:
                stab = np.asarray(gres["eigenvector_stability"])
                bar_colors = [
                    "#2CA02C" if s >= 0.9 else ("#FFAA00" if s >= 0.7 else "#D62728")
                    for s in stab
                ]
                ax_stab.bar(ranks, stab, color=bar_colors, width=0.8)
                ax_stab.axhline(0.9, color="#333333", linestyle="--",
                                linewidth=1.0, label="Stability threshold (0.9)")
                ax_stab.set_ylim(0, 1.05)
                ax_stab.set_xlabel("Component rank")
                ax_stab.set_ylabel("|cos theta| (mean across bootstrap)")
                ax_stab.set_title(f"{title} — eigenvector stability")
                ax_stab.set_xlim(0.5, n_comp + 0.5)
                ax_stab.legend(fontsize=7)
            else:
                ax_stab.text(
                    0.5, 0.5,
                    "Eigenvector stability\nnot available for cell_gram",
                    ha="center", va="center",
                    transform=ax_stab.transAxes, fontsize=9,
                )

        if groupby is not None:
            items = list(result["groups"].items())
            ngroups = len(items)
            fig, axes = plt.subplots(ngroups, 2, figsize=(12, 4 * ngroups), squeeze=False)
            for i, (label, gres) in enumerate(items):
                _draw_group(axes[i], gres, str(label))
        else:
            fig, axes = plt.subplots(1, 2, figsize=(12, 4))
            _draw_group(axes, result, "All cells")

        fig.suptitle(
            f"Bootstrap RMT spectrum — {matrix}"
            f"  (n_bootstrap={n_bootstrap}, frac={frac})",
            fontsize=10,
        )
        fig.tight_layout()
        return result, fig

    @staticmethod
    def _build_comparison_pairs(order, comparisons):
        """Return list of (group_a, group_b) pairs to test for significance.

        - `None`: no pairs (skip significance brackets).
        - `"all"`: every unordered pair.
        - `"adjacent"`: neighbours in `order`.
        - List of `(a, b)` tuples: explicit pairs.
        """
        if comparisons is None:
            return []
        if comparisons == "all":
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
    def _adjust_pvalues(pvalues, method="bh"):
        """Apply BH-FDR or Bonferroni correction across pairwise p-values."""
        pvals = np.asarray(pvalues, dtype=float)
        valid_mask = ~np.isnan(pvals)
        n_valid = int(valid_mask.sum())
        if n_valid == 0:
            return pvals.tolist()

        if method == "bonferroni":
            adjusted = pvals.copy()
            adjusted[valid_mask] = np.minimum(pvals[valid_mask] * n_valid, 1.0)
            return adjusted.tolist()

        if method == "bh":
            adjusted = pvals.copy()
            valid_idx = np.where(valid_mask)[0]
            order_local = np.argsort(pvals[valid_idx])
            sorted_idx = valid_idx[order_local]
            ranks = np.arange(1, n_valid + 1, dtype=float)
            sorted_p = pvals[sorted_idx]
            bh_vals = sorted_p * n_valid / ranks
            bh_vals = np.minimum(bh_vals, 1.0)
            # Enforce monotonicity from largest to smallest
            bh_vals = np.minimum.accumulate(bh_vals[::-1])[::-1]
            adjusted[sorted_idx] = bh_vals
            return adjusted.tolist()

        raise ValueError(f"Unknown method: {method}")

    @staticmethod
    def _default_significance_thresholds() -> list[tuple[float, str]]:
        """Default (p_threshold, label) ordering used to convert p to stars."""
        return [
            (1e-4, "****"),
            (1e-3, "***"),
            (1e-2, "**"),
            (5e-2, "*"),
        ]

    @staticmethod
    def _pvalue_to_stars(
        p,
        ns_label: str = "ns",
        thresholds: list[tuple[float, str]] | None = None,
    ) -> str:
        """Convert a p-value to a star annotation using configurable thresholds.

        Thresholds are a list of `(p_threshold, label)` tuples. The first
        tuple whose `p_threshold >= p` wins; if none match, `ns_label` is
        returned. Use stricter thresholds first.
        """
        if p is None:
            return ns_label
        try:
            p_f = float(p)
        except (TypeError, ValueError):
            return ns_label
        if np.isnan(p_f):
            return ns_label

        if thresholds is None:
            thresholds = Run._default_significance_thresholds()

        for threshold, label in thresholds:
            if p_f <= float(threshold):
                return label
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
        significance_thresholds: list[tuple[float, str]] | None = None,
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
            label = Run._pvalue_to_stars(
                pvalues[orig_idx],
                ns_label=ns_label,
                thresholds=significance_thresholds,
            )
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
        self._adata = filtered
        self._try_save_inplace(filtered)
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
            self._try_save_inplace(adata)

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

    @staticmethod
    def _balanced_group_indices(
        groups: pd.Series,
        n_per_group: int | None = None,
        random_state: int = 0,
        replace: bool = False,
    ) -> tuple[np.ndarray, dict[str, int], int]:
        """Return balanced per-group indices with optional target size."""
        group_series = groups.astype(str)
        counts = group_series.value_counts().sort_index()
        if counts.empty:
            raise ValueError("group_series is empty")

        count_dict = {str(group): int(count) for group, count in counts.items()}
        min_count = int(counts.min())

        if n_per_group is None:
            target_n = min_count
        else:
            target_n = int(n_per_group)
            if target_n <= 0:
                raise ValueError("n_per_group must be > 0")
            if not replace and target_n > min_count:
                raise ValueError(
                    "n_per_group exceeds smallest group size. "
                    "Set replace=True or omit n_per_group to use the minimum."
                )

        rng = np.random.default_rng(random_state)
        selected: list[np.ndarray] = []

        for group in sorted(count_dict.keys()):
            group_idx = np.where(group_series.values == group)[0]
            if len(group_idx) == target_n and not replace:
                selected.append(group_idx)
            else:
                selected.append(
                    rng.choice(group_idx, size=target_n, replace=replace)
                )

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
            self._try_save_inplace(adata)

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

    def gate_cell_cycle(
        self,
        marker_map: dict[str, str] | None = None,
        layer: str | None = None,
        thresholds: dict[str, float] | None = None,
        quantile_thresholds: dict[str, float] | None = None,
        threshold_methods: dict[str, str] | None = None,
        inplace: bool = True,
    ) -> dict[str, Any]:
        """Assign cells to cell-cycle phases using hierarchical marker gating.

        Uses arcsinh-transformed data (``adata.X`` by default, or a named
        layer). Do **not** pass the ``zscore`` layer — gating thresholds are
        calibrated for arcsinh-scale values.

        Gating hierarchy (mutually exclusive, applied in order):
          IdU+      → S_phase
          pH3+      → M_phase
          CyclinB1+ → G2_phase
          pRb+      → Cycling_G1   (if pRb in marker_map)
          pRb-      → G0_or_quiescent
          remaining → G1_or_quiescent (when pRb absent)

        Args:
            marker_map: Dict mapping role → actual ``var_name`` column.
                Required roles: ``IdU``, ``pH3``, ``CyclinB1``.
                Optional: ``pRb``. If ``None``, auto-detection is attempted
                from ``adata.var_names`` using known aliases.
            layer: AnnData layer to read expression from. ``None`` uses
                ``adata.X`` (typically arcsinh-transformed after ingestion).
                Avoid ``"zscore"`` — thresholds are not calibrated for it.
            thresholds: User-supplied thresholds by role. Missing roles use
                the automatic method.
            quantile_thresholds: Override the default quantile per role
                (used when strategy is ``"quantile"``).
            threshold_methods: Per-role auto-threshold strategy.
                ``"otsu"`` (default for IdU/pH3/CyclinB1) finds the valley in
                bimodal distributions. ``"quantile"`` (default for pRb) uses
                a fixed percentile.
            inplace: Persist updated AnnData to disk after gating.

        Returns:
            Dict with:
            - ``gated_df``: cell-level DataFrame with gate columns and
              ``cell_cycle_phase``.
            - ``summary``: phase counts and fractions.
            - ``thresholds``: thresholds used (role → value).

        Raises:
            ValueError: If required markers cannot be found in ``adata.var_names``.
        """
        from cytofstandard.cell_cycle import (
            auto_detect_cell_cycle_markers,
            REQUIRED_ROLES,
            gate_cell_cycle as _gate_cc,
        )

        self.require_ingested()
        adata = self.read_adata()

        if marker_map is None:
            detected = auto_detect_cell_cycle_markers(list(adata.var_names))
            missing_required = [r for r in REQUIRED_ROLES if r not in detected]
            if missing_required:
                raise ValueError(
                    f"Could not auto-detect required cell-cycle markers. "
                    f"Missing: {missing_required}. "
                    f"Pass marker_map explicitly, e.g. "
                    f"marker_map={{'IdU': 'IdU', 'pH3': 'pH3', 'CyclinB1': 'CyclinB1', ...}}"
                )
            marker_map = detected

        result = _gate_cc(
            adata,
            marker_map=marker_map,
            layer=layer,
            thresholds=thresholds,
            quantile_thresholds=quantile_thresholds,
            threshold_methods=threshold_methods,
            apply_arcsinh=False,
            return_adata=True,
        )

        adata = result["adata"]
        self._adata = adata

        if "cell_cycle_gating" not in adata.uns:
            adata.uns["cell_cycle_gating"] = {}
        from datetime import datetime as _dt
        adata.uns["cell_cycle_gating"]["latest"] = {
            "timestamp": _dt.utcnow().isoformat(),
            "marker_map": marker_map,
            "thresholds": result["thresholds"],
            "layer": layer,
            "n_cells": int(adata.n_obs),
        }

        if inplace:
            self._try_save_inplace(adata)

        self._log_run_event(
            "cell_cycle_gated",
            {
                "marker_map": marker_map,
                "thresholds": result["thresholds"],
                "layer": layer,
                "n_cells": int(adata.n_obs),
                "phase_counts": {
                    row["cell_cycle_phase"]: int(row["n_cells"])
                    for _, row in result["summary"].iterrows()
                },
            },
        )

        return {
            "gated_df": result["gated_df"],
            "summary": result["summary"],
            "thresholds": result["thresholds"],
        }

    def add_cell_cycle_pseudotime(
        self,
        marker_cols: dict[str, str] | None = None,
        phase_col: str = "cell_cycle_phase",
        output_col: str = "cell_cycle_pseudotime",
        angle_col: str = "cell_cycle_angle",
        on_cycle_col: str = "cell_cycle_on_cycle",
        phase_order: list[str] | None = None,
        phase_widths: dict[str, float] | None = None,
        overwrite: bool = False,
        inplace: bool = True,
    ) -> "anndata.AnnData":
        """Add a continuous cell-cycle pseudotime coordinate to this run.

        Requires that :meth:`gate_cell_cycle` has already been run so that
        ``adata.obs["cell_cycle_phase"]`` exists.

        Within each phase, ordering is based on marker intensity ranks:

        * **G0/G1 phases** — ``pRb`` (increases as CDK4/6 phosphorylate Rb)
        * **S phase** — ``DNA`` (content increases through replication)
        * **G2 phase** — ``CyclinB1`` (accumulates before M entry)
        * **M phase** — ``pH3`` (peaks during mitosis)

        Args:
            marker_cols: Dict mapping role → actual column name, e.g.
                ``{"pRb": "pRb", "DNA": "DNA", "CyclinB1": "CyclinB1",
                "pH3": "pH3", "IdU": "IdU"}``.
                If ``None``, the gating marker map stored in
                ``adata.uns["cell_cycle_gating"]`` is used for the gating
                markers (IdU/pH3/CyclinB1/pRb); **DNA must still be added
                manually** if it is needed for within-S ordering.
            phase_col: Column in ``adata.obs`` with categorical phase labels.
            output_col: Output column name for pseudotime.
            angle_col: Output column name for the 2π angle.
            on_cycle_col: Output column for the boolean on-cycle flag
                (``True`` for G1/S/G2/M, ``False`` for G0/quiescent).
            phase_order: Custom biological ordering of phase labels.
            phase_widths: Dict of phase → fractional arc width (summed to 1).
            overwrite: Overwrite existing pseudotime columns. Default False.
            inplace: Persist updated AnnData to disk. Default True.

        Returns:
            Updated AnnData object (with pseudotime columns in ``obs``).

        Example::

            run.add_cell_cycle_pseudotime(
                marker_cols={
                    "pRb": "pRb", "IdU": "IdU",
                    "CyclinB1": "CyclinB1", "pH3": "pH3",
                    "DNA": "DNA",
                },
                overwrite=True,
            )
        """
        from cytofstandard.cell_cycle import (
            add_cell_cycle_pseudotime as _add_pt,
        )

        self.require_ingested()
        adata = self.read_adata()

        # Auto-fill from stored gating marker_map if not provided
        if marker_cols is None:
            gating = adata.uns.get("cell_cycle_gating", {}).get("latest", {})
            stored = gating.get("marker_map", {})
            if not stored:
                raise ValueError(
                    "marker_cols not provided and no stored cell_cycle_gating "
                    "marker_map found in adata.uns. Run gate_cell_cycle first, "
                    "or pass marker_cols explicitly."
                )
            marker_cols = dict(stored)

        adata = _add_pt(
            adata,
            phase_col=phase_col,
            marker_cols=marker_cols,
            output_col=output_col,
            angle_col=angle_col,
            on_cycle_col=on_cycle_col,
            phase_order=phase_order,
            phase_widths=phase_widths,
            overwrite=overwrite,
            copy=False,
        )

        self._adata = adata

        if inplace:
            self._try_save_inplace(adata)

        self._log_run_event(
            "cell_cycle_pseudotime",
            {
                "phase_col": phase_col,
                "output_col": output_col,
                "marker_cols": marker_cols,
                "n_cells": int(adata.n_obs),
            },
        )

        return adata

    def open_cell_cycle_app(self, port: int = 8502) -> "subprocess.Popen":
        """Launch the cell-cycle gating Streamlit app pre-loaded with this run.

        Opens http://localhost:<port> in the browser automatically.
        The app shows Auto-threshold and Interactive-slider tabs — click
        "Apply Gating" to write the chosen thresholds back to the stored adata.

        Afterwards, call the Section 9 cell in the notebook to pick up the
        saved thresholds from ``adata.uns["cell_cycle_gating"]["latest"]``.

        Args:
            port: Local port to serve on (default 8502 to avoid clashing with
                  the full CyTOF app on 8501).

        Returns:
            The ``subprocess.Popen`` handle — call ``proc.terminate()`` to stop.
        """
        import os as _os
        import subprocess
        import sys
        from pathlib import Path as _Path

        app_script = _Path(__file__).parent.parent / "app" / "cell_cycle_app.py"
        if not app_script.exists():
            raise FileNotFoundError(f"cell_cycle_app.py not found at {app_script}")

        env = _os.environ.copy()
        env["CYTOF_PROJECT_PATH"] = str(self.project.path)
        env["CYTOF_RUN_ID"]       = self.run_id

        proc = subprocess.Popen(
            [
                sys.executable, "-m", "streamlit", "run", str(app_script),
                "--server.port", str(port),
                "--server.headless", "false",
                "--browser.gatherUsageStats", "false",
            ],
            env=env,
        )
        print(f"Cell-cycle app launched (PID {proc.pid})  →  http://localhost:{port}")
        print("Call proc.terminate() to stop it.")
        return proc

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

    def save_as_new_run(
        self,
        adata: anndata.AnnData,
        new_run_id: str,
        run_name: str | None = None,
        notes: str | None = None,
    ) -> "Run":
        """Persist an arbitrary AnnData as a brand-new run in the same project.

        Unlike ``save_adata(adata)``, which overwrites the *current* run, this
        method registers a fresh run and writes the provided AnnData to it.
        Panel metadata (panel_id, instrument, …) is copied from the source run.

        Typical use case — save a balanced subsample without destroying the full
        run's data::

            adata_sub = run.subsample_by_group("sample_id", n_per_group=5000)
            sub_run = run.save_as_new_run(adata_sub, new_run_id="CyTOF1_sub")

        Args:
            adata: AnnData to persist (e.g. from ``subsample_by_group``).
            new_run_id: ID for the new run.  Must not already exist in the project.
            run_name: Optional display name; defaults to
                ``"<source run name> (subset)"``.
            notes: Optional free-text notes attached to the new run's metadata.

        Returns:
            The newly created and saved ``Run`` object.
        """
        new_run = self.project.add_run(
            run_id=new_run_id,
            run_name=run_name
            or f"{self.run_config.get('run_name', self.run_id)} (subset)",
            panel_id=self.run_config.get("panel_id"),
            acquisition_date=self.run_config.get("acquisition_date"),
            instrument=self.run_config.get("instrument"),
            operator=self.run_config.get("operator"),
            notes=notes,
        )
        adata = adata.copy()
        adata.uns["derived_from"] = {
            "project_id": self.project.project_id,
            "source_run_id": self.run_id,
            "created_at": datetime.utcnow().isoformat(),
        }
        new_run._adata = adata
        new_run.save()
        self._log_run_event(
            "saved_as_new_run",
            {
                "new_run_id": new_run_id,
                "n_cells_source": int(self.read_adata().n_obs),
                "n_cells_new": int(adata.n_obs),
            },
        )
        new_run._log_run_event(
            "run_initialized_from_adata",
            {"source_run_id": self.run_id, "n_cells": int(adata.n_obs)},
        )
        return new_run

    def zarr_path(self) -> Path:
        """Get the path to the Zarr file.

        Returns:
            Path to Zarr file
        """
        return self.path / "processed" / f"{self.run_id}.zarr"

    @staticmethod
    def _strip_hidden_files(path: Path) -> None:
        """Remove dot-files (e.g. .DS_Store) that block rmdir on macOS/Dropbox."""
        import os as _os
        for root, _dirs, files in _os.walk(str(path)):
            for fname in files:
                if fname.startswith("."):
                    try:
                        _os.unlink(_os.path.join(root, fname))
                    except OSError:
                        pass

    def _write_zarr(
        self, adata: anndata.AnnData, zarr_path: Path | None = None
    ) -> None:
        """Write *adata* to zarr, converting PermissionErrors to ZarrLockedError.

        On macOS/Dropbox, ``shutil.rmtree`` inside zarr can fail with
        ``OSError: ENOTEMPTY`` because Finder or Dropbox inserts dot-files
        (e.g. ``.DS_Store``) between the recursive scan and the final
        ``os.rmdir`` call.  We proactively strip those files before writing and
        also catch-and-retry if they reappear during the write.
        """
        import errno as _errno

        target = zarr_path or self.zarr_path()
        # Proactively remove dot-files so zarr's shutil.rmtree succeeds
        if target.exists():
            self._strip_hidden_files(target)
        try:
            adata.write_zarr(target)
        except PermissionError as exc:
            locked: list[str] = []
            if target.exists():
                try:
                    locked = get_locked_zarr_parts(target)
                except Exception:
                    pass
            parts_str = (
                ", ".join(f'"{p}"' for p in locked) if locked else "unknown parts"
            )
            raise ZarrLockedError(
                f"Cannot write to run '{self.run_id}': zarr store is read-only "
                f"(locked parts: {parts_str}). "
                f"Call run.unlock_zarr_parts() to restore write access, "
                f"or run.locked_zarr_parts() to inspect what is locked."
            ) from exc
        except OSError as exc:
            if exc.errno != _errno.ENOTEMPTY:
                raise
            # Dot-files reappeared during the write — strip again and retry once.
            if target.exists():
                self._strip_hidden_files(target)
            adata.write_zarr(target)

    def _try_save_inplace(self, adata: anndata.AnnData | None = None) -> bool:
        """Persist *adata* to disk; on ZarrLockedError warn and keep in memory.

        Used by computation methods (``compute_umap``, ``cluster_leiden``, …)
        when ``inplace=True``.  Unlike a direct ``save()`` call, this does not
        crash: the result stays accessible via ``run.read_adata()`` and the
        user is told how to persist it once the store is unlocked.

        Returns:
            ``True`` if the save succeeded, ``False`` if blocked by a lock.
        """
        import warnings as _warnings
        if adata is not None:
            self._adata = adata
        try:
            self.save()
            return True
        except ZarrLockedError as exc:
            _warnings.warn(
                f"Computation complete but result NOT saved to disk: {exc}  "
                "Data is available in memory via run.read_adata(). "
                "To persist: run.unlock_zarr_parts() then run.save().",
                UserWarning,
                stacklevel=3,
            )
            return False

    def locked_zarr_parts(self) -> list[str]:
        """Return which zarr store parts are currently locked (read-only).

        Scans the run's zarr store and returns the logical parts — such as
        ``"obs"``, ``"layers/raw"`` — that contain any read-only files.
        An empty list means the store is fully writable.

        Returns:
            Sorted list of locked part paths.  ``["."]`` means the entire
            store is locked.

        Raises:
            RunNotIngestedError: If the run has not been ingested yet.
        """
        self.require_ingested()
        return get_locked_zarr_parts(self.zarr_path())

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
