"""Project class for cytofstandard."""

import shutil
import yaml
import pandas as pd
from pathlib import Path
from datetime import datetime

from cytofstandard.exceptions import (
    CytofStandardError,
    ProjectExistsError,
    ProjectValidationError,
    RunExistsError,
    RunNotFoundError,
    RunValidationError,
)
from cytofstandard.storage import (
    ensure_directory,
    write_yaml,
    write_parquet,
    read_parquet,
    compute_sha256,
)
from cytofstandard.markers import MarkerRegistry
from cytofstandard.provenance import (
    ProvenanceLogger,
    log_project_created,
    log_project_loaded,
    log_run_removed,
    log_run_renamed,
)
from cytofstandard.run import Run


class Project:
    """CytOF project container for managing runs."""

    def __init__(
        self,
        path: str,
        project_id: str,
        project_name: str,
        standard_marker_file: str,
        marker_alias_file: str,
    ):
        """Initialize a project (internal use - use create() or load()).

        Args:
            path: Path to project directory
            project_id: Project identifier
            project_name: Project name
            standard_marker_file: Path to standard marker Parquet file
            marker_alias_file: Path to marker alias YAML file
        """
        self.path = Path(path)
        self.project_id = project_id
        self.project_name = project_name
        self.standard_marker_file = standard_marker_file
        self.marker_alias_file = marker_alias_file

        # Initialize paths
        self.metadata_dir = self.path / "metadata"
        self.runs_dir = self.path / "runs"
        self.integrations_dir = self.path / "integrations"
        self.project_tables_dir = self.path / "project_tables"
        self.logs_dir = self.path / "logs"

        # Load marker registry
        self._marker_registry = MarkerRegistry.from_project(self)

        # Load run registry
        self._run_registry = self._load_run_registry()

        # Initialize provenance logger
        self._provenance_logger = ProvenanceLogger(self.logs_dir / "provenance.jsonl")

        # Log project loaded event
        log_project_loaded(self._provenance_logger, self.project_id)

    @classmethod
    def create(
        cls,
        path: str,
        project_id: str,
        project_name: str,
        standard_marker_file: str,
        marker_alias_file: str | None = None,
        overwrite: bool = False,
    ) -> "Project":
        """Create a new project.

        Args:
            path: Path to create project directory
            project_id: Project identifier
            project_name: Project name
            standard_marker_file: Path to standard marker CSV/Parquet file
            marker_alias_file: Path to marker alias YAML file (optional)
            overwrite: Whether to overwrite existing project

        Returns:
            New Project instance

        Raises:
            ProjectExistsError if project already exists and overwrite=False
        """
        project_path = Path(path)

        if project_path.exists() and not overwrite:
            raise ProjectExistsError(f"Project already exists at: {path}")

        # Create directory structure
        ensure_directory(project_path)
        ensure_directory(project_path / "metadata")
        ensure_directory(project_path / "runs")
        ensure_directory(project_path / "integrations")
        ensure_directory(project_path / "project_tables")
        ensure_directory(project_path / "logs")

        # Convert standard markers to Parquet if needed
        standard_markers_path = project_path / "metadata" / "standard_markers.parquet"
        cls._convert_markers_to_parquet(standard_marker_file, standard_markers_path)

        # Copy or create marker alias file
        if marker_alias_file:
            import shutil

            alias_path = project_path / "metadata" / "marker_aliases.yaml"
            shutil.copy2(marker_alias_file, alias_path)
        else:
            # Create empty aliases file
            alias_path = project_path / "metadata" / "marker_aliases.yaml"
            with open(alias_path, "w") as f:
                f.write("# Marker aliases\n")

        # Initialize empty run registry
        runs_path = project_path / "project_tables" / "runs.parquet"
        empty_df = pd.DataFrame(
            columns=[
                "run_id",
                "run_name",
                "panel_id",
                "acquisition_date",
                "instrument",
                "operator",
                "created_at",
                "status",
                "path",
            ]
        )
        write_parquet(runs_path, empty_df)

        # Initialize provenance log
        provenance_path = project_path / "logs" / "provenance.jsonl"
        if not provenance_path.exists():
            provenance_path.touch()

        # Write project.yaml
        project_yaml_path = project_path / "project.yaml"
        project_config = {
            "project_id": project_id,
            "project_name": project_name,
            "created_at": datetime.utcnow().isoformat(),
            "package_version": "0.1.0",
            "storage_format": "anndata_zarr",
            "cell_uuid_strategy": "uuid5_project_run_filehash_eventindex",
            "standard_marker_file": "metadata/standard_markers.parquet",
            "marker_alias_file": "metadata/marker_aliases.yaml",
        }
        write_yaml(project_yaml_path, project_config)

        # Create project instance
        project = cls(
            path=str(project_path),
            project_id=project_id,
            project_name=project_name,
            standard_marker_file="metadata/standard_markers.parquet",
            marker_alias_file="metadata/marker_aliases.yaml",
        )

        # Log project creation
        log_project_created(project._provenance_logger, project_id, project_name)

        return project

    @classmethod
    def load(cls, path: str, validate: bool = True) -> "Project":
        """Load an existing project.

        Args:
            path: Path to project directory
            validate: Whether to validate project structure

        Returns:
            Project instance

        Raises:
            ProjectValidationError if validation fails and validate=True
        """
        project_path = Path(path)

        if validate:
            cls._validate_project_structure(project_path)

        # Load project.yaml
        project_yaml_path = project_path / "project.yaml"
        with open(project_yaml_path, "r") as f:
            project_config = yaml.safe_load(f)

        project_id = project_config.get("project_id", project_path.name)
        project_name = project_config.get("project_name", "Unknown")
        standard_marker_file = project_config.get(
            "standard_marker_file", "metadata/standard_markers.parquet"
        )
        marker_alias_file = project_config.get(
            "marker_alias_file", "metadata/marker_aliases.yaml"
        )

        # Create project instance
        project = cls(
            path=str(project_path),
            project_id=project_id,
            project_name=project_name,
            standard_marker_file=standard_marker_file,
            marker_alias_file=marker_alias_file,
        )

        return project

    @staticmethod
    def _convert_markers_to_parquet(input_file: str, output_path: str):
        """Convert marker file to Parquet format.

        Args:
            input_file: Input CSV or Parquet file
            output_path: Output Parquet file path
        """
        input_path = Path(input_file)

        if input_path.suffix == ".csv":
            df = pd.read_csv(input_file)
        elif input_path.suffix == ".parquet":
            df = pd.read_parquet(input_file)
        else:
            # Try CSV first
            df = pd.read_csv(input_file)

        # Ensure required columns exist
        required = ["standard_marker_name", "marker_class"]
        for col in required:
            if col not in df.columns:
                df[col] = "unknown" if col == "marker_class" else f"unknown_{col}"

        # Write to Parquet
        import os

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_parquet(output_path, index=False)

    @staticmethod
    def _validate_project_structure(path: Path):
        """Validate that project directory has required structure.

        Args:
            path: Path to project directory

        Raises:
            ProjectValidationError if validation fails
        """
        required = [
            "project.yaml",
            "metadata/standard_markers.parquet",
            "project_tables/runs.parquet",
            "runs/",
            "logs/provenance.jsonl",
        ]

        missing = []
        for item in required:
            if not (path / item).exists():
                missing.append(item)

        if missing:
            raise ProjectValidationError(
                f"Project validation failed. Missing: {missing}"
            )

    def _load_run_registry(self) -> pd.DataFrame:
        """Load the run registry from disk.

        Returns:
            DataFrame with run registry
        """
        runs_path = self.project_tables_dir / "runs.parquet"
        if runs_path.exists():
            return read_parquet(str(runs_path))
        return pd.DataFrame()

    def _save_run_registry(self):
        """Save the run registry to disk."""
        write_parquet(self.project_tables_dir / "runs.parquet", self._run_registry)

    def add_run(
        self,
        run_id: str,
        run_name: str | None = None,
        panel_id: str | None = None,
        acquisition_date: str | None = None,
        instrument: str | None = None,
        operator: str | None = None,
        notes: str | None = None,
    ) -> Run:
        """Register a new run in the project.

        Args:
            run_id: Unique run identifier
            run_name: Human-readable run name
            panel_id: Panel identifier
            acquisition_date: Date of acquisition (YYYY-MM-DD)
            instrument: Instrument name
            operator: Operator name
            notes: Additional notes

        Returns:
            New Run instance

        Raises:
            RunExistsError if run_id already exists
        """
        if self.has_run(run_id):
            raise RunExistsError(f"Run '{run_id}' already exists in project")

        # Create run directory structure
        run_path = self.runs_dir / run_id
        ensure_directory(run_path)
        ensure_directory(run_path / "raw")
        ensure_directory(run_path / "metadata")
        ensure_directory(run_path / "processed")
        ensure_directory(run_path / "logs")

        # Write run.yaml
        run_config = {
            "run_id": run_id,
            "run_name": run_name or f"Run {run_id}",
            "project_id": self.project_id,
            "panel_id": panel_id,
            "acquisition_date": acquisition_date,
            "instrument": instrument,
            "operator": operator,
            "notes": notes,
            "status": "registered",
            "created_at": datetime.utcnow().isoformat(),
        }
        write_yaml(run_path / "run.yaml", run_config)

        # Update run registry
        new_run = pd.DataFrame(
            [
                {
                    "run_id": run_id,
                    "run_name": run_config["run_name"],
                    "panel_id": panel_id,
                    "acquisition_date": acquisition_date,
                    "instrument": instrument,
                    "operator": operator,
                    "created_at": run_config["created_at"],
                    "status": "registered",
                    "path": str(run_path),
                }
            ]
        )
        self._run_registry = pd.concat([self._run_registry, new_run], ignore_index=True)
        self._save_run_registry()

        # Log event
        from cytofstandard.provenance import log_run_registered

        log_run_registered(
            self._provenance_logger,
            self.project_id,
            run_id,
            run_config["run_name"],
        )

        return Run(self, run_id, run_config)

    def get_run(self, run_id: str, validate: bool = True) -> Run:
        """Get a run by ID.

        Args:
            run_id: Run identifier
            validate: Whether to validate run exists

        Returns:
            Run instance

        Raises:
            RunNotFoundError if run not found
            RunValidationError if validation fails
        """
        if not self.has_run(run_id):
            raise RunNotFoundError(f"Run '{run_id}' not found in project")

        run_info = self._run_registry[self._run_registry["run_id"] == run_id].iloc[0]
        run_path = Path(run_info["path"])

        if validate:
            self._validate_run_structure(run_path)

        # Load run.yaml
        run_yaml_path = run_path / "run.yaml"
        with open(run_yaml_path, "r") as f:
            run_config = yaml.safe_load(f)

        from cytofstandard.provenance import log_run_loaded

        log_run_loaded(self._provenance_logger, self.project_id, run_id)

        return Run(self, run_id, run_config)

    def remove_run(self, run_id: str) -> None:
        """Remove a run from the project and delete all run files.

        Args:
            run_id: Run identifier

        Raises:
            RunNotFoundError if run is not found
        """
        if not self.has_run(run_id):
            raise RunNotFoundError(f"Run '{run_id}' not found in project")

        run_info = self._run_registry[self._run_registry["run_id"] == run_id].iloc[0]
        run_path = Path(run_info["path"])
        had_files = run_path.exists()

        if run_path.exists():
            shutil.rmtree(run_path)

        self._run_registry = (
            self._run_registry[self._run_registry["run_id"] != run_id]
            .copy()
            .reset_index(drop=True)
        )
        self._save_run_registry()

        log_run_removed(
            self._provenance_logger,
            self.project_id,
            run_id,
            str(run_info.get("run_name", "")),
            had_files,
        )

    def rename_run(self, run_id: str, new_run_name: str) -> None:
        """Rename an existing run (run_name) and persist metadata updates.

        Args:
            run_id: Existing run identifier
            new_run_name: New human-readable run name

        Raises:
            RunNotFoundError: If run_id does not exist
            ValueError: If new_run_name is empty
        """
        if not self.has_run(run_id):
            raise RunNotFoundError(f"Run '{run_id}' not found in project")

        new_name = str(new_run_name).strip()
        if not new_name:
            raise ValueError("new_run_name must be a non-empty string")

        run_idx = self._run_registry["run_id"] == run_id
        run_info = self._run_registry[run_idx].iloc[0]
        old_name = str(run_info.get("run_name", ""))
        if old_name == new_name:
            return

        run_path = Path(run_info["path"])
        run_yaml_path = run_path / "run.yaml"
        if not run_yaml_path.exists():
            raise RunValidationError(
                f"Cannot rename run '{run_id}' because run.yaml is missing at: {run_yaml_path}"
            )

        with open(run_yaml_path, "r") as f:
            run_config = yaml.safe_load(f)
        run_config["run_name"] = new_name
        write_yaml(run_yaml_path, run_config)

        self._run_registry.loc[run_idx, "run_name"] = new_name
        self._save_run_registry()

        # Project-level provenance
        log_run_renamed(
            self._provenance_logger,
            self.project_id,
            run_id,
            old_name,
            new_name,
        )

        # Run-level provenance
        run_logger = ProvenanceLogger(run_path / "logs" / "provenance.jsonl")
        run_logger.log(
            event_type="run_renamed",
            data={
                "project_id": self.project_id,
                "run_id": run_id,
                "old_run_name": old_name,
                "new_run_name": new_name,
            },
        )

    def has_run(self, run_id: str) -> bool:
        """Check if a run exists in the project.

        Args:
            run_id: Run identifier

        Returns:
            True if run exists
        """
        return run_id in self._run_registry["run_id"].values

    def list_runs(self) -> pd.DataFrame:
        """List all runs in the project.

        Returns:
            DataFrame with run information
        """
        return self._run_registry.copy()

    def get_marker_registry(self) -> MarkerRegistry:
        """Get the marker registry for this project.

        Returns:
            MarkerRegistry instance
        """
        return self._marker_registry

    @staticmethod
    def _validate_run_structure(path: Path):
        """Validate run directory structure.

        Args:
            path: Path to run directory

        Raises:
            RunValidationError if validation fails
        """
        required = [
            "run.yaml",
            "raw/",
            "metadata/",
            "processed/",
            "logs/",
        ]

        missing = []
        for item in required:
            if not (path / item).exists():
                missing.append(item)

        if missing:
            raise RunValidationError(f"Run validation failed. Missing: {missing}")


__all__ = ["Project"]
