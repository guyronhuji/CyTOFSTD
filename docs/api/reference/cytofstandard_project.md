# `cytofstandard.project`

- Source: `cytofstandard/project.py`

Project class for cytofstandard.

## Public Exports (`__all__`)

- `Project`

## Top-level Functions

No public top-level functions.

## Classes

### `Project`

CytOF project container for managing runs.

#### Methods

##### `add_run(self, run_id: str, run_name: str | None = None, panel_id: str | None = None, acquisition_date: str | None = None, instrument: str | None = None, operator: str | None = None, notes: str | None = None) -> Run`

Register a new run in the project.

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

##### `create(cls, path: str, project_id: str, project_name: str, standard_marker_file: str, marker_alias_file: str | None = None, overwrite: bool = False) -> 'Project'`

- Decorators: `classmethod`

Create a new project.

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

##### `get_marker_registry(self) -> MarkerRegistry`

Get the marker registry for this project.

Returns:
    MarkerRegistry instance

##### `get_run(self, run_id: str, validate: bool = True) -> Run`

Get a run by ID.

Args:
    run_id: Run identifier
    validate: Whether to validate run exists

Returns:
    Run instance

Raises:
    RunNotFoundError if run not found
    RunValidationError if validation fails

##### `has_run(self, run_id: str) -> bool`

Check if a run exists in the project.

Args:
    run_id: Run identifier

Returns:
    True if run exists

##### `list_runs(self) -> pd.DataFrame`

List all runs in the project.

Returns:
    DataFrame with run information

##### `load(cls, path: str, validate: bool = True) -> 'Project'`

- Decorators: `classmethod`

Load an existing project.

Args:
    path: Path to project directory
    validate: Whether to validate project structure

Returns:
    Project instance

Raises:
    ProjectValidationError if validation fails and validate=True

##### `remove_run(self, run_id: str) -> None`

Remove a run from the project and delete all run files.

Args:
    run_id: Run identifier

Raises:
    RunNotFoundError if run is not found

##### `rename_run(self, run_id: str, new_run_name: str) -> None`

Rename an existing run (run_name) and persist metadata updates.

Args:
    run_id: Existing run identifier
    new_run_name: New human-readable run name

Raises:
    RunNotFoundError: If run_id does not exist
    ValueError: If new_run_name is empty
