# `cytofstandard.provenance`

- Source: `cytofstandard/provenance.py`

Provenance logging utilities for cytofstandard.

## Top-level Functions

### `get_timestamp() -> str`

Return current timestamp in ISO format.

### `log_ingestion_failed(logger: ProvenanceLogger, project_id: str, run_id: str, error: str)`

Log ingestion failure event.

### `log_ingestion_started(logger: ProvenanceLogger, project_id: str, run_id: str, n_files: int)`

Log ingestion start event.

### `log_project_created(logger: ProvenanceLogger, project_id: str, project_name: str)`

Log project creation event.

### `log_project_loaded(logger: ProvenanceLogger, project_id: str)`

Log project load event.

### `log_run_ingested(logger: ProvenanceLogger, project_id: str, run_id: str, n_cells: int, n_markers: int, n_files: int)`

Log successful ingestion event.

### `log_run_loaded(logger: ProvenanceLogger, project_id: str, run_id: str)`

Log run load event.

### `log_run_registered(logger: ProvenanceLogger, project_id: str, run_id: str, run_name: str)`

Log run registration event.

### `log_run_removed(logger: ProvenanceLogger, project_id: str, run_id: str, run_name: str | None, had_files: bool)`

Log run removal event.

### `log_run_renamed(logger: ProvenanceLogger, project_id: str, run_id: str, old_run_name: str | None, new_run_name: str)`

Log run rename event.

## Classes

### `ProvenanceLogger`

Logger for provenance events.

#### Methods

##### `log(self, event_type: str, data: dict)`

Log an event.

Args:
    event_type: Type of event (e.g., 'project_created', 'run_ingested')
    data: Event data dictionary
