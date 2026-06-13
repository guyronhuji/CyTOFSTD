# `cytofstandard.uuid_utils`

- Source: `cytofstandard/uuid_utils.py`

UUID generation utilities for cytofstandard.

## Public Exports (`__all__`)

- `generate_cell_uuid`
- `generate_file_uuid`

## Top-level Functions

### `generate_cell_uuid(project_id: str, run_id: str, file_hash: str, event_index: int) -> str`

Generate deterministic UUID for a cell/event.

Uses UUIDv5 with a namespace derived from project_id.

Args:
    project_id: Project identifier
    run_id: Run identifier
    file_hash: SHA256 hash of the source file
    event_index: Zero-based event index in the file
    
Returns:
    UUID string for the cell

### `generate_file_uuid(project_id: str, file_path: str, file_hash: str) -> str`

Generate deterministic UUID for a file.

Args:
    project_id: Project identifier
    file_path: Original file path
    file_hash: File hash
    
Returns:
    UUID string for the file

## Classes

No public classes.
