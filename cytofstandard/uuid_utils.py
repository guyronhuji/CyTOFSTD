"""UUID generation utilities for cytofstandard."""

import uuid


def generate_cell_uuid(
    project_id: str,
    run_id: str,
    file_hash: str,
    event_index: int,
) -> str:
    """Generate deterministic UUID for a cell/event.
    
    Uses UUIDv5 with a namespace derived from project_id.
    
    Args:
        project_id: Project identifier
        run_id: Run identifier
        file_hash: SHA256 hash of the source file
        event_index: Zero-based event index in the file
        
    Returns:
        UUID string for the cell
    """
    namespace = uuid.uuid5(uuid.NAMESPACE_DNS, project_id)
    uuid_input = f"{project_id}|{run_id}|{file_hash}|{event_index}"
    return str(uuid.uuid5(namespace, uuid_input))


def generate_file_uuid(project_id: str, file_path: str, file_hash: str) -> str:
    """Generate deterministic UUID for a file.
    
    Args:
        project_id: Project identifier
        file_path: Original file path
        file_hash: File hash
        
    Returns:
        UUID string for the file
    """
    namespace = uuid.uuid5(uuid.NAMESPACE_DNS, project_id)
    uuid_input = f"{project_id}|file|{file_path}|{file_hash}"
    return str(uuid.uuid5(namespace, uuid_input))


__all__ = [
    "generate_cell_uuid",
    "generate_file_uuid",
]
