"""Tests for UUID generation."""

import pytest
from cytofstandard.uuid_utils import generate_cell_uuid


def test_deterministic_uuid():
    """Test that same inputs produce same UUID."""
    project_id = "TEST_PROJECT"
    run_id = "run_001"
    file_hash = "abc123def456"
    event_index = 42
    
    uuid1 = generate_cell_uuid(project_id, run_id, file_hash, event_index)
    uuid2 = generate_cell_uuid(project_id, run_id, file_hash, event_index)
    
    assert uuid1 == uuid2


def test_different_event_index_different_uuid():
    """Test that different event indices produce different UUIDs."""
    project_id = "TEST_PROJECT"
    run_id = "run_001"
    file_hash = "abc123def456"
    
    uuid1 = generate_cell_uuid(project_id, run_id, file_hash, 0)
    uuid2 = generate_cell_uuid(project_id, run_id, file_hash, 1)
    
    assert uuid1 != uuid2


def test_different_file_hash_different_uuid():
    """Test that different file hashes produce different UUIDs."""
    project_id = "TEST_PROJECT"
    run_id = "run_001"
    event_index = 42
    
    uuid1 = generate_cell_uuid(project_id, run_id, "hash1", event_index)
    uuid2 = generate_cell_uuid(project_id, run_id, "hash2", event_index)
    
    assert uuid1 != uuid2


def test_different_run_id_different_uuid():
    """Test that different run IDs produce different UUIDs."""
    project_id = "TEST_PROJECT"
    file_hash = "abc123def456"
    event_index = 42
    
    uuid1 = generate_cell_uuid(project_id, "run_001", file_hash, event_index)
    uuid2 = generate_cell_uuid(project_id, "run_002", file_hash, event_index)
    
    assert uuid1 != uuid2


def test_different_project_id_different_uuid():
    """Test that different project IDs produce different UUIDs."""
    file_hash = "abc123def456"
    event_index = 42
    
    uuid1 = generate_cell_uuid("PROJECT_1", "run_001", file_hash, event_index)
    uuid2 = generate_cell_uuid("PROJECT_2", "run_001", file_hash, event_index)
    
    assert uuid1 != uuid2


def test_uuid_format():
    """Test that UUID has correct format."""
    uuid = generate_cell_uuid("TEST", "run_001", "hash", 0)
    
    # UUID should be 36 chars with 4 hyphens
    assert len(uuid) == 36
    assert uuid.count("-") == 4


def test_uuid_v5():
    """Test that UUIDs are version 5 (SHA1)."""
    import uuid
    uuid_str = generate_cell_uuid("TEST", "run_001", "hash", 0)
    uuid_obj = uuid.UUID(uuid_str)
    
    assert uuid_obj.version == 5
