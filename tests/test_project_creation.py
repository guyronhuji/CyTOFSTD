"""Tests for project creation functionality."""

import pytest
from pathlib import Path
from cytofstandard import Project
from cytofstandard.exceptions import ProjectExistsError


def test_project_create(temp_dir, example_standard_markers_path, example_marker_aliases_path):
    """Test creating a new project."""
    project_path = temp_dir / "test_project"
    
    project = Project.create(
        path=str(project_path),
        project_id="TEST_001",
        project_name="Test Project",
        standard_marker_file=str(example_standard_markers_path),
        marker_alias_file=str(example_marker_aliases_path),
    )
    
    assert project.project_id == "TEST_001"
    assert project.project_name == "Test Project"
    assert project.path == project_path


def test_project_structure(temp_dir, example_standard_markers_path, example_marker_aliases_path):
    """Test that project directory structure is created correctly."""
    project_path = temp_dir / "test_project"
    
    Project.create(
        path=str(project_path),
        project_id="TEST_002",
        project_name="Test Project",
        standard_marker_file=str(example_standard_markers_path),
        marker_alias_file=str(example_marker_aliases_path),
    )
    
    # Check required files exist
    assert (project_path / "project.yaml").exists()
    assert (project_path / "metadata" / "standard_markers.parquet").exists()
    assert (project_path / "metadata" / "marker_aliases.yaml").exists()
    assert (project_path / "runs").exists()
    assert (project_path / "integrations").exists()
    assert (project_path / "project_tables" / "runs.parquet").exists()
    assert (project_path / "logs" / "provenance.jsonl").exists()


def test_project_exists_error(temp_dir, example_standard_markers_path, example_marker_aliases_path):
    """Test that creating existing project raises error."""
    project_path = temp_dir / "test_project"
    
    # Create project
    Project.create(
        path=str(project_path),
        project_id="TEST_003",
        project_name="Test Project",
        standard_marker_file=str(example_standard_markers_path),
        marker_alias_file=str(example_marker_aliases_path),
    )
    
    # Try to create again
    with pytest.raises(ProjectExistsError):
        Project.create(
            path=str(project_path),
            project_id="TEST_003",
            project_name="Test Project",
            standard_marker_file=str(example_standard_markers_path),
            marker_alias_file=str(example_marker_aliases_path),
        )


def test_project_overwrite(temp_dir, example_standard_markers_path, example_marker_aliases_path):
    """Test that overwrite=True allows recreating project."""
    project_path = temp_dir / "test_project"
    
    # Create project
    Project.create(
        path=str(project_path),
        project_id="TEST_004",
        project_name="Test Project",
        standard_marker_file=str(example_standard_markers_path),
        marker_alias_file=str(example_marker_aliases_path),
    )
    
    # Recreate with overwrite
    project = Project.create(
        path=str(project_path),
        project_id="TEST_004",
        project_name="Test Project Overwritten",
        standard_marker_file=str(example_standard_markers_path),
        marker_alias_file=str(example_marker_aliases_path),
        overwrite=True,
    )
    
    assert project.project_name == "Test Project Overwritten"
