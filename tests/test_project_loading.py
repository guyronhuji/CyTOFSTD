"""Tests for project loading functionality."""

import pytest
from pathlib import Path
from cytofstandard import Project
from cytofstandard.exceptions import ProjectValidationError
import pandas as pd


def test_project_load(example_project):
    """Test loading an existing project."""
    project = example_project
    
    assert project.project_id == "TEST_PROJ_001"
    assert project.project_name == "Test Project"


def test_project_load_validation_failure(temp_dir):
    """Test that loading invalid project raises error."""
    project_path = temp_dir / "invalid_project"
    project_path.mkdir()
    
    (project_path / "project.yaml").write_text("project_id: TEST\n")
    
    with pytest.raises(ProjectValidationError):
        Project.load(str(project_path), validate=True)


def test_project_list_runs_empty(simple_project):
    """Test listing runs in empty project."""
    runs = simple_project.list_runs()
    
    assert len(runs) == 0
    assert "run_id" in runs.columns


def test_project_has_run(simple_project):
    """Test checking if run exists."""
    assert not simple_project.has_run("run_001")
    
    run = simple_project.add_run(
        run_id="run_001",
        run_name="Test Run",
    )
    
    assert simple_project.has_run("run_001")


def test_project_load_with_validation_disabled(temp_dir):
    """Test loading project with validation disabled."""
    project_path = temp_dir / "test_project"
    
    project_path.mkdir(parents=True)
    (project_path / "project.yaml").write_text("""
project_id: TEST_005
project_name: Test Project
created_at: "2026-05-15T12:00:00"
package_version: "0.1.1"
storage_format: "anndata_zarr"
cell_uuid_strategy: "uuid5_project_run_filehash_eventindex"
standard_marker_file: "metadata/standard_markers.parquet"
marker_alias_file: "metadata/marker_aliases.yaml"
    """)
    (project_path / "metadata").mkdir()
    (project_path / "runs").mkdir()
    (project_path / "integrations").mkdir()
    (project_path / "project_tables").mkdir()
    (project_path / "logs").mkdir()
    (project_path / "logs" / "provenance.jsonl").touch()
    
    # Create a valid standard markers parquet file
    markers_df = pd.DataFrame({
        "standard_marker_name": ["H3", "ECad"],
        "marker_class": ["histone", "identity"],
    })
    markers_df.to_parquet(project_path / "metadata" / "standard_markers.parquet", index=False)
    
    # Create empty alias file
    (project_path / "metadata" / "marker_aliases.yaml").write_text("# Aliases\n")
    
    # Load with validation disabled
    project = Project.load(str(project_path), validate=False)
    
    assert project.project_id == "TEST_005"
    assert project.project_name == "Test Project"
