"""Tests for run registration functionality."""

import pytest
from pathlib import Path
from cytofstandard import Project
from cytofstandard.exceptions import RunExistsError


def test_add_run(simple_project):
    """Test adding a run to a project."""
    run = simple_project.add_run(
        run_id="run_001",
        run_name="Test Run",
        panel_id="test_panel",
        acquisition_date="2026-05-15",
        instrument="Helios",
        operator="GR",
        notes="Test run",
    )
    
    assert run.run_id == "run_001"
    assert run.run_config["run_name"] == "Test Run"
    assert run.run_config["status"] == "registered"


def test_run_directory_structure(simple_project):
    """Test that run directory is created with correct structure."""
    run = simple_project.add_run(run_id="run_001")
    
    run_path = simple_project.path / "runs" / "run_001"
    
    assert run_path.exists()
    assert (run_path / "run.yaml").exists()
    assert (run_path / "raw").exists()
    assert (run_path / "metadata").exists()
    assert (run_path / "processed").exists()
    assert (run_path / "logs").exists()


def test_run_yaml_created(simple_project):
    """Test that run.yaml is created with correct content."""
    run = simple_project.add_run(
        run_id="run_001",
        run_name="Test Run",
    )
    
    run_yaml = simple_project.path / "runs" / "run_001" / "run.yaml"
    
    with open(run_yaml, "r") as f:
        import yaml
        config = yaml.safe_load(f)
    
    assert config["run_id"] == "run_001"
    assert config["run_name"] == "Test Run"
    assert config["status"] == "registered"


def test_run_in_list_runs(simple_project):
    """Test that added run appears in list_runs()."""
    run = simple_project.add_run(run_id="run_001")
    
    runs = simple_project.list_runs()
    
    assert len(runs) == 1
    assert runs.iloc[0]["run_id"] == "run_001"
    assert runs.iloc[0]["status"] == "registered"


def test_add_multiple_runs(simple_project):
    """Test adding multiple runs."""
    run1 = simple_project.add_run(run_id="run_001")
    run2 = simple_project.add_run(run_id="run_002")
    run3 = simple_project.add_run(run_id="run_003")
    
    runs = simple_project.list_runs()
    
    assert len(runs) == 3
    run_ids = set(runs["run_id"].tolist())
    assert run_ids == {"run_001", "run_002", "run_003"}


def test_run_exists_error(simple_project):
    """Test that adding duplicate run raises error."""
    simple_project.add_run(run_id="run_001")
    
    with pytest.raises(RunExistsError):
        simple_project.add_run(run_id="run_001")


def test_add_run_after_loading(simple_project):
    """Test adding run after reloading project."""
    project_path = simple_project.path
    
    # Reload project
    reloaded = Project.load(str(project_path))
    
    # Add run
    run = reloaded.add_run(run_id="run_001")
    
    assert run.run_id == "run_001"
    
    # Verify run appears in reloaded project
    runs = reloaded.list_runs()
    assert len(runs) == 1


def test_run_status_values(simple_project):
    """Test that status values are correct."""
    run = simple_project.add_run(run_id="run_001")
    
    assert run.status == "registered"
    assert run.status in ["registered", "ingested", "failed_ingestion"]
