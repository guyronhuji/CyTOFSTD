"""Tests for run rename functionality."""

import json

import pytest

from cytofstandard.exceptions import RunNotFoundError


def _read_jsonl(path):
    lines = [line for line in path.read_text().splitlines() if line.strip()]
    return [json.loads(line) for line in lines]


def test_project_rename_run_updates_metadata(simple_project):
    """Project.rename_run should update run.yaml and run registry."""
    simple_project.add_run(run_id="run_001", run_name="Old Name")

    simple_project.rename_run("run_001", "New Name")

    loaded = simple_project.get_run("run_001")
    assert loaded.run_config["run_name"] == "New Name"

    runs_df = simple_project.list_runs()
    row = runs_df[runs_df["run_id"] == "run_001"].iloc[0]
    assert row["run_name"] == "New Name"


def test_run_rename_convenience_method(simple_project):
    """Run.rename should call project rename and update this object."""
    run = simple_project.add_run(run_id="run_001", run_name="Old Name")

    run.rename("Renamed From Run")

    assert run.run_config["run_name"] == "Renamed From Run"
    loaded = simple_project.get_run("run_001")
    assert loaded.run_config["run_name"] == "Renamed From Run"


def test_rename_run_not_found_raises(simple_project):
    """Renaming unknown run should raise RunNotFoundError."""
    with pytest.raises(RunNotFoundError):
        simple_project.rename_run("missing_run", "Anything")


def test_rename_run_writes_provenance(simple_project):
    """Renaming a run should be recorded in both project and run provenance."""
    simple_project.add_run(run_id="run_001", run_name="Old Name")
    simple_project.rename_run("run_001", "New Name")

    project_log = simple_project.path / "logs" / "provenance.jsonl"
    run_log = simple_project.path / "runs" / "run_001" / "logs" / "provenance.jsonl"

    project_events = [e["event_type"] for e in _read_jsonl(project_log)]
    run_events = [e["event_type"] for e in _read_jsonl(run_log)]

    assert "run_renamed" in project_events
    assert "run_renamed" in run_events
