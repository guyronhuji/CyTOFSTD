"""Tests for run removal functionality."""

import shutil

import pytest

from cytofstandard.exceptions import RunNotFoundError


def test_remove_registered_run_deletes_files_and_registry(simple_project):
    """Removing a registered run deletes its directory and registry row."""
    run = simple_project.add_run(run_id="run_001", run_name="To delete")
    run_path = simple_project.path / "runs" / "run_001"
    assert run_path.exists()

    simple_project.remove_run("run_001")

    assert not run_path.exists()
    assert not simple_project.has_run("run_001")
    assert len(simple_project.list_runs()) == 0
    with pytest.raises(RunNotFoundError):
        simple_project.get_run("run_001")


def test_remove_run_not_found_raises(simple_project):
    """Removing an unknown run raises RunNotFoundError."""
    with pytest.raises(RunNotFoundError):
        simple_project.remove_run("does_not_exist")


def test_remove_ingested_run_removes_zarr(
    simple_project,
    simple_csv_file,
    sample_metadata_file,
):
    """Removing an ingested run removes processed artifacts too."""
    run = simple_project.add_run(run_id="run_001")

    raw_dir = simple_project.path / "runs" / "run_001" / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy(simple_csv_file, raw_dir / simple_csv_file.name)

    run.ingest(
        files=[str(raw_dir / simple_csv_file.name)],
        sample_metadata=str(sample_metadata_file),
        copy_raw=False,
        strict_markers=True,
    )

    zarr_path = run.zarr_path()
    assert zarr_path.exists()

    simple_project.remove_run("run_001")
    assert not zarr_path.exists()
    assert not simple_project.has_run("run_001")


def test_remove_run_handles_missing_directory(simple_project):
    """Registry cleanup still works if run directory is already missing."""
    simple_project.add_run(run_id="run_001")
    run_path = simple_project.path / "runs" / "run_001"
    shutil.rmtree(run_path)

    simple_project.remove_run("run_001")

    assert not simple_project.has_run("run_001")
    assert len(simple_project.list_runs()) == 0
