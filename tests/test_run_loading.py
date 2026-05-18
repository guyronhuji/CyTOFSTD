"""Tests for run loading functionality."""

import pytest
import shutil
import stat
import pandas as pd
from pathlib import Path
from cytofstandard import Project
from cytofstandard.exceptions import RunNotFoundError, RunNotIngestedError


def test_get_run(simple_project):
    """Test getting a registered run."""
    run = simple_project.add_run(run_id="run_001")

    loaded_run = simple_project.get_run("run_001")

    assert loaded_run.run_id == "run_001"
    assert loaded_run.run_config["run_id"] == "run_001"


def test_get_run_not_found(simple_project):
    """Test getting non-existent run raises error."""
    with pytest.raises(RunNotFoundError):
        simple_project.get_run("non_existent_run")


def test_get_run_validation(simple_project):
    """Test run validation in get_run."""
    run = simple_project.add_run(run_id="run_001")

    # Disable validation
    loaded = simple_project.get_run("run_001", validate=False)
    assert loaded.run_id == "run_001"


def test_read_adata_not_ingested(simple_project):
    """Test reading non-ingested run raises error."""
    run = simple_project.add_run(run_id="run_001")

    with pytest.raises(RunNotIngestedError):
        run.read_adata()


def test_is_ingested_false(simple_project):
    """Test is_ingested returns False for un-ingested run."""
    run = simple_project.add_run(run_id="run_001")

    assert not run.is_ingested()


def test_require_ingested_not_ingested(simple_project):
    """Test require_ingested raises error for non-ingested run."""
    run = simple_project.add_run(run_id="run_001")

    with pytest.raises(RunNotIngestedError):
        run.require_ingested()


def test_load_run_from_disk(simple_project):
    """Test reloading a run from disk."""
    # Add and save run
    run1 = simple_project.add_run(run_id="run_001")
    run1.run_config["status"] = "ingested"
    run1.run_config["created_at"] = "2026-05-15T12:00:00"

    # Reload project
    reloaded = Project.load(str(simple_project.path))

    # Get run
    run2 = reloaded.get_run("run_001")

    assert run2.run_id == "run_001"


def test_save_persists_external_modifications(
    simple_project,
    simple_csv_file,
    sample_metadata_file,
):
    """Run.save() persists adata modifications done outside cytofstandard."""
    run = simple_project.add_run(run_id="run_001")

    raw_dir = simple_project.path / "runs" / "run_001" / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy(simple_csv_file, raw_dir / simple_csv_file.name)

    sample_df = pd.read_csv(sample_metadata_file)
    sample_df["file_name"] = simple_csv_file.name
    sample_df.to_csv(sample_metadata_file, index=False)

    run.ingest(
        files=[str(raw_dir / simple_csv_file.name)],
        sample_metadata=str(sample_metadata_file),
        copy_raw=False,
        strict_markers=True,
    )

    adata = run.read_adata()
    adata.obs["external_flag"] = "ok"
    adata.uns["external_block"] = {"version": 1, "note": "saved"}

    run.save()

    reloaded_project = Project.load(str(simple_project.path))
    reloaded_run = reloaded_project.get_run("run_001")
    reloaded_adata = reloaded_run.read_adata()

    assert "external_flag" in reloaded_adata.obs.columns
    assert (reloaded_adata.obs["external_flag"] == "ok").all()
    assert reloaded_adata.uns["external_block"]["note"] == "saved"


def test_save_without_adata_or_zarr_raises(simple_project):
    """Run.save() should fail for un-ingested runs with no provided adata."""
    run = simple_project.add_run(run_id="run_001")

    with pytest.raises(RunNotIngestedError):
        run.save()


def test_get_run_populates_marker_variables_without_read_adata(
    simple_project,
    simple_csv_file,
    sample_metadata_file,
):
    """Project.get_run should return marker categories for ingested runs."""
    run = simple_project.add_run(run_id="run_001")

    raw_dir = simple_project.path / "runs" / "run_001" / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy(simple_csv_file, raw_dir / simple_csv_file.name)

    sample_df = pd.read_csv(sample_metadata_file)
    sample_df["file_name"] = simple_csv_file.name
    sample_df.to_csv(sample_metadata_file, index=False)

    run.ingest(
        files=[str(raw_dir / simple_csv_file.name)],
        sample_metadata=str(sample_metadata_file),
        copy_raw=False,
        strict_markers=True,
    )

    reloaded_project = Project.load(str(simple_project.path))
    loaded_run = reloaded_project.get_run("run_001")

    # Access marker variables directly, before calling read_adata().
    assert "H3" in loaded_run.markers_all
    assert "H3" in loaded_run.markers_core_histones
    assert "ECad" in loaded_run.markers_identity


def test_lock_and_unlock_selected_zarr_part(
    simple_project,
    simple_csv_file,
    sample_metadata_file,
):
    """Run can lock and unlock selected Zarr subpaths."""
    run = simple_project.add_run(run_id="run_001")

    raw_dir = simple_project.path / "runs" / "run_001" / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy(simple_csv_file, raw_dir / simple_csv_file.name)

    sample_df = pd.read_csv(sample_metadata_file)
    sample_df["file_name"] = simple_csv_file.name
    sample_df.to_csv(sample_metadata_file, index=False)

    run.ingest(
        files=[str(raw_dir / simple_csv_file.name)],
        sample_metadata=str(sample_metadata_file),
        copy_raw=False,
        strict_markers=True,
    )

    raw_layer_path = run.zarr_path() / "layers" / "raw"
    assert raw_layer_path.exists()

    obs_path = run.zarr_path() / "obs"
    obs_mode_before = obs_path.stat().st_mode

    locked_parts = run.lock_zarr_parts(parts=["layers/raw"])
    assert locked_parts == ["layers/raw"]
    assert (raw_layer_path.stat().st_mode & stat.S_IWUSR) == 0
    assert obs_path.stat().st_mode == obs_mode_before

    unlocked_parts = run.unlock_zarr_parts(parts=["layers/raw"])
    assert unlocked_parts == ["layers/raw"]
    assert (raw_layer_path.stat().st_mode & stat.S_IWUSR) != 0


def test_lock_missing_zarr_part_raises(
    simple_project,
    simple_csv_file,
    sample_metadata_file,
):
    """Locking a missing Zarr part raises FileNotFoundError in strict mode."""
    run = simple_project.add_run(run_id="run_001")

    raw_dir = simple_project.path / "runs" / "run_001" / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy(simple_csv_file, raw_dir / simple_csv_file.name)

    sample_df = pd.read_csv(sample_metadata_file)
    sample_df["file_name"] = simple_csv_file.name
    sample_df.to_csv(sample_metadata_file, index=False)

    run.ingest(
        files=[str(raw_dir / simple_csv_file.name)],
        sample_metadata=str(sample_metadata_file),
        copy_raw=False,
        strict_markers=True,
    )

    with pytest.raises(FileNotFoundError):
        run.lock_zarr_parts(parts=["layers/not_real"], strict=True)
