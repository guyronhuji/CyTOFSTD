"""Tests for Run.annotate_clusters."""

import warnings

import numpy as np
import pytest

from cytofstandard import Project


def _build_run_with_clusters(simple_project, temp_dir, run_id="run_annot"):
    """Create an ingested run with a synthetic cluster obs column."""
    run = simple_project.add_run(run_id=run_id)

    file_a = temp_dir / f"{run_id}_A.csv"
    file_a.write_text("H3,ECad\n0,10\n2,12\n")
    file_b = temp_dir / f"{run_id}_B.csv"
    file_b.write_text("H3,ECad\n10,20\n12,22\n")

    meta = temp_dir / f"{run_id}_meta.csv"
    meta.write_text(
        "file_name,sample_id,line_id\n"
        f"{file_a.name},S001,MCF7\n"
        f"{file_b.name},S002,T47D\n"
    )
    run.ingest(
        files=[str(file_a), str(file_b)],
        sample_metadata=str(meta),
        copy_raw=False,
        strict_markers=True,
    )

    adata = run.read_adata()
    # Assign alternating cluster labels: cells 0,2 → cluster "0"; cells 1,3 → cluster "1"
    adata.obs["test_leiden"] = ["0", "1", "0", "1"]
    run.save(adata)
    return run


def test_annotate_creates_output_column(simple_project, temp_dir):
    run = _build_run_with_clusters(simple_project, temp_dir)
    run.annotate_clusters("test_leiden", {"0": "T cell", "1": "B cell"})

    adata = run.read_adata()
    assert "test_leiden_annotated" in adata.obs.columns
    assert set(adata.obs["test_leiden_annotated"]) == {"T cell", "B cell"}


def test_annotate_correct_values(simple_project, temp_dir):
    run = _build_run_with_clusters(simple_project, temp_dir)
    annotated = run.annotate_clusters("test_leiden", {"0": "T cell", "1": "B cell"})

    expected = ["T cell", "B cell", "T cell", "B cell"]
    assert list(annotated) == expected


def test_annotate_custom_output_key(simple_project, temp_dir):
    run = _build_run_with_clusters(simple_project, temp_dir)
    run.annotate_clusters(
        "test_leiden", {"0": "T cell"}, output_key="cell_type_manual"
    )
    adata = run.read_adata()
    assert "cell_type_manual" in adata.obs.columns


def test_annotate_partial_mapping_passthrough(simple_project, temp_dir):
    """Unmapped clusters should keep their original label."""
    run = _build_run_with_clusters(simple_project, temp_dir)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        annotated = run.annotate_clusters("test_leiden", {"0": "T cell"})

    # cluster "1" is unmapped — should stay as "1"
    assert "T cell" in list(annotated)
    assert "1" in list(annotated)


def test_annotate_unknown_cluster_key_raises(simple_project, temp_dir):
    run = _build_run_with_clusters(simple_project, temp_dir)
    with pytest.raises(ValueError, match="cluster_key"):
        run.annotate_clusters("nonexistent_col", {"0": "T cell"})


def test_annotate_unknown_map_key_warns(simple_project, temp_dir):
    run = _build_run_with_clusters(simple_project, temp_dir)
    with pytest.warns(UserWarning, match="not found"):
        run.annotate_clusters("test_leiden", {"99": "Ghost cell"})


def test_annotate_persists_to_zarr(simple_project, temp_dir):
    run = _build_run_with_clusters(simple_project, temp_dir)
    run.annotate_clusters("test_leiden", {"0": "T cell", "1": "B cell"}, inplace=True)

    # Re-read from disk to verify persistence
    fresh = run.read_adata()
    assert "test_leiden_annotated" in fresh.obs.columns


def test_annotate_inplace_false_does_not_persist(simple_project, temp_dir):
    run = _build_run_with_clusters(simple_project, temp_dir)
    run.annotate_clusters(
        "test_leiden", {"0": "T cell", "1": "B cell"}, inplace=False
    )

    # Force fresh read from zarr
    run._adata = None
    fresh = run.read_adata()
    assert "test_leiden_annotated" not in fresh.obs.columns


def test_annotate_return_is_series(simple_project, temp_dir):
    import pandas as pd

    run = _build_run_with_clusters(simple_project, temp_dir)
    result = run.annotate_clusters("test_leiden", {"0": "T cell", "1": "B cell"})
    assert isinstance(result, pd.Series)
    assert len(result) == 4
