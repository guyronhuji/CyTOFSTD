"""Tests for Run.ingestion_summary."""

import pandas as pd
import pytest

from cytofstandard import Project
from cytofstandard.exceptions import RunNotIngestedError


def _build_two_sample_run(simple_project, temp_dir, run_id="run_summary"):
    run = simple_project.add_run(run_id=run_id)

    file_a = temp_dir / f"{run_id}_A.csv"
    file_a.write_text("H3,ECad\n0,10\n2,12\n4,14\n")
    file_b = temp_dir / f"{run_id}_B.csv"
    file_b.write_text("H3,ECad\n10,20\n12,22\n14,24\n")

    meta = temp_dir / f"{run_id}_meta.csv"
    meta.write_text(
        "file_name,sample_id,line_id,condition\n"
        f"{file_a.name},S001,MCF7,control\n"
        f"{file_b.name},S002,T47D,treated\n"
    )
    run.ingest(
        files=[str(file_a), str(file_b)],
        sample_metadata=str(meta),
        copy_raw=False,
        strict_markers=True,
    )
    return run


def test_returns_dataframe(simple_project, temp_dir):
    run = _build_two_sample_run(simple_project, temp_dir)
    result = run.ingestion_summary()
    assert isinstance(result, pd.DataFrame)


def test_one_row_per_file(simple_project, temp_dir):
    run = _build_two_sample_run(simple_project, temp_dir)
    result = run.ingestion_summary()
    assert len(result) == 2


def test_required_columns_present(simple_project, temp_dir):
    run = _build_two_sample_run(simple_project, temp_dir)
    result = run.ingestion_summary()
    for col in ("file_name", "sample_id", "line_id", "n_events_file",
                "n_cells_obs", "file_type", "ingested_at"):
        assert col in result.columns, f"Missing column: {col}"


def test_sample_ids_correct(simple_project, temp_dir):
    run = _build_two_sample_run(simple_project, temp_dir)
    result = run.ingestion_summary()
    assert set(result["sample_id"]) == {"S001", "S002"}


def test_n_events_file_correct(simple_project, temp_dir):
    run = _build_two_sample_run(simple_project, temp_dir)
    result = run.ingestion_summary()
    assert result["n_events_file"].sum() == 6


def test_n_cells_obs_equals_n_events_before_qc(simple_project, temp_dir):
    run = _build_two_sample_run(simple_project, temp_dir)
    result = run.ingestion_summary()
    assert result["n_cells_obs"].sum() == result["n_events_file"].sum()


def test_extra_metadata_column_included(simple_project, temp_dir):
    """Sample-level metadata columns (e.g. condition) should appear."""
    run = _build_two_sample_run(simple_project, temp_dir)
    result = run.ingestion_summary()
    assert "condition" in result.columns
    assert set(result["condition"]) == {"control", "treated"}


def test_not_ingested_raises(simple_project):
    run = simple_project.add_run(run_id="run_nodata")
    with pytest.raises(Exception):
        run.ingestion_summary()


def test_file_type_is_csv(simple_project, temp_dir):
    run = _build_two_sample_run(simple_project, temp_dir)
    result = run.ingestion_summary()
    assert (result["file_type"] == "csv").all()
