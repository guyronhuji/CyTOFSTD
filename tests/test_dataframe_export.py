"""Tests for Run.to_dataframe export helper."""

import shutil

import pandas as pd
import pytest


def _ingest_simple_run(simple_project, simple_csv_file, sample_metadata_file):
    run = simple_project.add_run(run_id="run_df")

    raw_dir = simple_project.path / "runs" / "run_df" / "raw"
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
    return run


def test_to_dataframe_mixes_markers_and_obs(
    simple_project,
    simple_csv_file,
    sample_metadata_file,
):
    """to_dataframe should export marker and obs fields in requested order."""
    run = _ingest_simple_run(simple_project, simple_csv_file, sample_metadata_file)

    df = run.to_dataframe(["sample_id", "H3", "line_id", "ECad"], layer="raw")

    assert df.columns.tolist() == ["sample_id", "H3", "line_id", "ECad"]
    assert len(df) == 3
    assert (df["sample_id"] == "S001").all()
    assert (df["line_id"] == "MCF7").all()


def test_to_dataframe_unknown_field_raises(
    simple_project,
    simple_csv_file,
    sample_metadata_file,
):
    """Unknown requested fields should raise a clear error."""
    run = _ingest_simple_run(simple_project, simple_csv_file, sample_metadata_file)

    with pytest.raises(ValueError):
        run.to_dataframe(["H3", "NOT_A_FIELD"], layer="raw")
