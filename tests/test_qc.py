"""Tests for QC plotting and gating."""

import shutil

import matplotlib
import pytest


matplotlib.use("Agg")


def _ingest_simple_csv(run, project, simple_csv_file, sample_metadata_file):
    raw_dir = project.path / "runs" / run.run_id / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy(simple_csv_file, raw_dir / simple_csv_file.name)

    run.ingest(
        files=[str(raw_dir / simple_csv_file.name)],
        sample_metadata=str(sample_metadata_file),
        copy_raw=False,
        strict_markers=True,
    )


def test_plot_marker_histograms(simple_project, simple_csv_file, sample_metadata_file):
    """Histogram plotting creates one axis per marker and sample."""
    run = simple_project.add_run(run_id="run_qc_plot")
    _ingest_simple_csv(run, simple_project, simple_csv_file, sample_metadata_file)

    fig, axes = run.plot_marker_histograms(markers=["H3", "ECad"], layer="raw")

    assert axes.shape == (1, 2)
    fig.canvas.draw()


def test_qc_gate_numeric_bounds(simple_project, simple_csv_file, sample_metadata_file):
    """Numeric gates remove failing cells and persist filtered data."""
    run = simple_project.add_run(run_id="run_qc_numeric")
    _ingest_simple_csv(run, simple_project, simple_csv_file, sample_metadata_file)

    pass_mask = run.qc_gate(
        {
            "H3": {
                "lower": 110,
                "upper": 140,
            }
        },
        layer="raw",
    )

    adata = run.read_adata()
    assert int(pass_mask.sum()) == 1
    assert adata.n_obs == 1
    assert "qc" in adata.uns
    assert adata.uns["qc"]["latest"]["n_removed"] == 2


def test_qc_gate_percentile_bounds(
    simple_project, simple_csv_file, sample_metadata_file
):
    """Percentile gates are resolved per marker and applied."""
    run = simple_project.add_run(run_id="run_qc_percentile")
    _ingest_simple_csv(run, simple_project, simple_csv_file, sample_metadata_file)

    run.qc_gate(
        {
            "H3": {
                "lower": "p50",
                "upper": None,
            }
        },
        layer="raw",
    )

    adata = run.read_adata()
    assert adata.n_obs == 2


def test_qc_gate_inplace_false(simple_project, simple_csv_file, sample_metadata_file):
    """When inplace is False, gate mask is returned without mutating run zarr."""
    run = simple_project.add_run(run_id="run_qc_no_inplace")
    _ingest_simple_csv(run, simple_project, simple_csv_file, sample_metadata_file)

    pass_mask = run.qc_gate(
        {
            "H3": (110, 140),
        },
        layer="raw",
        inplace=False,
    )

    adata = run.read_adata()
    assert int(pass_mask.sum()) == 1
    assert adata.n_obs == 3


def test_qc_gate_unknown_marker_raises(
    simple_project, simple_csv_file, sample_metadata_file
):
    """Unknown marker in gates raises a clear error."""
    run = simple_project.add_run(run_id="run_qc_unknown")
    _ingest_simple_csv(run, simple_project, simple_csv_file, sample_metadata_file)

    with pytest.raises(ValueError):
        run.qc_gate({"NOT_A_MARKER": {"lower": 1, "upper": 2}}, layer="raw")
