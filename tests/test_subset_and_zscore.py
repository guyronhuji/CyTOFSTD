"""Tests for subset-run creation and balanced z-scoring."""

import pandas as pd


def _ingest_uneven_two_sample_run(simple_project, temp_dir, run_id="run_001"):
    run = simple_project.add_run(run_id=run_id)

    file_a = temp_dir / f"{run_id}_sample_A.csv"
    file_a.write_text("H3,ECad\n0,10\n2,12\n")

    # Deliberately larger sample for bias test.
    file_b = temp_dir / f"{run_id}_sample_B.csv"
    file_b.write_text("H3,ECad\n10,20\n10,20\n10,20\n10,20\n")

    sample_meta = temp_dir / f"{run_id}_sample_metadata.csv"
    sample_meta.write_text(
        "file_name,sample_id,line_id\n"
        f"{file_a.name},S001,MCF7\n"
        f"{file_b.name},S002,T47D\n"
    )

    run.ingest(
        files=[str(file_a), str(file_b)],
        sample_metadata=str(sample_meta),
        copy_raw=False,
        strict_markers=True,
    )
    return run


def test_create_subset_run_by_sample_ids(simple_project, temp_dir):
    """Subset run keeps only selected sample IDs and writes a new run."""
    run = _ingest_uneven_two_sample_run(simple_project, temp_dir, run_id="run_source")

    subset = run.create_subset_run(
        new_run_id="run_subset_s001",
        sample_ids=["S001"],
    )

    subset_adata = subset.read_adata()
    assert subset_adata.n_obs == 2
    assert set(subset_adata.obs["sample_id"].astype(str).unique()) == {"S001"}
    assert "derived_from" in subset_adata.uns
    assert subset_adata.uns["derived_from"]["source_run_id"] == "run_source"


def test_create_subset_run_by_line_ids(simple_project, temp_dir):
    """Subset run can filter by line IDs."""
    run = _ingest_uneven_two_sample_run(
        simple_project, temp_dir, run_id="run_source_line"
    )

    subset = run.create_subset_run(
        new_run_id="run_subset_mcf7",
        line_ids=["MCF7"],
    )

    subset_adata = subset.read_adata()
    assert subset_adata.n_obs == 2
    assert set(subset_adata.obs["line_id"].astype(str).unique()) == {"MCF7"}


def test_zscore_markers_balanced_uses_equal_sample_contribution(
    simple_project, temp_dir
):
    """Z-score stats are computed from balanced per-sample reference cells."""
    run = _ingest_uneven_two_sample_run(simple_project, temp_dir, run_id="run_zscore")

    summary = run.zscore_markers_balanced(
        source_layer="raw",
        output_layer="z_balanced",
        groupby_col="sample_id",
        random_state=0,
    )

    adata = run.read_adata()
    assert "z_balanced" in adata.layers
    assert summary["balanced_target_per_group"] == 2
    assert summary["n_reference_cells"] == 4

    # Balanced H3 reference values are [0,2] from S001 and [10,10] from S002.
    # mean = 5.5 (different from unbalanced full-data mean 7.0).
    h3_mean = float(adata.var.loc["H3", "z_balanced_mean"])
    assert abs(h3_mean - 5.5) < 1e-6


def test_zscore_markers_balanced_single_sample_falls_back_to_all_cells(
    simple_project,
    simple_csv_file,
    sample_metadata_file,
):
    """Single-sample runs should use all cells for z-score estimation."""
    run = simple_project.add_run(run_id="run_single_sample_z")
    run.ingest(
        files=[str(simple_csv_file)],
        sample_metadata=str(sample_metadata_file),
        copy_raw=False,
        strict_markers=True,
    )

    summary = run.zscore_markers_balanced(source_layer="raw", output_layer="z_single")
    assert summary["n_reference_cells"] == run.read_adata().n_obs


def test_zscore_default_source_layer_uses_normlized_alias(simple_project, temp_dir):
    """Default source layer 'normlized' should resolve to 'normalized' when present."""
    run = _ingest_uneven_two_sample_run(
        simple_project, temp_dir, run_id="run_default_layer"
    )

    adata = run.read_adata()
    adata.layers["normalized"] = adata.layers["raw"].copy()
    run.save()

    summary = run.zscore_markers_balanced(output_layer="z_default")
    assert summary["source_layer"] == "normlized"
    assert "z_default" in run.read_adata().layers
