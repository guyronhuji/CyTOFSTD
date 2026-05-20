"""Tests for Run.differential_abundance."""

import numpy as np
import pandas as pd
import pytest

from cytofstandard import Project


def _build_run_with_abundance_data(simple_project, temp_dir, run_id="run_da"):
    """Ingested run with engineered cluster proportions across three groups.

    Group G1: 8 cells in cluster A, 2 cells in cluster B  (skewed toward A)
    Group G2: 2 cells in cluster A, 8 cells in cluster B  (skewed toward B)
    Group G3: 5 cells in cluster A, 5 cells in cluster B  (balanced)
    """
    run = simple_project.add_run(run_id=run_id)

    # Build CSV files — one per 'sample'; the actual marker values don't matter
    # for DA tests, but we need a valid ingestion.
    rows_g1a = "\n".join(f"0,10" for _ in range(8))
    rows_g1b = "\n".join(f"10,0" for _ in range(2))
    rows_g2a = "\n".join(f"0,10" for _ in range(2))
    rows_g2b = "\n".join(f"10,0" for _ in range(8))
    rows_g3 = "\n".join(f"5,5" for _ in range(10))

    file_g1 = temp_dir / f"{run_id}_G1.csv"
    file_g1.write_text("H3,ECad\n" + rows_g1a + "\n" + rows_g1b)
    file_g2 = temp_dir / f"{run_id}_G2.csv"
    file_g2.write_text("H3,ECad\n" + rows_g2a + "\n" + rows_g2b)
    file_g3 = temp_dir / f"{run_id}_G3.csv"
    file_g3.write_text("H3,ECad\n" + rows_g3)

    meta = temp_dir / f"{run_id}_meta.csv"
    meta.write_text(
        "file_name,sample_id,line_id\n"
        f"{file_g1.name},S_G1,G1\n"
        f"{file_g2.name},S_G2,G2\n"
        f"{file_g3.name},S_G3,G3\n"
    )
    run.ingest(
        files=[str(file_g1), str(file_g2), str(file_g3)],
        sample_metadata=str(meta),
        copy_raw=False,
        strict_markers=True,
    )

    adata = run.read_adata()
    # Assign cluster labels to match the engineered proportions above.
    # G1 file: first 8 rows → cluster "A", last 2 → cluster "B"
    # G2 file: first 2 rows → cluster "A", last 8 → cluster "B"
    # G3 file: 10 rows → alternating A/B (5 each)
    clusters = (
        ["A"] * 8 + ["B"] * 2      # G1
        + ["A"] * 2 + ["B"] * 8    # G2
        + ["A", "B"] * 5           # G3
    )
    adata.obs["test_leiden"] = clusters
    run.save(adata)
    return run


def test_da_returns_dataframe(simple_project, temp_dir):
    run = _build_run_with_abundance_data(simple_project, temp_dir)
    result = run.differential_abundance("test_leiden", groupby="line_id")
    assert isinstance(result, pd.DataFrame)


def test_da_columns(simple_project, temp_dir):
    run = _build_run_with_abundance_data(simple_project, temp_dir)
    result = run.differential_abundance("test_leiden", groupby="line_id")
    expected = {
        "cluster", "group_a", "group_b",
        "n_a", "n_b", "N_a", "N_b",
        "freq_a", "freq_b", "p_value", "p_adj",
    }
    assert expected.issubset(set(result.columns))


def test_da_all_comparisons_row_count(simple_project, temp_dir):
    """3 groups × 2 clusters = 2 clusters per pair; 3 pairs → 6 rows."""
    run = _build_run_with_abundance_data(simple_project, temp_dir)
    result = run.differential_abundance(
        "test_leiden", groupby="line_id", comparisons="all"
    )
    # 3 pairs (G1-G2, G1-G3, G2-G3) × 2 clusters = 6 rows
    assert len(result) == 6


def test_da_explicit_comparisons(simple_project, temp_dir):
    run = _build_run_with_abundance_data(simple_project, temp_dir)
    result = run.differential_abundance(
        "test_leiden",
        groupby="line_id",
        comparisons=[("G1", "G2")],
    )
    # 1 pair × 2 clusters = 2 rows
    assert len(result) == 2
    assert set(result["group_a"]) == {"G1"}
    assert set(result["group_b"]) == {"G2"}


def test_da_skewed_cluster_is_significant(simple_project, temp_dir):
    """G1 vs G2 should show a significant difference for both clusters."""
    run = _build_run_with_abundance_data(simple_project, temp_dir)
    result = run.differential_abundance(
        "test_leiden",
        groupby="line_id",
        comparisons=[("G1", "G2")],
        method="fisher",
        multitest="bh",
    )
    # Both clusters A and B are skewed between G1 and G2
    cluster_a = result[result["cluster"] == "A"].iloc[0]
    assert cluster_a["p_adj"] < 0.05


def test_da_method_chi2(simple_project, temp_dir):
    run = _build_run_with_abundance_data(simple_project, temp_dir)
    result = run.differential_abundance(
        "test_leiden", groupby="line_id", method="chi2"
    )
    assert result["p_value"].notna().all()
    assert (result["p_value"] >= 0).all()
    assert (result["p_value"] <= 1).all()


def test_da_no_multitest(simple_project, temp_dir):
    run = _build_run_with_abundance_data(simple_project, temp_dir)
    result = run.differential_abundance(
        "test_leiden", groupby="line_id", multitest=None
    )
    assert result["p_adj"].isna().all()


def test_da_invalid_cluster_key(simple_project, temp_dir):
    run = _build_run_with_abundance_data(simple_project, temp_dir)
    with pytest.raises(ValueError, match="cluster_key"):
        run.differential_abundance("nonexistent", groupby="line_id")


def test_da_invalid_groupby(simple_project, temp_dir):
    run = _build_run_with_abundance_data(simple_project, temp_dir)
    with pytest.raises(ValueError, match="groupby"):
        run.differential_abundance("test_leiden", groupby="nonexistent_col")


def test_da_invalid_method(simple_project, temp_dir):
    run = _build_run_with_abundance_data(simple_project, temp_dir)
    with pytest.raises(ValueError, match="method"):
        run.differential_abundance("test_leiden", groupby="line_id", method="wald")


def test_da_plot_returns_tuple(simple_project, temp_dir):
    import matplotlib
    matplotlib.use("Agg")
    run = _build_run_with_abundance_data(simple_project, temp_dir)
    result = run.differential_abundance(
        "test_leiden", groupby="line_id", plot=True
    )
    assert isinstance(result, tuple)
    df, (fig, ax) = result
    assert isinstance(df, pd.DataFrame)
    import matplotlib.pyplot as plt
    assert isinstance(fig, plt.Figure)
    plt.close("all")


def test_da_freq_values_sum_to_one_per_group(simple_project, temp_dir):
    """freq_a and freq_b should be proportions (each pair sums to 1.0 within a group)."""
    run = _build_run_with_abundance_data(simple_project, temp_dir)
    result = run.differential_abundance(
        "test_leiden",
        groupby="line_id",
        comparisons=[("G1", "G2")],
    )
    freq_a_sum = result.groupby(["group_a", "group_b"])["freq_a"].sum().values[0]
    freq_b_sum = result.groupby(["group_a", "group_b"])["freq_b"].sum().values[0]
    assert abs(freq_a_sum - 1.0) < 1e-9
    assert abs(freq_b_sum - 1.0) < 1e-9
