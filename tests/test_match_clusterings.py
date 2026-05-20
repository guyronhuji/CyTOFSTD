"""Tests for Run.match_clusterings."""

import numpy as np
import pandas as pd
import pytest

from cytofstandard import Project


def _build_run_with_two_clusterings(simple_project, temp_dir, run_id="run_match"):
    """Create an ingested run with two synthetic obs clustering columns."""
    run = simple_project.add_run(run_id=run_id)

    file_a = temp_dir / f"{run_id}_A.csv"
    file_a.write_text("H3,ECad\n0,10\n2,12\n4,14\n6,16\n")
    file_b = temp_dir / f"{run_id}_B.csv"
    file_b.write_text("H3,ECad\n10,20\n12,22\n14,24\n16,26\n")

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
    # 8 cells: 4 per file
    # CL:  0 0 1 1 0 0 1 1  (two clusters)
    # CL2: a a a a b b b b  (splits cleanly by file → enriched matches)
    adata.obs["CL"] = ["0", "0", "1", "1", "0", "0", "1", "1"]
    adata.obs["CL2"] = ["a", "a", "a", "a", "b", "b", "b", "b"]
    run.save(adata)
    return run


def test_returns_dict_with_expected_keys(simple_project, temp_dir):
    run = _build_run_with_two_clusterings(simple_project, temp_dir)
    report = run.match_clusterings("CL", "CL2", n_permutations=0)

    assert isinstance(report, dict)
    for key in ("A_to_B", "B_to_A", "one_to_one_matches", "contingency",
                "score_matrix", "ari", "nmi",
                "global_many_to_one_score_A_to_B",
                "global_many_to_one_score_B_to_A",
                "global_one_to_one_score"):
        assert key in report, f"Missing key: {key}"


def test_A_to_B_is_dataframe_with_correct_columns(simple_project, temp_dir):
    run = _build_run_with_two_clusterings(simple_project, temp_dir)
    report = run.match_clusterings("CL", "CL2", n_permutations=0)

    atb = report["A_to_B"]
    assert isinstance(atb, pd.DataFrame)
    for col in ("cluster_A", "matched_cluster_B", "score", "overlap",
                "size_A", "size_B", "frac_A_in_B", "frac_B_from_A",
                "hypergeom_p", "hypergeom_q"):
        assert col in atb.columns, f"A_to_B missing column: {col}"
    # One row per unique label in CL (2 clusters)
    assert len(atb) == 2


def test_B_to_A_is_dataframe_with_correct_columns(simple_project, temp_dir):
    run = _build_run_with_two_clusterings(simple_project, temp_dir)
    report = run.match_clusterings("CL", "CL2", n_permutations=0)

    bta = report["B_to_A"]
    assert isinstance(bta, pd.DataFrame)
    for col in ("cluster_B", "matched_cluster_A", "score", "overlap"):
        assert col in bta.columns, f"B_to_A missing column: {col}"
    # One row per unique label in CL2 (2 groups)
    assert len(bta) == 2


def test_one_to_one_matches_dataframe(simple_project, temp_dir):
    run = _build_run_with_two_clusterings(simple_project, temp_dir)
    report = run.match_clusterings("CL", "CL2", n_permutations=0)

    o2o = report["one_to_one_matches"]
    assert isinstance(o2o, pd.DataFrame)
    assert "cluster_A" in o2o.columns
    assert "cluster_B" in o2o.columns


def test_contingency_is_crosstab_dataframe(simple_project, temp_dir):
    run = _build_run_with_two_clusterings(simple_project, temp_dir)
    report = run.match_clusterings("CL", "CL2", n_permutations=0)

    ct = report["contingency"]
    assert isinstance(ct, pd.DataFrame)
    # 2×2: CL (0,1) × CL2 (a,b)
    assert ct.shape == (2, 2)
    # Total cells across contingency equals total cells
    assert ct.to_numpy().sum() == 8


def test_ari_and_nmi_are_floats(simple_project, temp_dir):
    run = _build_run_with_two_clusterings(simple_project, temp_dir)
    report = run.match_clusterings("CL", "CL2", n_permutations=0)

    assert isinstance(report["ari"], float)
    assert isinstance(report["nmi"], float)


def test_perfect_match_gives_high_scores(simple_project, temp_dir):
    """CL matches sample_id perfectly — expect high ARI and NMI."""
    run = _build_run_with_two_clusterings(simple_project, temp_dir)

    # Assign CL to directly mirror sample_id (S001 vs S002)
    adata = run.read_adata()
    adata.obs["CL_perfect"] = adata.obs["sample_id"].map({"S001": "0", "S002": "1"})
    run.save(adata)

    report = run.match_clusterings("CL_perfect", "sample_id", n_permutations=0)
    assert report["ari"] > 0.9
    assert report["nmi"] > 0.9


def test_invalid_key_a_raises_value_error(simple_project, temp_dir):
    run = _build_run_with_two_clusterings(simple_project, temp_dir)
    with pytest.raises(ValueError, match="nonexistent_key_A"):
        run.match_clusterings("nonexistent_key_A", "CL2", n_permutations=0)


def test_invalid_key_b_raises_value_error(simple_project, temp_dir):
    run = _build_run_with_two_clusterings(simple_project, temp_dir)
    with pytest.raises(ValueError, match="nonexistent_key_B"):
        run.match_clusterings("CL", "nonexistent_key_B", n_permutations=0)


def test_n_permutations_zero_yields_none_p_values(simple_project, temp_dir):
    run = _build_run_with_two_clusterings(simple_project, temp_dir)
    report = run.match_clusterings("CL", "CL2", n_permutations=0)

    assert report["permutation_p_A_to_B"] is None
    assert report["permutation_p_B_to_A"] is None
    assert report["permutation_p_one_to_one"] is None


def test_permutations_yield_float_p_values(simple_project, temp_dir):
    run = _build_run_with_two_clusterings(simple_project, temp_dir)
    report = run.match_clusterings("CL", "CL2", n_permutations=50, random_state=42)

    for key in ("permutation_p_A_to_B", "permutation_p_B_to_A",
                "permutation_p_one_to_one"):
        assert report[key] is not None
        assert 0.0 <= report[key] <= 1.0


def test_score_mode_overlap(simple_project, temp_dir):
    run = _build_run_with_two_clusterings(simple_project, temp_dir)
    report = run.match_clusterings("CL", "CL2", score_mode="overlap", n_permutations=0)
    # Overlap scores should equal raw counts — all values >= 0
    assert (report["score_matrix"].to_numpy() >= 0).all()


def test_score_mode_a(simple_project, temp_dir):
    run = _build_run_with_two_clusterings(simple_project, temp_dir)
    report = run.match_clusterings("CL", "CL2", score_mode="a", n_permutations=0)
    vals = report["score_matrix"].to_numpy()
    assert np.all(vals >= 0) and np.all(vals <= 1.0 + 1e-9)


def test_score_mode_b(simple_project, temp_dir):
    run = _build_run_with_two_clusterings(simple_project, temp_dir)
    report = run.match_clusterings("CL", "CL2", score_mode="b", n_permutations=0)
    vals = report["score_matrix"].to_numpy()
    assert np.all(vals >= 0) and np.all(vals <= 1.0 + 1e-9)


def test_invalid_score_mode_raises(simple_project, temp_dir):
    run = _build_run_with_two_clusterings(simple_project, temp_dir)
    with pytest.raises(ValueError, match="score_mode"):
        run.match_clusterings("CL", "CL2", score_mode="bad_mode", n_permutations=0)


def test_cluster_to_sample_matching(simple_project, temp_dir):
    """Matching a clustering against sample_id works (mix of obs types)."""
    run = _build_run_with_two_clusterings(simple_project, temp_dir)
    report = run.match_clusterings("CL", "sample_id", n_permutations=0)

    assert isinstance(report, dict)
    # sample_id has 2 values (S001, S002) → 2×2 contingency
    assert report["contingency"].shape == (2, 2)


def test_global_scores_are_non_negative(simple_project, temp_dir):
    run = _build_run_with_two_clusterings(simple_project, temp_dir)
    report = run.match_clusterings("CL", "CL2", n_permutations=0)

    assert report["global_many_to_one_score_A_to_B"] >= 0
    assert report["global_many_to_one_score_B_to_A"] >= 0
    assert report["global_one_to_one_score"] >= 0


def test_b_receives_from_a_and_a_receives_from_b_present(simple_project, temp_dir):
    run = _build_run_with_two_clusterings(simple_project, temp_dir)
    report = run.match_clusterings("CL", "CL2", n_permutations=0)

    assert "B_receives_from_A" in report
    assert "A_receives_from_B" in report
    assert isinstance(report["B_receives_from_A"], pd.DataFrame)
    assert isinstance(report["A_receives_from_B"], pd.DataFrame)


# --- plot= parameter tests ---

def test_plot_false_returns_dict(simple_project, temp_dir):
    run = _build_run_with_two_clusterings(simple_project, temp_dir)
    result = run.match_clusterings("CL", "CL2", n_permutations=0, plot=False)
    assert isinstance(result, dict)


def test_plot_heatmap_returns_report_and_figure(simple_project, temp_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    run = _build_run_with_two_clusterings(simple_project, temp_dir)
    result = run.match_clusterings("CL", "CL2", n_permutations=0, plot="heatmap")
    assert isinstance(result, tuple) and len(result) == 2
    report, fig = result
    assert isinstance(report, dict)
    assert hasattr(fig, "savefig")
    plt.close("all")


def test_plot_heatmap_matrix_overlap(simple_project, temp_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    run = _build_run_with_two_clusterings(simple_project, temp_dir)
    report, fig = run.match_clusterings(
        "CL", "CL2", n_permutations=0, plot="heatmap", matrix="overlap"
    )
    assert isinstance(report, dict)
    assert hasattr(fig, "savefig")
    plt.close("all")


def test_plot_clustermap_returns_report_and_clustergrid(simple_project, temp_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    run = _build_run_with_two_clusterings(simple_project, temp_dir)
    result = run.match_clusterings("CL", "CL2", n_permutations=0, plot="clustermap")
    assert isinstance(result, tuple) and len(result) == 2
    report, grid = result
    assert isinstance(report, dict)
    assert isinstance(grid, sns.matrix.ClusterGrid)
    plt.close("all")


def test_plot_invalid_value_raises(simple_project, temp_dir):
    run = _build_run_with_two_clusterings(simple_project, temp_dir)
    with pytest.raises(ValueError, match="plot"):
        run.match_clusterings("CL", "CL2", n_permutations=0, plot="bad")


def test_matrix_invalid_value_raises(simple_project, temp_dir):
    run = _build_run_with_two_clusterings(simple_project, temp_dir)
    with pytest.raises(ValueError, match="matrix"):
        run.match_clusterings("CL", "CL2", n_permutations=0, plot="heatmap", matrix="bad")
