"""Tests for Run.plot_cluster_composition."""

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd
import pytest

from cytofstandard import Project


def _build_run_with_clusters(simple_project, temp_dir, run_id="run_comp"):
    run = simple_project.add_run(run_id=run_id)

    file_a = temp_dir / f"{run_id}_A.csv"
    file_a.write_text("H3,ECad\n0,10\n2,12\n4,14\n")
    file_b = temp_dir / f"{run_id}_B.csv"
    file_b.write_text("H3,ECad\n10,20\n12,22\n14,24\n")

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
    # 6 cells: cluster 0 = cells 0,2,4; cluster 1 = cells 1,3,5
    adata.obs["CL"] = ["0", "1", "0", "1", "0", "1"]
    run.save(adata)
    return run


def test_returns_fig_ax_tuple(simple_project, temp_dir):
    import matplotlib.pyplot as plt
    run = _build_run_with_clusters(simple_project, temp_dir)
    result = run.plot_cluster_composition("CL", "sample_id")
    assert isinstance(result, tuple) and len(result) == 2
    fig, ax = result
    assert hasattr(fig, "savefig")
    assert hasattr(ax, "get_xlabel")
    plt.close("all")


def test_normalize_cluster_bars_equal_n_clusters(simple_project, temp_dir):
    import matplotlib.pyplot as plt
    run = _build_run_with_clusters(simple_project, temp_dir)
    fig, ax = run.plot_cluster_composition("CL", "sample_id", normalize="cluster")
    # x-axis should have 2 ticks (cluster 0 and cluster 1)
    assert len(ax.get_xticklabels()) == 2
    plt.close("all")


def test_normalize_group_bars_equal_n_groups(simple_project, temp_dir):
    import matplotlib.pyplot as plt
    run = _build_run_with_clusters(simple_project, temp_dir)
    fig, ax = run.plot_cluster_composition("CL", "sample_id", normalize="group")
    # x-axis should have 2 ticks (S001 and S002)
    assert len(ax.get_xticklabels()) == 2
    plt.close("all")


def test_ylim_is_0_to_1(simple_project, temp_dir):
    import matplotlib.pyplot as plt
    run = _build_run_with_clusters(simple_project, temp_dir)
    fig, ax = run.plot_cluster_composition("CL", "sample_id")
    ymin, ymax = ax.get_ylim()
    assert ymin == 0.0
    assert ymax == pytest.approx(1.0, abs=0.05)
    plt.close("all")


def test_invalid_normalize_raises(simple_project, temp_dir):
    run = _build_run_with_clusters(simple_project, temp_dir)
    with pytest.raises(ValueError, match="normalize"):
        run.plot_cluster_composition("CL", "sample_id", normalize="bad")


def test_invalid_cluster_key_raises(simple_project, temp_dir):
    run = _build_run_with_clusters(simple_project, temp_dir)
    with pytest.raises(ValueError, match="nonexistent_CL"):
        run.plot_cluster_composition("nonexistent_CL", "sample_id")


def test_invalid_groupby_raises(simple_project, temp_dir):
    run = _build_run_with_clusters(simple_project, temp_dir)
    with pytest.raises(ValueError, match="nonexistent_group"):
        run.plot_cluster_composition("CL", "nonexistent_group")


def test_stacked_bars_sum_to_one_cluster_mode(simple_project, temp_dir):
    """Each cluster bar must sum to 1.0 (fully stacked)."""
    import matplotlib.pyplot as plt
    run = _build_run_with_clusters(simple_project, temp_dir)
    fig, ax = run.plot_cluster_composition("CL", "sample_id", normalize="cluster")

    # Sum heights of all bar patches per x position
    bar_sums = {}
    for patch in ax.patches:
        x = round(patch.get_x() + patch.get_width() / 2, 2)
        bar_sums[x] = bar_sums.get(x, 0.0) + patch.get_height()

    for x, total in bar_sums.items():
        assert total == pytest.approx(1.0, abs=1e-6), f"Bar at x={x} sums to {total}"
    plt.close("all")


def test_stacked_bars_sum_to_one_group_mode(simple_project, temp_dir):
    """Each group bar must sum to 1.0 (fully stacked)."""
    import matplotlib.pyplot as plt
    run = _build_run_with_clusters(simple_project, temp_dir)
    fig, ax = run.plot_cluster_composition("CL", "sample_id", normalize="group")

    bar_sums = {}
    for patch in ax.patches:
        x = round(patch.get_x() + patch.get_width() / 2, 2)
        bar_sums[x] = bar_sums.get(x, 0.0) + patch.get_height()

    for x, total in bar_sums.items():
        assert total == pytest.approx(1.0, abs=1e-6), f"Bar at x={x} sums to {total}"
    plt.close("all")


def test_accepts_existing_ax(simple_project, temp_dir):
    """When ax= is provided, figure is not created internally."""
    import matplotlib.pyplot as plt
    run = _build_run_with_clusters(simple_project, temp_dir)
    fig_pre, ax_pre = plt.subplots()
    fig_out, ax_out = run.plot_cluster_composition("CL", "sample_id", ax=ax_pre)
    assert ax_out is ax_pre
    plt.close("all")


def test_custom_palette_string(simple_project, temp_dir):
    """String palette name should not raise."""
    import matplotlib.pyplot as plt
    run = _build_run_with_clusters(simple_project, temp_dir)
    fig, ax = run.plot_cluster_composition("CL", "sample_id", palette="Set2")
    assert fig is not None
    plt.close("all")
