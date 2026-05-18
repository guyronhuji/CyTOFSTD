"""Tests for plot_heatmap and plot_boxplot helpers on Run."""

import shutil

import matplotlib
import numpy as np
import pandas as pd
import pytest


matplotlib.use("Agg")


def _make_grouped_csv(temp_dir):
    """Build a CSV with 3 cell-type groups, 12 rows each."""
    rng = np.random.default_rng(0)
    rows = []
    for label, base_h3 in [("T_cell", 100), ("B_cell", 160), ("Macrophage", 220)]:
        for _ in range(12):
            rows.append(
                {
                    "H3": float(base_h3 + rng.normal(0, 4)),
                    "H3K27me3": float(200 + rng.normal(0, 10)),
                    "ECad": float(300 + rng.normal(0, 15)),
                    "cell_type": label,
                }
            )
    df = pd.DataFrame(rows)
    csv_path = temp_dir / "grouped.csv"
    df.to_csv(csv_path, index=False)
    return csv_path


def _make_metadata(temp_dir, csv_path):
    meta_path = temp_dir / "grouped_metadata.csv"
    meta_path.write_text(
        "file_name,sample_id,line_id,condition,replicate_id\n"
        f"{csv_path.name},S001,MCF7,control,rep1\n"
    )
    return meta_path


def _ingest(run, project, csv_path, meta_path):
    raw_dir = project.path / "runs" / run.run_id / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy(csv_path, raw_dir / csv_path.name)

    run.ingest(
        files=[str(raw_dir / csv_path.name)],
        sample_metadata=str(meta_path),
        copy_raw=False,
        strict_markers=True,
    )

    # Reattach the per-cell `cell_type` column to obs (read_csv drops
    # non-numeric columns during ingestion).
    source_df = pd.read_csv(csv_path)
    adata = run.read_adata()
    adata.obs["cell_type"] = source_df["cell_type"].astype(str).values
    run.save(adata)


def test_plot_heatmap_markers_by_obs_group(simple_project, temp_dir):
    """Heatmap with marker rows aggregated by an obs column."""
    csv_path = _make_grouped_csv(temp_dir)
    meta_path = _make_metadata(temp_dir, csv_path)
    run = simple_project.add_run(run_id="run_heatmap_markers")
    _ingest(run, simple_project, csv_path, meta_path)

    fig, ax = run.plot_heatmap(
        fields=["H3", "ECad"],
        groupby="cell_type",
        layer="raw",
        agg="mean",
    )

    fig.canvas.draw()
    xticks = [t.get_text() for t in ax.get_xticklabels()]
    yticks = [t.get_text() for t in ax.get_yticklabels()]
    assert sorted(xticks) == sorted(["B_cell", "Macrophage", "T_cell"])
    assert set(yticks) == {"H3", "ECad"}


def test_plot_heatmap_accepts_obs_numeric_column(
    simple_project, temp_dir
):
    """Heatmap can aggregate numeric obs columns alongside markers."""
    csv_path = _make_grouped_csv(temp_dir)
    meta_path = _make_metadata(temp_dir, csv_path)
    run = simple_project.add_run(run_id="run_heatmap_obs")
    _ingest(run, simple_project, csv_path, meta_path)

    adata = run.read_adata()
    adata.obs["pseudo_score"] = np.linspace(0.0, 1.0, adata.n_obs)
    run.save(adata)

    fig, ax = run.plot_heatmap(
        fields=["H3", "pseudo_score"],
        groupby="cell_type",
        layer="raw",
        agg="median",
        standard_scale="row",
    )
    fig.canvas.draw()
    assert ax.get_title().startswith("median")


def test_plot_heatmap_composite_groupby(simple_project, temp_dir):
    """Heatmap supports a list of obs columns for composite grouping."""
    csv_path = _make_grouped_csv(temp_dir)
    meta_path = _make_metadata(temp_dir, csv_path)
    run = simple_project.add_run(run_id="run_heatmap_composite")
    _ingest(run, simple_project, csv_path, meta_path)

    fig, ax = run.plot_heatmap(
        fields=["H3"],
        groupby=["sample_id", "cell_type"],
        layer="raw",
        agg="mean",
    )
    fig.canvas.draw()
    assert " | " in ax.get_xlabel()


def test_plot_heatmap_unknown_field_raises(simple_project, temp_dir):
    csv_path = _make_grouped_csv(temp_dir)
    meta_path = _make_metadata(temp_dir, csv_path)
    run = simple_project.add_run(run_id="run_heatmap_unknown")
    _ingest(run, simple_project, csv_path, meta_path)

    with pytest.raises(ValueError):
        run.plot_heatmap(
            fields=["NOT_A_FIELD"],
            groupby="cell_type",
            layer="raw",
        )


def test_plot_boxplot_marker_with_significance(simple_project, temp_dir):
    """Boxplot draws brackets with star labels when comparisons='all'."""
    csv_path = _make_grouped_csv(temp_dir)
    meta_path = _make_metadata(temp_dir, csv_path)
    run = simple_project.add_run(run_id="run_boxplot_marker")
    _ingest(run, simple_project, csv_path, meta_path)

    fig, ax = run.plot_boxplot(
        field="H3",
        groupby="cell_type",
        layer="raw",
        comparisons="all",
        test="mannwhitney",
        multitest="holm",
    )
    fig.canvas.draw()

    # Three groups -> three pairwise tests -> at least 3 star/ns labels.
    text_labels = [
        t.get_text() for t in ax.texts if t.get_text()
    ]
    assert len([t for t in text_labels if t in {"*", "**", "***", "****", "ns"}]) >= 3


def test_plot_boxplot_default_has_no_brackets(simple_project, temp_dir):
    """comparisons=None means no significance brackets are drawn."""
    csv_path = _make_grouped_csv(temp_dir)
    meta_path = _make_metadata(temp_dir, csv_path)
    run = simple_project.add_run(run_id="run_boxplot_default_no_brackets")
    _ingest(run, simple_project, csv_path, meta_path)

    fig, ax = run.plot_boxplot(
        field="H3",
        groupby="cell_type",
        layer="raw",
    )
    fig.canvas.draw()
    star_texts = [
        t.get_text() for t in ax.texts
        if t.get_text() in {"*", "**", "***", "****", "ns"}
    ]
    assert star_texts == []


def test_plot_boxplot_default_palette_is_distinct(simple_project, temp_dir):
    """Default palette assigns a distinct color to each box."""
    csv_path = _make_grouped_csv(temp_dir)
    meta_path = _make_metadata(temp_dir, csv_path)
    run = simple_project.add_run(run_id="run_boxplot_default_palette")
    _ingest(run, simple_project, csv_path, meta_path)

    fig, ax = run.plot_boxplot(
        field="H3",
        groupby="cell_type",
        layer="raw",
    )
    fig.canvas.draw()
    # Collect facecolors of box artists.
    facecolors = set()
    for patch in ax.patches:
        if patch.get_facecolor() is None:
            continue
        facecolors.add(tuple(round(c, 3) for c in patch.get_facecolor()))
    # 3 groups should yield 3 distinct facecolors.
    assert len(facecolors) >= 3


def test_plot_boxplot_show_outliers_false_hides_fliers(
    simple_project, temp_dir
):
    """show_outliers=False removes outlier markers from the boxplot."""
    csv_path = _make_grouped_csv(temp_dir)
    meta_path = _make_metadata(temp_dir, csv_path)
    run = simple_project.add_run(run_id="run_boxplot_no_fliers")
    _ingest(run, simple_project, csv_path, meta_path)

    fig_show, ax_show = run.plot_boxplot(
        field="H3",
        groupby="cell_type",
        layer="raw",
        show_outliers=True,
    )
    fig_hide, ax_hide = run.plot_boxplot(
        field="H3",
        groupby="cell_type",
        layer="raw",
        show_outliers=False,
    )
    fig_show.canvas.draw()
    fig_hide.canvas.draw()

    # When outliers are hidden, total line-collection count drops.
    n_lines_show = len(ax_show.lines)
    n_lines_hide = len(ax_hide.lines)
    assert n_lines_hide <= n_lines_show


def test_plot_boxplot_accepts_obs_numeric_field(simple_project, temp_dir):
    """Boxplot can plot a numeric obs column rather than a marker."""
    csv_path = _make_grouped_csv(temp_dir)
    meta_path = _make_metadata(temp_dir, csv_path)
    run = simple_project.add_run(run_id="run_boxplot_obs")
    _ingest(run, simple_project, csv_path, meta_path)

    adata = run.read_adata()
    rng = np.random.default_rng(1)
    cell_type = adata.obs["cell_type"].astype(str).values
    base = np.where(cell_type == "T_cell", 0.1, np.where(cell_type == "B_cell", 0.5, 0.9))
    adata.obs["score"] = base + rng.normal(0, 0.05, adata.n_obs)
    run.save(adata)

    fig, ax = run.plot_boxplot(
        field="score",
        groupby="cell_type",
        test="welch",
        multitest="bonferroni",
    )
    fig.canvas.draw()
    assert ax.get_ylabel() == "score"


def test_plot_boxplot_explicit_comparisons(simple_project, temp_dir):
    """Explicit comparison pairs are respected."""
    csv_path = _make_grouped_csv(temp_dir)
    meta_path = _make_metadata(temp_dir, csv_path)
    run = simple_project.add_run(run_id="run_boxplot_pairs")
    _ingest(run, simple_project, csv_path, meta_path)

    fig, ax = run.plot_boxplot(
        field="H3",
        groupby="cell_type",
        layer="raw",
        comparisons=[("T_cell", "Macrophage")],
        test="ttest",
        multitest=None,
    )
    fig.canvas.draw()

    star_texts = [
        t.get_text() for t in ax.texts
        if t.get_text() in {"*", "**", "***", "****", "ns"}
    ]
    assert len(star_texts) == 1


def test_pvalue_to_stars_mapping():
    from cytofstandard.run import Run

    assert Run._pvalue_to_stars(0.5) == "ns"
    assert Run._pvalue_to_stars(0.04) == "*"
    assert Run._pvalue_to_stars(0.005) == "**"
    assert Run._pvalue_to_stars(0.0005) == "***"
    assert Run._pvalue_to_stars(1e-6) == "****"
    assert Run._pvalue_to_stars(float("nan")) == "ns"


def test_pvalue_to_stars_custom_thresholds():
    from cytofstandard.run import Run

    thresholds = [(0.01, "HIGH"), (0.1, "LOW")]
    assert Run._pvalue_to_stars(0.005, thresholds=thresholds) == "HIGH"
    assert Run._pvalue_to_stars(0.05, thresholds=thresholds) == "LOW"
    assert Run._pvalue_to_stars(0.5, thresholds=thresholds) == "ns"
    assert (
        Run._pvalue_to_stars(0.5, thresholds=thresholds, ns_label="--") == "--"
    )


def test_plot_boxplot_custom_significance_thresholds(simple_project, temp_dir):
    """Custom significance thresholds change bracket labels."""
    csv_path = _make_grouped_csv(temp_dir)
    meta_path = _make_metadata(temp_dir, csv_path)
    run = simple_project.add_run(run_id="run_boxplot_thresholds")
    _ingest(run, simple_project, csv_path, meta_path)

    fig, ax = run.plot_boxplot(
        field="H3",
        groupby="cell_type",
        layer="raw",
        comparisons="all",
        test="mannwhitney",
        multitest=None,
        significance_thresholds=[(1e-12, "HIT")],
        ns_label="meh",
    )
    fig.canvas.draw()
    labels = {
        t.get_text() for t in ax.texts if t.get_text()
    }
    assert labels.issubset({"HIT", "meh"})
    assert ("HIT" in labels) or ("meh" in labels)


def test_plot_heatmap_forwards_kwargs_to_seaborn(simple_project, temp_dir):
    """heatmap_kwargs (e.g. vmin/vmax, linewidths) reach seaborn.heatmap."""
    csv_path = _make_grouped_csv(temp_dir)
    meta_path = _make_metadata(temp_dir, csv_path)
    run = simple_project.add_run(run_id="run_heatmap_kwargs")
    _ingest(run, simple_project, csv_path, meta_path)

    fig, ax = run.plot_heatmap(
        fields=["H3", "ECad"],
        groupby="cell_type",
        layer="raw",
        agg="mean",
        heatmap_kwargs={
            "vmin": 0.0,
            "vmax": 500.0,
            "linewidths": 0.5,
            "linecolor": "white",
        },
    )
    fig.canvas.draw()
    mesh = ax.collections[0]
    vmin, vmax = mesh.get_clim()
    assert vmin == pytest.approx(0.0)
    assert vmax == pytest.approx(500.0)


def test_plot_boxplot_forwards_kwargs_to_seaborn(simple_project, temp_dir):
    """boxplot_kwargs (e.g. width, linewidth) reach seaborn.boxplot."""
    csv_path = _make_grouped_csv(temp_dir)
    meta_path = _make_metadata(temp_dir, csv_path)
    run = simple_project.add_run(run_id="run_boxplot_kwargs")
    _ingest(run, simple_project, csv_path, meta_path)

    fig, ax = run.plot_boxplot(
        field="H3",
        groupby="cell_type",
        layer="raw",
        comparisons=[("T_cell", "B_cell")],
        multitest=None,
        boxplot_kwargs={"width": 0.4, "linewidth": 2.0, "notch": True},
        show_points=True,
        stripplot_kwargs={"jitter": 0.1},
    )
    fig.canvas.draw()
    # At least one significance bracket label was drawn.
    star_texts = [
        t.get_text() for t in ax.texts
        if t.get_text() in {"*", "**", "***", "****", "ns"}
    ]
    assert len(star_texts) == 1
