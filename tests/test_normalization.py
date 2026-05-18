"""Tests for normalization integration with external cytof_transform module."""

import types

import numpy as np
import pandas as pd


def _install_fake_cytof_transform(monkeypatch):
    """Install a fake cytof_transform module in sys.modules."""
    fake_module = types.ModuleType("fake_cytof_transform")

    class FakeConfig:
        def __init__(
            self,
            control_markers,
            markers_to_correct,
            use_compartments=False,
            n_pcs_for_T=1,
            anchor_to_median=True,
            zscore=True,
            line_col=None,
        ):
            self.control_markers = control_markers
            self.markers_to_correct = markers_to_correct
            self.use_compartments = use_compartments
            self.n_pcs_for_T = n_pcs_for_T
            self.anchor_to_median = anchor_to_median
            self.zscore = zscore
            self.line_col = line_col

    def fake_transform_global(asinh_data, config):
        corrected = asinh_data.copy()
        for marker in config.markers_to_correct:
            corrected[marker] = corrected[marker] + 1.0
        residuals_z = corrected.copy()
        tech = pd.Series(
            np.linspace(0.0, 1.0, corrected.shape[0]),
            index=corrected.index,
            name="tech1",
        )
        gamma = {marker: 0.1 for marker in config.markers_to_correct}
        alpha = {marker: 0.0 for marker in config.markers_to_correct}
        return types.SimpleNamespace(
            corrected=corrected,
            residuals_z=residuals_z,
            tech_factor=tech,
            gamma=gamma,
            alpha=alpha,
        )

    def fake_plot_tech_factor_qc(**kwargs):
        return "tech_factor_plot"

    def fake_plot_marker_correlations_qc(**kwargs):
        return pd.DataFrame({"marker": ["H3"], "corr_pre": [0.1], "corr_post": [0.01]})

    def fake_plot_gamma_qc(**kwargs):
        return "gamma_plot"

    fake_module.CytofTransformConfig = FakeConfig
    fake_module.cytof_transform_global = fake_transform_global
    fake_module.plot_tech_factor_qc = fake_plot_tech_factor_qc
    fake_module.plot_marker_correlations_qc = fake_plot_marker_correlations_qc
    fake_module.plot_gamma_qc = fake_plot_gamma_qc

    monkeypatch.setitem(__import__("sys").modules, "fake_cytof_transform", fake_module)


def _ingest_two_sample_run(simple_project, temp_dir):
    run = simple_project.add_run(run_id="run_norm")

    file_a = temp_dir / "sample_A.csv"
    file_a.write_text("H3,H3K27me3,ECad\n100,200,300\n150,250,350\n")
    file_b = temp_dir / "sample_B.csv"
    file_b.write_text("H3,H3K27me3,ECad\n110,210,310\n160,260,360\n")

    sample_meta = temp_dir / "sample_metadata.csv"
    sample_meta.write_text(
        "file_name,sample_id,line_id\nsample_A.csv,S001,MCF7\nsample_B.csv,S002,T47D\n"
    )

    run.ingest(
        files=[str(file_a), str(file_b)],
        sample_metadata=str(sample_meta),
        copy_raw=False,
        strict_markers=True,
    )
    return run


def test_normalize_with_cytof_transform_per_sample(
    simple_project, temp_dir, monkeypatch
):
    """Normalization runs per sample and stores outputs in layers/metadata."""
    _install_fake_cytof_transform(monkeypatch)
    run = _ingest_two_sample_run(simple_project, temp_dir)

    summary = run.normalize_with_cytof_transform(
        control_markers=["H3"],
        markers_to_correct=["ECad", "H3K27me3"],
        source_layer="raw",
        groupby_col="sample_id",
        module_name="fake_cytof_transform",
    )

    adata = run.read_adata()
    assert "normalized" in adata.layers
    assert "normalized_z" in adata.layers
    assert "norm_tech_factor" in adata.obs.columns
    assert not adata.obs["norm_tech_factor"].isna().any()
    assert summary["groupby_col"] == "sample_id"
    assert sorted(summary["groups"]) == ["S001", "S002"]
    assert set(summary["gamma_by_group"].keys()) == {"S001", "S002"}


def test_normalization_plot_wrappers(simple_project, temp_dir, monkeypatch):
    """QC plotting wrappers call the external cytof_transform plot helpers."""
    _install_fake_cytof_transform(monkeypatch)
    run = _ingest_two_sample_run(simple_project, temp_dir)

    run.normalize_with_cytof_transform(
        control_markers=["H3"],
        markers_to_correct=["ECad"],
        source_layer="raw",
        groupby_col="sample_id",
        module_name="fake_cytof_transform",
    )

    tech_plot = run.plot_normalization_tech_factor_qc(
        control_markers=["H3"],
        layer="raw",
        group_value="S001",
        module_name="fake_cytof_transform",
    )
    corr_df = run.plot_normalization_marker_correlations_qc(
        pre_layer="raw",
        post_layer="normalized",
        group_value="S001",
        module_name="fake_cytof_transform",
    )
    gamma_plot = run.plot_normalization_gamma_qc(
        group_value="S001",
        module_name="fake_cytof_transform",
    )

    assert tech_plot == "tech_factor_plot"
    assert "corr_pre" in corr_df.columns
    assert gamma_plot == "gamma_plot"
