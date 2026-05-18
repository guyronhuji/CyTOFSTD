"""Tests for per-run provenance logging."""

import json
import types

import numpy as np
import pandas as pd


def _read_jsonl(path):
    if not path.exists():
        return []
    lines = [line for line in path.read_text().splitlines() if line.strip()]
    return [json.loads(line) for line in lines]


def _install_fake_cytof_transform(monkeypatch):
    fake_module = types.ModuleType("fake_cytof_transform_prov")

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

    def fake_transform_global(asinh_data, config):
        corrected = asinh_data.copy()
        residuals_z = asinh_data.copy()
        tech = pd.Series(
            np.linspace(0.0, 1.0, corrected.shape[0]),
            index=corrected.index,
            name="tech1",
        )
        gamma = {m: 0.1 for m in config.markers_to_correct}
        alpha = {m: 0.0 for m in config.markers_to_correct}
        return types.SimpleNamespace(
            corrected=corrected,
            residuals_z=residuals_z,
            tech_factor=tech,
            gamma=gamma,
            alpha=alpha,
        )

    fake_module.CytofTransformConfig = FakeConfig
    fake_module.cytof_transform_global = fake_transform_global
    monkeypatch.setitem(
        __import__("sys").modules, "fake_cytof_transform_prov", fake_module
    )


def test_run_logs_include_zscore_event(simple_project, temp_dir):
    """Per-run provenance should log z-scoring operations."""
    run = simple_project.add_run(run_id="run_prov_z")

    file_a = temp_dir / "prov_a.csv"
    file_a.write_text("H3,ECad\n0,10\n2,12\n")
    sample_meta = temp_dir / "prov_meta.csv"
    sample_meta.write_text("file_name,sample_id,line_id\nprov_a.csv,S001,MCF7\n")

    run.ingest(
        files=[str(file_a)],
        sample_metadata=str(sample_meta),
        copy_raw=False,
        strict_markers=True,
    )
    run.zscore_markers_balanced(source_layer="raw", output_layer="z_test")

    log_path = run.path / "logs" / "provenance.jsonl"
    events = _read_jsonl(log_path)
    event_types = [entry.get("event_type") for entry in events]
    assert "run_zscored" in event_types


def test_run_logs_include_normalization_event(simple_project, temp_dir, monkeypatch):
    """Per-run provenance should log normalization operations."""
    _install_fake_cytof_transform(monkeypatch)

    run = simple_project.add_run(run_id="run_prov_norm")
    file_a = temp_dir / "prov_norm_a.csv"
    file_a.write_text("H3,H3K27me3,ECad\n100,200,300\n150,250,350\n")
    sample_meta = temp_dir / "prov_norm_meta.csv"
    sample_meta.write_text("file_name,sample_id,line_id\nprov_norm_a.csv,S001,MCF7\n")

    run.ingest(
        files=[str(file_a)],
        sample_metadata=str(sample_meta),
        copy_raw=False,
        strict_markers=True,
    )

    run.normalize_with_cytof_transform(
        control_markers=["H3"],
        markers_to_correct=["ECad"],
        source_layer="raw",
        groupby_col="sample_id",
        module_name="fake_cytof_transform_prov",
    )

    log_path = run.path / "logs" / "provenance.jsonl"
    events = _read_jsonl(log_path)
    event_types = [entry.get("event_type") for entry in events]
    assert "run_normalized" in event_types
