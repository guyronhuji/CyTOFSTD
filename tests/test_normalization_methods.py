"""Integration tests for normalization against the real cytof_transform (>=0.2.0).

These exercise every normalization family and option the external module exposes.
The module is not a hard dependency of cytofstandard, so the whole file is
skipped when it is not importable.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pytest

# cytof-transform is usually checked out beside this repo rather than installed.
_SIBLING = Path(__file__).resolve().parent.parent.parent / "cytof-transform"
if _SIBLING.is_dir() and str(_SIBLING) not in sys.path:
    sys.path.insert(0, str(_SIBLING))

ct = pytest.importorskip("cytof_transform", reason="cytof_transform not importable")

pytestmark = pytest.mark.skipif(
    not hasattr(ct, "cytof_normalize"),
    reason="requires cytof_transform >=0.2.0 (cytof_normalize)",
)

# Only markers in the example standard registry survive ingestion.
MARKERS = ["H3", "H3K27me3", "H3K9me2", "ECad", "EpCAM", "Ki67", "DAPI"]
CONTROLS = ["H3", "H3K27me3", "DAPI"]
TARGETS = ["ECad", "EpCAM"]
COVARIATE = "Ki67"


@pytest.fixture
def norm_run(simple_project, temp_dir):
    """A two-sample run with enough cells for regression and shrinkage."""
    run = simple_project.add_run(run_id="run_ct")
    rng = np.random.default_rng(0)

    files = []
    for name in ("A", "B"):
        path = temp_dir / f"sample_{name}.csv"
        data = np.abs(rng.normal(200, 50, size=(120, len(MARKERS))))
        body = "\n".join(",".join(f"{v:.2f}" for v in row) for row in data)
        path.write_text(",".join(MARKERS) + "\n" + body + "\n")
        files.append(str(path))

    meta = temp_dir / "sample_metadata.csv"
    meta.write_text(
        "file_name,sample_id,line_id,condition\n"
        "sample_A.csv,S001,MCF7,ctrl\n"
        "sample_B.csv,S002,T47D,treat\n"
    )
    run.ingest(
        files=files,
        sample_metadata=str(meta),
        copy_raw=False,
        strict_markers=False,
    )
    return run


def _assert_normalized(summary, run, expected_method):
    """Common post-conditions for any successful normalization."""
    adata = run.read_adata()
    assert "normalized" in adata.layers
    assert "normalized_z" in adata.layers
    assert not adata.obs["norm_tech_factor"].isna().any()
    assert adata.layers["normalized"].shape == (adata.n_obs, adata.n_vars)
    assert np.isfinite(adata.layers["normalized"]).all()

    assert summary["method"] == expected_method
    assert summary["entry_point"] == "cytof_normalize"
    # The summary is persisted into uns as JSON, so it must serialize.
    json.dumps(summary, sort_keys=True)
    assert adata.uns["normalization"]["latest"]["method"] == expected_method
    return adata


def test_regress_is_the_default(norm_run):
    summary = norm_run.normalize_with_cytof_transform(
        control_markers=CONTROLS,
        markers_to_correct=TARGETS,
    )
    _assert_normalized(summary, norm_run, "regress")
    assert summary["gamma_mode"] == "per_marker"
    assert summary["tech_factor_kind"] == "pc1"
    assert set(summary["gamma_by_group"]) == {"S001", "S002"}


def test_divide_method_returns_arcsinh_scale(norm_run):
    summary = norm_run.normalize_with_cytof_transform(
        control_markers=CONTROLS,
        markers_to_correct=TARGETS,
        method="divide",
    )
    adata = _assert_normalized(summary, norm_run, "divide")
    assert summary["tech_factor_kind"] == "permeabilization_factor"
    # raw -> divide -> arcsinh, so values stay on the arcsinh scale.
    assert adata.layers["normalized"].max() < 20


def test_divide_rejects_arcsinh_input(norm_run):
    with pytest.raises(ValueError, match="raw"):
        norm_run.normalize_with_cytof_transform(
            control_markers=CONTROLS,
            markers_to_correct=TARGETS,
            method="divide",
            input_is_arcsinh=True,
        )


def test_rejects_unknown_method(norm_run):
    with pytest.raises(ValueError, match="Unknown method"):
        norm_run.normalize_with_cytof_transform(
            control_markers=CONTROLS,
            markers_to_correct=TARGETS,
            method="bogus",
        )


@pytest.mark.parametrize("gamma_mode", ["per_marker", "single", "shrink"])
def test_gamma_modes(norm_run, gamma_mode):
    summary = norm_run.normalize_with_cytof_transform(
        control_markers=CONTROLS,
        markers_to_correct=TARGETS,
        gamma_mode=gamma_mode,
    )
    _assert_normalized(summary, norm_run, "regress")
    assert summary["gamma_mode"] == gamma_mode
    assert summary["gamma_shrink_by_group"], "gamma diagnostics should be recorded"


@pytest.mark.parametrize("shrink_target", ["control", "global"])
def test_shrink_targets(norm_run, shrink_target):
    summary = norm_run.normalize_with_cytof_transform(
        control_markers=CONTROLS,
        markers_to_correct=TARGETS,
        gamma_mode="shrink",
        shrink_target=shrink_target,
    )
    _assert_normalized(summary, norm_run, "regress")
    assert summary["shrink_target"] == shrink_target


def test_shrink_stability_uses_obs_column(norm_run):
    summary = norm_run.normalize_with_cytof_transform(
        control_markers=CONTROLS,
        markers_to_correct=TARGETS,
        gamma_mode="shrink_stability",
        stability_group_col="condition",
        min_group_cells=10,
    )
    adata = _assert_normalized(summary, norm_run, "regress")
    assert summary["stability_group_col"] == "condition"
    assert summary["min_group_cells"] == 10
    # The helper column handed to cytof_transform must not leak into the layers.
    assert adata.layers["normalized"].shape[1] == adata.n_vars


def test_protect_covariates(norm_run):
    summary = norm_run.normalize_with_cytof_transform(
        control_markers=CONTROLS,
        markers_to_correct=TARGETS,
        protect_covariates=[COVARIATE],
    )
    _assert_normalized(summary, norm_run, "regress")
    assert summary["protect_covariates"] == [COVARIATE]


def test_compartments(norm_run):
    summary = norm_run.normalize_with_cytof_transform(
        control_markers=CONTROLS,
        markers_to_correct=TARGETS,
        groupby_col="line_id",
        compartment_col="condition",
    )
    _assert_normalized(summary, norm_run, "regress")
    assert summary["compartment_col"] == "condition"
    assert summary["compartments"] == ["ctrl", "treat"]


def test_missing_obs_column_is_reported(norm_run):
    with pytest.raises(ValueError, match="compartment_col"):
        norm_run.normalize_with_cytof_transform(
            control_markers=CONTROLS,
            markers_to_correct=TARGETS,
            compartment_col="not_a_column",
        )


def test_parameters_are_logged_to_provenance(norm_run):
    norm_run.normalize_with_cytof_transform(
        control_markers=CONTROLS,
        markers_to_correct=["ECad"],
        gamma_mode="shrink",
        shrink_target="global",
        protect_covariates=[COVARIATE],
    )

    log_path = norm_run.path / "logs" / "provenance.jsonl"
    events = [json.loads(line) for line in log_path.read_text().splitlines() if line]
    normalized = [e for e in events if e["event_type"] == "run_normalized"]
    assert normalized, "run_normalized event was not logged"

    event = normalized[-1]
    for key in (
        "method",
        "entry_point",
        "module",
        "module_version",
        "gamma_mode",
        "shrink_target",
        "protect_covariates",
        "stability_group_col",
        "min_group_cells",
        "compartment_col",
        "control_markers",
        "markers_to_correct",
        "arcsinh_cofactor",
        "anchor_to_median",
        "zscore",
    ):
        assert key in event, f"{key} missing from provenance event"

    assert event["gamma_mode"] == "shrink"
    assert event["shrink_target"] == "global"
    assert event["protect_covariates"] == [COVARIATE]


def test_evaluate_marker_intensity_regime(norm_run):
    table = norm_run.evaluate_marker_intensity_regime(
        candidate_markers=TARGETS,
        layer="raw",
    )
    assert len(table) == len(TARGETS)


def test_compute_marker_tech_correlations_from_controls(norm_run):
    corr, tech_factor = norm_run.compute_marker_tech_correlations(
        control_markers=CONTROLS,
        layer="raw",
    )
    assert set(TARGETS).issubset(set(corr.index))
    assert len(tech_factor) == norm_run.read_adata().n_obs


def test_compute_marker_tech_correlations_from_stored_factor(norm_run):
    norm_run.normalize_with_cytof_transform(
        control_markers=CONTROLS, markers_to_correct=TARGETS
    )
    corr, tech_factor = norm_run.compute_marker_tech_correlations(
        layer="raw",
        use_stored_tech_factor=True,
    )
    assert set(TARGETS).issubset(set(corr.index))
    # The stored per-cell factor is reused verbatim, not recomputed.
    stored = norm_run.read_adata().obs["norm_tech_factor"].astype(float)
    assert np.allclose(tech_factor.values, stored.values)


def test_compute_marker_tech_correlations_requires_a_factor(norm_run):
    with pytest.raises(ValueError, match="control_markers is required"):
        norm_run.compute_marker_tech_correlations(layer="raw")


def test_qc_wrappers_reject_empty_group(norm_run):
    with pytest.raises(ValueError, match="No cells found"):
        norm_run.evaluate_marker_intensity_regime(
            candidate_markers=TARGETS,
            group_value="not_a_sample",
        )


def test_history_accumulates_each_call(norm_run):
    norm_run.normalize_with_cytof_transform(
        control_markers=CONTROLS, markers_to_correct=TARGETS
    )
    norm_run.normalize_with_cytof_transform(
        control_markers=CONTROLS, markers_to_correct=TARGETS, gamma_mode="single"
    )

    history = norm_run.read_adata().uns["normalization"]["history"]
    assert len(history) == 2
    assert [json.loads(h)["gamma_mode"] for h in history] == ["per_marker", "single"]
