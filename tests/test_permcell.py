"""Tests for PermCell integration on Run."""

import types

import numpy as np
import pandas as pd

from cytofstandard import Project


def _install_fake_permcell(monkeypatch):
    """Install fake PermCell module with required APIs."""
    module = types.ModuleType("fake_permcell")

    def gaussian_smooth_all_torch(
        data,
        positions,
        bandwidth,
        k=64,
        radius=None,
        chunk_size=2048,
        device=None,
    ):
        data = np.asarray(data, dtype=np.float32)
        return data + 1.0

    def sipsic_like_scores_v3(
        adata,
        marker_sets,
        layer="scaled",
        n_perm=2000,
        seed=0,
        exclude_set=True,
        two_sided=False,
        abs_variant=True,
        exact_max_combinations=50_000,
        progress=True,
        normalize_set_weights=None,
        use_sparse_W=False,
        prefer_permutation=True,
        perm_batch=1024,
    ):
        idx = adata.obs_names
        cols = list(marker_sets.keys())
        z = pd.DataFrame(1.0, index=idx, columns=cols)
        p = pd.DataFrame(0.05, index=idx, columns=cols)
        zdir = pd.DataFrame(1.0, index=idx, columns=cols)
        if abs_variant:
            zabs = pd.DataFrame(2.0, index=idx, columns=cols)
            pabs = pd.DataFrame(0.01, index=idx, columns=cols)
            return z, p, zabs, pabs, zdir
        return z, p, zdir

    module.gaussian_smooth_all_torch = gaussian_smooth_all_torch
    module.sipsic_like_scores_v3 = sipsic_like_scores_v3

    monkeypatch.setitem(__import__("sys").modules, "fake_permcell", module)
    return module


def _build_ingested_run(
    temp_dir, example_standard_markers_path, example_marker_aliases_path
):
    project_path = temp_dir / "permcell_project"
    project = Project.create(
        path=str(project_path),
        project_id="PERMCELL_TEST",
        project_name="PermCell Test",
        standard_marker_file=str(example_standard_markers_path),
        marker_alias_file=str(example_marker_aliases_path),
    )

    run = project.add_run(run_id="run_001")

    file_a = temp_dir / "sample_A.csv"
    file_a.write_text("H3,H3K27me3,ECad\n0,10,5\n2,12,7\n")
    sample_meta = temp_dir / "sample_meta.csv"
    sample_meta.write_text("file_name,sample_id,line_id\nsample_A.csv,S001,MCF7\n")

    run.ingest(
        files=[str(file_a)],
        sample_metadata=str(sample_meta),
        copy_raw=False,
        strict_markers=True,
    )

    adata = run.read_adata()
    adata.obsm["X_umap"] = np.array([[0.0, 0.0], [1.0, 1.0]], dtype=np.float32)
    run.save()
    return run


def test_run_permcell_stores_results(
    temp_dir,
    example_standard_markers_path,
    example_marker_aliases_path,
    monkeypatch,
):
    """run_permcell should store smoothed matrix + score outputs in obsm."""
    _install_fake_permcell(monkeypatch)
    run = _build_ingested_run(
        temp_dir,
        example_standard_markers_path,
        example_marker_aliases_path,
    )

    signatures = {
        "sigA": {"up": ["H3", "ECad"], "down": ["H3K27me3"]},
        "sigB": ["H3"],
    }

    result = run.run_permcell(
        signatures=signatures,
        source_layer="raw",
        positions_key="X_umap",
        result_prefix="pc",
        module_name="fake_permcell",
    )

    adata = run.read_adata()
    assert "pc_smoothed" in adata.obsm
    assert "pc_z" in adata.obsm
    assert "pc_p" in adata.obsm
    assert "pc_zdir" in adata.obsm
    assert "pc_zabs" in adata.obsm
    assert "pc_pabs" in adata.obsm
    assert "pc_raw_z" in adata.obsm
    assert "pc_raw_p" in adata.obsm
    assert "pc_raw_zdir" in adata.obsm
    assert "permcell" in adata.uns
    assert "pc" in adata.uns["permcell"]
    assert result["z"].shape[1] == 2
    assert result["raw_z"].shape[1] == 2


def test_run_permcell_overwrites_same_prefix(
    temp_dir,
    example_standard_markers_path,
    example_marker_aliases_path,
    monkeypatch,
):
    """Calling run_permcell twice with same prefix should overwrite metadata/results."""
    _install_fake_permcell(monkeypatch)
    run = _build_ingested_run(
        temp_dir,
        example_standard_markers_path,
        example_marker_aliases_path,
    )

    run.run_permcell(
        signatures={"sigA": ["H3"]},
        source_layer="raw",
        positions_key="X_umap",
        result_prefix="pc",
        abs_variant=True,
        module_name="fake_permcell",
    )
    run.run_permcell(
        signatures={"sigB": ["ECad"]},
        source_layer="raw",
        positions_key="X_umap",
        result_prefix="pc",
        abs_variant=False,
        module_name="fake_permcell",
    )

    adata = run.read_adata()
    meta = adata.uns["permcell"]["pc"]
    assert meta["params"]["abs_variant"] is False
    assert meta["signature_names"] == ["sigB"]
    assert meta["result_keys_smoothed"]["zabs"] is None
    assert "pc_zabs" not in adata.obsm
    assert "pc_pabs" not in adata.obsm


def test_run_permcell_accepts_explicit_module_object(
    temp_dir,
    example_standard_markers_path,
    example_marker_aliases_path,
    monkeypatch,
):
    """run_permcell should work when module object is passed directly."""
    mod = _install_fake_permcell(monkeypatch)
    run = _build_ingested_run(
        temp_dir,
        example_standard_markers_path,
        example_marker_aliases_path,
    )

    # Remove named module to ensure explicit object is used.
    __import__("sys").modules.pop("fake_permcell", None)

    run.run_permcell(
        signatures={"sigA": ["H3"]},
        source_layer="raw",
        positions_key="X_umap",
        result_prefix="pc_obj",
        permcell_module=mod,
        module_name="fake_permcell",
    )

    adata = run.read_adata()
    assert "pc_obj_z" in adata.obsm


def test_permcell_to_adata_smoothed_and_raw(
    temp_dir,
    example_standard_markers_path,
    example_marker_aliases_path,
    monkeypatch,
):
    """permcell_to_adata should export selected score matrix with signature var_names."""
    _install_fake_permcell(monkeypatch)
    run = _build_ingested_run(
        temp_dir,
        example_standard_markers_path,
        example_marker_aliases_path,
    )

    run.run_permcell(
        signatures={"sigA": ["H3"], "sigB": ["ECad"]},
        source_layer="raw",
        positions_key="X_umap",
        result_prefix="pc",
        module_name="fake_permcell",
    )

    ad_smoothed = run.permcell_to_adata(
        result_prefix="pc",
        score="z",
        smoothed=True,
    )
    ad_raw = run.permcell_to_adata(
        result_prefix="pc",
        score="z",
        smoothed=False,
    )

    assert ad_smoothed.n_vars == 2
    assert ad_smoothed.var_names.tolist() == ["sigA", "sigB"]
    assert ad_smoothed.X.shape[1] == 2
    assert "X_umap" in ad_smoothed.obsm

    assert ad_raw.n_vars == 2
    assert ad_raw.var_names.tolist() == ["sigA", "sigB"]
