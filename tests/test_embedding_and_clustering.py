"""Tests for X layer selection, UMAP, and Leiden clustering."""

import types

import numpy as np
import pandas as pd

from cytofstandard import Project


def _install_fake_mlx_umap(monkeypatch):
    """Install fake mlx_umap and mlx_umap._knn modules for tests."""
    fake_module = types.ModuleType("fake_mlx_umap")
    fake_knn = types.ModuleType("fake_mlx_umap._knn")

    class FakeUMAP:
        def __init__(self, **kwargs):
            self.kwargs = kwargs

        def fit_transform(self, x):
            x = np.asarray(x, dtype=np.float32)
            if x.shape[1] == 1:
                return np.concatenate([x, x], axis=1)
            return x[:, :2]

    def compute_knn(x, k, method="brute", random_state=42, verbose=False):
        x = np.asarray(x, dtype=np.float32)
        n = x.shape[0]
        dmat = np.sqrt(((x[:, None, :] - x[None, :, :]) ** 2).sum(axis=2))
        order = np.argsort(dmat, axis=1)[:, :k]
        dists = np.take_along_axis(dmat, order, axis=1)
        return order.astype(np.int32), dists.astype(np.float32)

    fake_module.UMAP = FakeUMAP
    fake_knn.compute_knn = compute_knn

    monkeypatch.setitem(__import__("sys").modules, "fake_mlx_umap", fake_module)
    monkeypatch.setitem(__import__("sys").modules, "fake_mlx_umap._knn", fake_knn)


def _install_fake_igraph(monkeypatch):
    """Install fake igraph module with community_leiden API."""
    fake_igraph = types.ModuleType("igraph")

    class FakePart:
        def __init__(self, membership):
            self.membership = membership

    class FakeEdges(dict):
        pass

    class FakeGraph:
        def __init__(self, n, edges, directed=False):
            self.n = n
            self.edges = edges
            self.directed = directed
            self.es = FakeEdges()

        def community_leiden(
            self,
            objective_function="modularity",
            weights="weight",
            resolution=1.0,
            n_iterations=2,
            beta=0.01,
        ):
            # Deterministic 2-cluster assignment for testing.
            labels = [0 if i % 2 == 0 else 1 for i in range(self.n)]
            return FakePart(labels)

    fake_igraph.Graph = FakeGraph
    monkeypatch.setitem(__import__("sys").modules, "igraph", fake_igraph)


def _build_ingested_run(
    temp_dir, example_standard_markers_path, example_marker_aliases_path
):
    project_path = temp_dir / "embed_project"
    project = Project.create(
        path=str(project_path),
        project_id="EMBED_TEST",
        project_name="Embed Test",
        standard_marker_file=str(example_standard_markers_path),
        marker_alias_file=str(example_marker_aliases_path),
    )

    run = project.add_run(run_id="run_001")

    file_a = temp_dir / "sample_A.csv"
    file_a.write_text("H3,H3K27me3,ECad\n0,10,5\n2,12,7\n")
    file_b = temp_dir / "sample_B.csv"
    file_b.write_text("H3,H3K27me3,ECad\n10,20,15\n12,22,17\n")

    sample_meta = temp_dir / "sample_meta.csv"
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


def test_set_x_from_layer_overwrites_x(
    temp_dir, example_standard_markers_path, example_marker_aliases_path
):
    """set_x_from_layer should directly overwrite adata.X."""
    run = _build_ingested_run(
        temp_dir, example_standard_markers_path, example_marker_aliases_path
    )
    adata = run.read_adata()

    shifted = np.asarray(adata.layers["raw"], dtype=np.float32) + 3.0
    adata.layers["shifted"] = shifted
    run.save()

    run.set_x_from_layer("shifted")
    updated = run.read_adata()
    assert np.allclose(np.asarray(updated.X), shifted)


def test_compute_umap_stores_embedding_and_graph(
    temp_dir,
    example_standard_markers_path,
    example_marker_aliases_path,
    monkeypatch,
):
    """compute_umap stores embedding metadata and graph artifacts."""
    _install_fake_mlx_umap(monkeypatch)
    run = _build_ingested_run(
        temp_dir, example_standard_markers_path, example_marker_aliases_path
    )

    meta = run.compute_umap(
        markers=["H3", "ECad"],
        source_layer="raw",
        embedding_name="my_umap",
        verbose=True,
        module_name="fake_mlx_umap",
    )

    adata = run.read_adata()
    assert "my_umap" in adata.obsm
    assert "my_umap_connectivities" in adata.obsp
    assert "my_umap_distances" in adata.obsp
    assert "my_umap_knn_indices" in adata.obsm
    assert "my_umap_knn_distances" in adata.obsm
    assert adata.obsm["my_umap"].shape[1] == 2
    assert meta["embedding_key"] == "my_umap"
    assert meta["verbose"] is True
    assert "umap_sec" in meta
    assert "knn_sec" in meta


def test_compute_umap_overwrites_same_name(
    temp_dir,
    example_standard_markers_path,
    example_marker_aliases_path,
    monkeypatch,
):
    """Recomputing same embedding name should overwrite metadata/artifacts."""
    _install_fake_mlx_umap(monkeypatch)
    run = _build_ingested_run(
        temp_dir, example_standard_markers_path, example_marker_aliases_path
    )

    run.compute_umap(
        markers=["H3", "ECad"],
        source_layer="raw",
        embedding_name="my_umap",
        n_neighbors=5,
        module_name="fake_mlx_umap",
    )
    meta2 = run.compute_umap(
        markers=["H3", "H3K27me3"],
        source_layer="raw",
        embedding_name="my_umap",
        n_neighbors=7,
        module_name="fake_mlx_umap",
    )

    adata = run.read_adata()
    assert adata.uns["embeddings"]["my_umap"]["n_neighbors"] == 7
    assert adata.uns["embeddings"]["my_umap"]["markers"] == ["H3", "H3K27me3"]
    assert meta2["n_neighbors"] == 7


def test_cluster_leiden_from_embedding(
    temp_dir,
    example_standard_markers_path,
    example_marker_aliases_path,
    monkeypatch,
):
    """cluster_leiden should cluster from stored embedding graph and persist metadata."""
    _install_fake_mlx_umap(monkeypatch)
    _install_fake_igraph(monkeypatch)
    run = _build_ingested_run(
        temp_dir, example_standard_markers_path, example_marker_aliases_path
    )

    run.compute_umap(
        markers=["H3", "ECad"],
        source_layer="raw",
        embedding_name="my_umap",
        module_name="fake_mlx_umap",
    )

    meta = run.cluster_leiden(
        embedding_name="my_umap",
        resolution=1.2,
        n_iterations=3,
        verbose=True,
    )

    adata = run.read_adata()
    assert "my_umap_leiden" in adata.obs.columns
    assert "my_umap_leiden" in adata.uns["clusterings"]
    assert meta["resolution"] == 1.2
    assert meta["verbose"] is True
    assert "leiden_sec" in meta


def test_cluster_leiden_overwrites_same_cluster_key(
    temp_dir,
    example_standard_markers_path,
    example_marker_aliases_path,
    monkeypatch,
):
    """Repeated clustering with same key should overwrite previous metadata."""
    _install_fake_mlx_umap(monkeypatch)
    _install_fake_igraph(monkeypatch)
    run = _build_ingested_run(
        temp_dir, example_standard_markers_path, example_marker_aliases_path
    )

    run.compute_umap(
        markers=["H3", "ECad"],
        source_layer="raw",
        embedding_name="my_umap",
        module_name="fake_mlx_umap",
    )
    run.cluster_leiden(
        embedding_name="my_umap", cluster_key="leiden_custom", resolution=0.5
    )
    meta2 = run.cluster_leiden(
        embedding_name="my_umap",
        cluster_key="leiden_custom",
        resolution=2.0,
    )

    adata = run.read_adata()
    assert adata.uns["clusterings"]["leiden_custom"]["resolution"] == 2.0
    assert meta2["resolution"] == 2.0
