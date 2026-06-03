"""Compute — UMAP embedding and cell clustering from within the app."""

import sys
from pathlib import Path

import streamlit as st

sys.path.insert(0, str(Path(__file__).parent.parent))

from _shared import bump_adata_version, get_run, page_header, render_sidebar  # noqa: E402

st.set_page_config(page_title="Compute | CyTOF Standard", layout="wide")
render_sidebar()
page_header("Compute", subtitle="Run UMAP embedding and cell clustering", icon="⚙️")

run = get_run()

if run.status != "ingested":
    st.warning(f"This run has not been ingested yet (status: **{run.status}**).")
    st.stop()

adata = run.read_adata()
all_markers = run.markers_all or list(adata.var_names)

# ── Current state ─────────────────────────────────────────────────────────────
existing_embeddings = [k for k in adata.obsm if k.startswith("X_")]
cluster_cols = [
    c for c in adata.obs.columns
    if "leiden" in c.lower() or "cluster" in c.lower() or "flowsom" in c.lower()
]

c1, c2, c3 = st.columns(3)
c1.metric("Cells", f"{adata.n_obs:,}")
c2.metric("Embeddings", len(existing_embeddings))
c3.metric("Clusterings", len(cluster_cols))

if existing_embeddings:
    st.caption("Embeddings: " + "  ·  ".join(existing_embeddings))
if cluster_cols:
    st.caption("Clusterings: " + "  ·  ".join(cluster_cols))

st.divider()

tab_umap, tab_leiden, tab_flowsom = st.tabs(["UMAP", "Leiden", "FlowSOM"])

# ── UMAP ──────────────────────────────────────────────────────────────────────
with tab_umap:
    st.subheader("Compute UMAP embedding")
    st.write(
        "Computes a 2-D embedding using `mlx-umap` and stores the KNN graph required for Leiden clustering."
    )

    umap_markers = st.multiselect(
        "Markers",
        all_markers,
        default=all_markers,
        help="Markers used to compute the embedding. Typically all functional/extracellular markers.",
        key="umap_markers",
    )

    col_a, col_b, col_c = st.columns(3)
    umap_n_neighbors = col_a.slider(
        "n_neighbors", min_value=5, max_value=50, value=15,
        help="Size of the local neighbourhood for graph construction. Higher = more global structure.",
    )
    umap_min_dist = col_b.slider(
        "min_dist", min_value=0.01, max_value=0.99, value=0.10, step=0.01,
        help="Minimum distance between points in the 2-D layout. Lower = tighter clusters.",
    )
    umap_name = col_c.text_input(
        "Embedding name", value="X_umap",
        help="Key stored in adata.obsm. Use a new name to keep existing embeddings.",
    )

    if umap_name in existing_embeddings:
        st.warning(f"`{umap_name}` already exists and will be overwritten.")

    if st.button("Compute UMAP", type="primary", disabled=not umap_markers):
        with st.spinner(
            f"Computing UMAP — {len(umap_markers)} markers × {adata.n_obs:,} cells. "
            "This may take a few minutes…"
        ):
            try:
                meta = run.compute_umap(
                    markers=umap_markers,
                    embedding_name=umap_name,
                    n_neighbors=umap_n_neighbors,
                    min_dist=umap_min_dist,
                    inplace=True,
                )
                bump_adata_version()
                umap_sec = meta.get("umap_sec", 0)
                knn_sec  = meta.get("knn_sec", 0)
                total    = umap_sec + knn_sec
                st.success(
                    f"Done — `{umap_name}` saved. "
                    f"KNN {knn_sec:.1f}s · UMAP {umap_sec:.1f}s · total {total:.1f}s"
                )
                st.info("Go to **Explore** to visualise the embedding.")
                st.rerun()
            except Exception as exc:
                st.error(f"UMAP failed: {exc}")

# ── Leiden clustering ──────────────────────────────────────────────────────────
with tab_leiden:
    st.subheader("Leiden clustering")
    st.write(
        "Clusters cells using the KNN graph produced by UMAP. "
        "Requires a UMAP embedding to have been computed first."
    )

    adata_now = run.read_adata()
    embed_keys = [k for k in adata_now.obsm if k.startswith("X_") and not any(
        s in k for s in ("_knn_", "_flowsom", "_node")
    )]

    if not embed_keys:
        st.info("No embeddings found. Compute UMAP first (UMAP tab).")
    else:
        col_a, col_b = st.columns(2)
        leiden_embedding = col_a.selectbox(
            "Embedding",
            embed_keys,
            help="KNN graph from this embedding is used for clustering.",
        )
        leiden_resolution = col_b.slider(
            "Resolution", min_value=0.1, max_value=5.0, value=1.0, step=0.1,
            help="Higher resolution = more, smaller clusters.",
        )
        leiden_key_input = st.text_input(
            "Cluster key",
            value="",
            placeholder=f"auto → {leiden_embedding}_leiden",
            help="Column name written to adata.obs. Leave blank to use the default.",
        )
        leiden_key = leiden_key_input.strip() or None

        if st.button("Run Leiden", type="primary"):
            with st.spinner("Running Leiden clustering…"):
                try:
                    meta = run.cluster_leiden(
                        embedding_name=leiden_embedding,
                        cluster_key=leiden_key,
                        resolution=leiden_resolution,
                        inplace=True,
                    )
                    bump_adata_version()
                    final_key  = meta.get("cluster_key", leiden_key or f"{leiden_embedding}_leiden")
                    n_clusters = meta.get("n_clusters", "?")
                    leiden_sec = meta.get("leiden_sec", 0)
                    st.success(
                        f"Done — `{final_key}` with **{n_clusters}** clusters "
                        f"({leiden_sec:.1f}s)."
                    )
                    st.info("Go to **Explore** to colour the UMAP by cluster, or **Analysis** for composition plots.")
                    st.rerun()
                except Exception as exc:
                    st.error(f"Leiden clustering failed: {exc}")

# ── FlowSOM ────────────────────────────────────────────────────────────────────
with tab_flowsom:
    st.subheader("FlowSOM clustering")

    _flowsom_ok = False
    try:
        import flowsom as _fs_check  # noqa: F401
        _flowsom_ok = True
    except ImportError:
        pass

    if not _flowsom_ok:
        st.info(
            "The `flowsom` package is not installed.\n\n"
            "```\npip install flowsom\n```"
        )
    else:
        st.write(
            "Trains a Self-Organising Map on the selected markers, "
            "then derives metaclusters via FlowSOM's built-in metaclustering. "
            "Results are stored in `adata.obs['X_flowsom']`."
        )

        fs_markers = st.multiselect(
            "Markers",
            all_markers,
            default=all_markers,
            key="fs_markers",
        )

        col_a, col_b = st.columns(2)
        fs_grid = col_a.slider(
            "SOM grid size", min_value=5, max_value=20, value=10,
            help="Grid will be grid_size × grid_size nodes.",
        )
        fs_n_meta_raw = col_b.number_input(
            "Meta-clusters", min_value=0, max_value=50, value=0,
            help="Number of FlowSOM metaclusters. Set 0 for automatic.",
        )
        fs_n_meta = int(fs_n_meta_raw) if int(fs_n_meta_raw) > 0 else None

        if st.button("Run FlowSOM", type="primary", disabled=not fs_markers):
            with st.spinner(
                f"Running FlowSOM — {len(fs_markers)} markers × {adata.n_obs:,} cells…"
            ):
                try:
                    meta = run.compute_flowsom(
                        markers=fs_markers,
                        grid_size=fs_grid,
                        n_meta_clusters=fs_n_meta,
                        inplace=True,
                    )
                    bump_adata_version()
                    n_clusters = meta.get("n_clusters", "?")
                    n_nodes    = meta.get("n_nodes", "?")
                    fit_sec    = meta.get("fit_sec", 0)
                    device     = meta.get("device_used", "cpu")
                    st.success(
                        f"Done — **{n_clusters}** metaclusters, {n_nodes} SOM nodes "
                        f"({fit_sec:.1f}s on {device})."
                    )
                    st.info("Go to **Explore** to colour the UMAP by FlowSOM cluster.")
                    st.rerun()
                except Exception as exc:
                    st.error(f"FlowSOM failed: {exc}")
