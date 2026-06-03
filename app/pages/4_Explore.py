"""UMAP explorer — interactive scatter with Plotly, cluster annotation."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).parent.parent))

from _shared import bump_adata_version, get_run, page_header, render_sidebar  # noqa: E402

st.set_page_config(page_title="Explore | CyTOF Standard", layout="wide")
render_sidebar()
page_header("UMAP Explorer", subtitle="Interactive embedding coloured by marker or cluster", icon="🗺️")

run = get_run()

if run.status != "ingested":
    st.warning(f"This run has not been ingested yet (status: **{run.status}**).")
    st.stop()

adata = run.read_adata()

embedding_keys = [k for k in adata.obsm if k.startswith("X_")]
if not embedding_keys:
    st.info(
        "No UMAP or other embeddings found. "
        "Compute UMAP first via the Python API (`run.compute_umap()`)."
    )
    st.stop()

# ── Controls + plot ───────────────────────────────────────────────────────────
ctrl_col, _, plot_col = st.columns([1, 0.05, 3])

with ctrl_col:
    st.subheader("Controls")

    umap_key = st.selectbox("Embedding", embedding_keys, index=0)

    obs_cols     = [c for c in adata.obs.columns if c not in {"cell_uuid", "source_file", "source_file_hash"}]
    marker_names = list(adata.var_names)
    color_opts   = obs_cols + marker_names

    color_by = st.selectbox("Colour by", color_opts, index=0)

    subsample_n = st.slider(
        "Subsample cells",
        min_value=5_000,
        max_value=min(adata.n_obs, 200_000),
        value=min(adata.n_obs, 80_000),
        step=5_000,
        help="Subsample for browser performance. 50–100k is usually fast.",
    )

    point_size = st.slider("Point size", 1, 8, 2)

# ── Build plot data ───────────────────────────────────────────────────────────
rng = np.random.default_rng(42)
if adata.n_obs > subsample_n:
    idx = rng.choice(adata.n_obs, subsample_n, replace=False)
    idx.sort()
else:
    idx = np.arange(adata.n_obs)

coords = adata.obsm[umap_key][idx]

if color_by in adata.obs.columns:
    color_vals    = adata.obs[color_by].iloc[idx].astype(str).values
    is_categorical = True
elif color_by in adata.var_names:
    marker_idx    = list(adata.var_names).index(color_by)
    raw           = adata.X[idx, marker_idx]
    color_vals    = np.asarray(raw).flatten().astype(float)
    is_categorical = False
else:
    color_vals    = ["unknown"] * len(idx)
    is_categorical = True

# ── Plotly scatter ────────────────────────────────────────────────────────────
import plotly.express as px

with plot_col:
    plot_df      = pd.DataFrame({"UMAP 1": coords[:, 0], "UMAP 2": coords[:, 1], color_by: color_vals})
    color_scale  = "Viridis" if not is_categorical else None
    color_seq    = px.colors.qualitative.Set3 if is_categorical else None

    fig = px.scatter(
        plot_df,
        x="UMAP 1", y="UMAP 2",
        color=color_by,
        color_continuous_scale=color_scale,
        color_discrete_sequence=color_seq,
        render_mode="webgl",
        opacity=0.7,
    )
    fig.update_traces(marker_size=point_size)
    fig.update_layout(
        height=600,
        margin=dict(l=0, r=0, t=30, b=0),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(22,27,34,0.8)",
        font_family="IBM Plex Mono",
    )
    st.plotly_chart(fig, use_container_width=True)

    if adata.n_obs > subsample_n:
        st.caption(f"Showing {subsample_n:,} of {adata.n_obs:,} cells (random subsample).")

# ── Cluster annotation ────────────────────────────────────────────────────────
st.divider()
st.subheader("Annotate clusters")

cluster_cols = [
    c for c in adata.obs.columns
    if "leiden" in c.lower() or "cluster" in c.lower() or "label" in c.lower()
]

if not cluster_cols:
    st.info(
        "No cluster columns found. "
        "Run Leiden clustering first via the Python API (`run.cluster_leiden()`)."
    )
else:
    cluster_key = st.selectbox("Cluster column", cluster_cols)
    clusters    = sorted(adata.obs[cluster_key].astype(str).unique())

    st.caption(f"{len(clusters)} clusters in **{cluster_key}**")

    with st.expander("Set cluster annotations", expanded=False):
        annotation_map: dict[str, str] = {}
        ncols = 3
        rows  = [clusters[i : i + ncols] for i in range(0, len(clusters), ncols)]
        for row_items in rows:
            cols = st.columns(ncols)
            for col, cl in zip(cols, row_items):
                val = col.text_input(f"Cluster {cl}", key=f"ann_{cluster_key}_{cl}", placeholder="e.g. T cell")
                if val.strip():
                    annotation_map[cl] = val.strip()

        if st.button("Apply annotations", type="primary", disabled=not annotation_map):
            with st.spinner("Annotating clusters…"):
                try:
                    run.annotate_clusters(cluster_key, annotation_map)
                    bump_adata_version()
                    st.success(f"Annotations applied to {len(annotation_map)} clusters.")
                    st.rerun()
                except Exception as exc:
                    st.error(f"Annotation failed: {exc}")
