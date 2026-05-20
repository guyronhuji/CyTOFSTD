"""Analysis — differential abundance, cluster composition, marker heatmap."""

import io
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import streamlit as st

sys.path.insert(0, str(Path(__file__).parent.parent))

from _shared import get_run, render_sidebar  # noqa: E402

st.set_page_config(page_title="Analysis | CyTOF Standard", layout="wide")
render_sidebar()

st.title("📊 Analysis")

run = get_run()

if run.status != "ingested":
    st.warning(f"This run has not been ingested yet (status: **{run.status}**).")
    st.stop()

adata = run.read_adata()

# Reusable column sets
all_obs_cols = [c for c in adata.obs.columns if c not in {"cell_uuid", "source_file", "source_file_hash", "event_index"}]
cluster_cols = [
    c for c in adata.obs.columns
    if "leiden" in c.lower() or "cluster" in c.lower() or "annotated" in c.lower()
]
group_cols = [c for c in all_obs_cols if c not in cluster_cols]

# -----------------------------------------------------------------------
tab1, tab2, tab3 = st.tabs(["Differential Abundance", "Cluster Composition", "Marker Heatmap"])

# -----------------------------------------------------------------------
# Tab 1 — Differential abundance
# -----------------------------------------------------------------------
with tab1:
    st.subheader("Differential cluster abundance")
    st.write("Tests whether cluster proportions differ significantly between groups.")

    if not cluster_cols:
        st.info("No cluster columns found. Run clustering first via the Python API.")
    elif not group_cols:
        st.info("No group columns found in obs.")
    else:
        col_a, col_b = st.columns(2)
        da_cluster_key = col_a.selectbox("Cluster column", cluster_cols, key="da_ck")
        da_groupby     = col_b.selectbox("Group by", group_cols, key="da_gb")

        da_method   = st.radio("Test", ["fisher", "chi2"], horizontal=True)
        da_multitest = st.radio("Multiple testing", ["bh", "bonferroni", "none"], horizontal=True)
        da_multitest_val = None if da_multitest == "none" else da_multitest

        if st.button("Run differential abundance", type="primary"):
            with st.spinner("Running…"):
                try:
                    result_df, (fig, ax) = run.differential_abundance(
                        da_cluster_key,
                        da_groupby,
                        method=da_method,
                        multitest=da_multitest_val,
                        plot=True,
                    )
                    st.pyplot(fig)
                    buf = io.BytesIO()
                    fig.savefig(buf, bbox_inches="tight", dpi=300)
                    buf.seek(0)
                    st.download_button(
                        "Download plot",
                        buf,
                        "differential_abundance.png",
                        "image/png",
                        key="da_download",
                    )
                    st.dataframe(result_df, use_container_width=True, hide_index=True)
                except Exception as exc:
                    st.error(f"Differential abundance failed: {exc}")

# -----------------------------------------------------------------------
# Tab 2 — Cluster composition
# -----------------------------------------------------------------------
with tab2:
    st.subheader("Cluster composition")
    st.write("Stacked bar chart showing how each cluster is composed by group, or vice versa.")

    if not cluster_cols:
        st.info("No cluster columns found. Run clustering first via the Python API.")
    elif not group_cols:
        st.info("No group columns found in obs.")
    else:
        col_a, col_b = st.columns(2)
        cc_cluster_key = col_a.selectbox("Cluster column", cluster_cols, key="cc_ck")
        cc_groupby     = col_b.selectbox("Group by", group_cols, key="cc_gb")

        normalize = st.radio(
            "Normalize",
            ["cluster", "group"],
            horizontal=True,
            help=(
                "**cluster**: each bar = one cluster, filled by group fraction.  "
                "**group**: each bar = one group, filled by cluster fraction."
            ),
        )

        if st.button("Plot composition", type="primary"):
            with st.spinner("Plotting…"):
                try:
                    fig, ax = run.plot_cluster_composition(cc_cluster_key, cc_groupby, normalize=normalize)
                    st.pyplot(fig)
                    buf = io.BytesIO()
                    fig.savefig(buf, bbox_inches="tight", dpi=300)
                    buf.seek(0)
                    st.download_button(
                        "Download plot",
                        buf,
                        "cluster_composition.png",
                        "image/png",
                        key="cc_download",
                    )
                except Exception as exc:
                    st.error(f"Cluster composition failed: {exc}")

# -----------------------------------------------------------------------
# Tab 3 — Marker heatmap
# -----------------------------------------------------------------------
with tab3:
    st.subheader("Marker heatmap")
    st.write("Mean (or median) marker expression per group, shown as a heatmap.")

    all_markers = run.markers_all or list(adata.var_names)
    if not all_markers:
        st.info("No markers found.")
    elif not group_cols:
        st.info("No group columns found in obs.")
    else:
        hm_markers = st.multiselect(
            "Markers",
            all_markers,
            default=all_markers[:min(20, len(all_markers))],
        )
        hm_groupby = st.selectbox("Group by", group_cols + cluster_cols, key="hm_gb")
        hm_agg = st.radio("Aggregation", ["mean", "median"], horizontal=True)
        hm_scale = st.radio("Standard scale", ["none", "row", "column"], horizontal=True,
                            help="Apply z-scoring across rows (markers) or columns (groups).")
        hm_scale_val = None if hm_scale == "none" else hm_scale

        if st.button("Plot heatmap", type="primary", disabled=not hm_markers):
            with st.spinner("Plotting…"):
                try:
                    fig, ax = run.plot_heatmap(
                        hm_markers,
                        groupby=hm_groupby,
                        agg=hm_agg,
                        standard_scale=hm_scale_val,
                    )
                    st.pyplot(fig)
                    buf = io.BytesIO()
                    fig.savefig(buf, bbox_inches="tight", dpi=300)
                    buf.seek(0)
                    st.download_button(
                        "Download plot",
                        buf,
                        "marker_heatmap.png",
                        "image/png",
                        key="hm_download",
                    )
                except Exception as exc:
                    st.error(f"Heatmap failed: {exc}")
