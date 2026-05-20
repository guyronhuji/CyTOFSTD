"""Ingestion summary — files, samples, and cell counts."""

import sys
from pathlib import Path

import streamlit as st

sys.path.insert(0, str(Path(__file__).parent.parent))

from _shared import get_run, render_sidebar  # noqa: E402

st.set_page_config(page_title="Ingestion | CyTOF Standard", layout="wide")
render_sidebar()

st.title("📂 Ingestion Summary")

run = get_run()

if run.status != "ingested":
    st.warning(f"This run has not been ingested yet (status: **{run.status}**). Ingest data first via the Python API.")
    st.stop()

try:
    summary_df = run.ingestion_summary()
except Exception as exc:
    st.error(f"Could not load ingestion summary: {exc}")
    st.stop()

# Top metrics
adata = run.read_adata()
col1, col2, col3, col4 = st.columns(4)
col1.metric("Files", len(summary_df))
col2.metric("Total cells (after QC)", f"{adata.n_obs:,}")
col3.metric("Markers", len(run.markers_all))
col4.metric("Samples", summary_df["sample_id"].nunique() if "sample_id" in summary_df.columns else "—")

st.divider()

st.subheader("Per-file summary")

show_cols = [c for c in summary_df.columns if c not in {"file_hash_sha256", "path"}]
st.dataframe(summary_df[show_cols], use_container_width=True, hide_index=True)

st.divider()

st.subheader("Marker list")
tab_all, tab_intra, tab_extra = st.tabs(["All", "Intracellular", "Extracellular"])

with tab_all:
    st.write(", ".join(run.markers_all) if run.markers_all else "—")

with tab_intra:
    st.write(", ".join(run.markers_intracellular) if run.markers_intracellular else "—")

with tab_extra:
    st.write(", ".join(run.markers_extracellular) if run.markers_extracellular else "—")
