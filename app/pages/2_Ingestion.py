"""Ingestion summary — files, samples, and cell counts."""

import sys
from pathlib import Path

import streamlit as st

sys.path.insert(0, str(Path(__file__).parent.parent))

from _shared import get_run, page_header, render_sidebar  # noqa: E402

st.set_page_config(page_title="Ingestion | CyTOF Standard", layout="wide")
render_sidebar()
page_header("Ingestion Summary", subtitle="Source files, samples, and cell counts", icon="📂")

run = get_run()

if run.status != "ingested":
    st.warning(f"This run has not been ingested yet (status: **{run.status}**). Ingest via the Python API first.")
    st.stop()

try:
    summary_df = run.ingestion_summary()
except Exception as exc:
    st.error(f"Could not load ingestion summary: {exc}")
    st.stop()

adata = run.read_adata()

c1, c2, c3, c4 = st.columns(4)
c1.metric("Files",                len(summary_df))
c2.metric("Cells (after QC)",     f"{adata.n_obs:,}")
c3.metric("Markers",              len(run.markers_all))
c4.metric("Samples",              summary_df["sample_id"].nunique() if "sample_id" in summary_df.columns else "—")

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
