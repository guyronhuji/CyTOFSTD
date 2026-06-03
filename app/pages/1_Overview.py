"""Run list and status overview."""

import sys
from pathlib import Path

import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).parent.parent))

from _shared import get_project, page_header, render_sidebar, status_badge  # noqa: E402

st.set_page_config(page_title="Overview | CyTOF Standard", layout="wide")
render_sidebar()
page_header("Project Overview", subtitle="All registered runs and their status", icon="📋")

proj = get_project()
runs_df = proj.list_runs()

if runs_df.empty:
    st.info("No runs registered in this project yet.")
    st.stop()

# Summary metrics
n_total      = len(runs_df)
n_ingested   = int((runs_df["status"] == "ingested").sum())
n_registered = int((runs_df["status"] == "registered").sum())
n_failed     = int((runs_df["status"] == "failed_ingestion").sum())

c1, c2, c3, c4 = st.columns(4)
c1.metric("Total runs",  n_total)
c2.metric("Ingested",    n_ingested)
c3.metric("Registered",  n_registered)
c4.metric("Failed",      n_failed)

st.divider()

# Runs table
show_cols = [c for c in ["run_name", "status", "panel_id", "acquisition_date", "instrument", "operator", "created_at"]
             if c in runs_df.columns]
display_df = runs_df[show_cols].copy()
status_labels = {
    "ingested":         "ingested",
    "registered":       "registered",
    "failed_ingestion": "failed",
}
display_df["status"] = display_df["status"].map(lambda s: status_labels.get(s, s))

st.dataframe(display_df, use_container_width=True, hide_index=True)

st.divider()
st.subheader("Run details")

run_names = runs_df["run_name"].tolist()
run_ids   = runs_df["run_id"].tolist()

selected_name = st.selectbox("Select a run to inspect", run_names)
selected_idx  = run_names.index(selected_name)
selected_id   = run_ids[selected_idx]
selected_row  = runs_df.iloc[selected_idx]

st.session_state["run_id"] = selected_id

d1, d2, d3 = st.columns(3)
d1.metric("Status",     status_labels.get(selected_row.get("status", ""), selected_row.get("status", "")))
d2.metric("Panel",      selected_row.get("panel_id") or "—")
d3.metric("Instrument", selected_row.get("instrument") or "—")

status_val = selected_row.get("status", "unknown")
st.markdown(status_badge(status_val), unsafe_allow_html=True)
st.markdown("")

if selected_row.get("acquisition_date"):
    st.caption(f"Acquisition date: {selected_row['acquisition_date']}")
if selected_row.get("operator"):
    st.caption(f"Operator: {selected_row['operator']}")
st.caption(f"Run ID: `{selected_id}`")
st.caption(f"Path: `{selected_row.get('path', '')}`")
