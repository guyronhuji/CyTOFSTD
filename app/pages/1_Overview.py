"""Run list and status overview."""

import sys
from pathlib import Path

import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).parent.parent))

from _shared import get_project, render_sidebar  # noqa: E402

st.set_page_config(page_title="Overview | CyTOF Standard", layout="wide")
render_sidebar()

st.title("📋 Project Overview")

proj = get_project()
runs_df = proj.list_runs()

if runs_df.empty:
    st.info("No runs registered in this project yet.")
    st.stop()

# Status badge column
status_map = {
    "ingested":           "🟢 ingested",
    "registered":         "🟡 registered",
    "failed_ingestion":   "🔴 failed",
}
display_df = runs_df.copy()
display_df["status"] = display_df["status"].map(lambda s: status_map.get(s, s))

show_cols = [c for c in ["run_name", "status", "panel_id", "acquisition_date", "instrument", "operator", "created_at"] if c in display_df.columns]
st.dataframe(display_df[show_cols], use_container_width=True, hide_index=True)

st.divider()
st.subheader("Run details")

run_names = runs_df["run_name"].tolist()
run_ids = runs_df["run_id"].tolist()

selected_name = st.selectbox("Select a run to inspect", run_names)
selected_idx = run_names.index(selected_name)
selected_id = run_ids[selected_idx]
selected_row = runs_df.iloc[selected_idx]

st.session_state["run_id"] = selected_id

col1, col2, col3 = st.columns(3)
col1.metric("Status", status_map.get(selected_row.get("status", ""), selected_row.get("status", "")))
col2.metric("Panel", selected_row.get("panel_id") or "—")
col3.metric("Instrument", selected_row.get("instrument") or "—")

if selected_row.get("acquisition_date"):
    st.caption(f"Acquisition date: {selected_row['acquisition_date']}")
if selected_row.get("operator"):
    st.caption(f"Operator: {selected_row['operator']}")

st.caption(f"Run ID: `{selected_id}`")
st.caption(f"Path: `{selected_row.get('path', '')}`")
