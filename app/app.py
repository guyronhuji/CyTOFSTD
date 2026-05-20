"""Home page — project folder picker + navigation."""

import sys
from pathlib import Path

import streamlit as st

sys.path.insert(0, str(Path(__file__).parent.parent))

from _shared import _load_project, render_sidebar  # noqa: E402

st.set_page_config(
    page_title="CyTOF Standard",
    page_icon="🔬",
    layout="wide",
)

render_sidebar()

st.title("🔬 CyTOF Standard")
st.write("A guided interface for CyTOF mass cytometry analysis.")

st.divider()


def _pick_folder() -> str | None:
    """Open a native OS folder-picker dialog via subprocess (macOS main-thread safe)."""
    import subprocess

    script = (
        "import tkinter as tk, sys\n"
        "from tkinter import filedialog\n"
        "root = tk.Tk(); root.withdraw()\n"
        "root.wm_attributes('-topmost', 1)\n"
        "path = filedialog.askdirectory(title='Select CyTOF Standard project folder')\n"
        "root.destroy()\n"
        "print(path)\n"
    )
    try:
        result = subprocess.run(
            [sys.executable, "-c", script],
            capture_output=True,
            text=True,
            timeout=120,
        )
        return result.stdout.strip() or None
    except Exception:
        return None


col_browse, col_text = st.columns([1, 4])

with col_browse:
    st.write("")  # vertical alignment spacer
    if st.button("📁 Browse…", use_container_width=True):
        chosen = _pick_folder()
        if chosen:
            st.session_state["project_path_input"] = chosen

with col_text:
    path_input = st.text_input(
        "Project folder",
        value=st.session_state.get("project_path_input", st.session_state.get("project_path", "")),
        placeholder="/path/to/my_cytof_project",
        help="Path to the CyTOF Standard project directory. Use Browse… or type directly.",
        key="project_path_input",
        label_visibility="collapsed",
    )

st.write("")

if st.button("Load project", type="primary"):
    path_val = (path_input or "").strip()
    if not path_val:
        st.error("Please select or type a project folder path.")
    elif not Path(path_val).exists():
        st.error(f"Path does not exist: `{path_val}`")
    else:
        try:
            proj = _load_project(path_val)
            st.session_state["project_path"] = path_val
            st.session_state.pop("run_id", None)
            st.session_state["adata_version"] = 0
            st.rerun()
        except Exception as exc:
            st.error(f"Could not load project: {exc}")

# -----------------------------------------------------------------------
# Loaded project info + navigation
# -----------------------------------------------------------------------
if st.session_state.get("project_path"):
    path = st.session_state["project_path"]
    try:
        proj = _load_project(path)
        runs_df = proj.list_runs()
        n_runs = len(runs_df)

        st.success(f"Project loaded: **{Path(path).name}** — {n_runs} run{'s' if n_runs != 1 else ''}")

        st.markdown("### Navigation")
        cols = st.columns(4)
        cols[0].page_link("pages/1_Overview.py",  label="Overview",  icon="📋")
        cols[1].page_link("pages/3_QC.py",         label="QC",        icon="🔍")
        cols[2].page_link("pages/4_Explore.py",    label="Explore",   icon="🗺️")
        cols[3].page_link("pages/5_Analysis.py",   label="Analysis",  icon="📊")

    except Exception as exc:
        st.error(f"Error reading project: {exc}")
