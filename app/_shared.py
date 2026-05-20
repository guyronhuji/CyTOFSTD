"""Shared helpers: project/run loading, sidebar, session state."""

from __future__ import annotations

import os
import sys
from pathlib import Path

import streamlit as st

# Allow importing cytofstandard from the repo root when running locally
_repo_root = Path(__file__).parent.parent
if str(_repo_root) not in sys.path:
    sys.path.insert(0, str(_repo_root))

from cytofstandard import Project  # noqa: E402


# ---------------------------------------------------------------------------
# Caching
# ---------------------------------------------------------------------------

@st.cache_resource(show_spinner="Loading project…")
def _load_project(path: str) -> Project:
    return Project.load(path)


# ---------------------------------------------------------------------------
# State accessors
# ---------------------------------------------------------------------------

def get_project() -> Project:
    path = st.session_state.get("project_path", "").strip()
    if not path:
        st.error("No project loaded. Go to the Home page and enter a project path.")
        st.stop()
    try:
        return _load_project(path)
    except Exception as exc:
        st.error(f"Could not load project: {exc}")
        st.stop()


def get_run():
    proj = get_project()
    run_id = st.session_state.get("run_id")
    if not run_id:
        st.warning("Select a run from the sidebar.")
        st.stop()
    # Run paths in the registry may be relative to the project's parent directory
    # (recorded at creation time). chdir there so Path resolution works regardless
    # of where `streamlit run` was launched from.
    prev_cwd = os.getcwd()
    try:
        os.chdir(proj.path.parent)
        return proj.get_run(run_id, validate=False)
    except Exception as exc:
        st.error(f"Could not load run: {exc}")
        st.stop()
    finally:
        os.chdir(prev_cwd)


def bump_adata_version():
    """Increment cache-busting counter after any write that modifies adata."""
    st.session_state["adata_version"] = st.session_state.get("adata_version", 0) + 1


# ---------------------------------------------------------------------------
# Sidebar
# ---------------------------------------------------------------------------

def render_sidebar():
    """Render project + run selector in the sidebar. Call from every page."""
    with st.sidebar:
        st.header("CyTOF Standard")

        path = st.session_state.get("project_path", "")
        if not path:
            st.info("Enter a project path on the Home page.")
            return

        st.caption(f"Project: `{Path(path).name}`")
        st.divider()

        try:
            proj = _load_project(path)
            runs_df = proj.list_runs()
        except Exception as exc:
            st.error(f"Cannot list runs: {exc}")
            return

        if runs_df.empty:
            st.info("No runs registered yet.")
            return

        names = runs_df["run_name"].tolist()
        ids = runs_df["run_id"].tolist()

        current_id = st.session_state.get("run_id")
        default_idx = ids.index(current_id) if current_id in ids else 0

        selected_idx = st.selectbox(
            "Run",
            range(len(names)),
            index=default_idx,
            format_func=lambda i: names[i],
            key="_sidebar_run_select",
        )
        st.session_state["run_id"] = ids[selected_idx]

        row = runs_df.iloc[selected_idx]
        status = row.get("status", "unknown")
        badge = {"ingested": "🟢", "registered": "🟡", "failed_ingestion": "🔴"}.get(status, "⚪")
        st.caption(f"{badge} {status}")

        if row.get("acquisition_date"):
            st.caption(f"Date: {row['acquisition_date']}")
        if row.get("operator"):
            st.caption(f"Operator: {row['operator']}")
