"""Shared helpers: project/run loading, sidebar, session state, design system."""

from __future__ import annotations

import os
import sys
from pathlib import Path

import streamlit as st

_repo_root = Path(__file__).parent.parent
if str(_repo_root) not in sys.path:
    sys.path.insert(0, str(_repo_root))

from cytofstandard import Project  # noqa: E402


# ── Caching ──────────────────────────────────────────────────────────────────

@st.cache_resource(show_spinner="Loading project…")
def _load_project(path: str) -> Project:
    return Project.load(path)


# ── State accessors ──────────────────────────────────────────────────────────

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


# ── Design system ────────────────────────────────────────────────────────────

_CSS = """<style>
@import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500;600&family=IBM+Plex+Sans:ital,wght@0,300;0,400;0,500;0,600;1,300&display=swap');

/* ── Typography ── */
h1 {
    font-family: 'IBM Plex Mono', monospace !important;
    font-weight: 600 !important;
    font-size: 1.45rem !important;
    letter-spacing: -0.01em !important;
}
h2 {
    font-family: 'IBM Plex Sans', sans-serif !important;
    font-weight: 500 !important;
    font-size: 1.1rem !important;
}
h3 {
    font-family: 'IBM Plex Sans', sans-serif !important;
    font-weight: 400 !important;
    font-size: 0.8rem !important;
    text-transform: uppercase !important;
    letter-spacing: 0.09em !important;
    color: #6E7681 !important;
}
p, li { font-family: 'IBM Plex Sans', sans-serif !important; }

/* ── Metric cards ── */
div[data-testid="stMetric"] {
    background: #161B22 !important;
    border: 1px solid rgba(0,212,177,0.18) !important;
    border-radius: 6px !important;
    padding: 0.85rem 1.1rem !important;
}
div[data-testid="stMetricLabel"] p {
    font-family: 'IBM Plex Mono', monospace !important;
    font-size: 0.63rem !important;
    text-transform: uppercase !important;
    letter-spacing: 0.1em !important;
    color: #6E7681 !important;
}
div[data-testid="stMetricValue"] > div {
    font-family: 'IBM Plex Mono', monospace !important;
    font-size: 1.65rem !important;
    font-weight: 600 !important;
    color: #00D4B1 !important;
}

/* ── Primary buttons ── */
button[data-testid="baseButton-primary"] {
    background: #00D4B1 !important;
    border: none !important;
    color: #0E1117 !important;
    font-family: 'IBM Plex Mono', monospace !important;
    font-weight: 600 !important;
    font-size: 0.74rem !important;
    letter-spacing: 0.08em !important;
    text-transform: uppercase !important;
    border-radius: 4px !important;
}
button[data-testid="baseButton-primary"]:hover { opacity: 0.82 !important; }

/* ── Secondary buttons ── */
button[data-testid="baseButton-secondary"] {
    border: 1px solid rgba(0,212,177,0.3) !important;
    color: #00D4B1 !important;
    background: transparent !important;
    font-family: 'IBM Plex Mono', monospace !important;
    font-size: 0.72rem !important;
    letter-spacing: 0.04em !important;
    border-radius: 4px !important;
}

/* ── Download buttons ── */
.stDownloadButton button {
    border: 1px solid rgba(139,148,158,0.25) !important;
    color: #8B949E !important;
    background: transparent !important;
    font-family: 'IBM Plex Mono', monospace !important;
    font-size: 0.7rem !important;
    letter-spacing: 0.04em !important;
    text-transform: uppercase !important;
    border-radius: 4px !important;
}

/* ── Tabs ── */
.stTabs [data-testid="stTab"] p {
    font-family: 'IBM Plex Mono', monospace !important;
    font-size: 0.7rem !important;
    font-weight: 500 !important;
    text-transform: uppercase !important;
    letter-spacing: 0.1em !important;
}

/* ── Page links (nav) ── */
.stPageLink a {
    font-family: 'IBM Plex Mono', monospace !important;
    font-size: 0.75rem !important;
    letter-spacing: 0.05em !important;
    text-transform: uppercase !important;
    border: 1px solid rgba(0,212,177,0.2) !important;
    border-radius: 4px !important;
    color: #00D4B1 !important;
}
.stPageLink a:hover { background: rgba(0,212,177,0.07) !important; }

/* ── Sidebar ── */
section[data-testid="stSidebar"] {
    border-right: 1px solid rgba(0,212,177,0.1) !important;
}

/* ── Dividers ── */
hr { border-color: rgba(0,212,177,0.1) !important; margin: 1.1rem 0 !important; }

/* ── Inline code ── */
code {
    font-family: 'IBM Plex Mono', monospace !important;
    background: rgba(0,212,177,0.08) !important;
    color: #00D4B1 !important;
    border-radius: 3px !important;
    padding: 0.1em 0.35em !important;
}

/* ── Captions ── */
.stCaption p {
    font-family: 'IBM Plex Mono', monospace !important;
    font-size: 0.67rem !important;
    color: #6E7681 !important;
    letter-spacing: 0.02em !important;
}

/* ── Form labels ── */
.stSelectbox label p, .stMultiSelect label p, .stSlider label p,
.stNumberInput label p, .stTextInput label p, .stRadio label p {
    font-family: 'IBM Plex Mono', monospace !important;
    font-size: 0.7rem !important;
    text-transform: uppercase !important;
    letter-spacing: 0.08em !important;
    color: #8B949E !important;
}
.stTextInput input {
    font-family: 'IBM Plex Mono', monospace !important;
    font-size: 0.84rem !important;
}

/* ── Expander ── */
details summary p {
    font-family: 'IBM Plex Mono', monospace !important;
    font-size: 0.76rem !important;
    text-transform: uppercase !important;
    letter-spacing: 0.07em !important;
}

/* ── Page header ── */
.cytof-header {
    padding-bottom: 0.65rem;
    border-bottom: 1px solid rgba(0,212,177,0.2);
    margin-bottom: 1.2rem;
}
.cytof-subtitle {
    font-family: 'IBM Plex Sans', sans-serif;
    font-size: 0.84rem;
    color: #6E7681;
    margin-top: 0.2rem;
    margin-bottom: 0;
}

/* ── Status badges ── */
.cytof-badge {
    display: inline-block;
    padding: 2px 9px;
    border-radius: 10px;
    font-family: 'IBM Plex Mono', monospace;
    font-size: 0.63rem;
    font-weight: 600;
    letter-spacing: 0.05em;
    text-transform: uppercase;
}
.cytof-badge-ingested   { background:rgba(63,185,80,0.12);  color:#3FB950; border:1px solid rgba(63,185,80,0.3);   }
.cytof-badge-registered { background:rgba(227,179,65,0.12); color:#E3B341; border:1px solid rgba(227,179,65,0.3);  }
.cytof-badge-failed     { background:rgba(248,81,73,0.12);  color:#F85149; border:1px solid rgba(248,81,73,0.3);   }
.cytof-badge-unknown    { background:rgba(139,148,158,0.12);color:#8B949E; border:1px solid rgba(139,148,158,0.3); }

/* ── Sidebar helpers ── */
.sb-brand { font-family:'IBM Plex Mono',monospace; font-size:0.68rem; color:#00D4B1; letter-spacing:0.16em; text-transform:uppercase; font-weight:600; margin:0; }
.sb-proj  { font-family:'IBM Plex Mono',monospace; font-size:0.68rem; color:#8B949E; margin:0.4rem 0 0.1rem; }
.sb-label { font-family:'IBM Plex Mono',monospace; font-size:0.62rem; text-transform:uppercase; letter-spacing:0.1em; color:#6E7681; margin:0.5rem 0 0.1rem; }
.sb-meta  { font-family:'IBM Plex Mono',monospace; font-size:0.65rem; color:#6E7681; margin:0.1rem 0; }
</style>"""


def inject_css() -> None:
    st.markdown(_CSS, unsafe_allow_html=True)


def page_header(title: str, subtitle: str = "", icon: str = "") -> None:
    prefix = f"{icon}&nbsp;" if icon else ""
    sub = f'<p class="cytof-subtitle">{subtitle}</p>' if subtitle else ""
    st.markdown(
        f'<div class="cytof-header"><h1>{prefix}{title}</h1>{sub}</div>',
        unsafe_allow_html=True,
    )


def status_badge(status: str) -> str:
    cls = {
        "ingested":         "cytof-badge-ingested",
        "registered":       "cytof-badge-registered",
        "failed_ingestion": "cytof-badge-failed",
    }.get(status, "cytof-badge-unknown")
    return f'<span class="cytof-badge {cls}">{status.replace("_", " ")}</span>'


# ── Sidebar ──────────────────────────────────────────────────────────────────

def render_sidebar():
    inject_css()
    with st.sidebar:
        st.markdown('<p class="sb-brand">CyTOF Standard</p>', unsafe_allow_html=True)
        st.divider()

        path = st.session_state.get("project_path", "")
        if not path:
            st.markdown('<p class="sb-meta">Open a project on the Home page.</p>', unsafe_allow_html=True)
            return

        st.markdown(f'<p class="sb-proj">📁 {Path(path).name}</p>', unsafe_allow_html=True)

        try:
            proj = _load_project(path)
            runs_df = proj.list_runs()
        except Exception as exc:
            st.error(f"Cannot list runs: {exc}")
            return

        if runs_df.empty:
            st.markdown('<p class="sb-meta">No runs registered yet.</p>', unsafe_allow_html=True)
            return

        names = runs_df["run_name"].tolist()
        ids   = runs_df["run_id"].tolist()

        current_id  = st.session_state.get("run_id")
        default_idx = ids.index(current_id) if current_id in ids else 0

        st.markdown('<p class="sb-label">Run</p>', unsafe_allow_html=True)
        selected_idx = st.selectbox(
            "Run",
            range(len(names)),
            index=default_idx,
            format_func=lambda i: names[i],
            key="_sidebar_run_select",
            label_visibility="collapsed",
        )
        st.session_state["run_id"] = ids[selected_idx]

        row    = runs_df.iloc[selected_idx]
        status = row.get("status", "unknown")
        st.markdown(status_badge(status), unsafe_allow_html=True)

        if row.get("acquisition_date"):
            st.markdown(f'<p class="sb-meta">date &nbsp;{row["acquisition_date"]}</p>', unsafe_allow_html=True)
        if row.get("operator"):
            st.markdown(f'<p class="sb-meta">by &nbsp;&nbsp;&nbsp;{row["operator"]}</p>', unsafe_allow_html=True)
