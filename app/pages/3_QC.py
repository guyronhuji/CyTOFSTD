"""QC page — marker histograms and QC gating."""

import sys
from pathlib import Path

import streamlit as st

sys.path.insert(0, str(Path(__file__).parent.parent))

from _shared import bump_adata_version, get_run, render_sidebar  # noqa: E402

st.set_page_config(page_title="QC | CyTOF Standard", layout="wide")
render_sidebar()

st.title("🔍 Quality Control")

run = get_run()

if run.status != "ingested":
    st.warning(f"This run has not been ingested yet (status: **{run.status}**).")
    st.stop()

adata = run.read_adata()

# -----------------------------------------------------------------------
# Marker histograms
# -----------------------------------------------------------------------
st.subheader("Marker histograms")

all_markers = run.markers_all or list(adata.var_names)

selected_markers = st.multiselect(
    "Choose markers to plot",
    all_markers,
    default=all_markers[:min(6, len(all_markers))],
    key="qc_hist_markers",
)

if selected_markers:
    with st.spinner("Plotting histograms…"):
        try:
            fig, axes = run.plot_marker_histograms(selected_markers)
            st.pyplot(fig)
        except Exception as exc:
            st.error(f"Could not plot histograms: {exc}")
else:
    st.info("Select at least one marker above.")

st.divider()

# -----------------------------------------------------------------------
# QC gating
# -----------------------------------------------------------------------
st.subheader("Apply QC gates")
st.write(
    "Set lower and/or upper bounds for each marker. "
    "Cells outside any gate will be removed. "
    "Leave both sliders at their extremes to skip a marker."
)

n_before = adata.n_obs

gate_markers = st.multiselect(
    "Markers to gate",
    all_markers,
    default=[m for m in all_markers if "DNA" in m or "dna" in m.lower()][:2],
    key="qc_gate_markers",
)

gates: dict = {}
if gate_markers:
    # Read min/max values for each marker to set slider ranges
    for marker in gate_markers:
        if marker not in adata.var_names:
            continue

        marker_idx = list(adata.var_names).index(marker)
        vals = adata.X[:, marker_idx]
        try:
            import numpy as np
            arr = np.asarray(vals).flatten()
            v_min = float(arr.min())
            v_max = float(arr.max())
        except Exception:
            v_min, v_max = 0.0, 10.0

        st.write(f"**{marker}**")
        col_lo, col_hi = st.columns(2)
        lo = col_lo.number_input(
            f"Lower bound ({marker})",
            min_value=v_min,
            max_value=v_max,
            value=v_min,
            key=f"gate_lo_{marker}",
        )
        hi = col_hi.number_input(
            f"Upper bound ({marker})",
            min_value=v_min,
            max_value=v_max,
            value=v_max,
            key=f"gate_hi_{marker}",
        )
        # Only add gate if not at the extremes (user actually set it)
        if lo > v_min or hi < v_max:
            gates[marker] = {"lower": lo, "upper": hi}

if gate_markers:
    if gates:
        st.info(f"Active gates: {list(gates.keys())}")
    else:
        st.caption("No active gates — move sliders to set bounds.")

if st.button("Apply QC Gate", type="primary", disabled=not gates):
    with st.spinner("Applying QC gate…"):
        try:
            run.qc_gate(gates, inplace=True)
            bump_adata_version()
            adata_after = run.read_adata()
            n_after = adata_after.n_obs
            n_removed = n_before - n_after
            st.success(
                f"Gate applied. {n_before:,} → {n_after:,} cells "
                f"({n_removed:,} removed, {100 * n_removed / n_before:.1f}%)."
            )
            st.rerun()
        except Exception as exc:
            st.error(f"QC gate failed: {exc}")

# Cell count summary
st.divider()
col1, col2 = st.columns(2)
col1.metric("Cells in current data", f"{n_before:,}")
if run.ingestion_summary is not None:
    try:
        summ = run.ingestion_summary()
        n_events = int(summ["n_events_file"].sum()) if "n_events_file" in summ.columns else None
        if n_events:
            col2.metric("Events in raw files", f"{n_events:,}")
    except Exception:
        pass
