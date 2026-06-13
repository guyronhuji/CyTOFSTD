"""Standalone cell-cycle gating app — launched by run.open_cell_cycle_app().

Reads CYTOF_PROJECT_PATH and CYTOF_RUN_ID from the environment so the
notebook can spawn it pre-loaded with the correct run. No project/run
picker, no multi-page navigation — just the gating UI.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import streamlit as st

sys.path.insert(0, str(Path(__file__).parent.parent))

from cytofstandard import Project
from cytofstandard.cell_cycle import (
    CELL_CYCLE_MARKER_ALIASES,
    DEFAULT_QUANTILE_THRESHOLDS,
    DEFAULT_THRESHOLD_METHODS,
    PHASE_COLORS,
    PHASE_ORDER,
    REQUIRED_ROLES,
    OPTIONAL_ROLES,
    auto_detect_cell_cycle_markers,
    assign_cell_cycle_phase,
    calculate_thresholds,
    extract_marker_dataframe,
    gate_cell_cycle,
    summarize_cell_cycle,
    plot_cell_cycle_phase_fractions,
)

# ── Bootstrap ──────────────────────────────────────────────────────────────────

_proj_path = os.environ.get("CYTOF_PROJECT_PATH", "")
_run_id    = os.environ.get("CYTOF_RUN_ID", "")

st.set_page_config(page_title="Cell-Cycle Gating", layout="wide")

if not _proj_path or not _run_id:
    st.error(
        "Missing environment variables. Launch this app via "
        "`run.open_cell_cycle_app()` from your notebook."
    )
    st.stop()

@st.cache_resource(show_spinner="Loading run…")
def _load_run(proj_path: str, run_id: str):
    proj = Project.load(proj_path)
    return proj.get_run(run_id, validate=False)

try:
    run = _load_run(_proj_path, _run_id)
except Exception as exc:
    st.error(f"Could not load run '{_run_id}': {exc}")
    st.stop()

adata = run.read_adata()
all_var_names  = list(adata.var_names)
available_layers = ["X"] + list(adata.layers.keys())

# ── Header ────────────────────────────────────────────────────────────────────

st.markdown(
    f"### Cell-Cycle Gating &nbsp;·&nbsp; "
    f"`{_run_id}` &nbsp;·&nbsp; "
    f"<span style='color:#6E7681;font-size:0.85em'>{Path(_proj_path).name}</span>",
    unsafe_allow_html=True,
)

# ── Layer ─────────────────────────────────────────────────────────────────────

st.subheader("Data layer")
st.caption(
    "Use arcsinh-transformed data (typically `cc_asinh` or `X`). "
    "Avoid `zscore` — thresholds are calibrated for arcsinh-scale values."
)
layer_choice = st.selectbox("Layer", available_layers,
                             index=available_layers.index("cc_asinh")
                             if "cc_asinh" in available_layers else 0,
                             key="cc_layer")
layer_arg = None if layer_choice == "X" else layer_choice

if layer_arg == "zscore":
    st.warning("**Warning:** z-scored data is not recommended for cell-cycle gating.")

st.divider()

# ── Marker mapping ────────────────────────────────────────────────────────────

st.subheader("Marker mapping")
detected = auto_detect_cell_cycle_markers(all_var_names)

marker_map: dict[str, str] = {}
detection_status: dict[str, bool] = {}

st.caption(
    "**Required:** IdU, pH3, CyclinB1.  "
    "**Optional:** Ki67, pRb — used to distinguish Cycling G1 from G0/quiescent."
)
cols_detect = st.columns(5)
for i, role in enumerate(["IdU", "pH3", "CyclinB1", "Ki67", "pRb"]):
    is_required  = role in REQUIRED_ROLES
    default_val  = detected.get(role, "")
    options      = ["— not mapped —"] + all_var_names
    raw_idx      = (all_var_names.index(default_val) + 1) if default_val and default_val in all_var_names else 0
    label        = role if is_required else f"{role} (optional)"
    with cols_detect[i]:
        chosen = st.selectbox(label, options, index=raw_idx, key=f"cc_marker_{role}",
                              help=f"Aliases: {', '.join(CELL_CYCLE_MARKER_ALIASES[role][:3])}…")
        if chosen != "— not mapped —":
            marker_map[role] = chosen
            detection_status[role] = True
        else:
            detection_status[role] = False

missing_required = [r for r in REQUIRED_ROLES if not detection_status.get(r)]
if missing_required:
    st.error(f"Missing required markers: **{', '.join(missing_required)}**")
    st.stop()

with st.expander("Gating options", expanded=False):
    ambiguous_ki67_prb = st.checkbox(
        "Split Ki67+ / pRb− cells into Early_G1_or_ambiguous", value=True,
    )

st.divider()

# ── Tabs ──────────────────────────────────────────────────────────────────────

tab_auto, tab_slider = st.tabs(["Auto thresholds", "Interactive sliders"])


# ════════════════════════════════════════════════════════════════════════════════
# TAB 1 — AUTO
# ════════════════════════════════════════════════════════════════════════════════

with tab_auto:
    st.markdown(
        "**Otsu** finds the valley between negative and positive populations.  \n"
        "**Quantile** uses a fixed percentile."
    )
    active_roles = list(marker_map.keys())
    method_inputs: dict[str, str] = {}
    quant_inputs:  dict[str, float] = {}

    cols_method = st.columns(len(active_roles))
    for idx, role in enumerate(active_roles):
        default_method = DEFAULT_THRESHOLD_METHODS.get(role, "quantile")
        with cols_method[idx]:
            method_inputs[role] = st.radio(
                f"{role} method", ["otsu", "quantile"],
                index=0 if default_method == "otsu" else 1,
                key=f"cc_method_{role}",
            )

    quantile_roles = [r for r in active_roles if method_inputs[r] == "quantile"]
    if quantile_roles:
        st.caption("Quantile settings:")
        col_q1, col_q2 = st.columns(2)
        for idx, role in enumerate(quantile_roles):
            default_q = DEFAULT_QUANTILE_THRESHOLDS.get(role, 0.90)
            col = col_q1 if idx % 2 == 0 else col_q2
            with col:
                quant_inputs[role] = st.slider(
                    f"{role} quantile", 0.50, 0.999, default_q, step=0.005,
                    format="%.3f", key=f"cc_auto_q_{role}",
                )

    @st.cache_data(show_spinner=False)
    def _get_mdf_auto(run_id: str, layer_key: str | None, _mm: str):
        return extract_marker_dataframe(run.read_adata(), marker_map, layer=layer_key)

    try:
        mdf = _get_mdf_auto(run.run_id, layer_arg, str(sorted(marker_map.items())))
    except Exception as exc:
        st.error(f"Could not extract marker data: {exc}")
        st.stop()

    preview_thresholds = calculate_thresholds(mdf, marker_map,
                                              quantile_thresholds=quant_inputs,
                                              threshold_methods=method_inputs)

    thr_preview_df = pd.DataFrame([
        {"Role": role, "Method": method_inputs.get(role, "—"), "Column": col,
         "Threshold": f"{preview_thresholds[role]:.4f}",
         "% positive": f"{100.0 * (mdf[col].values > preview_thresholds[role]).mean():.1f}%"}
        for role, col in marker_map.items()
    ])
    st.dataframe(thr_preview_df, use_container_width=True, hide_index=True)

    with st.expander("Marker histograms with thresholds", expanded=True):
        fig_cols = st.columns(len(marker_map))
        for idx, (role, col) in enumerate(marker_map.items()):
            vals = pd.to_numeric(mdf[col], errors="coerce").dropna().values
            thr  = preview_thresholds[role]
            pct  = 100.0 * (vals > thr).sum() / len(vals)
            fig, ax = plt.subplots(figsize=(3.5, 2.8))
            ax.hist(vals, bins=80, color="#4a90d9", alpha=0.75, edgecolor="none")
            ax.axvline(thr, color="#e05151", lw=1.5, ls="--", label=f"{thr:.3f} ({pct:.1f}%)")
            ax.legend(fontsize=7, loc="upper right")
            ax.set_title(f"{role}  [{method_inputs.get(role, '?')}]", fontsize=9, fontweight="bold")
            ax.set_xlabel(col, fontsize=8)
            ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
            fig.tight_layout()
            with fig_cols[idx]:
                st.pyplot(fig, use_container_width=True)
            plt.close(fig)

    if st.button("Run Cell-Cycle Gating", type="primary", key="cc_auto_run"):
        with st.spinner("Gating cells…"):
            try:
                result = run.gate_cell_cycle(
                    marker_map=marker_map, layer=layer_arg,
                    quantile_thresholds=quant_inputs, threshold_methods=method_inputs,
                    ambiguous_ki67_prb=ambiguous_ki67_prb, inplace=True,
                )
                st.session_state["cc_result_auto"] = result
                st.success(f"Gating complete — {adata.n_obs:,} cells assigned.")
            except Exception as exc:
                st.error(f"Gating failed: {exc}")

    if "cc_result_auto" in st.session_state:
        result  = st.session_state["cc_result_auto"]
        summary = result["summary"]
        st.subheader("Results")
        metric_cols = st.columns(min(len(summary), 6))
        for i, row in summary.iterrows():
            metric_cols[i % len(metric_cols)].metric(
                row["cell_cycle_phase"].replace("_", " "),
                f"{100 * row['fraction']:.1f}%",
                help=f"{row['n_cells']:,} cells",
            )
        fig = plot_cell_cycle_phase_fractions(summary)
        st.pyplot(fig, use_container_width=True); plt.close(fig)
        st.dataframe(summary.style.format({"fraction": "{:.3f}"}),
                     use_container_width=True, hide_index=True)
        st.download_button("Download summary CSV",
                           summary.to_csv(index=False).encode(),
                           file_name="cell_cycle_summary.csv", mime="text/csv")


# ════════════════════════════════════════════════════════════════════════════════
# TAB 2 — INTERACTIVE 2D SCATTER GATING
# ════════════════════════════════════════════════════════════════════════════════

with tab_slider:
    st.markdown(
        "Drag each slider to reposition the gate line. "
        "Cells are coloured by their current phase assignment. "
        "Click **Apply Gating** when satisfied."
    )

    N_DISPLAY = 20_000

    @st.cache_data(show_spinner=False)
    def _load_mdf_2d(run_id: str, layer_key: str | None, _mm: str) -> pd.DataFrame:
        return extract_marker_dataframe(run.read_adata(), marker_map, layer=layer_key)

    @st.cache_data(show_spinner=False)
    def _init_thr_2d(run_id: str, layer_key: str | None, _mm: str) -> dict:
        return calculate_thresholds(_load_mdf_2d(run_id, layer_key, _mm), marker_map)

    try:
        _mm_key    = str(sorted(marker_map.items()))
        mdf_all_2d = _load_mdf_2d(run.run_id, layer_arg, _mm_key)
        init_thr_2d = _init_thr_2d(run.run_id, layer_arg, _mm_key)
    except Exception as exc:
        st.error(f"Could not load marker data: {exc}")
        st.stop()

    _rng2d   = np.random.default_rng(42)
    _disp_idx = np.sort(_rng2d.choice(len(mdf_all_2d), min(N_DISPLAY, len(mdf_all_2d)), replace=False))
    mdf_disp_2d = mdf_all_2d.iloc[_disp_idx]

    def _slider_range(role: str) -> tuple[float, float]:
        v = mdf_all_2d[marker_map[role]].values
        return float(np.percentile(v, 0.2)), float(np.percentile(v, 99.8))

    def _scatter_gate(ax, mdf_sub, gated_sub, x_role, y_role,
                      x_thr, y_thr, ref_thr=None):
        x_col, y_col = marker_map[x_role], marker_map[y_role]
        for phase in PHASE_ORDER:
            mask = (gated_sub["cell_cycle_phase"] == phase).values
            if not mask.any():
                continue
            ax.scatter(mdf_sub.loc[mask, x_col].values, mdf_sub.loc[mask, y_col].values,
                       c=PHASE_COLORS.get(phase, "#aaaaaa"), s=4, alpha=0.55,
                       linewidths=0, label=f"{phase.replace('_', ' ')} ({mask.sum():,})",
                       rasterized=True)
        if x_thr is not None:
            ax.axvline(x_thr, color="black", lw=1.6, ls="--", alpha=0.85)
        if y_thr is not None:
            ax.axhline(y_thr, color="black", lw=1.6, ls="--", alpha=0.85)
        if ref_thr:
            for _, rv_thr in ref_thr.items():
                if rv_thr is not None:
                    ax.axvline(rv_thr, color="#888", lw=1.0, ls=":", alpha=0.5)
        ax.set_xlabel(x_col, fontsize=9); ax.set_ylabel(y_col, fontsize=9)
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

    slider_thr_2d: dict[str, float] = {}

    # Gate 1: CyclinB1 vs IdU
    gate1_roles = [r for r in ["IdU", "CyclinB1"] if r in marker_map]
    if gate1_roles:
        st.markdown("#### Gate 1 — S-phase & G2-phase: CyclinB1 vs IdU")
        col_ctrl, col_plot = st.columns([1, 3])
        with col_ctrl:
            for role in ["IdU", "CyclinB1"]:
                if role not in marker_map:
                    continue
                lo, hi = _slider_range(role)
                v0 = float(np.clip(init_thr_2d.get(role, (lo + hi) / 2), lo, hi))
                slider_thr_2d[role] = st.slider(
                    f"{role} threshold", lo, hi, v0,
                    step=max(0.001, (hi - lo) / 500), format="%.3f", key=f"cc2d_{role.lower()}")
                n_pos = int((mdf_all_2d[marker_map[role]].values > slider_thr_2d[role]).sum())
                st.caption(f"{role}+: **{100*n_pos/len(mdf_all_2d):.1f}%** ({n_pos:,})")
        cur_thr_1 = {**init_thr_2d, **slider_thr_2d}
        gated_g1  = assign_cell_cycle_phase(mdf_disp_2d, marker_map, cur_thr_1, ambiguous_ki67_prb)
        with col_plot:
            fig, ax = plt.subplots(figsize=(6, 4.5))
            _scatter_gate(ax, mdf_disp_2d, gated_g1,
                          x_role="IdU" if "IdU" in marker_map else gate1_roles[0],
                          y_role="CyclinB1" if "CyclinB1" in marker_map else gate1_roles[-1],
                          x_thr=slider_thr_2d.get("IdU"), y_thr=slider_thr_2d.get("CyclinB1"))
            ax.set_title("CyclinB1 vs IdU  (S: right | G2: top-left)", fontsize=9)
            ax.legend(markerscale=3, fontsize=7, loc="upper left", framealpha=0.7)
            fig.tight_layout(); st.pyplot(fig, use_container_width=True); plt.close(fig)
        st.divider()

    # Gate 2: pH3 vs CyclinB1
    if "pH3" in marker_map:
        st.markdown("#### Gate 2 — M-phase: pH3 vs CyclinB1")
        col_ctrl2, col_plot2 = st.columns([1, 3])
        with col_ctrl2:
            lo, hi = _slider_range("pH3")
            v0 = float(np.clip(init_thr_2d.get("pH3", (lo + hi) / 2), lo, hi))
            slider_thr_2d["pH3"] = st.slider(
                "pH3 threshold", lo, hi, v0,
                step=max(0.001, (hi - lo) / 500), format="%.3f", key="cc2d_ph3")
            n_pos = int((mdf_all_2d[marker_map["pH3"]].values > slider_thr_2d["pH3"]).sum())
            st.caption(f"pH3+: **{100*n_pos/len(mdf_all_2d):.1f}%** ({n_pos:,})")
        cur_thr_2 = {**init_thr_2d, **slider_thr_2d}
        x_role2   = next((r for r in ["CyclinB1", "IdU"] if r in marker_map), None)
        gated_g2  = assign_cell_cycle_phase(mdf_disp_2d, marker_map, cur_thr_2, ambiguous_ki67_prb)
        with col_plot2:
            if x_role2:
                fig, ax = plt.subplots(figsize=(6, 4.5))
                _scatter_gate(ax, mdf_disp_2d, gated_g2, x_role=x_role2, y_role="pH3",
                              x_thr=None, y_thr=slider_thr_2d["pH3"],
                              ref_thr={x_role2: slider_thr_2d.get(x_role2)} if x_role2 in slider_thr_2d else None)
                ax.set_title("pH3 vs CyclinB1  (M-phase: above line)", fontsize=9)
                ax.legend(markerscale=3, fontsize=7, loc="upper left", framealpha=0.7)
                fig.tight_layout(); st.pyplot(fig, use_container_width=True); plt.close(fig)
        st.divider()

    # Gate 3: Ki67 vs pRb
    g3_roles = [r for r in ["Ki67", "pRb"] if r in marker_map]
    if g3_roles:
        x_role3 = "pRb"  if "pRb"  in marker_map else None
        y_role3 = "Ki67" if "Ki67" in marker_map else None
        st.markdown(f"#### Gate 3 — G1 vs G0: {' vs '.join(filter(None, [y_role3, x_role3]))}")
        col_ctrl3, col_plot3 = st.columns([1, 3])
        with col_ctrl3:
            for role in ["Ki67", "pRb"]:
                if role not in marker_map:
                    continue
                lo, hi = _slider_range(role)
                v0 = float(np.clip(init_thr_2d.get(role, (lo + hi) / 2), lo, hi))
                slider_thr_2d[role] = st.slider(
                    f"{role} threshold", lo, hi, v0,
                    step=max(0.001, (hi - lo) / 500), format="%.3f", key=f"cc2d_{role.lower()}")
                n_pos = int((mdf_all_2d[marker_map[role]].values > slider_thr_2d[role]).sum())
                st.caption(f"{role}+: **{100*n_pos/len(mdf_all_2d):.1f}%** ({n_pos:,})")
        cur_thr_3 = {**init_thr_2d, **slider_thr_2d}
        gated_g3  = assign_cell_cycle_phase(mdf_disp_2d, marker_map, cur_thr_3, ambiguous_ki67_prb)
        with col_plot3:
            if x_role3 and y_role3:
                fig, ax = plt.subplots(figsize=(6, 4.5))
                _scatter_gate(ax, mdf_disp_2d, gated_g3, x_role=x_role3, y_role=y_role3,
                              x_thr=slider_thr_2d.get("pRb"), y_thr=slider_thr_2d.get("Ki67"))
                ax.set_title("Ki67 vs pRb  (Cycling G1: top-right | G0: bottom-left)", fontsize=9)
                ax.legend(markerscale=3, fontsize=7, loc="upper left", framealpha=0.7)
                fig.tight_layout(); st.pyplot(fig, use_container_width=True); plt.close(fig)
        st.divider()

    # Live summary
    final_thr_2d = {**init_thr_2d, **slider_thr_2d}
    gated_final  = assign_cell_cycle_phase(mdf_disp_2d, marker_map, final_thr_2d, ambiguous_ki67_prb)
    summary_2d   = summarize_cell_cycle(gated_final)

    st.markdown("#### Live preview")
    st.caption(f"Showing {len(mdf_disp_2d):,} of {len(mdf_all_2d):,} cells.")
    prev_cols = st.columns(min(len(summary_2d), 6))
    for i, row in summary_2d.iterrows():
        prev_cols[i % len(prev_cols)].metric(
            row["cell_cycle_phase"].replace("_", " "), f"{100 * row['fraction']:.1f}%",
            help=f"{row['n_cells']:,} of {len(mdf_disp_2d):,} display cells")
    fig_sum = plot_cell_cycle_phase_fractions(summary_2d, figsize=(7, 3))
    st.pyplot(fig_sum, use_container_width=True); plt.close(fig_sum)

    if st.button("Apply Gating", type="primary", key="cc_2d_run"):
        with st.spinner("Applying gating to all cells…"):
            try:
                result_2d = run.gate_cell_cycle(
                    marker_map=marker_map, layer=layer_arg,
                    thresholds=final_thr_2d,
                    ambiguous_ki67_prb=ambiguous_ki67_prb, inplace=True,
                )
                st.session_state["cc_result_2d"] = result_2d
                st.success(f"Gating applied — {adata.n_obs:,} cells assigned.")
            except Exception as exc:
                st.error(f"Gating failed: {exc}")

    if "cc_result_2d" in st.session_state:
        res2d = st.session_state["cc_result_2d"]
        sum2d = res2d["summary"]
        st.subheader("Saved results (all cells)")
        st.dataframe(sum2d.style.format({"fraction": "{:.3f}"}),
                     use_container_width=True, hide_index=True)
        st.download_button("Download summary CSV",
                           sum2d.to_csv(index=False).encode(),
                           file_name="cell_cycle_summary.csv", mime="text/csv",
                           key="cc_2d_csv")
