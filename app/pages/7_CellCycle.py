"""Cell-cycle gating — auto quantile or interactive slider thresholds."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import streamlit as st

sys.path.insert(0, str(Path(__file__).parent.parent))
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from _shared import bump_adata_version, get_run, page_header, render_sidebar  # noqa: E402
from cytofstandard.cell_cycle import (  # noqa: E402
    CELL_CYCLE_MARKER_ALIASES,
    DEFAULT_QUANTILE_THRESHOLDS,
    DEFAULT_THRESHOLD_METHODS,
    PHASE_COLORS,
    PHASE_ORDER,
    REQUIRED_ROLES,
    auto_detect_cell_cycle_markers,
    assign_cell_cycle_phase,
    calculate_thresholds,
    extract_marker_dataframe,
    gate_cell_cycle,
    summarize_cell_cycle,
    plot_cell_cycle_phase_fractions,
)

st.set_page_config(page_title="Cell Cycle | CyTOF Standard", layout="wide")
render_sidebar()
page_header(
    "Cell-Cycle Gating",
    subtitle="Assign cells to G0/G1/S/G2/M phases using IdU, pH3, CyclinB1, and pRb",
    icon="🔄",
)

run = get_run()

if run.status != "ingested":
    st.warning(f"This run has not been ingested yet (status: **{run.status}**).")
    st.stop()

adata = run.read_adata()
all_var_names = list(adata.var_names)
available_layers = ["X"] + list(adata.layers.keys())

# ── Layer selection ────────────────────────────────────────────────────────────

st.subheader("Data layer")
st.caption(
    "Use arcsinh-transformed data (typically `X` or `normalized`). "
    "Avoid `zscore` — thresholds are calibrated for arcsinh-scale values."
)

layer_choice = st.selectbox(
    "Layer",
    available_layers,
    index=0,
    help="Which expression layer to use for gating.",
    key="cc_layer",
)
layer_arg = None if layer_choice == "X" else layer_choice

if layer_arg == "zscore":
    st.warning(
        "**Warning:** z-scored data is not recommended for cell-cycle gating. "
        "Thresholds set on arcsinh values will not be meaningful. "
        "Select `X` or `normalized` instead."
    )

st.divider()

# ── Marker detection / selection ───────────────────────────────────────────────

st.subheader("Marker mapping")

detected = auto_detect_cell_cycle_markers(all_var_names)

marker_map: dict[str, str] = {}
detection_status: dict[str, bool] = {}

st.caption(
    "**Required:** IdU, pH3, CyclinB1 — needed for S/M/G2 assignment.  "
    "**Optional:** pRb — distinguishes Cycling G1 from G0/quiescent."
)
cols_detect = st.columns(4)
for i, role in enumerate(["IdU", "pH3", "CyclinB1", "pRb"]):
    is_required = role in REQUIRED_ROLES
    default_val = detected.get(role, "")
    default_idx = (
        all_var_names.index(default_val)
        if default_val and default_val in all_var_names
        else 0
    )
    options = ["— not mapped —"] + all_var_names
    raw_idx = default_idx + 1 if default_val else 0
    label = role if is_required else f"{role} (optional)"

    with cols_detect[i]:
        chosen = st.selectbox(
            label,
            options,
            index=raw_idx,
            key=f"cc_marker_{role}",
            help=f"Aliases: {', '.join(CELL_CYCLE_MARKER_ALIASES.get(role, [role])[:3])}…",
        )
        if chosen != "— not mapped —":
            marker_map[role] = chosen
            detection_status[role] = True
        else:
            detection_status[role] = False

missing_required = [r for r in REQUIRED_ROLES if not detection_status.get(r)]

if missing_required:
    st.error(
        f"Missing **required** marker mappings: **{', '.join(missing_required)}**. "
        "Select a column for each required role above."
    )
    st.stop()

mapped_str = "  ·  ".join(f"**{role}** → `{col}`" for role, col in marker_map.items())
missing_optional = [r for r in ["pRb"] if r not in marker_map]
opt_note = f"  (omitted optional: {', '.join(missing_optional)})" if missing_optional else ""
st.caption(f"Mapped: {mapped_str}{opt_note}")

st.divider()

# ── Tabs ───────────────────────────────────────────────────────────────────────

tab_auto, tab_slider = st.tabs(["Auto thresholds", "Interactive sliders"])


# ════════════════════════════════════════════════════════════════════════════════
# TAB 1 — AUTO MODE
# ════════════════════════════════════════════════════════════════════════════════

with tab_auto:
    st.markdown(
        "**Otsu** (default for IdU / pH3 / CyclinB1) finds the valley between the "
        "negative and S/M/G2 populations — ideal for bimodal CyTOF distributions.  \n"
        "**Quantile** (default for pRb) uses a fixed percentile."
    )

    active_roles = list(marker_map.keys())
    method_inputs: dict[str, str] = {}
    quant_inputs: dict[str, float] = {}

    cols_method = st.columns(len(active_roles))
    for idx, role in enumerate(active_roles):
        default_method = DEFAULT_THRESHOLD_METHODS.get(role, "quantile")
        with cols_method[idx]:
            method_inputs[role] = st.radio(
                f"{role} method",
                ["otsu", "quantile"],
                index=0 if default_method == "otsu" else 1,
                key=f"cc_method_{role}",
                horizontal=False,
            )

    # Show quantile slider only for roles using quantile method
    quantile_roles = [r for r in active_roles if method_inputs[r] == "quantile"]
    if quantile_roles:
        st.caption("Quantile settings (applies only to roles using the quantile method):")
        col_q1, col_q2 = st.columns(2)
        for idx, role in enumerate(quantile_roles):
            default_q = DEFAULT_QUANTILE_THRESHOLDS.get(role, 0.90)
            col = col_q1 if idx % 2 == 0 else col_q2
            with col:
                quant_inputs[role] = st.slider(
                    f"{role} quantile",
                    min_value=0.50,
                    max_value=0.999,
                    value=default_q,
                    step=0.005,
                    format="%.3f",
                    key=f"cc_auto_q_{role}",
                )

    # Compute and display preview thresholds ──────────────────────────────────
    @st.cache_data(show_spinner=False)
    def _get_marker_df_auto(run_id: str, layer_key: str | None):
        ad = run.read_adata()
        return extract_marker_dataframe(ad, marker_map, layer=layer_key)

    try:
        mdf = _get_marker_df_auto(run.run_id, layer_arg)
    except Exception as exc:
        st.error(f"Could not extract marker data: {exc}")
        st.stop()

    preview_thresholds = calculate_thresholds(
        mdf, marker_map,
        quantile_thresholds=quant_inputs,
        threshold_methods=method_inputs,
    )

    # Show threshold table
    thr_preview_df = pd.DataFrame(
        [
            {
                "Role": role,
                "Method": method_inputs.get(role, "—"),
                "Column": col,
                "Threshold": f"{preview_thresholds[role]:.4f}",
                "% positive": f"{100.0 * (mdf[col].values > preview_thresholds[role]).mean():.1f}%",
            }
            for role, col in marker_map.items()
        ]
    )
    st.dataframe(thr_preview_df, width="stretch", hide_index=True)

    # Quick histograms row ─────────────────────────────────────────────────────
    with st.expander("Marker histograms with thresholds", expanded=True):
        n_markers = len(marker_map)
        fig_cols = st.columns(n_markers)
        for idx, (role, col) in enumerate(marker_map.items()):
            vals = pd.to_numeric(mdf[col], errors="coerce").dropna().values
            thr = preview_thresholds[role]
            n_pos = int((vals > thr).sum())
            pct = 100.0 * n_pos / len(vals)

            fig, ax = plt.subplots(figsize=(3.5, 2.8))
            ax.hist(vals, bins=80, color="#4a90d9", alpha=0.75, edgecolor="none")
            ax.axvline(thr, color="#e05151", linewidth=1.5, linestyle="--",
                       label=f"{thr:.3f} ({pct:.1f}%)")
            ax.legend(fontsize=7, loc="upper right")
            ax.set_title(
                f"{role}  [{method_inputs.get(role, '?')}]",
                fontsize=9, fontweight="bold",
            )
            ax.set_xlabel(col, fontsize=8)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            fig.tight_layout()

            with fig_cols[idx]:
                st.pyplot(fig, width="stretch")
            plt.close(fig)

    # Apply gating ─────────────────────────────────────────────────────────────
    if st.button("Run Cell-Cycle Gating", type="primary", key="cc_auto_run"):
        with st.spinner("Gating cells…"):
            try:
                result = run.gate_cell_cycle(
                    marker_map=marker_map,
                    layer=layer_arg,
                    quantile_thresholds=quant_inputs,
                    threshold_methods=method_inputs,
                    inplace=True,
                )
                bump_adata_version()
                st.session_state["cc_result_auto"] = result
                st.success(
                    f"Cell-cycle gating complete — {adata.n_obs:,} cells assigned."
                )
            except Exception as exc:
                st.error(f"Gating failed: {exc}")

    # Results ──────────────────────────────────────────────────────────────────
    if "cc_result_auto" in st.session_state:
        result = st.session_state["cc_result_auto"]
        summary = result["summary"]

        st.subheader("Results")
        metric_cols = st.columns(min(len(summary), 6))
        for i, row in summary.iterrows():
            col_idx = i % len(metric_cols)
            metric_cols[col_idx].metric(
                row["cell_cycle_phase"].replace("_", " "),
                f"{100 * row['fraction']:.1f}%",
                help=f"{row['n_cells']:,} cells",
            )

        fig = plot_cell_cycle_phase_fractions(summary)
        st.pyplot(fig, width="stretch")
        plt.close(fig)

        st.dataframe(summary.style.format({"fraction": "{:.3f}"}),
                     width="stretch", hide_index=True)

        csv_bytes = summary.to_csv(index=False).encode()
        st.download_button(
            "Download summary CSV",
            csv_bytes,
            file_name="cell_cycle_summary.csv",
            mime="text/csv",
        )


# ════════════════════════════════════════════════════════════════════════════════
# TAB 2 — INTERACTIVE 2D SCATTER GATING
# ════════════════════════════════════════════════════════════════════════════════

with tab_slider:
    st.markdown(
        "Drag each threshold slider to reposition the gate line on the scatter plot. "
        "Cells are coloured by their current phase assignment. "
        "Click **Apply Gating** when satisfied."
    )

    N_DISPLAY = 20_000

    @st.cache_data(show_spinner=False)
    def _load_mdf_2d(run_id: str, layer_key: str | None, _mm: str) -> pd.DataFrame:
        ad = run.read_adata()
        return extract_marker_dataframe(ad, marker_map, layer=layer_key)

    @st.cache_data(show_spinner=False)
    def _init_thr_2d(run_id: str, layer_key: str | None, _mm: str) -> dict:
        mdf = _load_mdf_2d(run_id, layer_key, _mm)
        return calculate_thresholds(mdf, marker_map)

    try:
        _mm_key = str(sorted(marker_map.items()))
        mdf_all_2d = _load_mdf_2d(run.run_id, layer_arg, _mm_key)
        init_thr_2d = _init_thr_2d(run.run_id, layer_arg, _mm_key)
    except Exception as exc:
        st.error(f"Could not load marker data: {exc}")
        st.stop()

    # Fixed subsample for display (reproducible)
    _rng2d = np.random.default_rng(42)
    _disp_idx = np.sort(_rng2d.choice(len(mdf_all_2d), min(N_DISPLAY, len(mdf_all_2d)), replace=False))
    mdf_disp_2d = mdf_all_2d.iloc[_disp_idx]

    def _slider_range(role: str) -> tuple[float, float]:
        col = marker_map[role]
        v = mdf_all_2d[col].values
        return float(np.percentile(v, 0.2)), float(np.percentile(v, 99.8))

    def _scatter_gate(
        ax,
        mdf_sub: pd.DataFrame,
        gated_sub: pd.DataFrame,
        x_role: str,
        y_role: str,
        x_thr: float | None,
        y_thr: float | None,
    ) -> None:
        """Draw a coloured phase scatter with threshold lines. y_role is always IdU."""
        x_col = marker_map[x_role]
        y_col = marker_map[y_role]
        for phase in PHASE_ORDER:
            mask = (gated_sub["cell_cycle_phase"] == phase).values
            if not mask.any():
                continue
            ax.scatter(
                mdf_sub.loc[mask, x_col].values,
                mdf_sub.loc[mask, y_col].values,
                c=PHASE_COLORS.get(phase, "#aaaaaa"),
                s=4, alpha=0.55, linewidths=0,
                label=f"{phase.replace('_', ' ')} ({mask.sum():,})",
                rasterized=True,
            )
        if x_thr is not None:
            ax.axvline(x_thr, color="#ff7f0e", lw=1.6, ls="--", alpha=0.9)
        if y_thr is not None:
            ax.axhline(y_thr, color="#d62728", lw=1.6, ls="--", alpha=0.9)
        ax.set_xlabel(x_col, fontsize=9)
        ax.set_ylabel(y_col, fontsize=9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    slider_thr_2d: dict[str, float] = {}

    # ── Gate 1: IdU (Y) vs CyclinB1 (X)  ────────────────────────────────────
    # Defines S-phase (IdU+) and G2-phase (CyclinB1+ & IdU-)
    gate1_roles = [r for r in ["IdU", "CyclinB1"] if r in marker_map]

    if gate1_roles:
        st.markdown("#### Gate 1 — S-phase & G2-phase: IdU vs CyclinB1")
        col_ctrl, col_plot = st.columns([1, 3])

        with col_ctrl:
            if "IdU" in marker_map:
                lo, hi = _slider_range("IdU")
                step = max(0.001, (hi - lo) / 500)
                v0 = float(np.clip(init_thr_2d.get("IdU", (lo + hi) / 2), lo, hi))
                slider_thr_2d["IdU"] = st.slider(
                    "IdU threshold", lo, hi, v0, step=step, format="%.3f", key="cc2d_idu"
                )
            if "CyclinB1" in marker_map:
                lo, hi = _slider_range("CyclinB1")
                step = max(0.001, (hi - lo) / 500)
                v0 = float(np.clip(init_thr_2d.get("CyclinB1", (lo + hi) / 2), lo, hi))
                slider_thr_2d["CyclinB1"] = st.slider(
                    "CyclinB1 threshold", lo, hi, v0, step=step, format="%.3f", key="cc2d_cycb"
                )

            for role in gate1_roles:
                col_name = marker_map[role]
                n_pos = int((mdf_all_2d[col_name].values > slider_thr_2d[role]).sum())
                pct = 100.0 * n_pos / len(mdf_all_2d)
                st.caption(f"{role}+: **{pct:.1f}%** ({n_pos:,} cells)")

        cur_thr_1 = {**init_thr_2d, **slider_thr_2d}
        gated_g1 = assign_cell_cycle_phase(mdf_disp_2d, marker_map, cur_thr_1)

        with col_plot:
            fig, ax = plt.subplots(figsize=(6, 4.5))
            _scatter_gate(
                ax, mdf_disp_2d, gated_g1,
                x_role="CyclinB1", y_role="IdU",
                x_thr=slider_thr_2d.get("CyclinB1"),
                y_thr=slider_thr_2d.get("IdU"),
            )
            ax.set_title("S / G1 / G2 gate", fontsize=9, fontweight="bold")
            ax.legend(markerscale=3, fontsize=7, loc="upper left",
                      framealpha=0.7, borderpad=0.5)
            fig.tight_layout()
            st.pyplot(fig, width="stretch")
            plt.close(fig)

        st.divider()

    else:
        cur_thr_1 = init_thr_2d.copy()

    # ── Gate 2: IdU (Y) vs pH3 (X)  ─────────────────────────────────────────
    # Defines M-phase (pH3+); IdU threshold from Gate 1 shown as reference
    if "pH3" in marker_map:
        st.markdown("#### Gate 2 — M-phase: IdU vs pH3")
        col_ctrl2, col_plot2 = st.columns([1, 3])

        with col_ctrl2:
            lo, hi = _slider_range("pH3")
            step = max(0.001, (hi - lo) / 500)
            v0 = float(np.clip(init_thr_2d.get("pH3", (lo + hi) / 2), lo, hi))
            slider_thr_2d["pH3"] = st.slider(
                "pH3 threshold", lo, hi, v0, step=step, format="%.3f", key="cc2d_ph3"
            )
            n_pos = int((mdf_all_2d[marker_map["pH3"]].values > slider_thr_2d["pH3"]).sum())
            st.caption(f"pH3+: **{100*n_pos/len(mdf_all_2d):.1f}%** ({n_pos:,} cells)")

        cur_thr_2 = {**init_thr_2d, **slider_thr_2d}
        gated_g2 = assign_cell_cycle_phase(mdf_disp_2d, marker_map, cur_thr_2)

        with col_plot2:
            fig, ax = plt.subplots(figsize=(6, 4.5))
            _scatter_gate(
                ax, mdf_disp_2d, gated_g2,
                x_role="pH3", y_role="IdU",
                x_thr=slider_thr_2d["pH3"],
                y_thr=slider_thr_2d.get("IdU"),
            )
            ax.set_title("M gate (pH3)", fontsize=9, fontweight="bold")
            ax.legend(markerscale=3, fontsize=7, loc="upper left",
                      framealpha=0.7, borderpad=0.5)
            fig.tight_layout()
            st.pyplot(fig, width="stretch")
            plt.close(fig)

        st.divider()

    # ── Gate 3: IdU (Y) vs pRb (X)  ─────────────────────────────────────────
    # Defines Cycling G1 (pRb+) vs G0/quiescent (pRb-)
    if "pRb" in marker_map:
        st.markdown("#### Gate 3 — G0 vs Cycling G1: IdU vs pRb")
        col_ctrl3, col_plot3 = st.columns([1, 3])

        with col_ctrl3:
            lo, hi = _slider_range("pRb")
            step = max(0.001, (hi - lo) / 500)
            v0 = float(np.clip(init_thr_2d.get("pRb", (lo + hi) / 2), lo, hi))
            slider_thr_2d["pRb"] = st.slider(
                "pRb threshold", lo, hi, v0, step=step, format="%.3f", key="cc2d_prb"
            )
            n_pos = int((mdf_all_2d[marker_map["pRb"]].values > slider_thr_2d["pRb"]).sum())
            st.caption(f"pRb+: **{100*n_pos/len(mdf_all_2d):.1f}%** ({n_pos:,} cells)")

        cur_thr_3 = {**init_thr_2d, **slider_thr_2d}
        gated_g3 = assign_cell_cycle_phase(mdf_disp_2d, marker_map, cur_thr_3)

        with col_plot3:
            fig, ax = plt.subplots(figsize=(6, 4.5))
            _scatter_gate(
                ax, mdf_disp_2d, gated_g3,
                x_role="pRb", y_role="IdU",
                x_thr=slider_thr_2d.get("pRb"),
                y_thr=slider_thr_2d.get("IdU"),
            )
            ax.set_title("G0 gate (pRb)", fontsize=9, fontweight="bold")
            ax.legend(markerscale=3, fontsize=7, loc="upper left",
                      framealpha=0.7, borderpad=0.5)
            fig.tight_layout()
            st.pyplot(fig, width="stretch")
            plt.close(fig)

        st.divider()

    # ── Live summary ───────────────────────────────────────────────────────────
    final_thr_2d = {**init_thr_2d, **slider_thr_2d}
    gated_final_2d = assign_cell_cycle_phase(mdf_disp_2d, marker_map, final_thr_2d)
    summary_2d = summarize_cell_cycle(gated_final_2d)

    st.markdown("#### Live preview")
    st.caption(
        f"Showing {len(mdf_disp_2d):,} of {len(mdf_all_2d):,} cells. "
        "Gating will be applied to all cells when you click Apply."
    )
    prev_cols = st.columns(min(len(summary_2d), 6))
    for i, row in summary_2d.iterrows():
        prev_cols[i % len(prev_cols)].metric(
            row["cell_cycle_phase"].replace("_", " "),
            f"{100 * row['fraction']:.1f}%",
            help=f"{row['n_cells']:,} of {len(mdf_disp_2d):,} display cells",
        )

    fig_sum2d = plot_cell_cycle_phase_fractions(summary_2d, figsize=(7, 3))
    st.pyplot(fig_sum2d, width="stretch")
    plt.close(fig_sum2d)

    # ── Apply gating to full dataset ───────────────────────────────────────────
    if st.button("Apply Gating", type="primary", key="cc_2d_run"):
        with st.spinner("Applying cell-cycle gating to all cells…"):
            try:
                result_2d = run.gate_cell_cycle(
                    marker_map=marker_map,
                    layer=layer_arg,
                    thresholds=final_thr_2d,
                    inplace=True,
                )
                bump_adata_version()
                st.session_state["cc_result_2d"] = result_2d
                st.success(
                    f"Cell-cycle gating applied — {adata.n_obs:,} cells assigned."
                )
            except Exception as exc:
                st.error(f"Gating failed: {exc}")

    # ── Saved results ──────────────────────────────────────────────────────────
    if "cc_result_2d" in st.session_state:
        res2d = st.session_state["cc_result_2d"]
        sum2d = res2d["summary"]
        st.subheader("Saved results (all cells)")
        st.dataframe(sum2d.style.format({"fraction": "{:.3f}"}),
                     width="stretch", hide_index=True)
        st.download_button(
            "Download summary CSV",
            sum2d.to_csv(index=False).encode(),
            file_name="cell_cycle_summary.csv",
            mime="text/csv",
            key="cc_2d_csv",
        )
