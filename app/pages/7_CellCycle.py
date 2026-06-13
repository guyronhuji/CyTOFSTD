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
    OPTIONAL_ROLES,
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
    subtitle="Assign cells to G0/G1/S/G2/M phases using IdU, pH3, CyclinB1, Ki67, and pRb",
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
    "**Optional:** Ki67, pRb — used to distinguish Cycling G1 from G0/quiescent."
)
cols_detect = st.columns(5)
for i, role in enumerate(["IdU", "pH3", "CyclinB1", "Ki67", "pRb"]):
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
            help=f"Aliases: {', '.join(CELL_CYCLE_MARKER_ALIASES[role][:3])}…",
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
missing_optional = [r for r in OPTIONAL_ROLES if r not in marker_map]
opt_note = f"  (omitted optional: {', '.join(missing_optional)})" if missing_optional else ""
st.caption(f"Mapped: {mapped_str}{opt_note}")

# ── Options ────────────────────────────────────────────────────────────────────

with st.expander("Gating options", expanded=False):
    ambiguous_ki67_prb = st.checkbox(
        "Split Ki67+ / pRb− cells into Early_G1_or_ambiguous",
        value=True,
        help=(
            "When enabled, Ki67-positive cells that are pRb-negative are assigned "
            "Early_G1_or_ambiguous instead of Unclassified."
        ),
    )

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
        "**Quantile** (default for Ki67 / pRb) uses a fixed percentile."
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
    st.dataframe(thr_preview_df, use_container_width=True, hide_index=True)

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
                st.pyplot(fig, use_container_width=True)
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
                    ambiguous_ki67_prb=ambiguous_ki67_prb,
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
        st.pyplot(fig, use_container_width=True)
        plt.close(fig)

        st.dataframe(summary.style.format({"fraction": "{:.3f}"}),
                     use_container_width=True, hide_index=True)

        csv_bytes = summary.to_csv(index=False).encode()
        st.download_button(
            "Download summary CSV",
            csv_bytes,
            file_name="cell_cycle_summary.csv",
            mime="text/csv",
        )


# ════════════════════════════════════════════════════════════════════════════════
# TAB 2 — INTERACTIVE SLIDERS
# ════════════════════════════════════════════════════════════════════════════════

with tab_slider:
    st.markdown(
        "Adjust each threshold with a slider and watch the histogram update in real time. "
        "Click **Apply Gating** when satisfied."
    )

    # Load and cache marker values for fast histogram re-renders
    @st.cache_data(show_spinner=False)
    def _get_marker_values_slider(run_id: str, layer_key: str | None, _mm: str):
        ad = run.read_adata()
        mdf = extract_marker_dataframe(ad, marker_map, layer=layer_key)
        return {role: pd.to_numeric(mdf[col], errors="coerce").dropna().values
                for role, col in marker_map.items()}

    try:
        marker_values = _get_marker_values_slider(
            run.run_id, layer_arg, str(sorted(marker_map.items()))
        )
    except Exception as exc:
        st.error(f"Could not extract marker data: {exc}")
        st.stop()

    # Compute sensible slider ranges from the data
    slider_thresholds: dict[str, float] = {}

    for role in [r for r in ["IdU", "pH3", "CyclinB1", "Ki67", "pRb"] if r in marker_map]:
        vals = marker_values[role]
        v_min = float(np.percentile(vals, 0.5))
        v_max = float(np.percentile(vals, 99.9))
        default_q = DEFAULT_QUANTILE_THRESHOLDS.get(role, 0.90)
        default_thr = float(np.quantile(vals, default_q))
        # Clamp default to [v_min, v_max]
        default_thr = float(np.clip(default_thr, v_min, v_max))

        step = max(0.001, (v_max - v_min) / 500)

        col_hist, col_ctrl = st.columns([3, 1])

        with col_ctrl:
            st.markdown(f"**{role}**")
            st.caption(f"`{marker_map[role]}`")
            thr = st.slider(
                f"Threshold ({role})",
                min_value=v_min,
                max_value=v_max,
                value=default_thr,
                step=step,
                format="%.3f",
                label_visibility="collapsed",
                key=f"cc_slider_{role}",
            )
            n_pos = int((vals > thr).sum())
            pct_pos = 100.0 * n_pos / len(vals) if len(vals) > 0 else 0.0
            st.metric(
                "Positive",
                f"{pct_pos:.1f}%",
                help=f"{n_pos:,} / {len(vals):,} cells above threshold",
            )

        with col_hist:
            # Subsample for fast rendering (max 50k cells)
            display_vals = vals
            if len(vals) > 50_000:
                rng = np.random.default_rng(42)
                display_vals = rng.choice(vals, size=50_000, replace=False)

            fig, ax = plt.subplots(figsize=(6, 2.2))
            ax.hist(display_vals, bins=80, color="#4a90d9", alpha=0.72, edgecolor="none")
            ax.axvline(thr, color="#e05151", linewidth=1.8, linestyle="--",
                       label=f"{thr:.3f}")
            ax.legend(fontsize=8, loc="upper right")
            ax.set_xlabel(marker_map[role], fontsize=8)
            ax.set_ylabel("Cells", fontsize=8)
            ax.tick_params(labelsize=7)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            fig.tight_layout(pad=0.4)
            st.pyplot(fig, use_container_width=True)
            plt.close(fig)

        slider_thresholds[role] = thr

        st.divider()

    # Live preview of gating composition ──────────────────────────────────────
    st.markdown("#### Preview gating (not yet saved)")

    try:
        ad_preview = run.read_adata()
        mdf_slider = extract_marker_dataframe(ad_preview, marker_map, layer=layer_arg)
        gated_preview = assign_cell_cycle_phase(
            mdf_slider, marker_map, slider_thresholds,
            ambiguous_ki67_prb=ambiguous_ki67_prb,
        )
        summary_preview = summarize_cell_cycle(gated_preview)

        prev_cols = st.columns(min(len(summary_preview), 6))
        for i, row in summary_preview.iterrows():
            prev_cols[i % len(prev_cols)].metric(
                row["cell_cycle_phase"].replace("_", " "),
                f"{100 * row['fraction']:.1f}%",
                help=f"{row['n_cells']:,} cells",
            )

        fig_prev = plot_cell_cycle_phase_fractions(summary_preview, figsize=(7, 3.5))
        st.pyplot(fig_prev, use_container_width=True)
        plt.close(fig_prev)
    except Exception as exc:
        st.warning(f"Preview unavailable: {exc}")

    # Apply gating ─────────────────────────────────────────────────────────────
    if st.button("Apply Gating", type="primary", key="cc_slider_run"):
        with st.spinner("Applying cell-cycle gating…"):
            try:
                result_slider = run.gate_cell_cycle(
                    marker_map=marker_map,
                    layer=layer_arg,
                    thresholds=slider_thresholds,
                    ambiguous_ki67_prb=ambiguous_ki67_prb,
                    inplace=True,
                )
                bump_adata_version()
                st.session_state["cc_result_slider"] = result_slider
                st.success(
                    f"Cell-cycle gating applied — {adata.n_obs:,} cells assigned."
                )
            except Exception as exc:
                st.error(f"Gating failed: {exc}")

    # Saved results ────────────────────────────────────────────────────────────
    if "cc_result_slider" in st.session_state:
        result_s = st.session_state["cc_result_slider"]
        summary_s = result_s["summary"]

        st.subheader("Saved results")
        st.dataframe(
            summary_s.style.format({"fraction": "{:.3f}"}),
            use_container_width=True,
            hide_index=True,
        )
        csv_s = summary_s.to_csv(index=False).encode()
        st.download_button(
            "Download summary CSV",
            csv_s,
            file_name="cell_cycle_summary_slider.csv",
            mime="text/csv",
            key="cc_slider_csv",
        )
