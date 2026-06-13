"""Cell-cycle gating for CyTOF data using IdU, pH3, CyclinB1, and pRb.

Implements a transparent, rule-based exclusion hierarchy:
  IdU → S phase
  pH3 → M phase
  CyclinB1 → G2 phase
  pRb → Cycling G1 vs G0/quiescent
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ── Constants ──────────────────────────────────────────────────────────────────

CELL_CYCLE_MARKER_ALIASES: dict[str, list[str]] = {
    "IdU": ["IdU", "IDU", "BrdU", "IdU_BrdU", "5-iodo-2'-deoxyuridine"],
    "pH3": ["pH3", "pHH3", "H3S10ph", "H3S28ph", "H3pS28", "H3pS10", "pHH3_S28", "pHH3_S10"],
    "CyclinB1": ["CyclinB1", "Cyclin B1", "CCNB1", "Cyclin_B1", "CycB1", "cyclinB1"],
    "pRb": ["pRb", "phospho-Rb", "pRB", "Rb_phospho", "pRb_S807", "Rb_S807", "Rb_pS807"],
}

# Default auto-threshold strategy per role.
# Bimodal markers (IdU, pH3, CyclinB1) use Otsu's method to find the valley
# between negative and positive populations. Continuous markers (Ki67, pRb)
# use quantile, because their positive fraction varies biologically.
DEFAULT_THRESHOLD_METHODS: dict[str, str] = {
    "IdU": "otsu",
    "pH3": "otsu",
    "CyclinB1": "otsu",
    "pRb": "quantile",
}

# Quantile fallbacks used when strategy is "quantile" (or Otsu fails).
DEFAULT_QUANTILE_THRESHOLDS: dict[str, float] = {
    "IdU": 0.95,
    "pH3": 0.98,
    "CyclinB1": 0.85,
    "pRb": 0.60,
}

REQUIRED_ROLES: set[str] = {"IdU", "pH3", "CyclinB1"}
OPTIONAL_ROLES: set[str] = {"pRb"}

PHASE_ORDER: list[str] = [
    "G0_or_quiescent",
    "G1_or_quiescent",
    "Cycling_G1",
    "S_phase",
    "G2_phase",
    "M_phase",
    "Unclassified",
]

PHASE_COLORS: dict[str, str] = {
    "G0_or_quiescent": "#9e9e9e",  # = manual G0
    "G1_or_quiescent": "#bdbdbd",  # light grey (G0/G1 ambiguous, no pRb)
    "Cycling_G1":      "#1f77b4",  # = manual G1 (pRb+)
    "S_phase":         "#2ca02c",  # = manual S
    "G2_phase":        "#ff7f0e",  # = manual G2
    "M_phase":         "#d62728",  # = manual M
    "Unclassified":    "#aaaaaa",
}


# ── Marker auto-detection ──────────────────────────────────────────────────────

def auto_detect_cell_cycle_markers(
    var_names: list[str],
    aliases: dict[str, list[str]] | None = None,
) -> dict[str, str]:
    """Try to find cell-cycle marker columns in a list of variable names.

    Returns a (possibly partial) dict mapping role → actual column name.
    Missing roles are omitted — check the returned keys against
    CELL_CYCLE_MARKER_ALIASES.keys() to see what was found.
    """
    if aliases is None:
        aliases = CELL_CYCLE_MARKER_ALIASES

    detected: dict[str, str] = {}
    var_set = set(var_names)
    var_lower = {v.lower(): v for v in var_names}

    for role, alias_list in aliases.items():
        for alias in alias_list:
            # 1. Exact match
            if alias in var_set:
                detected[role] = alias
                break
            # 2. Case-insensitive exact match
            a_lower = alias.lower()
            if a_lower in var_lower:
                detected[role] = var_lower[a_lower]
                break
            # 3. Prefix match (alias_lower is a prefix of var, separated by _ or -)
            if len(a_lower) >= 3:
                for vn_lower, vn_orig in var_lower.items():
                    if (
                        vn_lower == a_lower
                        or vn_lower.startswith(a_lower + "_")
                        or vn_lower.startswith(a_lower + "-")
                    ):
                        detected[role] = vn_orig
                        break
            if role in detected:
                break

    return detected


# ── Data extraction ────────────────────────────────────────────────────────────

def extract_marker_dataframe(
    data: Any,
    marker_map: dict[str, str],
    layer: str | None = None,
) -> pd.DataFrame:
    """Extract cell-cycle marker columns from a DataFrame or AnnData.

    Args:
        data: pandas DataFrame or AnnData.
        marker_map: Dict mapping role (e.g. ``"IdU"``) → actual column name.
        layer: AnnData layer key. ``None`` uses ``adata.X`` via ``to_df()``.

    Returns:
        DataFrame indexed identically to the input, containing only the
        requested marker columns (renamed to the actual column names).
    """
    try:
        import anndata as _ad
        is_adata = isinstance(data, _ad.AnnData)
    except ImportError:
        is_adata = False

    if is_adata:
        if layer is None:
            base_df = data.to_df()
        else:
            if layer not in data.layers:
                raise ValueError(
                    f"Layer '{layer}' not found. "
                    f"Available: {list(data.layers.keys())}"
                )
            base_df = pd.DataFrame(
                data.layers[layer],
                index=data.obs_names,
                columns=data.var_names,
            )
    elif isinstance(data, pd.DataFrame):
        base_df = data
    else:
        raise TypeError(
            f"data must be a pandas DataFrame or AnnData, got {type(data).__name__}"
        )

    result: dict[str, pd.Series] = {}
    for role, col in marker_map.items():
        if col not in base_df.columns:
            raise ValueError(
                f"Column '{col}' (role: '{role}') not found. "
                f"Available: {list(base_df.columns[:15])}"
            )
        result[col] = pd.to_numeric(base_df[col], errors="coerce")

    return pd.DataFrame(result, index=base_df.index)


# ── Thresholding ───────────────────────────────────────────────────────────────

def _otsu_threshold(values: np.ndarray, n_bins: int = 512) -> float:
    """Otsu's method: find the threshold that maximises inter-class variance.

    Works well for bimodal distributions (e.g. IdU, pH3, CyclinB1) where a
    clear valley separates negative and positive populations.

    Returns the threshold value in the original data units.
    """
    hist, edges = np.histogram(values, bins=n_bins)
    hist = hist.astype(np.float64)
    total = hist.sum()
    if total == 0:
        return float(np.median(values))
    hist /= total
    centers = (edges[:-1] + edges[1:]) / 2.0

    # Cumulative weight and mean of the lower class
    w0 = np.cumsum(hist)
    mu_total = np.sum(hist * centers)
    mu0 = np.cumsum(hist * centers)

    w1 = 1.0 - w0
    # Avoid division by zero at the extremes
    with np.errstate(divide="ignore", invalid="ignore"):
        mu0_safe = np.where(w0 > 0, mu0 / w0, 0.0)
        mu1_safe = np.where(w1 > 0, (mu_total - mu0) / w1, 0.0)
    sigma_b2 = w0 * w1 * (mu0_safe - mu1_safe) ** 2

    # Pick the midpoint of any plateau of maxima
    max_val = sigma_b2.max()
    candidates = np.where(sigma_b2 >= max_val * 0.9999)[0]
    best_idx = int(candidates[len(candidates) // 2])
    return float(centers[best_idx])


def calculate_thresholds(
    marker_df: pd.DataFrame,
    marker_map: dict[str, str],
    thresholds: dict[str, float] | None = None,
    quantile_thresholds: dict[str, float] | None = None,
    threshold_methods: dict[str, str] | None = None,
) -> dict[str, float]:
    """Compute per-marker gating thresholds.

    User-supplied thresholds (``thresholds`` dict) always take precedence.
    For automatic estimation, each role has a default strategy:

    - ``"otsu"``     — Otsu's method; optimal for bimodal distributions
                       (IdU, pH3, CyclinB1). Finds the valley between
                       negative and positive populations.
    - ``"quantile"`` — Fixed quantile of the distribution; appropriate for
                       continuous markers (Ki67, pRb) where the positive
                       fraction is biologically variable.

    Args:
        marker_df: DataFrame with marker columns.
        marker_map: Role → column name mapping.
        thresholds: User-supplied thresholds by role. Missing roles use auto.
        quantile_thresholds: Per-role quantile overrides (used when strategy
            is ``"quantile"`` or as a fallback if Otsu produces an extreme
            value). Defaults to :data:`DEFAULT_QUANTILE_THRESHOLDS`.
        threshold_methods: Per-role strategy override (``"otsu"`` or
            ``"quantile"``). Defaults to :data:`DEFAULT_THRESHOLD_METHODS`.

    Returns:
        Dict mapping role → threshold value (float).
    """
    if quantile_thresholds is None:
        quantile_thresholds = DEFAULT_QUANTILE_THRESHOLDS.copy()
    if threshold_methods is None:
        threshold_methods = DEFAULT_THRESHOLD_METHODS.copy()

    result: dict[str, float] = {}
    for role, col in marker_map.items():
        if thresholds and role in thresholds and thresholds[role] is not None:
            result[role] = float(thresholds[role])
            continue

        values = pd.to_numeric(marker_df[col], errors="coerce").dropna().values
        if len(values) == 0:
            raise ValueError(
                f"Marker '{col}' (role: '{role}') has no valid numeric values."
            )

        method = threshold_methods.get(role, "quantile")

        if method == "otsu":
            thr = _otsu_threshold(values)
            # Sanity check: if Otsu returns a value outside the central 10–99 %ile
            # range it likely failed (all-negative or all-positive panel).
            # Fall back to quantile in that case.
            p10 = float(np.percentile(values, 10))
            p99 = float(np.percentile(values, 99))
            if thr <= p10 or thr >= p99:
                q = quantile_thresholds.get(role, 0.90)
                thr = float(np.quantile(values, q))
        else:
            q = quantile_thresholds.get(role, 0.90)
            thr = float(np.quantile(values, q))

        result[role] = thr

    return result


# ── Gating hierarchy ───────────────────────────────────────────────────────────

def assign_cell_cycle_phase(
    marker_df: pd.DataFrame,
    marker_map: dict[str, str],
    thresholds: dict[str, float],
) -> pd.DataFrame:
    """Assign cells to mutually exclusive cell-cycle phases.

    Required roles: ``IdU``, ``pH3``, ``CyclinB1``.
    Optional role: ``pRb``.

    Gating hierarchy (order matters):
      1. IdU+    → S_phase
      2. pH3+    → M_phase
      3. CyclinB1+ → G2_phase
      For remaining cells:
        pRb present: pRb+ → Cycling_G1, pRb- → G0_or_quiescent
        pRb absent:  all  → G1_or_quiescent (cannot distinguish)

    Returns:
        Copy of marker_df extended with ``cell_cycle_gate_*_pos`` boolean
        columns and a ``cell_cycle_phase`` string column.
    """
    n = len(marker_df)
    phases = np.full(n, "Unclassified", dtype=object)

    def _pos(role: str) -> np.ndarray:
        col = marker_map[role]
        vals = pd.to_numeric(marker_df[col], errors="coerce").fillna(0).values
        return vals > thresholds[role]

    has_prb = "pRb" in marker_map

    idu_pos = _pos("IdU")
    ph3_pos = _pos("pH3")
    cycb_pos = _pos("CyclinB1")
    prb_pos = _pos("pRb") if has_prb else None

    phases[idu_pos] = "S_phase"

    remaining = phases == "Unclassified"
    phases[remaining & ph3_pos] = "M_phase"

    remaining = phases == "Unclassified"
    phases[remaining & cycb_pos] = "G2_phase"

    remaining = phases == "Unclassified"

    if has_prb:
        phases[remaining & prb_pos] = "Cycling_G1"
        remaining = phases == "Unclassified"
        phases[remaining & ~prb_pos] = "G0_or_quiescent"
    else:
        phases[remaining] = "G1_or_quiescent"

    out = marker_df.copy()
    out["cell_cycle_gate_IdU_pos"] = idu_pos
    out["cell_cycle_gate_pH3_pos"] = ph3_pos
    out["cell_cycle_gate_CyclinB1_pos"] = cycb_pos
    if has_prb:
        out["cell_cycle_gate_pRb_pos"] = prb_pos
    out["cell_cycle_phase"] = phases

    return out


# ── Summary ────────────────────────────────────────────────────────────────────

def summarize_cell_cycle(gated_df: pd.DataFrame) -> pd.DataFrame:
    """Return counts and fractions for each cell-cycle phase.

    All phases in :data:`PHASE_ORDER` are included; phases with zero cells
    are omitted from the output.

    Returns:
        DataFrame with columns ``cell_cycle_phase``, ``n_cells``, ``fraction``.
    """
    counts = (
        gated_df["cell_cycle_phase"]
        .value_counts()
        .rename_axis("cell_cycle_phase")
        .reset_index(name="n_cells")
    )
    template = pd.DataFrame({"cell_cycle_phase": PHASE_ORDER})
    summary = template.merge(counts, on="cell_cycle_phase", how="left")
    summary["n_cells"] = summary["n_cells"].fillna(0).astype(int)
    total = summary["n_cells"].sum()
    summary["fraction"] = summary["n_cells"] / total if total > 0 else 0.0
    return summary[summary["n_cells"] > 0].reset_index(drop=True)


# ── QC plots ───────────────────────────────────────────────────────────────────

def plot_cell_cycle_marker_qc(
    marker_df: pd.DataFrame,
    marker_map: dict[str, str],
    thresholds: dict[str, float] | None = None,
    n_bins: int = 80,
    figsize: tuple[float, float] | None = None,
    output_dir: str | None = None,
) -> list:
    """Plot the distribution of each marker with an optional threshold line.

    Args:
        marker_df: DataFrame with marker columns.
        marker_map: Role → column name mapping.
        thresholds: Dict of role → threshold to draw as a vertical line.
        n_bins: Number of histogram bins.
        figsize: Per-figure size tuple ``(width, height)``.
        output_dir: If provided, save each figure as
            ``cell_cycle_qc_{role}.png`` in this directory.

    Returns:
        List of matplotlib ``Figure`` objects (one per marker).
    """
    figs: list = []

    for role, col in marker_map.items():
        if col not in marker_df.columns:
            continue

        values = pd.to_numeric(marker_df[col], errors="coerce").dropna().values

        fig, ax = plt.subplots(figsize=figsize or (6, 3.5))
        ax.hist(values, bins=n_bins, color="#4a90d9", alpha=0.78, edgecolor="none")

        if thresholds and role in thresholds:
            thr = thresholds[role]
            ax.axvline(thr, color="#e05151", linewidth=1.6, linestyle="--",
                       label=f"threshold = {thr:.3f}")
            n_pos = int((values > thr).sum())
            pct = 100.0 * n_pos / len(values) if len(values) > 0 else 0.0
            ax.legend(
                title=f"{n_pos:,} pos ({pct:.1f}%)",
                fontsize=8, title_fontsize=8,
            )

        ax.set_xlabel(col, fontsize=10)
        ax.set_ylabel("Cells", fontsize=10)
        ax.set_title(f"{role}  ·  {col}", fontsize=10, fontweight="bold")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.tight_layout()

        if output_dir is not None:
            from pathlib import Path as _Path
            out_path = _Path(output_dir) / f"cell_cycle_qc_{role}.png"
            fig.savefig(out_path, dpi=120, bbox_inches="tight")

        figs.append(fig)

    return figs


def plot_marker_thresholds(
    df: pd.DataFrame,
    marker_map: dict[str, str],
    thresholds: dict[str, float] | None = None,
    n_bins: int = 80,
    figsize: tuple[float, float] | None = None,
) -> list:
    """Alias for :func:`plot_cell_cycle_marker_qc` (spec compatibility)."""
    return plot_cell_cycle_marker_qc(
        df, marker_map, thresholds=thresholds,
        n_bins=n_bins, figsize=figsize,
    )


def plot_cell_cycle_phase_fractions(
    summary: pd.DataFrame,
    figsize: tuple[float, float] | None = None,
    output_path: str | None = None,
) -> "plt.Figure":
    """Bar chart of cell-cycle phase fractions."""
    fig, ax = plt.subplots(figsize=figsize or (7, 4))

    phases = summary["cell_cycle_phase"].tolist()
    fracs = summary["fraction"].tolist()
    colors = [PHASE_COLORS.get(p, "#aaaaaa") for p in phases]

    bars = ax.bar(phases, fracs, color=colors, edgecolor="none", width=0.6)
    for bar, frac in zip(bars, fracs):
        if frac > 0.005:
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.004,
                f"{100 * frac:.1f}%",
                ha="center", va="bottom", fontsize=9,
            )

    ax.set_ylabel("Fraction of cells", fontsize=10)
    ax.set_title("Cell-cycle phase distribution", fontsize=11)
    ax.set_ylim(0, (max(fracs) if fracs else 1) * 1.2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.xticks(rotation=30, ha="right", fontsize=9)
    fig.tight_layout()

    if output_path is not None:
        fig.savefig(output_path, dpi=120, bbox_inches="tight")

    return fig


def plot_cell_cycle_fractions_by_group(
    gated_df: pd.DataFrame,
    groupby: str,
    figsize: tuple[float, float] | None = None,
) -> "plt.Figure":
    """Stacked bar chart of phase fractions per group."""
    if groupby not in gated_df.columns:
        raise ValueError(f"Column '{groupby}' not found in dataframe.")

    pivot = (
        gated_df.groupby([groupby, "cell_cycle_phase"])
        .size()
        .unstack(fill_value=0)
    )
    cols = [p for p in PHASE_ORDER if p in pivot.columns]
    pivot = pivot[cols]
    pivot_norm = pivot.div(pivot.sum(axis=1), axis=0)
    colors = [PHASE_COLORS.get(c, "#aaaaaa") for c in pivot_norm.columns]

    n_groups = len(pivot_norm)
    fig, ax = plt.subplots(figsize=figsize or (max(6, n_groups * 0.8 + 2), 5))
    pivot_norm.plot(kind="bar", stacked=True, ax=ax, color=colors, edgecolor="none", width=0.7)

    ax.set_ylabel("Fraction of cells", fontsize=10)
    ax.set_xlabel(groupby, fontsize=10)
    ax.set_title(f"Cell-cycle distribution by {groupby}", fontsize=11)
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8, title="Phase")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.xticks(rotation=30, ha="right", fontsize=9)
    fig.tight_layout()

    return fig


# ── Main wrapper ───────────────────────────────────────────────────────────────

def gate_cell_cycle(
    data: Any,
    marker_map: dict[str, str],
    layer: str | None = None,
    thresholds: dict[str, float] | None = None,
    quantile_thresholds: dict[str, float] | None = None,
    threshold_methods: dict[str, str] | None = None,
    apply_arcsinh: bool = False,
    cofactor: float = 5.0,
    return_adata: bool = True,
) -> dict[str, Any]:
    """Full cell-cycle gating workflow for a DataFrame or AnnData.

    Args:
        data: pandas DataFrame or AnnData. When AnnData, results are written
            to ``adata.obs`` and ``adata.uns`` if ``return_adata=True``.
        marker_map: Dict mapping role → actual column/var name.
            Required roles: ``IdU``, ``pH3``, ``CyclinB1``.
            Optional role: ``pRb``.
        layer: AnnData layer key. ``None`` uses ``adata.X`` (recommended:
            use arcsinh-transformed data, not z-scored).
        thresholds: User-supplied thresholds by role. Missing roles fall back
            to quantile-based estimates.
        quantile_thresholds: Per-role quantile overrides for auto-thresholding.
        apply_arcsinh: Apply arcsinh transform before gating (use only when
            data contains raw CyTOF intensities, not already-transformed values).
        cofactor: Arcsinh cofactor (default 5.0).
        return_adata: Write results back to ``adata.obs`` / ``adata.uns`` when
            input is AnnData.

    Returns:
        Dict with keys:
        - ``gated_df``: DataFrame with gate columns and ``cell_cycle_phase``.
        - ``summary``: Phase counts and fractions.
        - ``thresholds``: Thresholds used (role → value).
        - ``adata``: Updated AnnData (only when input is AnnData and
          ``return_adata=True``).
    """
    try:
        import anndata as _ad
        is_adata = isinstance(data, _ad.AnnData)
    except ImportError:
        is_adata = False

    marker_df = extract_marker_dataframe(data, marker_map, layer=layer)

    if apply_arcsinh:
        if cofactor <= 0:
            raise ValueError("cofactor must be > 0")
        for col in marker_df.columns:
            marker_df[col] = np.arcsinh(marker_df[col] / cofactor)

    resolved = calculate_thresholds(
        marker_df, marker_map,
        thresholds=thresholds,
        quantile_thresholds=quantile_thresholds,
        threshold_methods=threshold_methods,
    )

    gated_df = assign_cell_cycle_phase(marker_df, marker_map, resolved)
    summary = summarize_cell_cycle(gated_df)

    result: dict[str, Any] = {
        "gated_df": gated_df,
        "summary": summary,
        "thresholds": resolved,
    }

    if is_adata and return_adata:
        adata = data
        adata.obs["cell_cycle_phase"] = gated_df["cell_cycle_phase"].values
        for role in marker_map:
            gate_col = f"cell_cycle_gate_{role}_pos"
            if gate_col in gated_df.columns:
                adata.obs[gate_col] = gated_df[gate_col].values
        adata.uns["cell_cycle_gating_thresholds"] = resolved
        adata.uns["cell_cycle_gating_marker_map"] = marker_map
        result["adata"] = adata

    return result


# ── Per-group gating ───────────────────────────────────────────────────────────

def gate_cell_cycle_by_group(
    data: Any,
    marker_map: dict[str, str],
    groupby: str,
    layer: str | None = None,
    thresholds: dict[str, float] | None = None,
    quantile_thresholds: dict[str, float] | None = None,
    threshold_methods: dict[str, str] | None = None,
    apply_arcsinh: bool = False,
    cofactor: float = 5.0,
) -> dict[str, Any]:
    """Cell-cycle gating with independent threshold calculation per group.

    Useful when CyTOF signal intensity varies across batches or samples.

    Args:
        groupby: Column in ``adata.obs`` or DataFrame to split on.

    Returns:
        Dict with ``gated_df``, ``summary``, ``per_group_summary``,
        ``per_group_thresholds``.
    """
    try:
        import anndata as _ad
        is_adata = isinstance(data, _ad.AnnData)
    except ImportError:
        is_adata = False

    obs = data.obs if is_adata else data
    if not isinstance(obs, pd.DataFrame):
        raise TypeError("data must be a DataFrame or AnnData")
    if groupby not in obs.columns:
        raise ValueError(
            f"groupby column '{groupby}' not found. "
            f"Available: {list(obs.columns[:15])}"
        )

    groups = obs[groupby].unique()
    parts: list[pd.DataFrame] = []
    per_group_thresholds: dict[str, dict[str, float]] = {}

    for grp in groups:
        mask = (obs[groupby] == grp).values
        grp_data = data[mask] if is_adata else data.loc[mask]
        res = gate_cell_cycle(
            grp_data, marker_map, layer=layer,
            thresholds=thresholds, quantile_thresholds=quantile_thresholds,
            threshold_methods=threshold_methods,
            apply_arcsinh=apply_arcsinh, cofactor=cofactor,
            return_adata=False,
        )
        grp_df = res["gated_df"].copy()
        grp_df[groupby] = str(grp)
        parts.append(grp_df)
        per_group_thresholds[str(grp)] = res["thresholds"]

    all_gated = pd.concat(parts, axis=0)
    # Restore original cell order
    all_gated = all_gated.loc[obs.index]

    per_group_summary = {
        str(grp): summarize_cell_cycle(all_gated[all_gated[groupby] == str(grp)])
        for grp in groups
    }

    return {
        "gated_df": all_gated,
        "summary": summarize_cell_cycle(all_gated),
        "per_group_summary": per_group_summary,
        "per_group_thresholds": per_group_thresholds,
    }
