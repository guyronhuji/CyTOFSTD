"""Cell-cycle gating for CyTOF data using IdU, pH3, CyclinB1, and pRb.

Implements a transparent, rule-based exclusion hierarchy:
  IdU → S phase
  pH3 → M phase
  CyclinB1 → G2 phase
  pRb → Cycling G1 vs G0/quiescent
"""

from __future__ import annotations

import warnings
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


# ── Cell-cycle pseudotime ──────────────────────────────────────────────────────

# Biological ordering of all recognised phase labels (superset).
# Only phases actually present in the data are used; this list controls
# their relative order on the pseudotime circle.
DEFAULT_PSEUDOTIME_PHASE_ORDER: list[str] = [
    "G0", "G0_or_quiescent",
    "G1", "G1_or_quiescent", "early_G1",
    "late_G1", "Cycling_G1",
    "G1S",
    "S", "S_phase",
    "G2", "G2_phase",
    "G2M",
    "M", "M_phase",
]

# Within-phase scoring marker role(s) per phase label.
# Biological rationale:
#   G0/G1 phases — pRb increases as CDK4/6 phosphorylate Rb, preparing for S
#   S phase      — DNA content increases through replication (not IdU intensity,
#                  which identifies S cells but isn't monotone through S)
#   G2 phase     — CyclinB1 accumulates before M entry
#   G2/M         — average CyclinB1 and pH3 ranks
#   M phase      — pH3 (H3-Ser28/Ser10) peaks during mitosis
PHASE_SCORE_ROLES: dict[str, str | list[str]] = {
    "G0":              "pRb",
    "G0_or_quiescent": "pRb",
    "G1":              "pRb",
    "G1_or_quiescent": "pRb",
    "early_G1":        "pRb",
    "late_G1":         "pRb",
    "Cycling_G1":      "pRb",
    "G1S":             "pRb",
    "S":               "DNA",
    "S_phase":         "DNA",
    "G2":              "CyclinB1",
    "G2_phase":        "CyclinB1",
    "G2M":             ["CyclinB1", "pH3"],
    "M":               "pH3",
    "M_phase":         "pH3",
}

# Fallback role priority for unknown phase labels
_FALLBACK_ROLE_PRIORITY: list[str] = ["pRb", "DNA", "CyclinB1", "pH3", "IdU"]


# ── Pseudotime helper functions ────────────────────────────────────────────────

def _percentile_rank_0_1(values: np.ndarray) -> np.ndarray:
    """Percentile ranks scaled to [0, 1), NaN-aware, average-tied.

    - NaN inputs → NaN outputs (ignored when ranking).
    - Single non-NaN cell → 0.5.
    - Never returns exactly 1.0 (next phase starts at its own boundary).
    """
    values = np.asarray(values, dtype=float)
    out = np.full(len(values), np.nan)
    valid_mask = ~np.isnan(values)
    valid_vals = values[valid_mask]
    n = int(valid_mask.sum())

    if n == 0:
        return out
    if n == 1:
        out[valid_mask] = 0.5
        return out

    # Build 0-based rank via double-argsort
    order = np.argsort(valid_vals, kind="stable")
    rank = np.empty(n, dtype=float)
    rank[order] = np.arange(n, dtype=float)

    # Resolve ties: replace each tied group with its average rank
    sorted_vals = valid_vals[order]
    i = 0
    while i < n:
        j = i + 1
        while j < n and sorted_vals[j] == sorted_vals[i]:
            j += 1
        avg = (i + j - 1) / 2.0
        rank[order[i:j]] = avg
        i = j

    out[valid_mask] = rank / n  # maps to [0, (n-1)/n] ⊂ [0, 1)
    return out


def _get_marker_df(data: Any, marker_cols: dict[str, str]) -> pd.DataFrame:
    """Extract marker columns from AnnData (obs → X) or DataFrame.

    For AnnData, each column is looked up in ``adata.obs`` first, then in
    ``adata.to_df()`` (expression matrix), so both obs metadata and marker
    channels are supported transparently.
    """
    is_anndata = hasattr(data, "obs") and hasattr(data, "to_df")
    result: dict[str, pd.Series] = {}

    if is_anndata:
        expr_df: pd.DataFrame | None = None
        index = data.obs_names
        for role, col in marker_cols.items():
            if col in data.obs.columns:
                result[role] = pd.to_numeric(data.obs[col], errors="coerce")
            else:
                if expr_df is None:
                    expr_df = data.to_df()
                if col in expr_df.columns:
                    result[role] = pd.to_numeric(expr_df[col], errors="coerce")
                else:
                    warnings.warn(
                        f"Marker column '{col}' (role '{role}') not found in "
                        f"adata.obs or adata.X — skipping.",
                        UserWarning, stacklevel=4,
                    )
    else:
        index = data.index
        for role, col in marker_cols.items():
            if col in data.columns:
                result[role] = pd.to_numeric(data[col], errors="coerce")
            else:
                warnings.warn(
                    f"Marker column '{col}' (role '{role}') not found in "
                    f"DataFrame — skipping.",
                    UserWarning, stacklevel=4,
                )

    return pd.DataFrame(result, index=index)


def _get_observed_phase_order(
    phases: pd.Series,
    phase_order: list[str] | None,
) -> list[str]:
    """Filter *phase_order* to phases present in the data.

    Phases absent from the data are silently dropped. Phase labels present
    in the data but missing from *phase_order* are appended at the end with
    a ``UserWarning``. ``"Unclassified"`` is always excluded.
    """
    present = set(phases.dropna().astype(str).unique())
    present.discard("Unclassified")

    if phase_order is None:
        phase_order = DEFAULT_PSEUDOTIME_PHASE_ORDER

    ordered = [p for p in phase_order if p in present]
    unknown = present - set(phase_order)
    if unknown:
        warnings.warn(
            f"Phase label(s) not in the known phase_order: "
            f"{sorted(unknown)}. Appending after known phases. "
            "Pass a custom phase_order list to control their position.",
            UserWarning, stacklevel=4,
        )
        ordered.extend(sorted(unknown))

    return ordered


def _get_phase_widths(
    observed_phases: list[str],
    phase_widths: dict[str, float] | None,
) -> dict[str, float]:
    """Return normalised fractional widths for each observed phase.

    If *phase_widths* is ``None``, equal widths are assigned. Phases present
    in the data but missing from *phase_widths* receive an equal share of the
    remaining space (after summing provided widths). The result always sums
    to 1.0 across *observed_phases*.
    """
    n = len(observed_phases)
    if n == 0:
        return {}

    if phase_widths is None:
        w = 1.0 / n
        return {p: w for p in observed_phases}

    provided = {p: float(phase_widths[p]) for p in observed_phases if p in phase_widths}
    missing = [p for p in observed_phases if p not in provided]

    if missing:
        remaining = max(0.0, 1.0 - sum(provided.values()))
        share = remaining / len(missing) if remaining > 0 else 1.0 / n
        for p in missing:
            provided[p] = share

    total = sum(provided[p] for p in observed_phases)
    if total <= 0:
        total = 1.0
    return {p: provided[p] / total for p in observed_phases}


def _compute_within_phase_scores(
    obs: pd.DataFrame,
    markers: pd.DataFrame,
    phase_col: str,
    observed_phase_order: list[str],
) -> tuple[pd.Series, pd.Series]:
    """Compute within-phase [0, 1) score and integer phase index per cell.

    Uses :data:`PHASE_SCORE_ROLES` to choose the ranking marker(s) for each
    phase. Unknown phases fall back to the first available role in
    :data:`_FALLBACK_ROLE_PRIORITY`.

    Returns
    -------
    within_score : pd.Series of float in [0, 1), NaN for Unclassified cells
    phase_idx    : pd.Series of int (0-based index into observed_phase_order,
                   -1 for cells not in any observed phase)
    """
    within_score = pd.Series(np.nan, index=obs.index, dtype=float)
    phase_idx = pd.Series(-1, index=obs.index, dtype=int)
    phases = obs[phase_col].astype(str)

    for idx, phase in enumerate(observed_phase_order):
        mask = phases == phase
        if not mask.any():
            continue
        phase_idx[mask] = idx

        roles = PHASE_SCORE_ROLES.get(phase)

        if roles is None:
            # Unknown phase label — use first available fallback marker
            roles = next(
                (r for r in _FALLBACK_ROLE_PRIORITY if r in markers.columns),
                None,
            )
            if roles is None:
                within_score[mask] = 0.5
                continue

        if isinstance(roles, str):
            roles = [roles]

        score_arrays: list[np.ndarray] = []
        for role in roles:
            if role not in markers.columns:
                warnings.warn(
                    f"Marker role '{role}' needed for phase '{phase}' "
                    f"within-phase scoring is not available — skipping.",
                    UserWarning, stacklevel=5,
                )
                continue
            vals = markers.loc[mask, role].values
            score_arrays.append(_percentile_rank_0_1(vals))

        if not score_arrays:
            within_score[mask] = 0.5
        elif len(score_arrays) == 1:
            within_score[mask] = score_arrays[0]
        else:
            # Average the per-marker ranks (NaN-safe mean across markers)
            stacked = np.vstack(score_arrays)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                within_score[mask] = np.nanmean(stacked, axis=0)

    return within_score, phase_idx


# ── Main pseudotime function ───────────────────────────────────────────────────

def add_cell_cycle_pseudotime(
    data: Any,
    phase_col: str = "cell_cycle_phase",
    marker_cols: dict[str, str] | None = None,
    output_col: str = "cell_cycle_pseudotime",
    angle_col: str = "cell_cycle_angle",
    method: str = "rank_within_phase",
    phase_order: list[str] | None = None,
    phase_widths: dict[str, float] | None = None,
    overwrite: bool = False,
    copy: bool = True,
) -> Any:
    """Add a continuous cell-cycle pseudotime coordinate to gated CyTOF data.

    Requires that cells already carry a categorical phase label produced by
    :func:`gate_cell_cycle` (or equivalent). The phase labels are used as
    anchors; this function adds only the continuous within-cycle position.

    The pseudotime circle runs G0/G1 → S → G2 → M → back to G0/G1, encoded
    as a value in **[0, 1)** (and as an angle in **[0, 2π)**).

    Within each phase the ordering is determined by marker intensity ranks:

    * **G0 / G1 phases** — ``pRb``: increases as CDK4/6 phosphorylate Rb
    * **S phase** — ``DNA``: content increases through replication
    * **G2 phase** — ``CyclinB1``: accumulates before M entry
    * **G2/M** — average of ``CyclinB1`` and ``pH3`` ranks
    * **M phase** — ``pH3``: peaks at mitosis

    Args:
        data: pandas DataFrame or AnnData. For AnnData, phase labels and
            pseudotime output are stored in ``adata.obs``; marker values are
            read from ``adata.obs`` first, then ``adata.X``.
        phase_col: Column containing categorical cell-cycle phase labels.
        marker_cols: Dict mapping role → actual column name, e.g.
            ``{"DNA": "DNA1", "pRb": "pRb_S807"}``. Roles used for
            within-phase ordering: ``pRb``, ``DNA``, ``CyclinB1``, ``pH3``,
            ``IdU``. If ``None``, role names are used directly as column
            names.
        output_col: Name for the pseudotime output column.
        angle_col: Name for the 2π-scaled angle output column.
        method: Within-phase ordering method. Currently only
            ``"rank_within_phase"`` is supported.
        phase_order: Ordered list of phase labels (biological order G0→M).
            Only labels present in the data are used. Unknown labels are
            appended after known phases with a warning.
        phase_widths: Dict mapping phase label → fractional arc width.
            Normalised to sum to 1. Missing phases receive equal share of
            the remainder. ``None`` → equal widths.
        overwrite: If ``False`` (default), raise an error if the output
            columns already exist.
        copy: If ``True`` (default), return a modified copy. If ``False``,
            modify in-place.

    Returns:
        Modified DataFrame or AnnData with four new columns:

        * ``cell_cycle_pseudotime`` — continuous position in [0, 1)
        * ``cell_cycle_angle`` — angle in [0, 2π)
        * ``cell_cycle_phase_index`` — 0-based phase index (−1 for
          Unclassified / unrecognised cells)
        * ``cell_cycle_within_phase_rank`` — within-phase percentile rank

        ``"Unclassified"`` cells receive ``NaN`` in all four columns.

    Raises:
        ValueError: If *phase_col* is not found.
        ValueError: If output columns already exist and *overwrite* is False.
        NotImplementedError: If an unsupported *method* is requested.

    Example::

        adata = add_cell_cycle_pseudotime(
            adata,
            phase_col="cell_cycle_phase",
            marker_cols={"pRb": "pRb", "IdU": "IdU",
                         "CyclinB1": "CyclinB1", "pH3": "pH3",
                         "DNA": "DNA"},
            overwrite=True,
        )
        plot_cell_cycle_pseudotime_markers(adata)
        plot_cell_cycle_phase_circle(adata)
    """
    if method != "rank_within_phase":
        raise NotImplementedError(
            f"method='{method}' is not supported. "
            "Only 'rank_within_phase' is currently implemented."
        )

    is_anndata = hasattr(data, "obs") and hasattr(data, "to_df")

    if copy:
        data = data.copy()

    obs: pd.DataFrame = data.obs if is_anndata else data

    if phase_col not in obs.columns:
        raise ValueError(
            f"Phase column '{phase_col}' not found. "
            f"Available columns: {list(obs.columns[:15])}"
        )

    # Guard against overwriting existing results
    out_cols = [output_col, angle_col, "cell_cycle_phase_index",
                "cell_cycle_within_phase_rank"]
    existing = [c for c in out_cols if c in obs.columns]
    if existing and not overwrite:
        raise ValueError(
            f"Output columns already exist: {existing}. "
            "Pass overwrite=True to replace them."
        )

    if marker_cols is None:
        marker_cols = {r: r for r in ["pRb", "IdU", "CyclinB1", "pH3", "DNA"]}

    markers = _get_marker_df(data, marker_cols)

    phases = obs[phase_col].astype(str)
    observed_phases = _get_observed_phase_order(phases, phase_order)

    if not observed_phases:
        raise ValueError(
            f"No recognised phase labels found in '{phase_col}'. "
            f"Unique values present: {list(phases.unique())}"
        )

    widths = _get_phase_widths(observed_phases, phase_widths)

    # Compute cumulative start position for each phase
    phase_starts: dict[str, float] = {}
    cumulative = 0.0
    for p in observed_phases:
        phase_starts[p] = cumulative
        cumulative += widths[p]

    within_score, phase_idx = _compute_within_phase_scores(
        obs, markers, phase_col, observed_phases
    )

    # pseudotime = phase_start + within_rank * phase_width
    pseudotime = pd.Series(np.nan, index=obs.index, dtype=float)
    for phase in observed_phases:
        mask = phases == phase
        if not mask.any():
            continue
        pseudotime[mask] = phase_starts[phase] + within_score[mask] * widths[phase]

    # Clamp to [0, 1) via modulo (handles any floating-point drift)
    valid = pseudotime.notna()
    pseudotime[valid] = pseudotime[valid] % 1.0

    angle = 2.0 * np.pi * pseudotime

    if is_anndata:
        data.obs[output_col] = pseudotime.values
        data.obs[angle_col] = angle.values
        data.obs["cell_cycle_phase_index"] = phase_idx.values
        data.obs["cell_cycle_within_phase_rank"] = within_score.values
    else:
        data[output_col] = pseudotime.values
        data[angle_col] = angle.values
        data["cell_cycle_phase_index"] = phase_idx.values
        data["cell_cycle_within_phase_rank"] = within_score.values

    return data


# ── Pseudotime validation plots ────────────────────────────────────────────────

def plot_cell_cycle_pseudotime_markers(
    data: Any,
    pseudotime_col: str = "cell_cycle_pseudotime",
    phase_col: str = "cell_cycle_phase",
    marker_cols: list[str] | dict[str, str] | None = None,
    bins: int = 50,
    use_hexbin: bool = True,
    figsize_per_marker: tuple[float, float] = (5, 3),
) -> "plt.Figure":
    """Plot marker intensity vs cell-cycle pseudotime for validation.

    Each marker gets a row. With ``use_hexbin=True`` (default), a density
    hexbin is drawn; with ``use_hexbin=False``, cells are coloured by their
    gated phase. A binned mean trend line is overlaid in both cases. Vertical
    dashed lines mark the start of each phase.

    Expected biological patterns:

    * **pRb**: rises before S phase
    * **IdU**: high in S phase only
    * **DNA**: increases through S, stays high in G2/M
    * **CyclinB1**: accumulates in G2
    * **pH3**: peaks in M

    Args:
        data: pandas DataFrame or AnnData.
        pseudotime_col: Column with pseudotime values (from
            :func:`add_cell_cycle_pseudotime`).
        phase_col: Column with categorical phase labels.
        marker_cols: Markers to plot. List of column names, dict of
            role → column, or ``None`` to use
            ``["pRb", "IdU", "DNA", "CyclinB1", "pH3"]``.
        bins: Number of hexbin / trend-line bins along the pseudotime axis.
        use_hexbin: If ``True``, plot density hexbin. If ``False``, plot
            phase-coloured scatter.
        figsize_per_marker: ``(width, height)`` for each marker row.

    Returns:
        matplotlib Figure with one row per marker.
    """
    is_anndata = hasattr(data, "obs") and hasattr(data, "to_df")
    obs: pd.DataFrame = data.obs.copy() if is_anndata else data.copy()
    expr_df: pd.DataFrame | None = None

    if pseudotime_col not in obs.columns:
        raise ValueError(
            f"Column '{pseudotime_col}' not found. "
            "Run add_cell_cycle_pseudotime first."
        )

    if marker_cols is None:
        marker_cols = ["pRb", "IdU", "DNA", "CyclinB1", "pH3"]

    if isinstance(marker_cols, dict):
        col_list = list(marker_cols.values())
    else:
        col_list = list(marker_cols)

    # Resolve columns — obs first, then expression matrix
    resolved: list[str] = []
    for col in col_list:
        if col in obs.columns:
            resolved.append(col)
        else:
            if is_anndata and expr_df is None:
                expr_df = data.to_df()
            if expr_df is not None and col in expr_df.columns:
                obs[col] = expr_df[col].values
                resolved.append(col)
            else:
                warnings.warn(
                    f"Marker column '{col}' not found — skipping.",
                    UserWarning, stacklevel=2,
                )

    if not resolved:
        raise ValueError("No marker columns available to plot.")

    pt = pd.to_numeric(obs[pseudotime_col], errors="coerce").values
    has_phase = phase_col in obs.columns
    phases_arr = obs[phase_col].values if has_phase else None

    n_markers = len(resolved)
    fig, axes = plt.subplots(
        n_markers, 1,
        figsize=(figsize_per_marker[0], figsize_per_marker[1] * n_markers),
        squeeze=False,
    )

    for ax_row, col in zip(axes, resolved):
        ax = ax_row[0]
        vals = pd.to_numeric(obs[col], errors="coerce").values
        valid = ~np.isnan(pt) & ~np.isnan(vals)
        pt_v, vals_v = pt[valid], vals[valid]

        if use_hexbin:
            hb = ax.hexbin(pt_v, vals_v, gridsize=bins, cmap="Blues",
                           mincnt=1, linewidths=0.1)
            fig.colorbar(hb, ax=ax, pad=0.01, fraction=0.03, label="cells")
        else:
            # Phase-coloured scatter
            if phases_arr is not None:
                ph_v = phases_arr[valid]
                for phase in PHASE_ORDER:
                    ph_mask = ph_v == phase
                    if not ph_mask.any():
                        continue
                    ax.scatter(
                        pt_v[ph_mask], vals_v[ph_mask],
                        c=PHASE_COLORS.get(phase, "#aaaaaa"),
                        s=2, alpha=0.4, linewidths=0, rasterized=True,
                        label=phase.replace("_", " "),
                    )
            else:
                ax.scatter(pt_v, vals_v, s=2, alpha=0.3,
                           color="#4a90d9", rasterized=True)

        # Binned mean trend line
        bin_edges = np.linspace(0, 1, bins + 1)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0
        bin_means = np.array([
            np.nanmean(vals_v[(pt_v >= lo) & (pt_v < hi)])
            if ((pt_v >= lo) & (pt_v < hi)).sum() >= 3 else np.nan
            for lo, hi in zip(bin_edges[:-1], bin_edges[1:])
        ])
        ok = ~np.isnan(bin_means)
        if ok.sum() > 1:
            ax.plot(bin_centers[ok], bin_means[ok],
                    color="#d62728", lw=2.0, zorder=6, label="mean")

        # Phase boundary lines
        if has_phase:
            for phase in PHASE_ORDER:
                ph_mask = (phases_arr[valid] == phase)
                if not ph_mask.any():
                    continue
                start = np.nanmin(pt_v[ph_mask])
                color = PHASE_COLORS.get(phase, "#888888")
                ax.axvline(start, color=color, lw=1.0, ls="--", alpha=0.7)
                ax.text(start + 0.005, ax.get_ylim()[1],
                        phase.replace("_or_", "\n").replace("_", " "),
                        fontsize=6, va="top", color=color, alpha=0.85)

        ax.set_xlim(0, 1)
        ax.set_ylabel(col, fontsize=9)
        ax.set_xlabel("Cell-cycle pseudotime", fontsize=9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.suptitle("Marker dynamics along cell-cycle pseudotime",
                 fontsize=11, y=1.01)
    fig.tight_layout()
    return fig


def plot_cell_cycle_phase_circle(
    data: Any,
    angle_col: str = "cell_cycle_angle",
    phase_col: str = "cell_cycle_phase",
    n_subsample: int = 5000,
    seed: int = 42,
) -> "plt.Figure":
    """Polar scatter plot of cell positions on the cell-cycle circle.

    Each cell is placed at its cell-cycle angle. The radial axis is jittered
    for visibility. Intended as a quick sanity check that cells are
    distributed reasonably around the cycle and that phase colours match
    expectations.

    Args:
        data: pandas DataFrame or AnnData.
        angle_col: Column containing cell-cycle angle in [0, 2π)
            (from :func:`add_cell_cycle_pseudotime`).
        phase_col: Column with categorical phase labels.
        n_subsample: Maximum cells to plot (random subsample).
        seed: Random seed for reproducible subsampling.

    Returns:
        matplotlib Figure with a single polar axis.
    """
    is_anndata = hasattr(data, "obs") and hasattr(data, "to_df")
    obs: pd.DataFrame = data.obs if is_anndata else data

    if angle_col not in obs.columns:
        raise ValueError(
            f"Angle column '{angle_col}' not found. "
            "Run add_cell_cycle_pseudotime first."
        )

    angles = pd.to_numeric(obs[angle_col], errors="coerce")
    valid = ~angles.isna()
    angles_v = angles[valid].values
    phases_v = obs.loc[valid, phase_col].values if phase_col in obs.columns else None

    rng = np.random.default_rng(seed)
    n = len(angles_v)
    if n > n_subsample:
        idx = rng.choice(n, n_subsample, replace=False)
        angles_v = angles_v[idx]
        if phases_v is not None:
            phases_v = phases_v[idx]

    fig, ax = plt.subplots(figsize=(5.5, 5.5), subplot_kw={"projection": "polar"})

    # Jitter radius so overlapping points are visible
    r = rng.uniform(0.4, 1.0, len(angles_v))

    if phases_v is not None:
        for phase in PHASE_ORDER:
            mask = phases_v == phase
            if not mask.any():
                continue
            ax.scatter(
                angles_v[mask], r[mask],
                c=PHASE_COLORS.get(phase, "#aaaaaa"),
                s=3, alpha=0.5, linewidths=0, rasterized=True,
                label=phase.replace("_", " "),
            )
    else:
        ax.scatter(angles_v, r, s=3, alpha=0.5,
                   color="#4a90d9", linewidths=0, rasterized=True)

    # Orient so G0/G1 starts at top, cycle runs clockwise
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    ax.set_rlim(0, 1.3)
    ax.set_rticks([])
    ax.set_title("Cell-cycle pseudotime", fontsize=11, pad=15)

    if phases_v is not None:
        handles = [
            plt.Line2D([0], [0], marker="o", color="w", markersize=8,
                       markerfacecolor=PHASE_COLORS.get(p, "#aaaaaa"),
                       label=p.replace("_", " "))
            for p in PHASE_ORDER if (phases_v == p).any()
        ]
        ax.legend(handles=handles, bbox_to_anchor=(1.35, 1.0),
                  loc="upper left", fontsize=8, framealpha=0.8)

    fig.tight_layout()
    return fig
