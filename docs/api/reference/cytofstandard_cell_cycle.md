# `cytofstandard.cell_cycle`

- Source: `cytofstandard/cell_cycle.py`

Cell-cycle gating for CyTOF data using IdU, pH3, CyclinB1, and pRb.

Implements a transparent, rule-based exclusion hierarchy:
  IdU → S phase
  pH3 → M phase
  CyclinB1 → G2 phase
  pRb → Cycling G1 vs G0/quiescent

## Top-level Functions

### `add_cell_cycle_pseudotime(data: Any, phase_col: str = 'cell_cycle_phase', marker_cols: dict[str, str] | None = None, output_col: str = 'cell_cycle_pseudotime', angle_col: str = 'cell_cycle_angle', method: str = 'rank_within_phase', phase_order: list[str] | None = None, phase_widths: dict[str, float] | None = None, overwrite: bool = False, copy: bool = True) -> Any`

Add a continuous cell-cycle pseudotime coordinate to gated CyTOF data.

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

### `assign_cell_cycle_phase(marker_df: pd.DataFrame, marker_map: dict[str, str], thresholds: dict[str, float]) -> pd.DataFrame`

Assign cells to mutually exclusive cell-cycle phases.

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

### `auto_detect_cell_cycle_markers(var_names: list[str], aliases: dict[str, list[str]] | None = None) -> dict[str, str]`

Try to find cell-cycle marker columns in a list of variable names.

Returns a (possibly partial) dict mapping role → actual column name.
Missing roles are omitted — check the returned keys against
CELL_CYCLE_MARKER_ALIASES.keys() to see what was found.

### `calculate_thresholds(marker_df: pd.DataFrame, marker_map: dict[str, str], thresholds: dict[str, float] | None = None, quantile_thresholds: dict[str, float] | None = None, threshold_methods: dict[str, str] | None = None) -> dict[str, float]`

Compute per-marker gating thresholds.

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

### `extract_marker_dataframe(data: Any, marker_map: dict[str, str], layer: str | None = None) -> pd.DataFrame`

Extract cell-cycle marker columns from a DataFrame or AnnData.

Args:
    data: pandas DataFrame or AnnData.
    marker_map: Dict mapping role (e.g. ``"IdU"``) → actual column name.
    layer: AnnData layer key. ``None`` uses ``adata.X`` via ``to_df()``.

Returns:
    DataFrame indexed identically to the input, containing only the
    requested marker columns (renamed to the actual column names).

### `gate_cell_cycle(data: Any, marker_map: dict[str, str], layer: str | None = None, thresholds: dict[str, float] | None = None, quantile_thresholds: dict[str, float] | None = None, threshold_methods: dict[str, str] | None = None, apply_arcsinh: bool = False, cofactor: float = 5.0, return_adata: bool = True) -> dict[str, Any]`

Full cell-cycle gating workflow for a DataFrame or AnnData.

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

### `gate_cell_cycle_by_group(data: Any, marker_map: dict[str, str], groupby: str, layer: str | None = None, thresholds: dict[str, float] | None = None, quantile_thresholds: dict[str, float] | None = None, threshold_methods: dict[str, str] | None = None, apply_arcsinh: bool = False, cofactor: float = 5.0) -> dict[str, Any]`

Cell-cycle gating with independent threshold calculation per group.

Useful when CyTOF signal intensity varies across batches or samples.

Args:
    groupby: Column in ``adata.obs`` or DataFrame to split on.

Returns:
    Dict with ``gated_df``, ``summary``, ``per_group_summary``,
    ``per_group_thresholds``.

### `plot_cell_cycle_fractions_by_group(gated_df: pd.DataFrame, groupby: str, figsize: tuple[float, float] | None = None) -> 'plt.Figure'`

Stacked bar chart of phase fractions per group.

### `plot_cell_cycle_marker_qc(marker_df: pd.DataFrame, marker_map: dict[str, str], thresholds: dict[str, float] | None = None, n_bins: int = 80, figsize: tuple[float, float] | None = None, output_dir: str | None = None) -> list`

Plot the distribution of each marker with an optional threshold line.

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

### `plot_cell_cycle_phase_circle(data: Any, angle_col: str = 'cell_cycle_angle', phase_col: str = 'cell_cycle_phase', n_subsample: int = 5000, seed: int = 42) -> 'plt.Figure'`

Polar scatter plot of cell positions on the cell-cycle circle.

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

### `plot_cell_cycle_phase_fractions(summary: pd.DataFrame, figsize: tuple[float, float] | None = None, output_path: str | None = None) -> 'plt.Figure'`

Bar chart of cell-cycle phase fractions.

### `plot_cell_cycle_pseudotime_markers(data: Any, pseudotime_col: str = 'cell_cycle_pseudotime', phase_col: str = 'cell_cycle_phase', marker_cols: list[str] | dict[str, str] | None = None, bins: int = 50, use_hexbin: bool = True, figsize_per_marker: tuple[float, float] = (5, 3)) -> 'plt.Figure'`

Plot marker intensity vs cell-cycle pseudotime for validation.

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

### `plot_marker_thresholds(df: pd.DataFrame, marker_map: dict[str, str], thresholds: dict[str, float] | None = None, n_bins: int = 80, figsize: tuple[float, float] | None = None) -> list`

Alias for :func:`plot_cell_cycle_marker_qc` (spec compatibility).

### `summarize_cell_cycle(gated_df: pd.DataFrame) -> pd.DataFrame`

Return counts and fractions for each cell-cycle phase.

All phases in :data:`PHASE_ORDER` are included; phases with zero cells
are omitted from the output.

Returns:
    DataFrame with columns ``cell_cycle_phase``, ``n_cells``, ``fraction``.

## Classes

No public classes.
