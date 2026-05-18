# `cytofstandard.run`

- Source: `cytofstandard/run.py`

Run class for cytofstandard.

## Public Exports (`__all__`)

- `Run`

## Top-level Functions

No public top-level functions.

## Classes

### `Run`

CytOF run for managing ingestion and data access.

#### Methods

##### `cluster_leiden(self, embedding_name: str, cluster_key: str | None = None, resolution: float = 1.0, n_iterations: int = 2, beta: float = 0.01, objective_function: str = 'modularity', seed: int = 42, verbose: bool = False, inplace: bool = True) -> dict[str, Any]`

Cluster cells with Leiden using graph artifacts from an embedding.

Follows the Leiden settings used in the MetalUMAP notebook.
Existing clusterings with the same key are overwritten.

##### `compute_umap(self, markers: list[str], source_layer: str = 'X', embedding_name: str = 'X_umap', n_neighbors: int = 15, n_components: int = 2, min_dist: float = 0.1, metric: str = 'euclidean', random_state: int = 42, knn_method: str = 'brute', verbose: bool = False, module_name: str = 'mlx_umap', inplace: bool = True) -> dict[str, Any]`

Compute UMAP embedding from selected markers using mlx-umap.

Also computes/stores KNN and graph artifacts required for downstream Leiden.
Existing artifacts with the same `embedding_name` are overwritten.

##### `create_subset_run(self, new_run_id: str, sample_ids: list[str] | None = None, line_ids: list[str] | None = None, run_name: str | None = None, notes: str | None = None) -> 'Run'`

Create a new run from a subset of samples and/or lines in this run.

Args:
    new_run_id: New run ID to create.
    sample_ids: Sample IDs to keep.
    line_ids: Line IDs to keep.
    run_name: Optional name for new run.
    notes: Optional notes for new run metadata.

Returns:
    Newly created and persisted subset Run.

##### `ingest(self, files: list[str], sample_metadata: str, copy_raw: bool = True, strict_markers: bool = True, allow_extra_markers: bool = False, drop_columns: list[str] | None = None) -> None`

Ingest files into this run.

Args:
    files: List of file paths to ingest
    sample_metadata: Path to sample metadata CSV/Parquet file
    copy_raw: Whether to copy raw files to project
    strict_markers: Whether to fail on unknown markers
    allow_extra_markers: Whether to allow extra markers
    drop_columns: Column names to remove before marker processing

Raises:
    MetadataValidationError if metadata validation fails
    MarkerValidationError if marker validation fails
    IngestionError if ingestion fails

##### `is_ingested(self) -> bool`

Check if the run has been ingested.

Returns:
    True if run is ingested

##### `lock_zarr_parts(self, parts: list[str] | None = None, strict: bool = True) -> list[str]`

Make selected Zarr store parts read-only.

Args:
    parts: Relative paths in the run Zarr store (for example,
        ``["layers/raw", "obs"]``). If None, lock the full store.
    strict: Whether to raise if a requested part does not exist.

Returns:
    Normalized part paths that were locked. ``"."`` means full store.

##### `normalize_with_cytof_transform(self, control_markers: list[str], markers_to_correct: list[str], source_layer: str = 'raw', corrected_layer: str = 'normalized', z_layer: str = 'normalized_z', groupby_col: str = 'sample_id', input_is_arcsinh: bool = False, arcsinh_cofactor: float = 5.0, anchor_to_median: bool = True, zscore: bool = True, module_name: str = 'cytof_transform', inplace: bool = True) -> dict[str, Any]`

Normalize markers using external cytof_transform, per sample/line group.

Args:
    control_markers: Core markers (e.g. histones) used to estimate technical factor.
    markers_to_correct: Markers to normalize.
    source_layer: Layer used as normalization input.
    corrected_layer: Output layer for corrected (asinh-space) values.
    z_layer: Output layer for z-scored corrected values.
    groupby_col: Obs column used for per-group normalization (e.g. sample_id/line_id).
    input_is_arcsinh: If True, source layer is already arcsinh-transformed.
    arcsinh_cofactor: Cofactor used when transforming source data with arcsinh.
    anchor_to_median: Passed to cytof_transform config.
    zscore: Passed to cytof_transform config.
    module_name: Module name for importing cytof_transform.
    inplace: If True, persist updates to run zarr.

Returns:
    Summary dictionary of normalization outputs and settings.

##### `permcell_to_adata(self, result_prefix: str = 'permcell', score: str = 'z', smoothed: bool = True, embedding_key: str | None = None) -> anndata.AnnData`

Build an AnnData view of PermCell results.

Args:
    result_prefix: Prefix used in `run_permcell`.
    score: One of: z, p, zdir, zabs, pabs.
    smoothed: If True use smoothed PermCell results, else raw results.
    embedding_key: Optional embedding to copy into returned `.obsm`.
        If None, uses the run's stored positions_key for that PermCell run.

Returns:
    AnnData with:
      - X: selected PermCell score matrix
      - var_names: signature names
      - obs: copied run obs
      - obsm: selected embedding (if available)

##### `plot_boxplot(self, field: str, groupby, layer: str = 'X', order: list[str] | None = None, comparisons = None, test: str = 'mannwhitney', multitest: str | None = 'holm', show_points: bool = False, show_outliers: bool = True, max_points: int = 2000, point_alpha: float = 0.4, point_size: float = 2.0, palette = None, figsize: tuple[float, float] | None = None, ax = None, bracket_color: str = 'black', bracket_linewidth: float = 1.0, bracket_fontsize: float = 11.0, ns_label: str = 'ns', significance_thresholds: list[tuple[float, str]] | None = None, random_state: int = 0, boxplot_kwargs: dict | None = None, stripplot_kwargs: dict | None = None)`

Boxplot for a single field grouped by an obs column, with significance brackets.

The field can be a marker (in `adata.var_names`) or a numeric obs
column. `groupby` can be a single obs column or a list of obs columns
for composite grouping.

Args:
    field: Marker name or numeric obs column to plot.
    groupby: Single obs column or list of obs columns.
    layer: Layer used when `field` is a marker (default `X`).
    order: Optional explicit group order along the x-axis.
    comparisons: Pairs to test. Accepts:
        - `None` (default): no significance brackets are drawn.
        - `"all"`: all unordered pairs.
        - `"adjacent"`: only neighbours in `order`.
        - List of `(group_a, group_b)` tuples for explicit pairs.
    test: `"mannwhitney"` (default), `"ttest"`, or `"welch"`.
    multitest: `None`, `"holm"`, or `"bonferroni"`.
    show_points: Overlay a stripplot of per-cell points (subsampled).
    show_outliers: Whether to render boxplot outlier markers
        (maps to seaborn's `showfliers`). Defaults to True. Set to
        False to clean up dense plots (or pass `"showfliers": False`
        via `boxplot_kwargs`, which takes precedence).
    max_points: Maximum number of stripplot points (subsampled).
    point_alpha: Alpha for overlaid points.
    point_size: Size for overlaid points.
    palette: Seaborn palette name, color list, or dict. When None,
        a distinct color per group is generated automatically using
        the `"tab10"` palette.
    figsize: Optional figure size.
    ax: Optional matplotlib axes to draw into.
    bracket_color: Color used for significance brackets.
    bracket_linewidth: Line width used for significance brackets.
    bracket_fontsize: Font size used for significance labels.
    ns_label: Label used for non-significant comparisons.
    significance_thresholds: Ordered list of `(p_threshold, label)`
        tuples used to convert each (adjusted) p-value into a star
        label. The first tuple whose `p_threshold` is `>= p` wins.
        Provide thresholds from strictest to loosest, for example
        `[(1e-4, "****"), (1e-3, "***"), (1e-2, "**"), (5e-2, "*")]`
        (the default). Any p-value larger than every threshold is
        labelled with `ns_label`.
    random_state: Seed used when subsampling points.
    boxplot_kwargs: Extra keyword arguments forwarded to
        `seaborn.boxplot` (e.g. `width`, `linewidth`, `notch`,
        `whis`, `saturation`, `showfliers`). Keys here override the
        explicit arguments above.
    stripplot_kwargs: Extra keyword arguments forwarded to
        `seaborn.stripplot` when `show_points=True` (e.g. `jitter`,
        `dodge`).

Returns:
    Tuple of (figure, axes).

##### `plot_heatmap(self, fields: list[str], groupby, layer: str = 'X', agg: str = 'mean', standard_scale: str | None = None, cmap: str | None = None, center: float | None = None, annot: bool = False, fmt: str = '.2f', order: list[str] | None = None, figsize: tuple[float, float] | None = None, ax = None, heatmap_kwargs: dict | None = None)`

Plot a heatmap of aggregated field values grouped by an obs column.

Fields can be markers (in `adata.var_names`) or numeric `adata.obs`
columns. `groupby` can be a single obs column or a list of obs columns
for composite grouping (e.g. `["line_id", "condition"]`).

Args:
    fields: Markers or numeric obs columns to include as heatmap rows.
    groupby: Single obs column name or a list of obs column names.
    layer: Layer used when a field is a marker (default `X`).
    agg: `"mean"` or `"median"`.
    standard_scale: `None`, `"row"`, or `"column"`. Applies z-score over
        rows (fields) or columns (groups) after aggregation.
    cmap: Matplotlib colormap. Defaults to `viridis` (raw) or `RdBu_r`
        (when `standard_scale` is set).
    center: Value at which to center the colormap. Defaults to 0 when
        `standard_scale` is set, otherwise None.
    annot: If True, write the aggregated value in each cell.
    fmt: Format string used when `annot=True`.
    order: Optional explicit group order along the x-axis.
    figsize: Optional figure size.
    ax: Optional matplotlib axes to draw into.
    heatmap_kwargs: Extra keyword arguments forwarded to
        `seaborn.heatmap` (e.g. `linewidths`, `linecolor`,
        `cbar_kws`, `vmin`, `vmax`). Explicit arguments above take
        precedence over keys in this dict.

Returns:
    Tuple of (figure, axes).

##### `plot_marker_histograms(self, markers: list[str], layer: str = 'X', cofactor: float = 5.0, fill: bool = False, stat: str = 'density', element: str = 'step', bins: int | str = 'auto')`

Plot per-sample histograms for selected markers.

Args:
    markers: Marker names to plot (must exist in `adata.var_names`).
    layer: Expression layer to plot from (`X` or a layer key).
    cofactor: Cofactor used in arcsinh transform.
    fill: Passed to seaborn.histplot.
    stat: Passed to seaborn.histplot.
    element: Passed to seaborn.histplot.
    bins: Passed to seaborn.histplot.

Returns:
    Tuple of (figure, axes) from matplotlib.

##### `plot_normalization_gamma_qc(self, group_value: str, marker_groups: dict[str, list[str]] | None = None, module_name: str = 'cytof_transform')`

Plot gamma values for a normalized group using cytof_transform.plot_gamma_qc.

##### `plot_normalization_marker_correlations_qc(self, pre_layer: str = 'raw', post_layer: str = 'normalized', group_value: str | None = None, groupby_col: str = 'sample_id', input_pre_is_arcsinh: bool = False, arcsinh_cofactor: float = 5.0, top_n: int = 25, module_name: str = 'cytof_transform')`

Plot marker-tech correlations pre/post normalization for a run or group.

##### `plot_normalization_tech_factor_qc(self, control_markers: list[str], layer: str = 'raw', group_value: str | None = None, groupby_col: str = 'sample_id', input_is_arcsinh: bool = False, arcsinh_cofactor: float = 5.0, module_name: str = 'cytof_transform')`

Plot technical-factor QC using cytof_transform.plot_tech_factor_qc.

##### `qc_gate(self, gates: dict[str, Any], layer: str = 'X', inplace: bool = True) -> pd.Series`

Apply marker QC gates and optionally persist filtered cells.

Args:
    gates: Dict mapping marker -> gate specification.
        Gate spec can be:
        - {'lower': value, 'upper': value}
        - (lower, upper)
        where each bound is numeric, None, or percentile string (e.g., 'p1').
    layer: Expression layer used for gating.
    inplace: If True, persist filtered AnnData back to run zarr path.

Returns:
    Boolean pass mask indexed by original `obs_names`.

##### `read_adata(self, backed: bool = False)`

Read the ingested AnnData object.

Args:
    backed: Whether to use backed mode (currently ignored)

Returns:
    AnnData object

Raises:
    RunNotIngestedError if run has not been ingested

##### `rename(self, new_run_name: str) -> None`

Rename this run via project metadata (run_name only).

##### `require_ingested(self) -> None`

Require that the run has been ingested.

Raises:
    RunNotIngestedError if run has not been ingested

##### `run_permcell(self, signatures: dict[str, object], source_layer: str = 'X', positions_key: str = 'X_umap', result_prefix: str = 'permcell', smoothed_key: str | None = None, compute_unsmoothed: bool = True, bandwidth: float = -1, k: int | None = 64, radius: float | None = None, chunk_size: int = 2048, device: str | None = None, n_perm: int = 2000, seed: int = 0, exclude_set: bool = True, two_sided: bool = False, abs_variant: bool = True, exact_max_combinations: int = 50000, progress: bool = True, normalize_set_weights: str | None = None, use_sparse_W: bool = False, prefer_permutation: bool = True, perm_batch: int = 1024, permcell_module: Any | None = None, module_name: str = 'PermCell_Smooth', module_path: str | None = None, inplace: bool = True) -> dict[str, pd.DataFrame]`

Run PermCell smoothing + scoring and store results in AnnData.

The permanent outputs are stored in obsm/uns (no layers are modified).
Existing outputs with the same `result_prefix` are overwritten.

Returns:
    Dict of result DataFrames keyed by: z, p, zdir (and zabs/pabs when enabled).

##### `save(self, adata: anndata.AnnData | None = None) -> None`

Persist current run AnnData to disk.

Use this after external modifications to `run.read_adata()` output.

Args:
    adata: Optional AnnData object to save. If omitted, saves the current
        in-memory `self._adata`.

Raises:
    RunNotIngestedError: If run has no zarr and no adata was provided.

##### `save_adata(self, adata: anndata.AnnData | None = None) -> None`

Alias for `save()` for explicit naming in notebooks.

##### `set_x_from_layer(self, layer: str, inplace: bool = True) -> None`

Overwrite `adata.X` with values from a selected layer.

This operation is intentionally direct: no backup/history is created.

Args:
    layer: Source layer name (or "X").
    inplace: If True, persist the updated AnnData to run zarr.

##### `status(self) -> str`

- Decorators: `property`

Get the run status.

Returns:
    Status string (registered, ingested, failed_ingestion)

##### `to_dataframe(self, fields: list[str], layer: str = 'X') -> pd.DataFrame`

Return a DataFrame from selected marker and obs fields.

Args:
    fields: Ordered list of field names to include. Each name must be
        either a marker in `adata.var_names` or an `adata.obs` column.
    layer: Expression layer used for marker values (`X` or layer key).

Returns:
    DataFrame indexed by `adata.obs_names` with columns in `fields` order.

##### `unlock_zarr_parts(self, parts: list[str] | None = None, strict: bool = True) -> list[str]`

Make selected Zarr store parts owner-writable again.

Args:
    parts: Relative paths in the run Zarr store (for example,
        ``["layers/raw", "obs"]``). If None, unlock the full store.
    strict: Whether to raise if a requested part does not exist.

Returns:
    Normalized part paths that were unlocked. ``"."`` means full store.

##### `zarr_path(self) -> Path`

Get the path to the Zarr file.

Returns:
    Path to Zarr file

##### `zscore_markers_balanced(self, source_layer: str = 'normlized', output_layer: str = 'zscore', groupby_col: str = 'sample_id', random_state: int = 0, inplace: bool = True) -> dict[str, Any]`

Z-score all markers using a balanced subsample across sample IDs.

If multiple groups exist in `groupby_col`, z-score parameters (mean/std)
are estimated from an equal-size subsample per group to avoid bias.

Args:
    source_layer: Layer to z-score (`X` or layer key).
    output_layer: Target layer name for z-scored values.
    groupby_col: Obs column used for balancing (default: sample_id).
    random_state: RNG seed for balanced subsampling.
    inplace: If True, persist updated AnnData to run zarr.

Returns:
    Summary dictionary with balancing and z-score metadata.
