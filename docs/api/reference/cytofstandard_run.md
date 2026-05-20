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

##### `annotate_clusters(self, cluster_key: str, annotation_map: dict[str, str], output_key: str | None = None, inplace: bool = True) -> pd.Series`

Map cluster labels to named cell types.

Adds a new obs column with the mapped names.  Unmapped labels are
passed through unchanged so partial mappings are allowed.

Args:
    cluster_key: Obs column containing cluster labels (e.g.
        ``"my_umap_leiden"``).
    annotation_map: Dict from cluster label (str) to cell-type name.
        Keys are compared against string-cast cluster values, so
        integer labels like ``0`` should be passed as ``"0"``.
    output_key: Name for the new obs column.  Defaults to
        ``"{cluster_key}_annotated"``.
    inplace: If True, persist updated AnnData to Zarr.

Returns:
    pd.Series of annotated labels indexed by cell names.

Raises:
    ValueError: If ``cluster_key`` is not in ``adata.obs``.

##### `cluster_leiden(self, embedding_name: str, cluster_key: str | None = None, resolution: float = 1.0, n_iterations: int = 2, beta: float = 0.01, objective_function: str = 'modularity', seed: int = 42, verbose: bool = False, inplace: bool = True) -> dict[str, Any]`

Cluster cells with Leiden using graph artifacts from an embedding.

Follows the Leiden settings used in the MetalUMAP notebook.
Existing clusterings with the same key are overwritten.

##### `cluster_leiden_jaccard(self, embedding_name: str, cluster_key: str | None = None, jaccard_connectivities_key: str | None = None, min_jaccard: float = 0.0, resolution: float = 1.0, n_iterations: int = 2, beta: float = 0.01, objective_function: str = 'modularity', seed: int = 42, verbose: bool = False, inplace: bool = True) -> dict[str, Any]`

Cluster cells using a PhenoGraph-style Jaccard graph with Leiden.

Equivalent to PhenoGraph but uses Leiden instead of Louvain:
    1. KNN graph (pre-computed by `compute_umap`)
    2. Edge weights replaced by Jaccard similarity of neighbor sets
    3. Leiden community detection on the Jaccard-weighted graph

The KNN indices are read from the `obsm` artifacts stored by
`compute_umap`.  No recomputation of neighbors is performed.

Args:
    embedding_name: Name of the embedding whose KNN artifacts to use
        (must have been computed via `compute_umap`).
    cluster_key: Key under which cluster labels are stored in
        `adata.obs`. Defaults to
        `f"{embedding_name}_jaccard_leiden"`.
    jaccard_connectivities_key: Key used to store the Jaccard
        connectivity matrix in `adata.obsp`. Defaults to
        `f"{embedding_name}_jaccard_connectivities"`.
    min_jaccard: Edges with Jaccard similarity below this threshold
        are pruned before clustering.  Useful to remove very weak
        connections (default 0.0 keeps all edges).
    resolution: Leiden resolution parameter.
    n_iterations: Number of Leiden iterations.
    beta: Leiden randomness parameter.
    objective_function: `"modularity"` or `"CPM"`.
    seed: Random seed passed to Leiden.
    verbose: If True, print progress to stdout.
    inplace: If True, persist updated AnnData to Zarr after
        clustering.

Returns:
    Metadata dict stored under `adata.uns["clusterings"][cluster_key]`.

##### `compare_groups(self, field: str, groupby, layer: str = 'X', comparisons: str | list[tuple[str, str]] | None = 'all', method: str = 'ttest', equal_var: bool = True, order: list[str] | None = None, multitest: str | None = 'bh') -> pd.DataFrame`

Compute pairwise group comparisons for a marker or numeric obs field.

Args:
    field: Marker name (in `adata.var_names`) or numeric obs column.
    groupby: Single obs column or list of obs columns for composite grouping.
    layer: Layer used when `field` is a marker (default `X`).
    comparisons: `"all"` (default), `"adjacent"`, explicit list of
        (group_a, group_b) pairs, or None to skip comparisons.
    method: `"ttest"` or `"wald"`.
    equal_var: When method is `"ttest"`, use pooled variance (True)
        or Welch correction (False).
    order: Optional explicit group order; otherwise sorted unique groups.
    multitest: Optional p-value correction (`"bh"` or `"bonferroni"`).

Returns:
    DataFrame with columns:
        group_a, group_b, n_a, n_b, mean_a, mean_b, var_a, var_b,
        stat, p_value, p_adj, method

##### `compute_umap(self, markers: list[str], source_layer: str = 'X', embedding_name: str = 'X_umap', n_neighbors: int = 15, n_components: int = 2, min_dist: float = 0.1, metric: str = 'euclidean', random_state: int = 42, knn_method: str = 'brute', verbose: bool = False, module_name: str = 'mlx_umap', inplace: bool = True) -> dict[str, Any]`

Compute UMAP embedding from selected markers using mlx-umap.

Also computes/stores KNN and graph artifacts required for downstream Leiden.
Existing artifacts with the same `embedding_name` are overwritten.

##### `compute_umap_balanced(self, markers: list[str], source_layer: str = 'X', embedding_name: str = 'X_umap', groupby_col: str = 'sample_id', n_per_group: int | None = None, replace: bool = False, n_neighbors: int = 15, n_components: int = 2, min_dist: float = 0.1, metric: str = 'euclidean', random_state: int = 42, knn_method: str = 'brute', verbose: bool = False, module_name: str = 'mlx_umap', inplace: bool = True) -> dict[str, Any]`

Compute UMAP with fit on balanced subsample, then transform all cells.

The balanced subsample is defined by `groupby_col`, using `n_per_group`
cells per group (or the smallest group size if None). The UMAP model
is fit on the subsample, then used to transform all cells.

KNN and graph artifacts are still computed for the full dataset to
support downstream Leiden clustering.

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

##### `differential_abundance(self, cluster_key: str, groupby: str, comparisons: str | list[tuple[str, str]] | None = 'all', method: str = 'fisher', multitest: str | None = 'bh', order: list[str] | None = None, plot: bool = False, figsize: tuple[float, float] | None = None, ax = None) -> pd.DataFrame | tuple[pd.DataFrame, tuple]`

Test whether cluster proportions differ between groups.

For each cluster and each pair of groups, the observed cell counts are
compared using Fisher's exact test or a chi-squared test.  P-values are
optionally corrected across all clusters × all pairs.

Args:
    cluster_key: Obs column with cluster labels (e.g.
        ``"my_umap_leiden"``).
    groupby: Obs column defining the groups to compare (e.g.
        ``"line_id"`` or ``"condition"``).
    comparisons: ``"all"`` (every pair), ``"adjacent"`` (neighbours in
        ``order``), explicit list of ``(group_a, group_b)`` tuples, or
        ``None`` (skip — returns proportions without p-values).
    method: ``"fisher"`` (Fisher's exact) or ``"chi2"``
        (chi-squared contingency).
    multitest: P-value correction: ``"bh"`` (Benjamini-Hochberg),
        ``"bonferroni"``, or ``None``.
    order: Explicit group order.  Defaults to sorted unique groups.
    plot: If True, also return a stacked-bar proportion plot.
    figsize: Figure size when ``plot=True``.
    ax: Existing axes to draw on (single axes, ``plot=True`` only).

Returns:
    DataFrame with columns:
        ``cluster``, ``group_a``, ``group_b``,
        ``n_a``, ``n_b``, ``N_a``, ``N_b``,
        ``freq_a``, ``freq_b``, ``p_value``, ``p_adj``

    When ``plot=True``: ``(DataFrame, (fig, ax))``.

Raises:
    ValueError: If ``cluster_key`` or ``groupby`` are not in
        ``adata.obs``.

##### `ingest(self, files: list[str], sample_metadata: str, copy_raw: bool = True, strict_markers: bool = True, allow_extra_markers: bool = False, common_markers_only: bool = False, drop_columns: list[str] | None = None) -> None`

Ingest files into this run.

Args:
    files: List of file paths to ingest
    sample_metadata: Path to sample metadata CSV/Parquet file
    copy_raw: Whether to copy raw files to project
    strict_markers: Whether to fail on unknown markers
    allow_extra_markers: Whether to allow extra markers
    common_markers_only: If True, drop markers that are not present in
        all files for this run.
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

##### `match_clusterings(self, key_a: str, key_b: str, score_mode: str = 'jaccard', n_permutations: int = 1000, random_state: int = 0) -> dict`

Compare two obs columns (clusterings or any categorical labels).

Performs both many-to-one matching (every A label gets its best B
label and vice versa) and one-to-one optimal matching via the
Hungarian algorithm.  Hypergeometric enrichment and permutation-based
significance are computed for each match.

Args:
    key_a: First obs column (e.g. ``"CL"``).
    key_b: Second obs column (e.g. ``"sample_id"`` or ``"CL_ID"``).
    score_mode: Matching criterion — ``"jaccard"`` (default),
        ``"a"`` (fraction of A in B), ``"b"`` (fraction of B from A),
        or ``"overlap"`` / ``None`` (raw count).
    n_permutations: Permutations for global significance test.
        Set to 0 to skip.
    random_state: Random seed.

Returns:
    dict with keys:

    - ``A_to_B`` — every A label matched to its best B label
    - ``B_to_A`` — every B label matched to its best A label
    - ``one_to_one_matches`` — Hungarian optimal assignment
    - ``contingency`` — raw overlap crosstab
    - ``score_matrix`` — scoring matrix used for matching
    - ``B_receives_from_A``, ``A_receives_from_B`` — split/merge structure
    - ``ari``, ``nmi`` — global agreement metrics
    - ``global_many_to_one_score_A_to_B``, ``_B_to_A``, ``_one_to_one``
    - ``permutation_p_A_to_B``, ``_B_to_A``, ``_one_to_one``

Raises:
    ValueError: If either key is not in ``adata.obs``.

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

##### `plot_boxplot(self, field: str, groupby, layer: str = 'X', order: list[str] | None = None, comparisons = None, test: str = 'mannwhitney', multitest: str | None = 'bh', show_points: bool = False, show_outliers: bool = True, max_points: int = 2000, point_alpha: float = 0.4, point_size: float = 2.0, palette = None, figsize: tuple[float, float] | None = None, ax = None, bracket_color: str = 'black', bracket_linewidth: float = 1.0, bracket_fontsize: float = 11.0, ns_label: str = 'ns', significance_thresholds: list[tuple[float, str]] | None = None, random_state: int = 0, boxplot_kwargs: dict | None = None, stripplot_kwargs: dict | None = None)`

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
    multitest: `None`, `"bh"`, or `"bonferroni"`.
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

##### `plot_cluster_composition(self, cluster_key: str, groupby: str, normalize: str = 'cluster', palette = None, figsize: tuple[float, float] | None = None, legend_kwargs: dict | None = None, ax = None) -> tuple`

Stacked bar chart of label composition.

Two modes controlled by ``normalize``:

- ``"cluster"`` *(default)* — one bar per cluster, showing the
  fraction of cells in that cluster that come from each group (e.g.
  sample).  Answers: *"what samples make up each cluster?"*
- ``"group"`` — one bar per group, showing the fraction of cells in
  that group assigned to each cluster.  Answers: *"how are each
  sample's cells distributed across clusters?"*

Args:
    cluster_key: Obs column with cluster labels.
    groupby: Obs column with group labels (e.g. ``"sample_id"``).
    normalize: ``"cluster"`` or ``"group"``.
    palette: Colour palette passed to seaborn.  ``None`` uses the
        default categorical palette.
    figsize: Figure size.  Defaults to a width proportional to the
        number of bars.
    legend_kwargs: Extra kwargs forwarded to ``ax.legend()``.
    ax: Existing axes to draw on.

Returns:
    ``(fig, ax)`` tuple.

Raises:
    ValueError: If ``cluster_key`` or ``groupby`` are not in obs,
        or ``normalize`` is not ``"cluster"`` or ``"group"``.

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

##### `subsample_by_group(self, groupby_col: str, n_per_group: int | None = None, random_state: int = 0, replace: bool = False) -> anndata.AnnData`

Return a balanced subsample AnnData by group.

Args:
    groupby_col: Obs column to balance on.
    n_per_group: Number of cells per group. If None, uses the
        smallest group size.
    random_state: RNG seed for sampling.
    replace: Sample with replacement if True.

Returns:
    Subsampled AnnData copy.

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
