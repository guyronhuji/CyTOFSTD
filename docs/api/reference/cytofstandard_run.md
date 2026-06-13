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

##### `add_cell_cycle_pseudotime(self, marker_cols: dict[str, str] | None = None, phase_col: str = 'cell_cycle_phase', output_col: str = 'cell_cycle_pseudotime', angle_col: str = 'cell_cycle_angle', phase_order: list[str] | None = None, phase_widths: dict[str, float] | None = None, overwrite: bool = False, inplace: bool = True) -> 'anndata.AnnData'`

Add a continuous cell-cycle pseudotime coordinate to this run.

Requires that :meth:`gate_cell_cycle` has already been run so that
``adata.obs["cell_cycle_phase"]`` exists.

Within each phase, ordering is based on marker intensity ranks:

* **G0/G1 phases** — ``pRb`` (increases as CDK4/6 phosphorylate Rb)
* **S phase** — ``DNA`` (content increases through replication)
* **G2 phase** — ``CyclinB1`` (accumulates before M entry)
* **M phase** — ``pH3`` (peaks during mitosis)

Args:
    marker_cols: Dict mapping role → actual column name, e.g.
        ``{"pRb": "pRb", "DNA": "DNA", "CyclinB1": "CyclinB1",
        "pH3": "pH3", "IdU": "IdU"}``.
        If ``None``, the gating marker map stored in
        ``adata.uns["cell_cycle_gating"]`` is used for the gating
        markers (IdU/pH3/CyclinB1/pRb); **DNA must still be added
        manually** if it is needed for within-S ordering.
    phase_col: Column in ``adata.obs`` with categorical phase labels.
    output_col: Output column name for pseudotime.
    angle_col: Output column name for the 2π angle.
    phase_order: Custom biological ordering of phase labels.
    phase_widths: Dict of phase → fractional arc width (summed to 1).
    overwrite: Overwrite existing pseudotime columns. Default False.
    inplace: Persist updated AnnData to disk. Default True.

Returns:
    Updated AnnData object (with pseudotime columns in ``obs``).

Example::

    run.add_cell_cycle_pseudotime(
        marker_cols={
            "pRb": "pRb", "IdU": "IdU",
            "CyclinB1": "CyclinB1", "pH3": "pH3",
            "DNA": "DNA",
        },
        overwrite=True,
    )

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

##### `bootstrap_rmt_spectrum(self, markers = None, layer = None, groupby = None, matrix = 'marker_cov', standardize = True, n_cells = 10000, frac = 0.8, n_bootstrap = 200, random_state = 0, uns_key = 'rmt_bootstrap', plot = False, inplace = True)`

Bootstrap stability of the RMT eigenvalue spectrum.

Repeatedly subsamples a fraction of cells and recomputes the
eigenvalue spectrum to quantify estimation uncertainty.

- **Eigenvalue distributions**: spread of each rank's eigenvalue
  across resamples — tight = well-estimated, wide = noisy.
- **Eigenvector stability** (``matrix="marker_cov"`` only): mean
  absolute cosine similarity between the reference eigenvector and
  each bootstrap replicate.  Near 1 = reproducible; near 0 = unstable.

These complement the MP boundary from :meth:`compute_rmt_spectrum`:
RMT says "above the noise floor"; bootstrap says "stably estimated".

Args:
    markers: Subset of marker names to use.  ``None`` uses all.
    layer: AnnData layer to use.  ``None`` uses ``adata.X``.
    groupby: ``adata.obs`` column to group by.
    matrix: ``"marker_cov"`` or ``"cell_gram"``.
    standardize: Centre and scale each marker before analysis.
    n_cells: Maximum cell pool size per group before bootstrapping.
    frac: Fraction of the pool drawn without replacement per replicate.
    n_bootstrap: Number of bootstrap resamples.
    random_state: Seed for the RNG.
    uns_key: Key under which results are stored in ``adata.uns``.
    plot: If ``True``, return ``(result, fig)`` with a two-panel
        figure: eigenvalue distributions and eigenvector stability.
    inplace: If ``True``, persist the updated AnnData to disk.

Returns:
    Result dict (or ``(dict, Figure)`` when ``plot=True``).  Contains
    ``eigenvalue_ref``, ``eigenvalue_matrix`` (n_bootstrap x p),
    ``eigenvalue_mean``, ``eigenvalue_std``, ``eigenvalue_ci_low``,
    ``eigenvalue_ci_high``, ``eigenvector_stability`` (marker_cov
    only), ``lambda_max_mp_ref``, ``lambda_max_mp_distribution``,
    ``n_signal_ref``, ``markers``, ``matrix``.
    With ``groupby``, nested under ``"groups"``.

Raises:
    ValueError: Invalid parameter values or missing data.
    RunNotIngestedError: Run has not been ingested.

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

##### `compute_flowsom(self, markers: list[str], source_layer: str = 'X', grid_size: int | tuple[int, int] = 10, n_iterations: int = 100, neighborhood_radius: int = 1, learning_rate: float = 0.5, cluster_method: str = 'leiden', cluster_resolution: float = 1.0, n_meta_clusters: int | None = None, random_state: int = 42, verbose: bool = False, inplace: bool = True, use_gpu: bool = True) -> dict[str, Any]`

Compute FlowSOM clustering.

Uses the external ``flowsom`` package (no local reimplementation).
Results are copied back into this run's AnnData object:

- ``adata.obs['X_flowsom']``: FlowSOM metaclusters
- ``adata.obsm['X_flowsom_node']``: SOM node assignment per cell
- ``adata.varm['flowsom_weights']``: SOM codebook (markers x nodes)

Args:
    markers: Markers to use for clustering (must exist in `adata.var_names`).
    source_layer: Layer used for expression values (`X` or layer key).
    grid_size: SOM grid size (10x10) or (rows, cols) tuple.
    n_iterations: FlowSOM training passes (`rlen`).
    neighborhood_radius: Kept for compatibility (stored in metadata).
    learning_rate: Initial FlowSOM learning rate (`alpha[0]`).
    cluster_method: Kept for compatibility. FlowSOM metaclustering is
        always used as final labels.
    cluster_resolution: Backward-compatible hint used to derive default
        ``n_meta_clusters`` when not provided.
    n_meta_clusters: Explicit number of FlowSOM metaclusters.
    random_state: RNG seed.
    verbose: If True, print progress to stdout.
    inplace: If True, persist results to run zarr.
    use_gpu: If True, request MPS when available (best-effort; only
        used when the installed FlowSOM backend exposes a ``device``
        parameter).

Returns:
    Metadata dict with FlowSOM metadata stored under
    `adata.uns["clusterings"]["X_flowsom"]` and `adata.uns["flowsom"]`.

##### `compute_rmt_spectrum(self, markers = None, layer = None, groupby = None, matrix = 'marker_cov', standardize = True, n_cells = 10000, sigma_sq = None, random_state = 0, uns_key = 'rmt_spectrum', plot = False, inplace = True)`

Compute the marker eigenvalue spectrum and compare to the Marchenko-Pastur null.

Uses Random Matrix Theory (RMT) to identify which principal components
of the marker covariance (or cell Gram) matrix represent genuine
biological signal vs. sampling noise.  Eigenvalues above the MP upper
bound ``lambda_max = sigma2 * (1 + sqrt(q))^2`` are classed as signal.

The Marchenko-Pastur law states that for a random n x p matrix with
i.i.d. entries of variance sigma2, the eigenvalues of the sample
covariance ``C = X.T @ X / n`` lie within
``[sigma2*(1-sqrt(q))^2, sigma2*(1+sqrt(q))^2]`` with ``q = p/n``.

Args:
    markers: Subset of marker names to use.  ``None`` uses all markers.
    layer: AnnData layer to use.  ``None`` uses ``adata.X``.
    groupby: ``adata.obs`` column to group by.  When set, the spectrum
        is computed separately per group.
    matrix: ``"marker_cov"`` (default) builds the p x p marker
        covariance matrix.  ``"cell_gram"`` builds the n x n cell Gram
        matrix ``X @ X.T / p`` on a subsampled cell set.
    standardize: If ``True``, centre and scale each marker to zero mean
        and unit variance before analysis.  For ``"marker_cov"`` this
        yields a correlation matrix (sigma2 = 1 under the null).
    n_cells: Maximum cells per group.  Randomly subsampled when larger.
    sigma_sq: Noise variance for the MP distribution.  ``None``
        auto-estimates as ``trace(C) / p``.
    random_state: Seed for subsampling.
    uns_key: Key under which results are stored in ``adata.uns``.
    plot: If ``True``, return ``(result, fig)`` with a scree plot.
    inplace: If ``True``, persist the updated AnnData to disk.

Returns:
    Result dict (or ``(dict, Figure)`` when ``plot=True``).  Contains
    ``eigenvalues``, ``lambda_max_mp``, ``lambda_min_mp``,
    ``n_signal``, ``q``, ``sigma_sq``, ``n``, ``p``, ``markers``,
    ``matrix``.  With ``groupby``, nested under ``"groups"``.

Raises:
    ValueError: Invalid parameter values or missing data.
    RunNotIngestedError: Run has not been ingested.

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

##### `gate_cell_cycle(self, marker_map: dict[str, str] | None = None, layer: str | None = None, thresholds: dict[str, float] | None = None, quantile_thresholds: dict[str, float] | None = None, threshold_methods: dict[str, str] | None = None, inplace: bool = True) -> dict[str, Any]`

Assign cells to cell-cycle phases using hierarchical marker gating.

Uses arcsinh-transformed data (``adata.X`` by default, or a named
layer). Do **not** pass the ``zscore`` layer — gating thresholds are
calibrated for arcsinh-scale values.

Gating hierarchy (mutually exclusive, applied in order):
  IdU+      → S_phase
  pH3+      → M_phase
  CyclinB1+ → G2_phase
  pRb+      → Cycling_G1   (if pRb in marker_map)
  pRb-      → G0_or_quiescent
  remaining → G1_or_quiescent (when pRb absent)

Args:
    marker_map: Dict mapping role → actual ``var_name`` column.
        Required roles: ``IdU``, ``pH3``, ``CyclinB1``.
        Optional: ``pRb``. If ``None``, auto-detection is attempted
        from ``adata.var_names`` using known aliases.
    layer: AnnData layer to read expression from. ``None`` uses
        ``adata.X`` (typically arcsinh-transformed after ingestion).
        Avoid ``"zscore"`` — thresholds are not calibrated for it.
    thresholds: User-supplied thresholds by role. Missing roles use
        the automatic method.
    quantile_thresholds: Override the default quantile per role
        (used when strategy is ``"quantile"``).
    threshold_methods: Per-role auto-threshold strategy.
        ``"otsu"`` (default for IdU/pH3/CyclinB1) finds the valley in
        bimodal distributions. ``"quantile"`` (default for pRb) uses
        a fixed percentile.
    inplace: Persist updated AnnData to disk after gating.

Returns:
    Dict with:
    - ``gated_df``: cell-level DataFrame with gate columns and
      ``cell_cycle_phase``.
    - ``summary``: phase counts and fractions.
    - ``thresholds``: thresholds used (role → value).

Raises:
    ValueError: If required markers cannot be found in ``adata.var_names``.

##### `ingest(self, files: list[str], sample_metadata: str, copy_raw: bool = True, strict_markers: bool = True, allow_extra_markers: bool = False, common_markers_only: bool = False, drop_columns: list[str] | None = None, show_marker_coverage: bool = False) -> None`

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
    show_marker_coverage: If True, display a heatmap after ingestion
        showing which standard markers are present or absent in each
        file.  Only markers that are missing from at least one file are
        shown; if all markers are present in all files the plot is
        skipped and a message is printed.

Raises:
    MetadataValidationError if metadata validation fails
    MarkerValidationError if marker validation fails
    IngestionError if ingestion fails

##### `ingestion_summary(self) -> pd.DataFrame`

Return a DataFrame summarising every file and sample ingested into this run.

Columns always present:

- ``file_name`` — original filename
- ``sample_id`` — sample identifier from the metadata sheet
- ``line_id`` — cell-line identifier
- ``n_events_file`` — cells read from the file at ingest time
- ``n_cells_obs`` — cells currently in the zarr (after any QC gating)
- ``file_type`` — ``"csv"``, ``"fcs"``, ``"parquet"``, …
- ``file_size_bytes`` — original file size
- ``ingested_at`` — UTC timestamp of ingestion

Additional metadata columns from the sample metadata sheet (e.g.
``"condition"``, ``"replicate_id"``) are included when available.

Raises:
    RunNotIngestedError: If the run has not been ingested yet.

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

##### `locked_zarr_parts(self) -> list[str]`

Return which zarr store parts are currently locked (read-only).

Scans the run's zarr store and returns the logical parts — such as
``"obs"``, ``"layers/raw"`` — that contain any read-only files.
An empty list means the store is fully writable.

Returns:
    Sorted list of locked part paths.  ``["."]`` means the entire
    store is locked.

Raises:
    RunNotIngestedError: If the run has not been ingested yet.

##### `match_clusterings(self, key_a: str, key_b: str, score_mode: str = 'jaccard', n_permutations: int = 1000, random_state: int = 0, plot: str | bool = False, matrix: str = 'score', annot: bool = True, figsize: tuple[float, float] | None = None) -> dict`

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
    plot: ``False`` (default, no plot), ``"heatmap"`` (seaborn heatmap,
        no dendrogram), or ``"clustermap"`` (seaborn clustermap with
        row/column dendrograms).
    matrix: Which values to display — ``"score"`` (default, the
        scoring matrix used for matching) or ``"overlap"`` (raw cell
        counts from the contingency table).
    annot: Annotate each cell with its value.  Default ``True``.
    figsize: Figure size passed to the plot.  When ``None`` a
        sensible default is chosen based on the matrix dimensions.

Returns:
    When ``plot=False``: the report ``dict``.
    When ``plot`` is set: ``(report_dict, figure)`` where *figure* is
    the ``matplotlib.figure.Figure`` for a heatmap or the
    ``seaborn.ClusterGrid`` for a clustermap.

    Report dict keys:

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
    ValueError: If either key is not in ``adata.obs``, or ``plot``/
        ``matrix`` are invalid values.

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

##### `open_cell_cycle_app(self, port: int = 8502) -> 'subprocess.Popen'`

Launch the cell-cycle gating Streamlit app pre-loaded with this run.

Opens http://localhost:<port> in the browser automatically.
The app shows Auto-threshold and Interactive-slider tabs — click
"Apply Gating" to write the chosen thresholds back to the stored adata.

Afterwards, call the Section 9 cell in the notebook to pick up the
saved thresholds from ``adata.uns["cell_cycle_gating"]["latest"]``.

Args:
    port: Local port to serve on (default 8502 to avoid clashing with
          the full CyTOF app on 8501).

Returns:
    The ``subprocess.Popen`` handle — call ``proc.terminate()`` to stop.

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

##### `read_adata(self, backed: bool = False, force: bool = False)`

Read the ingested AnnData object.

Args:
    backed: Whether to use backed mode (currently ignored)
    force: Re-read from disk even if an in-memory copy is cached.
        Use this after external modifications (e.g. gating via the
        Streamlit app) to pick up changes written by another process.

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

##### `save_as_new_run(self, adata: anndata.AnnData, new_run_id: str, run_name: str | None = None, notes: str | None = None) -> 'Run'`

Persist an arbitrary AnnData as a brand-new run in the same project.

Unlike ``save_adata(adata)``, which overwrites the *current* run, this
method registers a fresh run and writes the provided AnnData to it.
Panel metadata (panel_id, instrument, …) is copied from the source run.

Typical use case — save a balanced subsample without destroying the full
run's data::

    adata_sub = run.subsample_by_group("sample_id", n_per_group=5000)
    sub_run = run.save_as_new_run(adata_sub, new_run_id="CyTOF1_sub")

Args:
    adata: AnnData to persist (e.g. from ``subsample_by_group``).
    new_run_id: ID for the new run.  Must not already exist in the project.
    run_name: Optional display name; defaults to
        ``"<source run name> (subset)"``.
    notes: Optional free-text notes attached to the new run's metadata.

Returns:
    The newly created and saved ``Run`` object.

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

##### `vendi_score(self, k: int = 15, markers: list[str] | None = None, layer: str | None = None, use_rep: str | None = None, metric: str = 'cosine', obs_key: str = 'vendi_score', groupby: str | list[str] | None = None, sigma: float | None = None, n_bins: int = 10, n_reps: int = 1, m: int | None = None, random_state: int = 0, return_eigenvalues: bool = False, inplace: bool = True) -> np.ndarray | pd.DataFrame | tuple[pd.DataFrame, dict[str, np.ndarray]]`

Compute Vendi score over local k-NN neighborhoods or whole groups.

**Per-cell mode** (``groupby=None``): for each cell i, features of its
k-nearest-neighbor neighborhood N_k(i) are binned with
``KBinsDiscretizer`` (uniform) then scored with
``vendi_score.vendi.score_dual``.  With rarefaction (``n_reps > 1``),
results are bootstrap-averaged and 95 % CI bounds are also stored in
``adata.obs`` as ``{obs_key}_ci_low`` / ``{obs_key}_ci_high``.

**Per-group mode** (``groupby`` set): for each unique value (or
combination of values when ``groupby`` is a list) in the chosen
``adata.obs`` column(s), all cells belonging to that group are
rarefied, binned, and scored.  Returns a DataFrame indexed by group
label (MultiIndex when ``groupby`` is a list) with columns
``vendi_score``, ``ci_low``, ``ci_high``, stored in
``adata.uns["vendi"][obs_key]``.

Args:
    k: Number of nearest neighbors per cell (per-cell mode only).
    markers: Subset of marker names (``adata.var_names``) to use as
        features. ``None`` uses all markers. Ignored when ``use_rep``
        is set.
    layer: AnnData layer to use as the feature matrix. ``None`` uses
        ``adata.X``.
    use_rep: Key in ``adata.obsm`` to use instead of marker expression
        (e.g. ``"X_umap"``). When provided, ``layer`` is ignored.
    metric: ``"cosine"`` uses ``score_dual`` on binned features;
        ``"rbf"`` builds a Gaussian kernel matrix and calls
        ``score_K``.
    obs_key: Storage key: ``adata.obs`` column (per-cell) or
        ``adata.uns["vendi"]`` key (per-group).
    groupby: ``adata.obs`` column(s) to group by (e.g. ``"leiden"``,
        ``["condition", "cell_cycle"]``). When a list is given, one
        score is produced per unique combination of values.
    sigma: Bandwidth for the RBF kernel. Defaults to the median
        pairwise distance in the data.
    n_bins: Number of uniform bins for ``KBinsDiscretizer``.
    n_reps: Bootstrap repetitions. ``1`` scores once without
        bootstrapping.
    m: Subsample size per repetition. Defaults to pool size
        (``k + 1`` per-cell, group size per-group).
    random_state: Seed for the bootstrap RNG.
    return_eigenvalues: Per-group mode only. When ``True``, also
        compute the eigenvalue spectrum of the dual kernel for each
        group (using the full binned pool, no subsampling) and return
        ``(DataFrame, {group: eigenvalues_array})``. Eigenvalues are
        also stored in ``adata.uns["vendi"][obs_key + "_eigenvalues"]``.
    inplace: If ``True``, persist the updated AnnData to disk.

Returns:
    Per-cell: ``np.ndarray`` of shape ``(n_cells,)`` stored in
    ``adata.obs[obs_key]``.

    Per-group (``return_eigenvalues=False``): ``pd.DataFrame`` indexed
    by group label, stored in ``adata.uns["vendi"][obs_key]``.

    Per-group (``return_eigenvalues=True``): ``(DataFrame, dict)``
    where the dict maps each group label to a sorted
    ``np.ndarray`` of eigenvalues (ascending).

Raises:
    ValueError: If ``metric`` is invalid, ``groupby`` column is
        missing, or the requested ``layer``/``use_rep`` does not exist.
    RunNotIngestedError: If the run has not been ingested yet.

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
