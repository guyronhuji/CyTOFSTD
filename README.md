# CyTOF Standard Package

Standard CyTOF analysis package - Phase 1: Ingestion and storage.

## Installation

```bash
pip install anndata zarr numpy pandas pyyaml fcsparser
```

## Quick Start

```python
from cytofstandard import Project

# Create a new project
project = Project.create(
    path="my_project",
    project_id="BRCA_CYTOF_2026",
    project_name="BRCA CyTOF histone panel",
    standard_marker_file="standard_markers.csv",
    marker_alias_file="marker_aliases.yaml",
)

# Register a run
run = project.add_run(
    run_id="run_001",
    run_name="First BRCA CyTOF run",
    panel_id="breast_histone_panel_v1",
    acquisition_date="2026-05-15",
    instrument="Helios",
    operator="GR",
)

# Ingest data
run.ingest(
    files=["sample_A.fcs", "sample_B.fcs"],
    sample_metadata="sample_metadata.csv",
    copy_raw=True,
    strict_markers=True,
)

# Load and use the data
adata = run.read_adata()
print(adata)

# If you modify adata externally, persist back to run storage
adata.obs["my_flag"] = "external"
run.save()  # or run.save_adata(adata)

# Rename run (run_name only)
run.rename("My renamed run")

# Lock/unlock selected Zarr parts on disk
run.lock_zarr_parts(parts=["layers/raw", "obs"])
run.unlock_zarr_parts(parts=["layers/raw", "obs"])

# Lock/unlock the full Zarr store
run.lock_zarr_parts()
run.unlock_zarr_parts()
```

## Normalization (cytof_transform)

Normalization calls the external `cytof_transform` module (it is not vendored into this package).

```python
summary = run.normalize_with_cytof_transform(
    control_markers=["H3.3", "H3", "H4"],
    markers_to_correct=["H3K27ac", "H3K4me3", "H3K9ac", "ER", "KI67"],
    source_layer="raw",
    groupby_col="sample_id",  # or "line_id"
    input_is_arcsinh=False,
    arcsinh_cofactor=5.0,
)

adata = run.read_adata()
print(adata.layers.keys())  # includes 'normalized', 'normalized_z'
```

QC plotting wrappers for normalization:

```python
run.plot_normalization_tech_factor_qc(
    control_markers=["H3.3", "H3", "H4"],
    layer="raw",
    group_value="S001",
)

run.plot_normalization_marker_correlations_qc(
    pre_layer="raw",
    post_layer="normalized",
    group_value="S001",
)

run.plot_normalization_gamma_qc(group_value="S001")
```

## Embeddings and Clustering

```python
# Set X from a selected layer (direct overwrite, no backup)
run.set_x_from_layer("normalized")

# Compute UMAP with mlx-umap from selected markers and layer
run.compute_umap(
    markers=["H3", "H3K27me3", "ECad", "EpCAM"],
    source_layer="normalized",
    embedding_name="norm_umap",
    n_neighbors=15,
    min_dist=0.1,
    verbose=True,
)

# Cluster with Leiden from stored graph for that embedding
run.cluster_leiden(
    embedding_name="norm_umap",
    resolution=1.0,
    verbose=True,
)

# Labels are stored in obs as: "norm_umap_leiden"
```

## PermCell

Run PermCell smoothing + signature scoring on a chosen layer/embedding:

```python
signatures = {
    "EPI_like": {"up": ["EpCAM", "KRT8-18"], "down": ["Vimentin"]},
    "Stem_like": ["CD44", "BMI1"],
}

# Import PermCell module first (run_permcell does not import it for you)
import PermCell_Smooth as PCS

res = run.run_permcell(
    signatures=signatures,
    source_layer="raw",
    positions_key="X_umap",     # any 2D coordinate in obsm
    result_prefix="permcell_epi",
    permcell_module=PCS,
)

# Stored outputs:
# - obsm:    permcell_epi_smoothed
# - obsm:    permcell_epi_z, permcell_epi_p, permcell_epi_zdir (plus zabs/pabs if enabled)
# - obsm:    permcell_epi_raw_z, permcell_epi_raw_p, ... (if compute_unsmoothed=True)
# - uns:     uns["permcell"]["permcell_epi"] metadata

# Export selected PermCell result to a compact AnnData object
perm_adata = run.permcell_to_adata(
    result_prefix="permcell_epi",
    score="z",
    smoothed=True,
)
```

## Heatmaps and Boxplots

Both helpers accept marker names (from `adata.var_names`) **or** numeric obs
columns (e.g. per-cell scores, cluster labels cast to numeric, manual gates),
and `groupby` can be a single obs column or a list of obs columns for composite
grouping.

```python
# Mean expression heatmap, markers as rows, groups as columns
run.plot_heatmap(
    fields=["H3", "H3K27me3", "ECad", "EpCAM"],
    groupby="sample_id",
    layer="normalized",
    agg="mean",
    standard_scale="row",   # z-score per marker for readability
)

# Composite grouping (line + condition)
run.plot_heatmap(
    fields=["H3K27ac", "H3K27me3"],
    groupby=["line_id", "condition"],
    layer="normalized",
    agg="median",
)

# Mix markers and numeric obs columns (e.g. a PermCell score)
run.plot_heatmap(
    fields=["H3K27me3", "permcell_epi_z"],
    groupby="leiden_cluster",
    layer="normalized",
    agg="mean",
    standard_scale="row",
    heatmap_kwargs={"linewidths": 0.5, "linecolor": "white", "vmax": 2.0},
)
```

```python
# Boxplot of one marker per group with significance brackets
run.plot_boxplot(
    field="H3K27me3",
    groupby="condition",
    layer="normalized",
    comparisons="all",    # None=no brackets, "all", "adjacent", or pair list
    test="mannwhitney",   # 'mannwhitney' | 'ttest' | 'welch'
    multitest="bh",       # None | 'bh' | 'bonferroni'
    show_points=True,     # overlay subsampled stripplot
    show_outliers=False,  # hide boxplot fliers
    boxplot_kwargs={"width": 0.5, "notch": True},
    stripplot_kwargs={"jitter": 0.15},
)

# Boxplot of a numeric obs column (e.g. PermCell z-score) per cluster
run.plot_boxplot(
    field="permcell_epi_z",
    groupby="leiden_cluster",
    test="welch",
    comparisons=[("0", "1"), ("0", "2")],   # only test selected pairs
)

# Customize significance star thresholds (strictest first)
run.plot_boxplot(
    field="H3K27me3",
    groupby="condition",
    comparisons="all",
    significance_thresholds=[
        (1e-3, "***"),
        (1e-2, "**"),
        (0.05, "*"),
    ],
    ns_label="ns",
)
```

Significance brackets use:

| p-value     | Label  |
| ----------- | ------ |
| `> 0.05`    | `ns`   |
| `<= 0.05`   | `*`    |
| `<= 0.01`   | `**`   |
| `<= 0.001`  | `***`  |
| `<= 0.0001` | `****` |


## Random Matrix Theory Spectrum

Identify which marker PCs represent genuine biology vs. noise using the
Marchenko-Pastur (MP) null distribution.

```python
# Marker covariance mode (p x p, fast, works on any n)
result = run.compute_rmt_spectrum(
    matrix="marker_cov",   # or "cell_gram" for the n x n Gram matrix
    standardize=True,      # use correlation matrix (sigma2=1 under null)
    n_cells=10_000,        # subsample cap per group
    groupby="condition",   # None for whole dataset
    plot=True,             # returns (result, fig)
)
# result["n_signal"]      – number of eigenvalues above the MP bound
# result["lambda_max_mp"] – MP upper bound λ_max = σ²(1+√q)²
# result["eigenvalues"]   – descending list of eigenvalues

# Check stability of the spectrum via bootstrap resampling
boot, fig = run.bootstrap_rmt_spectrum(
    matrix="marker_cov",
    frac=0.8,          # fraction of cells per resample
    n_bootstrap=200,
    plot=True,
)
# boot["eigenvector_stability"] – mean |cos θ| per component (marker_cov only)
# boot["eigenvalue_ci_low"]     – 2.5th percentile per rank across resamples
```

Both results are stored in `adata.uns` (keys `"rmt_spectrum"` and
`"rmt_bootstrap"` by default) and persisted automatically when `inplace=True`.

## API Documentation

This repository includes generated API description files:

- `docs/api/README.md` - index of all modules
- `docs/api/reference/*.md` - per-module API references
- `docs/api/api_manifest.json` - machine-readable API manifest

Regenerate API docs after any public API change:

```bash
python3 scripts/generate_api_docs.py
```
