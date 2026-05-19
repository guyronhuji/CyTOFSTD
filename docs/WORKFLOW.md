# CyTOF Standard Student Workflow

This guide walks through a standard CyTOF analysis workflow using the
`cytofstandard` package. It is designed to be followed end-to-end in a
notebook or script.

## 1) Install dependencies

```bash
pip install anndata zarr numpy pandas pyyaml fcsparser matplotlib seaborn pyarrow
```

## 2) Prepare marker registry files

You need:

- `standard_markers.csv` with at least columns:
  - `standard_marker_name`
  - `marker_class`
- `marker_aliases.yaml` (optional, but recommended)

## 3) Create a project

```python
from cytofstandard import Project

project = Project.create(
    path="my_project",
    project_id="BRCA_CYTOF_2026",
    project_name="BRCA CyTOF histone panel",
    standard_marker_file="standard_markers.csv",
    marker_alias_file="marker_aliases.yaml",
)
```

## 4) Prepare sample metadata

Your sample metadata file must include:

- `file_name`
- `sample_id`
- `line_id`

Optional columns include `condition`, `replicate_id`, `batch_id`, `barcoding_id`,
`notes`.

Example CSV:

```csv
file_name,sample_id,line_id,condition,replicate_id
sample_A.fcs,S001,MCF7,control,rep1
sample_B.fcs,S002,T47D,treated,rep1
```

## 5) Register a run and ingest data

Supported input formats: **CSV**, **FCS**, **Parquet**.

```python
run = project.add_run(
    run_id="run_001",
    run_name="First BRCA CyTOF run",
    panel_id="breast_histone_panel_v1",
    acquisition_date="2026-05-15",
    instrument="Helios",
    operator="GR",
)

run.ingest(
    files=["sample_A.fcs", "sample_B.fcs"],
    sample_metadata="sample_metadata.csv",
    copy_raw=True,
    strict_markers=True,
    allow_extra_markers=False,
    drop_columns=["Time"],
)
```

Notes:
- `copy_raw=True` stores original files under `runs/<run_id>/raw`.
- `strict_markers=True` fails if unknown markers are present.
- `drop_columns` removes raw channels before marker standardization.

## 6) Load data

```python
adata = run.read_adata()
```

If you modify `adata` externally, persist it:

```python
run.save(adata)
```

## 7) QC

Histograms:

```python
run.plot_marker_histograms(markers=["H3", "ECad"], layer="raw")
```

Simple gating:

```python
run.qc_gate({"H3": {"lower": 100, "upper": 400}}, layer="raw")
```

## 8) Normalization (cytof_transform)

```python
summary = run.normalize_with_cytof_transform(
    control_markers=["H3.3", "H3", "H4"],
    markers_to_correct=["H3K27ac", "H3K4me3"],
    source_layer="raw",
    groupby_col="sample_id",
    input_is_arcsinh=False,
)
```

## 9) Embedding and clustering

```python
run.set_x_from_layer("normalized")

run.compute_umap(
    markers=["H3", "H3K27me3", "ECad", "EpCAM"],
    source_layer="normalized",
    embedding_name="norm_umap",
    n_neighbors=15,
    min_dist=0.1,
)

run.cluster_leiden(
    embedding_name="norm_umap",
    resolution=1.0,
)

# PhenoGraph-style Jaccard + Leiden
run.cluster_leiden_jaccard(
    embedding_name="norm_umap",
    resolution=1.0,
    min_jaccard=0.0,
)
```

## 10) PermCell

```python
import PermCell_Smooth as PCS

signatures = {
    "EPI_like": {"up": ["EpCAM", "KRT8-18"], "down": ["Vimentin"]},
    "Stem_like": ["CD44", "BMI1"],
}

run.run_permcell(
    signatures=signatures,
    source_layer="raw",
    positions_key="X_umap",
    result_prefix="pc",
    permcell_module=PCS,
)

# Z-scores are written to obs for easy plotting:
#   PC_<signature>   (smoothed)
#   PC_R_<signature> (raw)
```

## 11) Plotting and statistics

Heatmap of markers or obs fields by group:

```python
run.plot_heatmap(
    fields=["H3", "H3K27me3", "PC_EPI_like"],
    groupby="sample_id",
    layer="normalized",
    agg="mean",
    standard_scale="row",
)
```

Boxplot with significance brackets:

```python
run.plot_boxplot(
    field="PC_EPI_like",
    groupby="condition",
    comparisons="all",
    test="ttest",
    multitest="bh",
)
```

Tabulate comparisons (t-test or Wald):

```python
table = run.compare_groups(
    field="PC_EPI_like",
    groupby="condition",
    method="ttest",
    multitest="bh",
)
```

## 12) Export for downstream analysis

```python
df = run.to_dataframe(fields=["H3", "H3K27me3", "sample_id"], layer="raw")
```

## 13) Locking and reproducibility

```python
run.lock_zarr_parts(parts=["layers/raw", "obs"])
run.unlock_zarr_parts(parts=["layers/raw", "obs"])
```

## 14) Common issues

- **Missing markers**: set `allow_extra_markers=True` or update
  `marker_aliases.yaml` to map alternate names.
- **Sample metadata mismatch**: ensure every input file appears in
  `sample_metadata` and that `sample_id` maps to exactly one `line_id`.
- **Large datasets**: use `show_points=False` for dense boxplots and
  consider subsetting before exploratory plotting.
