# CyTOFSTD Data Structure

Every run is stored as a single [AnnData](https://anndata.readthedocs.io/) object serialised to Zarr.  
Path: `<project_dir>/runs/<run_id>/processed/<run_id>.zarr`

This document describes every field that CyTOFSTD writes, the method that writes it, and what the values mean.

---

## Dimensions

| Dimension | Axis | Description |
|-----------|------|-------------|
| `n_obs` | rows | cells (one row per event) |
| `n_vars` | columns | markers (panel columns) |

---

## `adata.X` — Working expression matrix

The matrix pointed to by `X` changes over the processing pipeline:

| After step | `X` contains |
|------------|--------------|
| Ingestion | `arcsinh(raw / 5)` transformed values (same as `layers["raw"]` passed through arcsinh) |
| `normalize()` | unchanged — `X` still holds arcsinh values; normalized values go in `layers["normalized"]` |
| Gating | unchanged |

> All gating (QC, cell-cycle) operates on `X` by default.  Pass `layer=` to override.

---

## `adata.layers` — Preserved expression matrices

| Key | Written by | Contents |
|-----|-----------|----------|
| `"raw"` | `Run.ingest()` | Raw ion counts as loaded from FCS/CSV/Parquet, **before** any transformation |
| `"normalized"` | `Run.normalize()` | Variance-stabilised, batch-corrected values (method stored in `uns["normalization"]`) |
| `"zscore"` | `Run.zscore()` | Per-marker z-scores computed across the run or a reference group |

---

## `adata.obs` — Per-cell metadata

### Core columns (always present after ingestion)

| Column | Type | Written by | Description |
|--------|------|-----------|-------------|
| `cell_uuid` | str | `ingest()` | UUID uniquely identifying this cell event across all runs in the project |
| `sample_id` | str | `ingest()` | Sample identifier from sample metadata |
| `line_id` | str | `ingest()` | Cell-line / biological-replicate identifier |
| `source_file_hash` | str | `ingest()` | SHA-256 of the source FCS/CSV/Parquet file for provenance |

Any additional columns in the sample metadata CSV are also copied here verbatim.

---

### QC columns

| Column | Type | Written by | Description |
|--------|------|-----------|-------------|
| `qc_pass` | bool | `Run.gate_qc()` | `True` if cell passed all QC gates |
| `qc_reason` | str | `Run.gate_qc()` | Human-readable reason for failing QC, or `""` for passing cells |

---

### Embedding & clustering columns

| Column | Type | Written by | Description |
|--------|------|-----------|-------------|
| `leiden` | Categorical str | `Run.cluster()` / `Run.cluster_flowsom()` | Default Leiden cluster label |
| `*` (any `cluster_key`) | Categorical str | `Run.cluster()` | User-specified clustering stored under the chosen key |
| `X_flowsom` | str | `Run.cluster_flowsom()` | FlowSOM metacluster label |
| `cluster_annotation` | str | `Run.annotate_clusters()` | Manual annotation string mapped from cluster labels |

---

### Normalization columns

| Column | Type | Written by | Description |
|--------|------|-----------|-------------|
| `norm_tech_factor` | float | `Run.normalize()` | Per-cell technical correction factor applied during normalization |

---

### Z-score columns

| Column | Type | Written by | Description |
|--------|------|-----------|-------------|
| `<marker>_zscore` | float | `Run.zscore()` | Z-score for each marker (one column per marker). Column names follow the pattern `<var_name>_zscore`. |

---

### Cell-cycle gating columns

Written by `Run.gate_cell_cycle()` / `cytofstandard.cell_cycle.gate_cell_cycle()`.

| Column | Type | Description |
|--------|------|-------------|
| `cell_cycle_phase` | Categorical str | Primary gated phase label. Values: see **Phase labels** table below. |

#### Phase labels

| Label | Meaning | Requires |
|-------|---------|---------|
| `S_phase` | DNA synthesis (IdU+) | IdU > threshold |
| `M_phase` | Mitosis (pH3+, IdU−) | pH3 > threshold, IdU negative |
| `G2_phase` | G2 (CyclinB1+, IdU−, pH3−) | CyclinB1 > threshold |
| `Cycling_G1` | Active G1 (pRb+, all others −) | pRb > threshold |
| `G0_or_quiescent` | Quiescent / off-cycle (pRb−) | pRb present but below threshold |
| `G1_or_quiescent` | G1 or quiescent, ambiguous | pRb marker not available in panel |
| `Unclassified` | Doesn't fit any gate | Fallback for edge cases |

The gating hierarchy is applied in exclusion order:  
**IdU → S** → **pH3 → M** → **CyclinB1 → G2** → **pRb → Cycling_G1 / G0_or_quiescent**

---

### Cell-cycle pseudotime columns

Written by `Run.add_cell_cycle_pseudotime()` / `cytofstandard.cell_cycle.add_cell_cycle_pseudotime()`.

> G0/quiescent cells are **excluded** from the cyclic coordinate. They receive `NaN` in all four numeric columns and `False` in `cell_cycle_on_cycle`.

| Column | Type | Range | Description |
|--------|------|-------|-------------|
| `cell_cycle_pseudotime` | float | [0, 1) or NaN | Continuous cyclic position. 0 = early G1, values increase through S → G2 → M and wrap back to 0. NaN for G0/quiescent and Unclassified cells. |
| `cell_cycle_angle` | float | [0, 2π) or NaN | Same as `cell_cycle_pseudotime × 2π`. Used for polar plots. |
| `cell_cycle_phase_index` | int | ≥ 0, or −1 | Zero-based index of this cell's phase in the observed cycling phase order. −1 for G0, Unclassified, or cells not matched to any phase. |
| `cell_cycle_within_phase_rank` | float | [0, 1) or NaN | Percentile rank within the cell's own phase, determined by the phase-specific marker (pRb for G1, DNA for S, CyclinB1 for G2, pH3 for M). |
| `cell_cycle_on_cycle` | bool | True / False | `True` for cycling cells (G1/S/G2/M). `False` for G0/quiescent and Unclassified cells. Use this column to separate quiescence fraction from cycling position. |

#### Within-phase ordering marker logic

| Phase | Ordering marker | Biological rationale |
|-------|----------------|---------------------|
| G1, Cycling_G1, early_G1, late_G1, G1S | pRb | CDK4/6 phosphorylate Rb as cells commit to S phase |
| S, S_phase | DNA (intercalator) | DNA content increases monotonically through replication |
| G2, G2_phase | CyclinB1 | CyclinB1 accumulates before mitotic entry |
| G2M | average(CyclinB1, pH3) | Transition between G2 accumulation and mitotic phosphorylation |
| M, M_phase | pH3 | H3 Ser10/Ser28 phosphorylation peaks during chromosome condensation |

---

### Diversity (Vendi) columns

| Column | Type | Written by | Description |
|--------|------|-----------|-------------|
| `<obs_key>` | float | `Run.vendi_score()` | Per-cell Vendi diversity score (default `"vendi_score"`) |
| `<obs_key>_ci_low` | float | `Run.vendi_score(bootstrap=True)` | Lower 95 % bootstrap CI |
| `<obs_key>_ci_high` | float | `Run.vendi_score(bootstrap=True)` | Upper 95 % bootstrap CI |

---

## `adata.var` — Per-marker metadata

| Column | Type | Written by | Description |
|--------|------|-----------|-------------|
| `standard_name` | str | `ingest()` | Canonical marker name from `standard_markers.csv` |
| `channel` | str | `ingest()` | Instrument channel (e.g. `"Ir191"`) |
| `panel_name` | str | `ingest()` | Original column name from the source file |
| `is_marker` | bool | `ingest()` | `True` for protein markers; `False` for DNA/live-dead/bead channels |

---

## `adata.obsm` — Per-cell embeddings

| Key | Shape | Written by | Description |
|-----|-------|-----------|-------------|
| `"X_umap"` | (n_obs, 2) | `Run.embed()` | 2-D UMAP coordinates |
| `"X_umap_knn_indices"` | (n_obs, k) | `Run.embed()` | K-nearest-neighbour indices used to build the UMAP graph |
| `"X_umap_knn_distances"` | (n_obs, k) | `Run.embed()` | Euclidean distances to KNN neighbours |
| `"X_flowsom_node"` | (n_obs, 1) | `Run.cluster_flowsom()` | SOM node assignment (0-based) |

Additional `obsm` keys may be written by `Run.embed()` if a custom `embedding_key` is supplied.

---

## `adata.obsp` — Per-cell-pair sparse matrices

| Key | Written by | Description |
|-----|-----------|-------------|
| `"X_umap_connectivities"` | `Run.embed()` | Weighted KNN connectivity graph (CSR sparse) |
| `"X_umap_distances"` | `Run.embed()` | KNN distance graph (CSR sparse) |
| `"X_umap_jaccard_connectivities"` | `Run.cluster()` | Jaccard-weighted graph used for Leiden clustering |

---

## `adata.varm` — Per-marker embeddings

| Key | Shape | Written by | Description |
|-----|-------|-----------|-------------|
| `"flowsom_weights"` | (n_vars, n_nodes) | `Run.cluster_flowsom()` | SOM codebook: learned prototype vector for each SOM node |

---

## `adata.uns` — Unstructured metadata

### `uns["project"]`

Set at ingestion. Never changes.

```python
{
    "project_id": str,    # UUID of the project
    "project_name": str,  # human-readable name
}
```

### `uns["run"]`

Set at ingestion. Never changes.

```python
{
    "run_id": str,             # UUID of this run
    "run_name": str,
    "panel_id": str,
    "acquisition_date": str,
    "instrument": str,
    "operator": str,
}
```

### `uns["ingestion"]`

Set at ingestion.

```python
{
    "created_at": str,           # ISO timestamp
    "package_version": str,
    "strict_markers": bool,
    "common_markers_only": bool,
    "drop_columns": list[str],
    "n_files": int,
    "n_cells": int,
    "n_markers": int,
}
```

### `uns["metadata_hashes"]`

SHA-256 hashes of input metadata files, for reproducibility auditing.

```python
{
    "sample_metadata_sha256": str,
    "standard_markers_sha256": str,
    "marker_aliases_sha256": str,
}
```

### `uns["qc"]`

Written by `Run.gate_qc()`. Has `"latest"` and `"history"` entries.

```python
{
    "latest": {
        "timestamp": str,         # ISO
        "gates": {                 # one entry per marker
            "<marker>": {
                "low": float | None,
                "high": float | None,
                "method": str,
            }
        },
        "n_pass": int,
        "n_fail": int,
    },
    "history": [...]              # list of prior "latest" snapshots
}
```

### `uns["normalization"]`

Written by `Run.normalize()`.

```python
{
    "latest": {
        "timestamp": str,
        "method": str,            # e.g. "quantile", "cyclical_loess"
        "groupby": str | None,
        "reference_group": str | None,
        "n_cells": int,
        "markers": list[str],
    },
    "history": [...]
}
```

### `uns["zscore"]`

Written by `Run.zscore()`.

```python
{
    "latest": {
        "timestamp": str,
        "groupby": str | None,
        "reference_group": str | None,
        "markers": list[str],
    },
    "history": [...]
}
```

### `uns["cell_cycle_gating"]`

Written by `Run.gate_cell_cycle()`.

```python
{
    "latest": {
        "timestamp": str,
        "marker_map": {           # role -> actual column name used
            "IdU": str,
            "pH3": str,
            "CyclinB1": str,
            "pRb": str,           # optional
        },
        "thresholds": {           # final thresholds applied (arcsinh space)
            "IdU": float,
            "pH3": float,
            "CyclinB1": float,
            "pRb": float,
        },
        "layer": str,             # which layer was gated ("X", "normalized", …)
        "n_cells": int,
    }
}
```

> **No `"history"` list** — each `gate_cell_cycle()` call overwrites `"latest"`.  
> The `"marker_map"` entry is consumed by `Run.add_cell_cycle_pseudotime()` to auto-fill its `marker_cols` argument.

### `uns["clusterings"]`

Written by `Run.cluster()`. One sub-dict per clustering run.

```python
{
    "<cluster_key>": {
        "timestamp": str,
        "method": str,            # "leiden" | "louvain"
        "resolution": float,
        "n_clusters": int,
        "embedding_key": str,
        "connectivities_key": str,
    }
}
```

### `uns["flowsom"]`

Written by `Run.cluster_flowsom()`.

```python
{
    "xdim": int,
    "ydim": int,
    "n_metaclusters": int,
    "markers": list[str],
    "metacluster_map": dict,      # node_id -> metacluster_id
}
```

### `uns["flowsom_grid"]`

Written by `Run.cluster_flowsom()`. Numpy array of shape `(n_obs, 2)` — (x, y) grid coordinates on the SOM for each cell.

### `uns["permcell"]`

Written by `Run.permcell()`. Contains PermCell differential-abundance results keyed by `result_prefix`.

### `uns["vendi"]`

Written by `Run.vendi_score()`. Keyed by `obs_key`.

```python
{
    "<obs_key>": {
        "<group_label>": float,   # Vendi score per group
        ...
    },
    "<obs_key>_eigenvalues": {
        "<group_label>": list[float],
    }
}
```

### `uns["derived_from"]`

Present on subset runs created by `Run.create_subset()`.

```python
{
    "parent_run_id": str,
    "created_at": str,
    "filter_description": str,
    "n_cells_original": int,
    "n_cells_subset": int,
}
```

---

## Typical processing order

```
ingest()
    → adata.X (arcsinh), layers["raw"], obs core columns, uns project/run/ingestion

gate_qc()
    → obs["qc_pass"], obs["qc_reason"], uns["qc"]

normalize()
    → layers["normalized"], obs["norm_tech_factor"], uns["normalization"]

embed()
    → obsm["X_umap"], obsp connectivities/distances

cluster()
    → obs["leiden"] (or custom key), uns["clusterings"]

zscore()
    → obs["<marker>_zscore"] columns, uns["zscore"]

gate_cell_cycle()
    → obs["cell_cycle_phase"], uns["cell_cycle_gating"]

add_cell_cycle_pseudotime()
    → obs["cell_cycle_pseudotime"], obs["cell_cycle_angle"],
      obs["cell_cycle_phase_index"], obs["cell_cycle_within_phase_rank"],
      obs["cell_cycle_on_cycle"]
```
