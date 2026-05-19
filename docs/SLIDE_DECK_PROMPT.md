# LLM Prompt: Generate CyTOF Standard Slide Deck

Use this prompt verbatim (or paste into Claude / GPT-4) to generate a
python-pptx script that builds the slide deck.

---

## PROMPT START

You are an expert in Python and python-pptx. Generate a complete, runnable Python script
that builds a PowerPoint slide deck called `CyTOF_Standard_Overview.pptx`.

The deck is for a lab presentation introducing the `cytofstandard` analysis package
to a mixed audience of wet-lab biologists and computational analysts.
Tone: clear, professional, scientific. No filler slides.

### Design specifications

- Slide size: widescreen 16:9 (33.867 cm × 19.05 cm)
- Background: white (#FFFFFF)
- Title slides: dark navy header bar (#1F3864), white text
- Section divider slides: solid navy background (#1F3864), large white text
- Content slides: white background, navy headings (#1F3864), dark gray body (#333333)
- Accent color for highlights and table headers: steel blue (#2E75B6)
- Code blocks: light gray background (#F2F2F2), monospace font (Courier New 10pt)
- Body font: Calibri. Heading font: Calibri Bold.
- Keep slides uncluttered — max 6 bullet points per slide, max 2 code blocks per slide.

### Slide outline

Generate exactly the slides listed below, in this order.
For each slide, the content is specified. Implement it faithfully.

---

**Slide 1 — Title slide**
- Title: "CyTOF Standard Analysis Package"
- Subtitle: "Reproducible CyTOF analysis from submission to statistics"
- Footer (small): "cytofstandard · Lab internal · 2026"
- Layout: navy full-width header bar, white lower area

---

**Slide 2 — Agenda**
- Title: "What this deck covers"
- Two-column layout:
  - Left column header "For wet-lab users":
    - Why we switched to standardized analysis
    - What to submit after each CyTOF run
    - Folder structure and metadata rules
    - Submission checklist
  - Right column header "For analysis users":
    - Package architecture and data model
    - Full analysis workflow (API walkthrough)
    - Embedding, clustering, statistics
    - Export and reproducibility

---

**Slide 3 — Section divider**
- Navy background
- Large text: "Part 1"
- Subtitle: "Why we use cytofstandard"

---

**Slide 4 — The problem with the old workflow**
- Title: "The old workflow had hidden costs"
- Bullet points (use ✗ prefix, red-ish color #C00000):
  - Marker names differed across runs (H3K27Ac vs H3K27ac vs H3K27ac_168Er)
  - No permanent cell IDs — reanalysis could not match cells across sessions
  - Run-specific scripts made comparisons fragile
  - Sample/line linkage stored in filenames or analyst memory
  - No standardized processing log

---

**Slide 5 — What cytofstandard solves**
- Title: "cytofstandard: five core guarantees"
- Five rows in a clean table (2 columns: "Guarantee" | "What it means"):
  1. Permanent cell IDs | Every cell gets a UUID fixed at ingestion. Reanalysis finds the same cells.
  2. Standardized marker names | Alternate names mapped to lab standard at ingestion. No silent mismatches.
  3. Run-first analysis | Each run analyzed independently before any cross-run integration.
  4. Explicit sample → line linking | Every cell is traceable to its sample and biological line.
  5. Full provenance log | Every processing step (normalization, z-score, clustering) is recorded with parameters and version.

---

**Slide 6 — Section divider**
- Navy background
- Large text: "Part 2"
- Subtitle: "Data submission (wet-lab users)"

---

**Slide 7 — What to submit**
- Title: "For every CyTOF run, submit two things"
- Two large numbered items (big bold numerals):
  1. Debarcoded data files (FCS, CSV, or Parquet) — after preliminary QC by the CyTOF unit
  2. sample_metadata.csv — one row per submitted file, with file_name, sample_id, line_id
- Note box (yellow/green): "Do not submit pre-debarcoding multiplexed files."

---

**Slide 8 — Folder structure**
- Title: "Required folder structure"
- Code block (monospace):
```
run_2026_05_15/
    data_files/
        S001_MCF7_control_rep1.fcs
        S002_MCF7_treated_rep1.fcs
        S003_T47D_control_rep1.fcs
    sample_metadata.csv
```
- Short note below: "One folder per run. Do not rename files after creating the metadata."

---

**Slide 9 — sample_metadata.csv — required columns**
- Title: "sample_metadata.csv — required columns"
- Table (3 columns: Column | Required? | Meaning):
  - file_name | Yes | Exact filename. Case-sensitive. Include extension.
  - sample_id | Yes | Unique sample identifier for this run.
  - line_id   | Yes | Biological line or group (e.g. MCF7, T47D).
- Below table: "Recommended: condition, replicate_id, batch_id, acquisition_order, barcoding_id, notes"
- Code block (small example):
```
file_name,sample_id,line_id,condition,replicate_id
S001_MCF7_control_rep1.fcs,S001,MCF7,control,rep1
S002_T47D_treated_rep1.fcs,S002,T47D,treated,rep1
```

---

**Slide 10 — Metadata rules (the critical ones)**
- Title: "Four rules you must follow"
- Four bullet points with bold rule name followed by brief explanation:
  - Every submitted file must appear exactly once — if a file is missing from the metadata, ingestion stops.
  - File names must match exactly — case-sensitive, with extension. sample_A.fcs ≠ Sample_A.fcs
  - Each file maps to one sample and one line — duplicate rows for the same file with different lines will fail.
  - Use consistent spelling — MCF7 and MCF-7 are treated as different lines.

---

**Slide 11 — Marker names**
- Title: "Marker name standardization"
- Left side (60%): explanation paragraph:
  "The package maps common alternate names to the lab standard at ingestion.
   Use standard names whenever possible. The package recognizes these alternates:"
- Right side (40%): small 2-column table (Alternate → Standard):
  H3K27Ac → H3K27ac
  H3K9Ac → H3K9ac
  E-cad → E-cadherin
  PanKRT → Pan-KRT
  DNA1 → DNA
  DNA2 → DNA
- Warning box (orange): "Two DNA channels (DNA1 + DNA2) in the same file may clash. Notify the analysis person before ingestion."

---

**Slide 12 — Submission checklist**
- Title: "Submission checklist"
- Two-column checklist layout, checkboxes (☐):
  Left column:
  - ☐ Debarcoded files, post-QC
  - ☐ All files in data_files/
  - ☐ Metadata named sample_metadata.csv
  - ☐ Every file appears once in metadata
  - ☐ Each file → one sample_id, one line_id
  Right column:
  - ☐ file_name matches exactly (case, extension)
  - ☐ sample_id and line_id filled for all rows
  - ☐ Consistent spelling of lines and conditions
  - ☐ Acquisition order included if known
  - ☐ Problems noted in notes column

---

**Slide 13 — Section divider**
- Navy background
- Large text: "Part 3"
- Subtitle: "How the data is stored: AnnData"

---

**Slide 14 — AnnData concept**
- Title: "AnnData: one self-describing object per run"
- Central ASCII-art diagram (monospace, light gray box):
```
                  markers (var)
              H3  H3K27ac  ECad  ...
           ┌─────────────────────┐
 cell_0001 │  1.23   0.45  3.21 │ ← obs: sample_id, line_id,
 cell_0002 │  1.11   0.39  3.45 │        cell_uuid, cluster ...
 ...       │  ...    ...   ...  │
           └─────────────────────┘
                  X  (active layer)

layers["raw"]         obsm["X_umap"]
layers["normalized"]  obsp["connectivities"]
layers["normalized_z"] uns["provenance"]
```
- Short caption: "All processing layers, embeddings, and metadata travel together."

---

**Slide 15 — obs: cell metadata columns**
- Title: "obs — cell metadata (one row per cell)"
- Table (3 columns: Column | Added when | Description):
  cell_uuid      | Ingestion  | Permanent UUID. Never changes.
  sample_id      | Ingestion  | From sample_metadata.csv
  line_id        | Ingestion  | From sample_metadata.csv
  source_file    | Ingestion  | Original file the cell came from
  event_index    | Ingestion  | Row index in source file (0-based)
  qc_pass        | QC gating  | True if cell passed all gates
  norm_umap_leiden | Clustering | Integer cluster label
  PC_EPI_like    | PermCell   | Smoothed signature z-score

---

**Slide 16 — layers: expression matrices**
- Title: "layers — named snapshots of the expression matrix"
- Visual pipeline diagram (left-to-right flow boxes connected by arrows):
  [raw] → [normalized] → [normalized_z]
  Labels below each box:
    raw: "Original debarcoded values. Never overwritten."
    normalized: "After cytof_transform correction per sample."
    normalized_z: "After balanced z-scoring across groups."
- Below diagram: note "X points to the active layer. run.set_x_from_layer('normalized') switches X without copying data."

---

**Slide 17 — uns and obsm**
- Title: "uns and obsm — metadata and embeddings"
- Two-column layout:
  Left (uns):
  Title: "uns — run-level metadata"
  Bullets:
  - cytofstandard_version
  - marker_mapping (original → standard)
  - normalization parameters
  - z-score mean/std per group
  - clustering settings (resolution, seed)
  - provenance log (all steps + timestamps)
  Right (obsm / obsp):
  Title: "obsm / obsp — per-cell arrays"
  Bullets:
  - obsm["X_norm_umap"] — UMAP coordinates
  - obsm["norm_umap_indices"] — KNN indices
  - obsp["norm_umap_connectivities"] — KNN graph
  - obsp["norm_umap_jaccard_connectivities"] — Jaccard graph
  - obsm["pc_smoothed"] — PermCell score matrix

---

**Slide 18 — Disk layout**
- Title: "Project layout on disk"
- Code block (monospace, small font):
```
my_project/
├── project.json
├── standard_markers.csv
└── runs/
    └── run_001/
        ├── run.json
        ├── provenance.jsonl
        ├── raw/  (original submitted files)
        └── processed/
            └── run_001.zarr/
                ├── X/
                ├── obs/
                ├── var/
                ├── layers/raw, normalized, normalized_z
                ├── obsm/  (embeddings)
                ├── obsp/  (graphs)
                └── uns/   (metadata)
```

---

**Slide 19 — Section divider**
- Navy background
- Large text: "Part 4"
- Subtitle: "Analysis workflow (API walkthrough)"

---

**Slide 20 — Workflow overview**
- Title: "End-to-end analysis workflow"
- Numbered vertical pipeline (single column, large step numbers):
  1. Install & prepare marker files
  2. Create / load project
  3. Register run + ingest data
  4. QC: histograms + gating
  5. Normalization (cytof_transform)
  6. Balanced z-score
  7. Embedding (UMAP)
  8. Clustering (Leiden)
  9. PermCell signature scoring
  10. Plotting & statistics
  11. Export / lock

---

**Slide 21 — Project and run setup**
- Title: "Step 1–3: Project and run setup"
- Code block:
```python
from cytofstandard import Project

project = Project.create(
    path="my_project",
    project_id="BRCA_CYTOF_2026",
    project_name="BRCA CyTOF histone panel",
    standard_marker_file="standard_markers.csv",
    marker_alias_file="marker_aliases.yaml",
)

run = project.add_run(
    run_id="run_001",
    acquisition_date="2026-05-15",
    instrument="Helios",
    operator="GR",
)

run.ingest(
    files=["data_files/sample_A.fcs", "data_files/sample_B.fcs"],
    sample_metadata="sample_metadata.csv",
    copy_raw=True,
    strict_markers=True,
)
```

---

**Slide 22 — QC and normalization**
- Title: "Step 4–6: QC, normalization, z-scoring"
- Code block (two sections separated by a comment line):
```python
# QC
run.plot_marker_histograms(markers=["H3", "ECad"], layer="raw")
run.qc_gate({"H3": {"lower": 100, "upper": 400}}, layer="raw")

# Normalization
run.normalize_with_cytof_transform(
    control_markers=["H3.3", "H3", "H4"],
    markers_to_correct=["H3K27ac", "H3K4me3"],
    source_layer="raw",
    groupby_col="sample_id",
)

# Balanced z-score
run.zscore_markers_balanced(
    source_layer="normalized",
    output_layer="normalized_z",
    groupby_col="sample_id",
)
```

---

**Slide 23 — Embedding and clustering**
- Title: "Step 7–8: UMAP and Leiden clustering"
- Code block:
```python
run.set_x_from_layer("normalized")

# Balanced UMAP: fit on equal-size subsample, transform all cells
run.compute_umap_balanced(
    markers=["H3", "H3K27me3", "ECad", "EpCAM"],
    source_layer="normalized",
    embedding_name="norm_umap",
    groupby_col="sample_id",
    n_per_group=2000,
)

# Leiden clustering
run.cluster_leiden(embedding_name="norm_umap", resolution=1.0)

# PhenoGraph-style: Jaccard similarity + Leiden
run.cluster_leiden_jaccard(embedding_name="norm_umap", resolution=1.0)
```

---

**Slide 24 — PermCell signature scoring**
- Title: "Step 9: PermCell signature scoring"
- Left panel (explanation):
  "PermCell scores each cell for biological signatures using spatial smoothing on the UMAP graph.
   Results are stored directly in obs for easy plotting."
- Code block:
```python
import PermCell_Smooth as PCS

signatures = {
    "EPI_like":  {"up": ["EpCAM", "KRT8-18"], "down": ["Vimentin"]},
    "Stem_like": ["CD44", "BMI1"],
}

run.run_permcell(
    signatures=signatures,
    source_layer="raw",
    positions_key="X_norm_umap",
    result_prefix="pc",
    permcell_module=PCS,
)
# obs["PC_EPI_like"]   — smoothed z-score
# obs["PC_R_EPI_like"] — raw z-score
```

---

**Slide 25 — Plotting and statistics**
- Title: "Step 10: Plotting and statistics"
- Code block (three sections):
```python
# Heatmap: mean marker values per group
run.plot_heatmap(
    fields=["H3", "H3K27me3", "PC_EPI_like"],
    groupby="sample_id", layer="normalized",
    agg="mean", standard_scale="row",
)

# Boxplot with significance brackets
run.plot_boxplot(
    field="PC_EPI_like", groupby="condition",
    comparisons="all", test="mannwhitney", multitest="bh",
)

# Pairwise comparison table
table = run.compare_groups(
    field="PC_EPI_like", groupby="condition",
    method="ttest", multitest="bh",
)
```

---

**Slide 26 — Export and locking**
- Title: "Step 11: Export and reproducibility locking"
- Code block:
```python
# Export selected fields to DataFrame / CSV
df = run.to_dataframe(
    fields=["H3", "H3K27me3", "sample_id", "condition"],
    layer="raw",
)
df.to_csv("my_export.csv")

# Lock raw layers to prevent accidental modification
run.lock_zarr_parts(parts=["layers/raw", "obs"])

# Unlock when re-processing is needed
run.unlock_zarr_parts(parts=["layers/raw"])
```
- Note box: "Locking writes file-system permissions directly on the Zarr subpaths."

---

**Slide 27 — Common issues**
- Title: "Common issues and fixes"
- Table (2 columns: Issue | Fix):
  MetadataValidationError: file not in metadata | Check file_name exact match (case + extension)
  MarkerValidationError: unknown markers | Use allow_extra_markers=True or add to marker_aliases.yaml
  Markers differ across files | Set common_markers_only=True to keep the intersection
  UMAP biased by group size imbalance | Use compute_umap_balanced() with n_per_group
  Dense boxplots / slow plotting | Use show_points=False; subset first with create_subset_run()

---

**Slide 28 — Summary**
- Title: "Summary"
- Two-column layout:
  Left — "For wet-lab users":
  - Submit debarcoded FCS / CSV / Parquet files
  - Include sample_metadata.csv with file_name, sample_id, line_id
  - Follow naming and spelling rules exactly
  - Use the submission checklist before handing off
  Right — "For analysis users":
  - Project → add_run → ingest → QC → normalize → zscore → UMAP → cluster → PermCell → plot
  - All data stored in AnnData/Zarr with full provenance
  - Every cell has a permanent UUID traceable to source file and event index
  - Raw data is never overwritten; all processing stages retained as named layers

---

### python-pptx implementation notes

- Use `python-pptx` only (no external dependencies beyond `Pillow` for any image work).
- Section divider slides: set slide background fill to navy #1F3864.
- Content slide layout: use blank layout (index 6), position all elements manually with `Inches()` / `Cm()`.
- For the code blocks: add a rectangle shape with light gray fill (#F2F2F2), no line, add a text frame inside using Courier New 10pt, preserve leading whitespace using spaces (not tabs).
- For tables: use `add_table`, shade header row with steel blue #2E75B6 fill, white bold text.
- For the note/warning boxes: rectangle shape with fill color, no border, italic text inside.
- For the ASCII diagram on slide 14: place in a text box with Courier New 9pt, light gray background box behind it.
- The script should be self-contained and runnable with `python3 build_slide_deck.py`.
- Save the output as `CyTOF_Standard_Overview.pptx` in the current directory.
- Print a summary at the end: number of slides written and output path.

## PROMPT END
