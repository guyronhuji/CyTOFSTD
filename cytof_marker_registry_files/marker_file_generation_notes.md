# Marker file generation notes

Generated from `Markers_Names.xlsx` and `Mapping.xlsx`.

- Unique canonical markers from `Markers_Names.xlsx`: 95
- Extra canonical markers added from `Mapping.xlsx` STDName targets: 3
  - KLF4, MHC, mCD45
- Final canonical markers in `standard_markers.csv`: 98
- Duplicate canonical rows ignored from `Markers_Names.xlsx`: 2
  - H3K27me2, H3K9ac
- Canonical alias keys in `marker_aliases.yaml`: 98
- Total alias strings, including canonical self-aliases: 130

## Added-from-mapping markers needing review

- `KLF4`: appears as a `STDName` target in `Mapping.xlsx` but not in `Markers_Names.xlsx`; added with `marker_class=Unknown` and `needs_review=true`.
- `MHC`: appears as a `STDName` target in `Mapping.xlsx` but not in `Markers_Names.xlsx`; added with `marker_class=Unknown` and `needs_review=true`.
- `mCD45`: appears as a `STDName` target in `Mapping.xlsx` but not in `Markers_Names.xlsx`; added with `marker_class=Unknown` and `needs_review=true`.

## Duplicate canonical rows ignored

- `H3K27me2`
- `H3K9ac`

## Flag derivation

- `is_core_histone` set true for H2A, H2B, H3, H3.3, H4.
- `is_qc_marker` set true for core histones plus Iridium and DNA.
- `is_identity_marker` derived from broad phenotype groups such as Cancer, Immune, CAFs, Human-mouse, EMT, Stemness.
- `is_functional_marker` derived from Chromatin, Cell-cycle, Meth groups plus selected signaling/regulatory markers.
