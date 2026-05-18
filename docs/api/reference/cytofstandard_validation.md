# `cytofstandard.validation`

- Source: `cytofstandard/validation.py`

Validation utilities for cytofstandard.

## Top-level Functions

### `normalize_marker_name(name: str) -> str`

Normalize a marker name for matching.

Args:
    name: Original marker name

Returns:
    Normalized marker name

### `validate_file_list(files: list[str], sample_metadata_df: pd.DataFrame) -> None`

Validate that every ingested file has matching sample metadata.

Args:
    files: List of file paths
    sample_metadata_df: Sample metadata DataFrame

Raises:
    MetadataValidationError if validation fails

### `validate_marker_consistency(marker_mappings: list[pd.DataFrame], strict: bool = True, allow_extra_markers: bool = False) -> list[str]`

Validate that all files have consistent marker sets.

Args:
    marker_mappings: List of marker mapping DataFrames (one per file)
    strict: Whether to fail on consistency violations
    allow_extra_markers: Whether to allow extra markers

Returns:
    List of standard markers that all files have in common

Raises:
    MarkerValidationError if validation fails

### `validate_sample_id_unique(sample_metadata_df: pd.DataFrame, allow_duplicate_samples: bool = False) -> None`

Validate sample IDs are unique within a run.

Args:
    sample_metadata_df: Sample metadata DataFrame
    allow_duplicate_samples: Whether to allow duplicate sample IDs

Raises:
    MetadataValidationError if validation fails

### `validate_sample_metadata(sample_metadata_df: pd.DataFrame) -> None`

Validate sample metadata has required columns.

Args:
    sample_metadata_df: Sample metadata DataFrame

Raises:
    MetadataValidationError if validation fails

## Classes

No public classes.
