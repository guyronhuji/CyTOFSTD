"""Validation utilities for cytofstandard."""

import pandas as pd
from pathlib import Path
from typing import Optional

from cytofstandard.exceptions import MetadataValidationError


def validate_file_list(
    files: list[str],
    sample_metadata_df: pd.DataFrame,
) -> None:
    """Validate that every ingested file has matching sample metadata.

    Args:
        files: List of file paths
        sample_metadata_df: Sample metadata DataFrame

    Raises:
        MetadataValidationError if validation fails
    """
    file_names = [Path(f).name for f in files]
    metadata_files = sample_metadata_df["file_name"].tolist()

    # Check for missing files in metadata
    missing_in_metadata = set(file_names) - set(metadata_files)
    if missing_in_metadata:
        raise MetadataValidationError(
            f"The following input files are missing from sample_metadata: "
            f"{sorted(missing_in_metadata)}"
        )


def validate_sample_metadata(sample_metadata_df: pd.DataFrame) -> None:
    """Validate sample metadata has required columns.

    Args:
        sample_metadata_df: Sample metadata DataFrame

    Raises:
        MetadataValidationError if validation fails
    """
    required_cols = ["file_name", "sample_id", "line_id"]
    missing_cols = [
        col for col in required_cols if col not in sample_metadata_df.columns
    ]

    if missing_cols:
        raise MetadataValidationError(
            f"Sample metadata is missing required columns: {missing_cols}"
        )

    # Check for missing values in required columns
    for col in required_cols:
        if sample_metadata_df[col].isna().any():
            raise MetadataValidationError(
                f"Sample metadata has missing values in column: {col}"
            )

    # Check that sample_id maps to exactly one line_id
    sample_to_line = sample_metadata_df.groupby("sample_id")["line_id"].nunique()
    if (sample_to_line > 1).any():
        bad_samples = sample_to_line[sample_to_line > 1].index.tolist()
        raise MetadataValidationError(
            f"Samples mapped to multiple lines: {bad_samples}"
        )


def validate_sample_id_unique(
    sample_metadata_df: pd.DataFrame,
    allow_duplicate_samples: bool = False,
) -> None:
    """Validate sample IDs are unique within a run.

    Args:
        sample_metadata_df: Sample metadata DataFrame
        allow_duplicate_samples: Whether to allow duplicate sample IDs

    Raises:
        MetadataValidationError if validation fails
    """
    if not allow_duplicate_samples:
        sample_counts = sample_metadata_df["sample_id"].value_counts()
        duplicates = sample_counts[sample_counts > 1].index.tolist()
        if duplicates:
            raise MetadataValidationError(
                f"Duplicate sample IDs found: {duplicates}. "
                "Set allow_duplicate_samples=True to allow."
            )


def validate_marker_consistency(
    marker_mappings: list[pd.DataFrame],
    strict: bool = True,
    allow_extra_markers: bool = False,
) -> list[str]:
    """Validate that all files have consistent marker sets.

    Args:
        marker_mappings: List of marker mapping DataFrames (one per file)
        strict: Whether to fail on consistency violations
        allow_extra_markers: Whether to allow extra markers

    Returns:
        List of standard markers that all files have in common

    Raises:
        MarkerValidationError if validation fails
    """
    if not marker_mappings:
        return []

    # Get standard markers for each file
    file_markers = []
    for mapping in marker_mappings:
        matched = mapping[
            mapping["mapping_status"].isin(["matched_standard", "matched_alias"])
        ]
        file_markers.append(set(matched["standard_marker_name"].tolist()))

    if not file_markers:
        return []

    # Find common markers
    common = file_markers[0]
    for markers in file_markers[1:]:
        common = common.intersection(markers)

    # Check for missing markers
    for i, markers in enumerate(file_markers):
        missing = markers - common
        if missing and strict:
            raise MarkerValidationError(
                f"File {i} is missing markers that other files have: {missing}"
            )

    # Check for extra markers
    for i, markers in enumerate(file_markers):
        extra = markers - common
        if extra and not allow_extra_markers and strict:
            raise MarkerValidationError(
                f"File {i} has extra markers not in other files: {extra}. "
                "Set allow_extra_markers=True to allow."
            )

    return sorted(list(common))


def normalize_marker_name(name: str) -> str:
    """Normalize a marker name for matching.

    Args:
        name: Original marker name

    Returns:
        Normalized marker name
    """
    # Strip whitespace
    normalized = name.strip()
    # Collapse whitespace
    import re

    normalized = re.sub(r"\s+", " ", normalized)
    return normalized
