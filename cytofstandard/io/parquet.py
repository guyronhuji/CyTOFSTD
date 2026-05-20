"""Parquet file reader for cytofstandard."""

import pandas as pd
from typing import Tuple


def read_parquet(file_path: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Read a Parquet file and return data and marker metadata.

    Assumes:
    - Rows are cells/events
    - Numeric columns are marker columns
    - Non-numeric columns are metadata columns (ignored, as with CSV)

    Args:
        file_path: Path to Parquet file

    Returns:
        Tuple of (data_matrix, marker_metadata)
        - data_matrix: DataFrame with events as rows, markers as columns
        - marker_metadata: DataFrame with marker information
    """
    df = pd.read_parquet(file_path)

    numeric_cols = df.select_dtypes(include=["number"]).columns.tolist()

    data_df = df[numeric_cols].copy()

    marker_metadata = []
    for col in numeric_cols:
        marker_metadata.append(
            {
                "original_channel_name": col,
                "original_marker_name": col,
                "column_type": "numeric",
            }
        )

    marker_metadata_df = pd.DataFrame(marker_metadata)

    return data_df, marker_metadata_df


__all__ = ["read_parquet"]
