"""CSV file reader for cytofstandard."""

import pandas as pd
from pathlib import Path
from typing import Tuple


def read_csv(file_path: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Read a CSV file and return data and marker metadata.
    
    Assumes:
    - Rows are cells/events
    - Numeric columns are marker columns
    - Non-numeric columns are metadata columns
    
    Args:
        file_path: Path to CSV file
        
    Returns:
        Tuple of (data_matrix, marker_metadata)
        - data_matrix: DataFrame with events as rows, markers as columns
        - marker_metadata: DataFrame with marker information
    """
    # Read CSV
    df = pd.read_csv(file_path)
    
    # Separate numeric (marker) and non-numeric (metadata) columns
    numeric_cols = df.select_dtypes(include=["number"]).columns.tolist()
    non_numeric_cols = df.select_dtypes(exclude=["number"]).columns.tolist()
    
    # Data matrix is just thenumeric columns
    data_df = df[numeric_cols].copy()
    
    # Build marker metadata
    marker_metadata = []
    for col in numeric_cols:
        marker_metadata.append({
            "original_channel_name": col,
            "original_marker_name": col,
            "column_type": "numeric",
        })
    
    marker_metadata_df = pd.DataFrame(marker_metadata)
    
    return data_df, marker_metadata_df
