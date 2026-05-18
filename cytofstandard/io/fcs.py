"""FCS file reader for cytofstandard."""

import fcsparser
import pandas as pd
from pathlib import Path
from typing import Tuple, Optional


def read_fcs(file_path: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Read an FCS file and return data and marker metadata.
    
    Args:
        file_path: Path to FCS file
        
    Returns:
        Tuple of (data_matrix, marker_metadata)
        - data_matrix: DataFrame with events as rows, channels as columns
        - marker_metadata: DataFrame with marker/channel information
    """
    # Parse FCS file
    meta, data = fcsparser.read(file_path)
    
    # Get channel names
    if "PnN" in meta:
        channel_names = [meta.get(f"PnN" % (i + 1)) for i in range(len(data.columns))]
    else:
        channel_names = data.columns.tolist()
    
    # Create dataframe
    data_df = pd.DataFrame(data.values, columns=channel_names)
    
    # Build marker metadata
    marker_metadata = []
    for i, col in enumerate(channel_names):
        entry = {
            "original_channel_name": col,
            "original_marker_name": meta.get(f"PnN" % (i + 1), col),
            "metal": meta.get(f"Pn" % (i + 1), None),
            "mass": meta.get(f"Pn" % (i + 1), None),
            "fcs_parameter": f"P{i + 1}",
        }
        marker_metadata.append(entry)
    
    marker_metadata_df = pd.DataFrame(marker_metadata)
    
    return data_df, marker_metadata_df
