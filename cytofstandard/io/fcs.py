"""FCS file reader for cytofstandard."""

import fcsparser
import pandas as pd
from pathlib import Path
from typing import Tuple


def read_fcs(file_path: str, channel_naming: str = "$PnS") -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Read an FCS file and return data and marker metadata.

    Args:
        file_path: Path to FCS file
        channel_naming: FCS field used as column names. ``"$PnS"`` (default)
            gives the stain/label names (e.g. ``140Ce_Cytokeratin_5``) which
            the marker sanitizer can parse. ``"$PnN"`` gives the short metal
            names (e.g. ``Ce140Di``) which are unique but carry no marker info.

    Returns:
        Tuple of (data_matrix, marker_metadata)
        - data_matrix: DataFrame with events as rows, channels as columns
        - marker_metadata: DataFrame with one row per channel
    """
    meta, data = fcsparser.parse(file_path, channel_naming=channel_naming)

    marker_metadata = []
    for i, channel_name in enumerate(data.columns):
        idx = i + 1
        stain = meta.get(f"$P{idx}S", "")
        marker_metadata.append(
            {
                "original_channel_name": channel_name,
                "original_marker_name": stain or channel_name,
                "fcs_parameter": f"P{idx}",
            }
        )

    return data, pd.DataFrame(marker_metadata)


__all__ = ["read_fcs"]
