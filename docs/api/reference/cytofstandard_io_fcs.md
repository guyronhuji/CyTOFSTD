# `cytofstandard.io.fcs`

- Source: `cytofstandard/io/fcs.py`

FCS file reader for cytofstandard.

## Top-level Functions

### `read_fcs(file_path: str, channel_naming: str = '$PnS') -> Tuple[pd.DataFrame, pd.DataFrame]`

Read an FCS file and return data and marker metadata.

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

## Classes

No public classes.
