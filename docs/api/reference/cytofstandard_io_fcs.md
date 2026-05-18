# `cytofstandard.io.fcs`

- Source: `cytofstandard/io/fcs.py`

FCS file reader for cytofstandard.

## Top-level Functions

### `read_fcs(file_path: str) -> Tuple[pd.DataFrame, pd.DataFrame]`

Read an FCS file and return data and marker metadata.

Args:
    file_path: Path to FCS file
    
Returns:
    Tuple of (data_matrix, marker_metadata)
    - data_matrix: DataFrame with events as rows, channels as columns
    - marker_metadata: DataFrame with marker/channel information

## Classes

No public classes.
