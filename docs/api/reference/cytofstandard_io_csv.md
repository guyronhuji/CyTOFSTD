# `cytofstandard.io.csv`

- Source: `cytofstandard/io/csv.py`

CSV file reader for cytofstandard.

## Top-level Functions

### `read_csv(file_path: str) -> Tuple[pd.DataFrame, pd.DataFrame]`

Read a CSV file and return data and marker metadata.

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

## Classes

No public classes.
