# `cytofstandard.storage`

- Source: `cytofstandard/storage.py`

Storage utilities for cytofstandard.

## Public Exports (`__all__`)

- `compute_sha256`
- `compute_file_hashes`
- `read_yaml`
- `write_yaml`
- `read_parquet`
- `write_parquet`
- `ensure_directory`
- `get_directory_size`
- `copy_file`
- `file_exists`
- `directory_exists`
- `copy_files_to_directory`
- `set_zarr_parts_writable`
- `get_locked_zarr_parts`

## Top-level Functions

### `compute_file_hashes(file_paths: list[str]) -> dict[str, str]`

Compute SHA256 hashes for multiple files.

Args:
    file_paths: List of file paths
    
Returns:
    Dictionary mapping file path to hash

### `compute_sha256(file_path: str) -> str`

Compute SHA256 hash of a file.

Args:
    file_path: Path to the file
    
Returns:
    Hex digest of SHA256 hash

### `copy_file(src: str, dst: str)`

Copy a file from source to destination.

Args:
    src: Source file path
    dst: Destination file path

### `copy_files_to_directory(source_files: list[str], dest_dir: str, preserve_name: bool = True) -> dict[str, str]`

Copy multiple files to a directory.

Args:
    source_files: List of source file paths
    dest_dir: Destination directory
    preserve_name: Whether to preserve original filenames
    
Returns:
    Dictionary mapping source path to destination path

### `directory_exists(dir_path: str) -> bool`

Check if a directory exists.

Args:
    dir_path: Path to directory
    
Returns:
    True if directory exists

### `ensure_directory(path: str)`

Ensure a directory exists, creating it if necessary.

Args:
    path: Path to directory

### `file_exists(file_path: str) -> bool`

Check if a file exists.

Args:
    file_path: Path to file
    
Returns:
    True if file exists

### `get_directory_size(path: str) -> int`

Get total size of all files in a directory.

Args:
    path: Path to directory
    
Returns:
    Total size in bytes

### `get_locked_zarr_parts(zarr_root: str | Path) -> list[str]`

Return logical zarr parts that contain any read-only files.

Scans the top-level directories of the zarr store (and one level deeper
inside ``layers/`` to expose individual layer names) and returns those
whose file trees contain at least one file without owner-write permission.

Args:
    zarr_root: Root directory of the Zarr store.

Returns:
    Sorted list of locked part paths (e.g. ``["layers/raw", "obs"]``).
    An empty list means the store is fully writable.
    ``["."]`` means the entire store is locked.

### `read_parquet(file_path: str) -> pd.DataFrame`

Read a Parquet file.

Args:
    file_path: Path to Parquet file
    
Returns:
    DataFrame with file contents

### `read_yaml(file_path: str) -> dict`

Read a YAML file.

Args:
    file_path: Path to YAML file
    
Returns:
    Dictionary with YAML contents

### `set_zarr_parts_writable(zarr_root: str | Path, parts: list[str] | None = None, writable: bool = True, strict: bool = True) -> list[str]`

Set writable/read-only permissions for selected Zarr store parts.

Args:
    zarr_root: Root directory for the Zarr store.
    parts: Relative paths inside the store (e.g. "layers/raw", "obs").
        If None, applies to the full store.
    writable: If True, add owner write permission. If False, remove all
        write permissions.
    strict: If True, raise when a requested part does not exist.

Returns:
    Normalized part paths that were processed. "." means full store.

### `write_parquet(file_path: str, df: pd.DataFrame)`

Write DataFrame to Parquet file.

Args:
    file_path: Path to output Parquet file
    df: DataFrame to write

### `write_yaml(file_path: str, data: dict)`

Write data to a YAML file.

Args:
    file_path: Path to output YAML file
    data: Dictionary to write

## Classes

No public classes.
