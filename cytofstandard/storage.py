"""Storage utilities for cytofstandard."""

import hashlib
import os
import shutil
import stat
import yaml
import pandas as pd
from pathlib import Path


def compute_sha256(file_path: str) -> str:
    """Compute SHA256 hash of a file.
    
    Args:
        file_path: Path to the file
        
    Returns:
        Hex digest of SHA256 hash
    """
    sha256_hash = hashlib.sha256()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            sha256_hash.update(chunk)
    return sha256_hash.hexdigest()


def compute_file_hashes(file_paths: list[str]) -> dict[str, str]:
    """Compute SHA256 hashes for multiple files.
    
    Args:
        file_paths: List of file paths
        
    Returns:
        Dictionary mapping file path to hash
    """
    return {fp: compute_sha256(fp) for fp in file_paths}


def read_yaml(file_path: str) -> dict:
    """Read a YAML file.
    
    Args:
        file_path: Path to YAML file
        
    Returns:
        Dictionary with YAML contents
    """
    with open(file_path, "r") as f:
        return yaml.safe_load(f)


def write_yaml(file_path: str, data: dict):
    """Write data to a YAML file.
    
    Args:
        file_path: Path to output YAML file
        data: Dictionary to write
    """
    with open(file_path, "w") as f:
        yaml.dump(data, f, default_flow_style=False)


def read_parquet(file_path: str) -> pd.DataFrame:
    """Read a Parquet file.
    
    Args:
        file_path: Path to Parquet file
        
    Returns:
        DataFrame with file contents
    """
    return pd.read_parquet(file_path)


def write_parquet(file_path: str, df: pd.DataFrame):
    """Write DataFrame to Parquet file.
    
    Args:
        file_path: Path to output Parquet file
        df: DataFrame to write
    """
    df.to_parquet(file_path, index=False)


def ensure_directory(path: str):
    """Ensure a directory exists, creating it if necessary.
    
    Args:
        path: Path to directory
    """
    Path(path).mkdir(parents=True, exist_ok=True)


def get_directory_size(path: str) -> int:
    """Get total size of all files in a directory.
    
    Args:
        path: Path to directory
        
    Returns:
        Total size in bytes
    """
    total = 0
    for file_path in Path(path).rglob("*"):
        if file_path.is_file():
            total += file_path.stat().st_size
    return total


def copy_file(src: str, dst: str):
    """Copy a file from source to destination.
    
    Args:
        src: Source file path
        dst: Destination file path
    """
    Path(dst).parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def file_exists(file_path: str) -> bool:
    """Check if a file exists.
    
    Args:
        file_path: Path to file
        
    Returns:
        True if file exists
    """
    return Path(file_path).exists()


def directory_exists(dir_path: str) -> bool:
    """Check if a directory exists.
    
    Args:
        dir_path: Path to directory
        
    Returns:
        True if directory exists
    """
    return Path(dir_path).is_dir()


def copy_files_to_directory(
    source_files: list[str],
    dest_dir: str,
    preserve_name: bool = True,
) -> dict[str, str]:
    """Copy multiple files to a directory.
    
    Args:
        source_files: List of source file paths
        dest_dir: Destination directory
        preserve_name: Whether to preserve original filenames
        
    Returns:
        Dictionary mapping source path to destination path
    """
    dest_path = Path(dest_dir)
    dest_path.mkdir(parents=True, exist_ok=True)
    
    mappings = {}
    for src in source_files:
        src_path = Path(src)
        if preserve_name:
            dst = dest_path / src_path.name
        else:
            dst = dest_path / src_path.stem + src_path.suffix
        copy_file(src, str(dst))
        mappings[src] = str(dst)
    
    return mappings


def set_zarr_parts_writable(
    zarr_root: str | Path,
    parts: list[str] | None = None,
    writable: bool = True,
    strict: bool = True,
) -> list[str]:
    """Set writable/read-only permissions for selected Zarr store parts.

    Args:
        zarr_root: Root directory for the Zarr store.
        parts: Relative paths inside the store (e.g. "layers/raw", "obs").
            If None, applies to the full store.
        writable: If True, add owner write permission. If False, remove all
            write permissions.
        strict: If True, raise when a requested part does not exist.

    Returns:
        Normalized part paths that were processed. "." means full store.
    """
    root = Path(zarr_root)
    if not root.exists():
        raise FileNotFoundError(f"Zarr root not found: {root}")

    normalized_parts = _normalize_zarr_parts(parts)
    targets: list[tuple[str, Path]] = []
    root_resolved = root.resolve()

    for part in normalized_parts:
        target = root if part == "." else root / part
        target_resolved = target.resolve()

        if target_resolved != root_resolved and root_resolved not in target_resolved.parents:
            raise ValueError(f"Invalid zarr part outside store root: {part}")

        if not target.exists():
            if strict:
                raise FileNotFoundError(f"Zarr part not found: {part}")
            continue

        targets.append((part, target))

    seen = set()
    processed_parts = []
    for part, target in targets:
        if part in seen:
            continue
        seen.add(part)
        for path in _iter_paths_recursive(target):
            _set_single_path_writable(path, writable=writable)
        processed_parts.append(part)

    return processed_parts


def _normalize_zarr_parts(parts: list[str] | None) -> list[str]:
    """Return normalized relative part paths.

    "." represents the full store root.
    """
    if parts is None:
        return ["."]
    if not parts:
        raise ValueError("parts cannot be empty")

    normalized = []
    for part in parts:
        value = str(part).strip().replace("\\", "/")
        if not value:
            raise ValueError("zarr part paths must be non-empty")
        value = value.strip("/")
        normalized.append(value or ".")
    return normalized


def _iter_paths_recursive(root: Path):
    """Yield root and all descendants for permission updates."""
    yield root
    if root.is_dir():
        for child in root.rglob("*"):
            yield child


def _set_single_path_writable(path: Path, writable: bool) -> None:
    """Adjust a single path's write permissions."""
    mode = path.stat().st_mode
    if writable:
        new_mode = mode | stat.S_IWUSR
    else:
        new_mode = mode & ~stat.S_IWUSR & ~stat.S_IWGRP & ~stat.S_IWOTH

    if new_mode != mode:
        os.chmod(path, new_mode)
