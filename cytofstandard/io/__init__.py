"""Input/output adapters for file formats."""

from cytofstandard.io.fcs import read_fcs
from cytofstandard.io.csv import read_csv
from cytofstandard.io.parquet import read_parquet

__all__ = ["read_fcs", "read_csv", "read_parquet"]
