"""Configuration for cytofstandard."""

# Package version
__version__ = "0.1.0"

# Default FCS reader to use
DEFAULT_FCS_READER = "fcsparser"

# Default batch size for processing
DEFAULT_BATCH_SIZE = 10000

# Marker validation defaults
DEFAULT_STRICT_MARKERS = True
DEFAULT_ALLOW_EXTRA_MARKERS = False

# Storage defaults
DEFAULT_STORAGE_FORMAT = "anndata_zarr"
DEFAULT_ZARR_COMPRESSOR = "blosc"
DEFAULT_ZARR_COMPRESSION_LEVEL = 1

# UUID strategy
DEFAULT_UUID_STRATEGY = "uuid5_project_run_filehash_eventindex"

# File naming patterns
RUN_DIR_PATTERN = "run_{id:03d}"
ZARR_FILE_PATTERN = "{run_id}.zarr"

# Status values
STATUS_REGISTERED = "registered"
STATUS_INGESTED = "ingested"
STATUS_FAILED_INGESTION = "failed_ingestion"
STATUS_VALUES = [STATUS_REGISTERED, STATUS_INGESTED, STATUS_FAILED_INGESTION]
