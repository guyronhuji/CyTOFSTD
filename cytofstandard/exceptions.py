"""Custom exception classes for cytofstandard package."""


class CytofStandardError(Exception):
    """Base exception for all cytofstandard errors."""


class ProjectExistsError(CytofStandardError):
    """Raised when trying to create a project that already exists."""


class ProjectValidationError(CytofStandardError):
    """Raised when project validation fails."""


class RunExistsError(CytofStandardError):
    """Raised when trying to add a run that already exists."""


class RunNotFoundError(CytofStandardError):
    """Raised when a run is not found in the project."""


class RunValidationError(CytofStandardError):
    """Raised when run validation fails."""


class RunNotIngestedError(CytofStandardError):
    """Raised when trying to read a run that has not been ingested."""


class MetadataValidationError(CytofStandardError):
    """Raised when metadata validation fails."""


class MarkerValidationError(CytofStandardError):
    """Raised when marker validation fails."""


class IngestionError(CytofStandardError):
    """Raised when ingestion fails."""


class FileFormatError(CytofStandardError):
    """Raised when file format is invalid."""
