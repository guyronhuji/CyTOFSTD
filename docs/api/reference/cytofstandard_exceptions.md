# `cytofstandard.exceptions`

- Source: `cytofstandard/exceptions.py`

Custom exception classes for cytofstandard package.

## Top-level Functions

No public top-level functions.

## Classes

### `CytofStandardError`

- Inherits: `Exception`

Base exception for all cytofstandard errors.

#### Methods

No public methods.

### `FileFormatError`

- Inherits: `CytofStandardError`

Raised when file format is invalid.

#### Methods

No public methods.

### `IngestionError`

- Inherits: `CytofStandardError`

Raised when ingestion fails.

#### Methods

No public methods.

### `MarkerValidationError`

- Inherits: `CytofStandardError`

Raised when marker validation fails.

#### Methods

No public methods.

### `MetadataValidationError`

- Inherits: `CytofStandardError`

Raised when metadata validation fails.

#### Methods

No public methods.

### `ProjectExistsError`

- Inherits: `CytofStandardError`

Raised when trying to create a project that already exists.

#### Methods

No public methods.

### `ProjectValidationError`

- Inherits: `CytofStandardError`

Raised when project validation fails.

#### Methods

No public methods.

### `RunExistsError`

- Inherits: `CytofStandardError`

Raised when trying to add a run that already exists.

#### Methods

No public methods.

### `RunNotFoundError`

- Inherits: `CytofStandardError`

Raised when a run is not found in the project.

#### Methods

No public methods.

### `RunNotIngestedError`

- Inherits: `CytofStandardError`

Raised when trying to read a run that has not been ingested.

#### Methods

No public methods.

### `RunValidationError`

- Inherits: `CytofStandardError`

Raised when run validation fails.

#### Methods

No public methods.
