# `cytofstandard.markers`

- Source: `cytofstandard/markers.py`

Marker registry for standardizing marker names.

## Top-level Functions

No public top-level functions.

## Classes

### `MarkerRegistry`

Registry for standard marker names and aliases.

#### Methods

##### `from_project(cls, project) -> 'MarkerRegistry'`

- Decorators: `classmethod`

Create a MarkerRegistry from a Project.

Args:
    project: Project instance
    
Returns:
    MarkerRegistry loaded from project files

##### `get_marker_info(self, marker_name: str) -> Optional[dict]`

Get info about a standard marker.

Args:
    marker_name: Standard marker name
    
Returns:
    Marker info dictionary or None if not found

##### `get_standard_markers(self) -> list[str]`

Get list of all standard marker names.

Returns:
    List of standard marker names

##### `standardize_marker_names(self, observed_names: list[str], strict: bool = True) -> pd.DataFrame`

Standardize a list of marker names.

Args:
    observed_names: List of observed marker names
    strict: Whether to fail on unknown markers
    
Returns:
    DataFrame with standardization results
