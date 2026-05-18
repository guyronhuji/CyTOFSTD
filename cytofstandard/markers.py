"""Marker registry for standardizing marker names."""

import yaml
import pandas as pd
from pathlib import Path
from typing import Optional

from cytofstandard.exceptions import MarkerValidationError


class MarkerRegistry:
    """Registry for standard marker names and aliases."""
    
    def __init__(
        self,
        standard_markers: pd.DataFrame,
        aliases: dict[str, list[str]],
    ):
        """Initialize marker registry.
        
        Args:
            standard_markers: DataFrame with standard marker info
            aliases: Dictionary mapping standard markers to their aliases
        """
        self.standard_markers = standard_markers
        self.aliases = aliases
        
        # Build reverse lookup: alias -> standard marker
        self._alias_to_standard = {}
        for standard, aliases_list in aliases.items():
            for alias in aliases_list:
                self._alias_to_standard[alias] = standard
    
    @classmethod
    def from_project(cls, project) -> "MarkerRegistry":
        """Create a MarkerRegistry from a Project.
        
        Args:
            project: Project instance
            
        Returns:
            MarkerRegistry loaded from project files
        """
        standard_file = project.metadata_dir / "standard_markers.parquet"
        alias_file = project.metadata_dir / "marker_aliases.yaml"
        
        standard_markers = pd.read_parquet(standard_file)
        
        with open(alias_file, "r") as f:
            aliases_content = f.read()
            if not aliases_content.strip():
                aliases = {}
            else:
                aliases = yaml.safe_load(aliases_content)
        
        if aliases is None:
            aliases = {}
        
        return cls(standard_markers, aliases)
    
    def standardize_marker_names(
        self,
        observed_names: list[str],
        strict: bool = True,
    ) -> pd.DataFrame:
        """Standardize a list of marker names.
        
        Args:
            observed_names: List of observed marker names
            strict: Whether to fail on unknown markers
            
        Returns:
            DataFrame with standardization results
        """
        results = []
        
        for name in observed_names:
            # Check for exact match first
            if name in self.standard_markers["standard_marker_name"].values:
                results.append({
                    "original_marker_name": name,
                    "standard_marker_name": name,
                    "mapping_status": "matched_standard",
                    "mapping_source": "exact",
                    "message": None,
                })
            # Check for alias match
            elif name in self._alias_to_standard:
                standard = self._alias_to_standard[name]
                results.append({
                    "original_marker_name": name,
                    "standard_marker_name": standard,
                    "mapping_status": "matched_alias",
                    "mapping_source": "alias",
                    "message": None,
                })
            else:
                results.append({
                    "original_marker_name": name,
                    "standard_marker_name": None,
                    "mapping_status": "unknown",
                    "mapping_source": None,
                    "message": f"Unknown marker: {name}",
                })
        
        result_df = pd.DataFrame(results)
        
        # Check for duplicates after mapping
        matched = result_df[
            result_df["mapping_status"].isin(["matched_standard", "matched_alias"])
        ]
        if len(matched) > 0:
            duplicates = matched[matched.duplicated(subset=["standard_marker_name"], keep=False)]
            if len(duplicates) > 0:
                for standard_name in duplicates["standard_marker_name"].unique():
                    dup_rows = duplicates[duplicates["standard_marker_name"] == standard_name]
                    original_names = dup_rows["original_marker_name"].tolist()
                    result_df.loc[
                        result_df["standard_marker_name"] == standard_name,
                        "mapping_status"
                    ] = "duplicate_after_mapping"
                    result_df.loc[
                        result_df["standard_marker_name"] == standard_name,
                        "message"
                    ] = f"Multiple columns map to {standard_name}: {original_names}"
        
        # Validation
        unknown = result_df[result_df["mapping_status"] == "unknown"]
        if strict and len(unknown) > 0:
            unknown_names = unknown["original_marker_name"].tolist()
            raise MarkerValidationError(
                f"Unknown marker names found: {unknown_names}. "
                "Add to marker_aliases.yaml or set strict_markers=False."
            )
        
        # Check for duplicates in strict mode
        if strict:
            duplicate = result_df[result_df["mapping_status"] == "duplicate_after_mapping"]
            if len(duplicate) > 0:
                raise MarkerValidationError(
                    f"Duplicate marker mappings found: {duplicate['standard_marker_name'].tolist()}. "
                    "This is ambiguous and ingestion was stopped."
                )
        
        return result_df
    
    def get_standard_markers(self) -> list[str]:
        """Get list of all standard marker names.
        
        Returns:
            List of standard marker names
        """
        return self.standard_markers["standard_marker_name"].tolist()
    
    def get_marker_info(self, marker_name: str) -> Optional[dict]:
        """Get info about a standard marker.
        
        Args:
            marker_name: Standard marker name
            
        Returns:
            Marker info dictionary or None if not found
        """
        marker = self.standard_markers[
            self.standard_markers["standard_marker_name"] == marker_name
        ]
        if len(marker) > 0:
            return marker.iloc[0].to_dict()
        return None
