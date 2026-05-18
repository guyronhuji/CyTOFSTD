"""Tests for marker standardization."""

import pytest
import pandas as pd
from cytofstandard import MarkerRegistry
from cytofstandard.exceptions import MarkerValidationError
from pathlib import Path


@pytest.fixture
def marker_registry():
    """Create a marker registry for testing."""
    standard_markers = pd.DataFrame({
        "standard_marker_name": ["H3", "H3K27me3", "ECad"],
        "marker_class": ["histone", "histone", "identity"],
    })
    
    aliases = {
        "H3": ["H3", "Histone H3", "Ir-H3"],
        "H3K27me3": ["H3K27me3", "H3K27me3_176"],
        "ECad": ["ECad", "E-Cadherin", "CDH1"],
    }
    
    return MarkerRegistry(standard_markers, aliases)


def test_exact_match(marker_registry):
    """Test exact match to standard marker."""
    result = marker_registry.standardize_marker_names(["H3"], strict=False)
    
    assert len(result) == 1
    assert result.iloc[0]["original_marker_name"] == "H3"
    assert result.iloc[0]["standard_marker_name"] == "H3"
    assert result.iloc[0]["mapping_status"] == "matched_standard"


def test_alias_match(marker_registry):
    """Test match via alias."""
    result = marker_registry.standardize_marker_names(["E-Cadherin"], strict=False)
    
    assert len(result) == 1
    assert result.iloc[0]["original_marker_name"] == "E-Cadherin"
    assert result.iloc[0]["standard_marker_name"] == "ECad"
    assert result.iloc[0]["mapping_status"] == "matched_alias"


def test_unknown_marker_strict(marker_registry):
    """Test unknown marker fails in strict mode."""
    with pytest.raises(MarkerValidationError):
        marker_registry.standardize_marker_names(["UnknownMarker"], strict=True)


def test_unknown_marker_non_strict(marker_registry):
    """Test unknown marker passes in non-strict mode."""
    result = marker_registry.standardize_marker_names(["UnknownMarker"], strict=False)
    
    assert len(result) == 1
    assert result.iloc[0]["mapping_status"] == "unknown"


def test_duplicate_mapping(marker_registry):
    """Test duplicate mapping raises error."""
    with pytest.raises(MarkerValidationError) as exc_info:
        marker_registry.standardize_marker_names(
            ["E-Cadherin", "ECad"],
            strict=True,
        )
    
    assert "duplicate" in str(exc_info.value).lower() or "ambiguous" in str(exc_info.value).lower()


def test_get_standard_markers(marker_registry):
    """Test getting list of standard markers."""
    markers = marker_registry.get_standard_markers()
    
    assert set(markers) == {"H3", "H3K27me3", "ECad"}


def test_get_marker_info(marker_registry):
    """Test getting marker info."""
    info = marker_registry.get_marker_info("H3")
    
    assert info is not None
    assert info["standard_marker_name"] == "H3"
    assert info["marker_class"] == "histone"


def test_get_marker_info_not_found(marker_registry):
    """Test getting info for non-existent marker."""
    info = marker_registry.get_marker_info("NonExistent")
    
    assert info is None
