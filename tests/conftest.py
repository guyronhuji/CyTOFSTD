"""Pytest fixtures for cytofstandard tests."""

import pytest
import shutil
import tempfile
from pathlib import Path
from cytofstandard import Project
from cytofstandard.exceptions import RunNotFoundError


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test artifacts."""
    tmpdir = tempfile.mkdtemp()
    yield Path(tmpdir)
    shutil.rmtree(tmpdir)


@pytest.fixture
def example_project_path(temp_dir):
    """Create an example project for testing."""
    project_path = temp_dir / "test_project"
    standard_file = Path(__file__).parent.parent / "example_standard_markers.csv"
    alias_file = Path(__file__).parent.parent / "example_marker_aliases.yaml"
    
    project = Project.create(
        path=str(project_path),
        project_id="TEST_PROJ_001",
        project_name="Test Project",
        standard_marker_file=str(standard_file),
        marker_alias_file=str(alias_file),
    )
    
    return project_path


@pytest.fixture
def example_project(example_project_path):
    """Load example project for testing."""
    return Project.load(str(example_project_path))


@pytest.fixture
def simple_csv_file(temp_dir):
    """Create a simple CSV file for testing."""
    csv_path = temp_dir / "simple_test.csv"
    csv_path.write_text(
        "H3,H3K27me3,ECad,cell_type\n"
        "100,200,300,T_cell\n"
        "150,250,350,B_cell\n"
        "120,220,320,Macrophage\n"
    )
    return csv_path


@pytest.fixture
def sample_metadata_file(temp_dir, simple_csv_file):
    """Create a sample metadata file."""
    meta_path = temp_dir / "sample_metadata.csv"
    meta_path.write_text(
        "file_name,sample_id,line_id,condition,replicate_id\n"
        f"{simple_csv_file.name},S001,MCF7,control,rep1\n"
    )
    return meta_path


@pytest.fixture
def simple_project_path(temp_dir):
    """Create a simple test project."""
    project_path = temp_dir / "simple_project"
    standard_file = Path(__file__).parent.parent / "example_standard_markers.csv"
    alias_file = Path(__file__).parent.parent / "example_marker_aliases.yaml"
    
    project = Project.create(
        path=str(project_path),
        project_id="SIMPLE_TEST",
        project_name="Simple Test",
        standard_marker_file=str(standard_file),
        marker_alias_file=str(alias_file),
    )
    
    return project_path


@pytest.fixture
def simple_project(simple_project_path):
    """Load simple test project."""
    return Project.load(str(simple_project_path))


@pytest.fixture
def example_standard_markers_path():
    """Get path to example standard markers file."""
    return Path(__file__).parent.parent / "example_standard_markers.csv"


@pytest.fixture
def example_marker_aliases_path():
    """Get path to example marker aliases file."""
    return Path(__file__).parent.parent / "example_marker_aliases.yaml"
