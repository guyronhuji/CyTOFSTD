"""Tests for metadata validation."""

import pytest
import pandas as pd
from pathlib import Path
from cytofstandard import Project
from cytofstandard.exceptions import (
    MetadataValidationError,
    MarkerValidationError,
)


def test_missing_file_in_metadata(
    temp_dir, example_standard_markers_path, example_marker_aliases_path
):
    """Test that missing file in metadata raises error."""
    project_path = temp_dir / "test_project"

    project = Project.create(
        path=str(project_path),
        project_id="META_TEST",
        project_name="Metadata Test",
        standard_marker_file=str(example_standard_markers_path),
        marker_alias_file=str(example_marker_aliases_path),
    )

    dummy_file = temp_dir / "dummy.csv"
    dummy_file.write_text("H3,H3K27me3,ECad\n1,2,3\n")

    meta_path = temp_dir / "metadata.csv"
    meta_path.write_text("file_name,sample_id,line_id\nsample_A.csv,S001,MCF7\n")

    run = project.add_run(run_id="run_001")

    with pytest.raises(MetadataValidationError) as exc_info:
        run.ingest(
            files=[str(dummy_file)],
            sample_metadata=str(meta_path),
            copy_raw=False,
            strict_markers=False,
        )

    assert "missing from sample_metadata" in str(exc_info.value)


def test_metadata_extra_rows_without_input_samples_are_allowed(
    temp_dir,
    example_standard_markers_path,
    example_marker_aliases_path,
):
    """Metadata can contain rows for files that are not included in this ingest."""
    project_path = temp_dir / "test_project"

    project = Project.create(
        path=str(project_path),
        project_id="META_TEST2",
        project_name="Metadata Test",
        standard_marker_file=str(example_standard_markers_path),
        marker_alias_file=str(example_marker_aliases_path),
    )

    file_a = temp_dir / "sample_A.csv"
    file_a.write_text("H3,H3K27me3,ECad\n1,2,3\n")

    meta_path = temp_dir / "metadata.csv"
    meta_path.write_text(
        "file_name,sample_id,line_id\nsample_A.csv,S001,MCF7\nsample_C.csv,S002,T47D\n"
    )

    run = project.add_run(run_id="run_001")

    run.ingest(
        files=[str(file_a)],
        sample_metadata=str(meta_path),
        copy_raw=False,
        strict_markers=False,
    )

    assert run.is_ingested()


def test_missing_sample_id(
    temp_dir, example_standard_markers_path, example_marker_aliases_path
):
    """Test that missing sample_id raises error."""
    project_path = temp_dir / "test_project"

    project = Project.create(
        path=str(project_path),
        project_id="META_TEST3",
        project_name="Metadata Test",
        standard_marker_file=str(example_standard_markers_path),
        marker_alias_file=str(example_marker_aliases_path),
    )

    file_path = temp_dir / "sample.csv"
    file_path.write_text("H3,H3K27me3,ECad\n1,2,3\n")

    meta_path = temp_dir / "metadata.csv"
    meta_path.write_text("file_name,sample_id,line_id\nsample.csv,,MCF7\n")

    run = project.add_run(run_id="run_001")

    with pytest.raises(MetadataValidationError):
        run.ingest(
            files=[str(file_path)],
            sample_metadata=str(meta_path),
            copy_raw=False,
            strict_markers=False,
        )


def test_missing_line_id(
    temp_dir, example_standard_markers_path, example_marker_aliases_path
):
    """Test that missing line_id raises error."""
    project_path = temp_dir / "test_project"

    project = Project.create(
        path=str(project_path),
        project_id="META_TEST4",
        project_name="Metadata Test",
        standard_marker_file=str(example_standard_markers_path),
        marker_alias_file=str(example_marker_aliases_path),
    )

    file_path = temp_dir / "sample.csv"
    file_path.write_text("H3,H3K27me3,ECad\n1,2,3\n")

    meta_path = temp_dir / "metadata.csv"
    meta_path.write_text("file_name,sample_id,line_id\nsample.csv,S001,\n")

    run = project.add_run(run_id="run_001")

    with pytest.raises(MetadataValidationError):
        run.ingest(
            files=[str(file_path)],
            sample_metadata=str(meta_path),
            copy_raw=False,
            strict_markers=False,
        )


def test_sample_maps_to_multiple_lines(
    temp_dir, example_standard_markers_path, example_marker_aliases_path
):
    """Test that sample mapping to multiple lines raises error."""
    project_path = temp_dir / "test_project"

    project = Project.create(
        path=str(project_path),
        project_id="META_TEST5",
        project_name="Metadata Test",
        standard_marker_file=str(example_standard_markers_path),
        marker_alias_file=str(example_marker_aliases_path),
    )

    # Create multiple CSV files with same sample but different lines
    file1 = temp_dir / "sample1.csv"
    file1.write_text("H3,H3K27me3,ECad\n1,2,3\n")
    file2 = temp_dir / "sample2.csv"
    file2.write_text("H3,H3K27me3,ECad\n4,5,6\n")

    # S001 maps to both MCF7 and T47D
    meta_path = temp_dir / "metadata.csv"
    meta_path.write_text(
        "file_name,sample_id,line_id\nsample1.csv,S001,MCF7\nsample2.csv,S001,T47D\n"
    )

    run = project.add_run(run_id="run_001")

    with pytest.raises(MetadataValidationError) as exc_info:
        run.ingest(
            files=[str(file1), str(file2)],
            sample_metadata=str(meta_path),
            copy_raw=False,
            strict_markers=False,
        )

    assert "mapped to multiple lines" in str(exc_info.value)


def test_unknown_marker_in_strict_mode(
    temp_dir, example_standard_markers_path, example_marker_aliases_path
):
    """Test that unknown marker fails in strict mode."""
    project_path = temp_dir / "test_project"

    project = Project.create(
        path=str(project_path),
        project_id="META_TEST6",
        project_name="Metadata Test",
        standard_marker_file=str(example_standard_markers_path),
        marker_alias_file=str(example_marker_aliases_path),
    )

    file_path = temp_dir / "sample.csv"
    file_path.write_text("H3,UnknownMarker,ECad\n1,2,3\n")

    meta_path = temp_dir / "metadata.csv"
    meta_path.write_text("file_name,sample_id,line_id\nsample.csv,S001,MCF7\n")

    run = project.add_run(run_id="run_001")

    with pytest.raises(MarkerValidationError):
        run.ingest(
            files=[str(file_path)],
            sample_metadata=str(meta_path),
            copy_raw=False,
            strict_markers=True,
        )


def test_unknown_marker_in_non_strict_mode(
    temp_dir, example_standard_markers_path, example_marker_aliases_path
):
    """Test that unknown marker passes in non-strict mode."""
    project_path = temp_dir / "test_project"

    project = Project.create(
        path=str(project_path),
        project_id="META_TEST7",
        project_name="Metadata Test",
        standard_marker_file=str(example_standard_markers_path),
        marker_alias_file=str(example_marker_aliases_path),
    )

    file_path = temp_dir / "sample.csv"
    file_path.write_text("H3,UnknownMarker,ECad\n1,2,3\n")

    meta_path = temp_dir / "metadata.csv"
    meta_path.write_text("file_name,sample_id,line_id\nsample.csv,S001,MCF7\n")

    run = project.add_run(run_id="run_001")

    run.ingest(
        files=[str(file_path)],
        sample_metadata=str(meta_path),
        copy_raw=False,
        strict_markers=False,
    )

    assert run.is_ingested()
