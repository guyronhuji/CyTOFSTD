"""Tests for CSV ingestion."""

import pytest
import pandas as pd
from pathlib import Path
from cytofstandard import Project
import tempfile
import shutil


@pytest.fixture
def csv_project(temp_dir, example_standard_markers_path, example_marker_aliases_path):
    """Create a project for CSV ingestion tests."""
    project_path = temp_dir / "csv_project"
    project = Project.create(
        path=str(project_path),
        project_id="CSV_TEST",
        project_name="CSV Test",
        standard_marker_file=str(example_standard_markers_path),
        marker_alias_file=str(example_marker_aliases_path),
    )
    return project


def test_csv_ingestion(temp_dir, csv_project, simple_csv_file, sample_metadata_file):
    """Test basic CSV ingestion."""
    run = csv_project.add_run(run_id="run_001")

    raw_dir = csv_project.path / "runs" / "run_001" / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy(simple_csv_file, raw_dir / simple_csv_file.name)

    sample_df = pd.read_csv(sample_metadata_file)
    sample_df["file_name"] = simple_csv_file.name
    sample_df.to_csv(sample_metadata_file, index=False)

    run.ingest(
        files=[str(raw_dir / simple_csv_file.name)],
        sample_metadata=str(sample_metadata_file),
        copy_raw=False,
        strict_markers=True,
    )

    assert run.is_ingested()


def test_ingested_data_structure(
    temp_dir, csv_project, simple_csv_file, sample_metadata_file
):
    """Test that ingested data has correct structure."""
    run = csv_project.add_run(run_id="run_001")

    raw_dir = csv_project.path / "runs" / "run_001" / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy(simple_csv_file, raw_dir / simple_csv_file.name)

    sample_df = pd.read_csv(sample_metadata_file)
    sample_df["file_name"] = simple_csv_file.name
    sample_df.to_csv(sample_metadata_file, index=False)

    run.ingest(
        files=[str(raw_dir / simple_csv_file.name)],
        sample_metadata=str(sample_metadata_file),
        copy_raw=False,
        strict_markers=True,
    )

    adata = run.read_adata()

    assert "cell_uuid" in adata.obs.columns
    assert "sample_id" in adata.obs.columns
    assert "line_id" in adata.obs.columns
    assert "source_file" in adata.obs.columns

    assert "standard_marker_name" in adata.var.columns

    assert "raw" in adata.layers

    assert adata.n_obs == 3
    assert adata.n_vars == 3


def test_ingested_data_has_run_metadata(
    temp_dir, csv_project, simple_csv_file, sample_metadata_file
):
    """Test that ingested data has run-level metadata."""
    run = csv_project.add_run(
        run_id="run_001",
        run_name="Test Run",
    )

    raw_dir = csv_project.path / "runs" / "run_001" / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy(simple_csv_file, raw_dir / simple_csv_file.name)

    sample_df = pd.read_csv(sample_metadata_file)
    sample_df["file_name"] = simple_csv_file.name
    sample_df.to_csv(sample_metadata_file, index=False)

    run.ingest(
        files=[str(raw_dir / simple_csv_file.name)],
        sample_metadata=str(sample_metadata_file),
        copy_raw=False,
        strict_markers=True,
    )

    adata = run.read_adata()

    assert "run" in adata.uns
    assert adata.uns["run"]["run_id"] == "run_001"
    assert adata.uns["run"]["run_name"] == "Test Run"


def test_ingested_data_has_project_metadata(
    temp_dir, csv_project, simple_csv_file, sample_metadata_file
):
    """Test that ingested data has project-level metadata."""
    run = csv_project.add_run(run_id="run_001")

    raw_dir = csv_project.path / "runs" / "run_001" / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy(simple_csv_file, raw_dir / simple_csv_file.name)

    sample_df = pd.read_csv(sample_metadata_file)
    sample_df["file_name"] = simple_csv_file.name
    sample_df.to_csv(sample_metadata_file, index=False)

    run.ingest(
        files=[str(raw_dir / simple_csv_file.name)],
        sample_metadata=str(sample_metadata_file),
        copy_raw=False,
        strict_markers=True,
    )

    adata = run.read_adata()

    assert "project" in adata.uns
    assert adata.uns["project"]["project_id"] == "CSV_TEST"


def test_csv_ingestion_drop_columns(
    temp_dir, csv_project, simple_csv_file, sample_metadata_file
):
    """Test ingestion with explicit dropped columns."""
    run = csv_project.add_run(run_id="run_001")

    raw_dir = csv_project.path / "runs" / "run_001" / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy(simple_csv_file, raw_dir / simple_csv_file.name)

    sample_df = pd.read_csv(sample_metadata_file)
    sample_df["file_name"] = simple_csv_file.name
    sample_df.to_csv(sample_metadata_file, index=False)

    run.ingest(
        files=[str(raw_dir / simple_csv_file.name)],
        sample_metadata=str(sample_metadata_file),
        copy_raw=False,
        strict_markers=True,
        drop_columns=["ECad"],
    )

    adata = run.read_adata()

    assert adata.n_vars == 2
    assert "ECad" not in adata.var_names
    assert adata.uns["ingestion"]["drop_columns"] == ["ECad"]


def test_csv_ingestion_with_alias_column_names(
    temp_dir, csv_project, sample_metadata_file
):
    """Test ingestion renames alias marker columns to standard names first."""
    run = csv_project.add_run(run_id="run_001")

    alias_csv = temp_dir / "alias_markers.csv"
    alias_csv.write_text(
        "H3,H3K27me3,E-Cadherin,cell_type\n"
        "100,200,300,T_cell\n"
        "150,250,350,B_cell\n"
        "120,220,320,Macrophage\n"
    )

    raw_dir = csv_project.path / "runs" / "run_001" / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy(alias_csv, raw_dir / alias_csv.name)

    sample_df = pd.read_csv(sample_metadata_file)
    sample_df["file_name"] = alias_csv.name
    sample_df.to_csv(sample_metadata_file, index=False)

    run.ingest(
        files=[str(raw_dir / alias_csv.name)],
        sample_metadata=str(sample_metadata_file),
        copy_raw=False,
        strict_markers=True,
    )

    adata = run.read_adata()

    assert "ECad" in adata.var_names
    assert "E-Cadherin" not in adata.var_names


def test_ingestion_handles_mixed_bool_var_columns(
    temp_dir, example_marker_aliases_path
):
    """Test ingestion when var metadata columns mix bool and missing values."""
    standard_file = temp_dir / "standard_markers_mixed_bool.csv"
    standard_file.write_text(
        "standard_marker_name,marker_class,is_core_histone\n"
        "H3,histone,True\n"
        "ECad,identity,\n"
    )

    project_path = temp_dir / "mixed_bool_project"
    project = Project.create(
        path=str(project_path),
        project_id="MIXED_BOOL_TEST",
        project_name="Mixed Bool Test",
        standard_marker_file=str(standard_file),
        marker_alias_file=str(example_marker_aliases_path),
    )

    run = project.add_run(run_id="run_001")

    csv_path = temp_dir / "sample.csv"
    csv_path.write_text("H3,ECad\n100,300\n150,350\n")

    sample_meta = temp_dir / "sample_metadata.csv"
    sample_meta.write_text("file_name,sample_id,line_id\nsample.csv,S001,MCF7\n")

    run.ingest(
        files=[str(csv_path)],
        sample_metadata=str(sample_meta),
        copy_raw=False,
        strict_markers=True,
    )

    adata = run.read_adata()
    assert "is_core_histone" in adata.var.columns


def test_csv_ingestion_strips_metal_prefixes_in_marker_names(
    temp_dir,
    csv_project,
    sample_metadata_file,
):
    """Marker names like 168Er_H3K27me3 are sanitized before registry matching."""
    run = csv_project.add_run(run_id="run_001")

    prefixed_csv = temp_dir / "prefixed_markers.csv"
    prefixed_csv.write_text(
        "168Er_H3K27me3,141Pr_EpCAM,163Dy_ER\n200,300,50\n250,350,55\n"
    )

    raw_dir = csv_project.path / "runs" / "run_001" / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy(prefixed_csv, raw_dir / prefixed_csv.name)

    sample_df = pd.read_csv(sample_metadata_file)
    sample_df["file_name"] = prefixed_csv.name
    sample_df.to_csv(sample_metadata_file, index=False)

    run.ingest(
        files=[str(raw_dir / prefixed_csv.name)],
        sample_metadata=str(sample_metadata_file),
        copy_raw=False,
        strict_markers=True,
    )

    adata = run.read_adata()
    assert "H3K27me3" in adata.var_names
    assert "EpCAM" in adata.var_names
    assert "ER" in adata.var_names


def test_run_marker_variables_populated_after_ingestion(
    temp_dir,
    csv_project,
    simple_csv_file,
    sample_metadata_file,
):
    """Run exposes marker category variables for markers present in run."""
    run = csv_project.add_run(run_id="run_001")

    raw_dir = csv_project.path / "runs" / "run_001" / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy(simple_csv_file, raw_dir / simple_csv_file.name)

    sample_df = pd.read_csv(sample_metadata_file)
    sample_df["file_name"] = simple_csv_file.name
    sample_df.to_csv(sample_metadata_file, index=False)

    run.ingest(
        files=[str(raw_dir / simple_csv_file.name)],
        sample_metadata=str(sample_metadata_file),
        copy_raw=False,
        strict_markers=True,
    )

    assert "H3" in run.markers_all
    assert "ECad" in run.markers_all
    assert "H3" in run.markers_core_histones
    assert "ECad" in run.markers_identity
    assert "H3K27me3" in run.markers_functional
    assert isinstance(run.markers_by_class, dict)
    assert "histone" in run.markers_by_class


def test_run_uses_intra_extra_column_for_marker_partition(
    temp_dir, example_marker_aliases_path
):
    """Run intra/extra marker lists should follow intra_extra registry column."""
    standard_file = temp_dir / "standard_markers_intra_extra.csv"
    standard_file.write_text(
        "standard_marker_name,intra_extra,marker_class,is_core_histone,is_identity_marker,is_functional_marker,is_qc_marker\n"
        "H3,Intra,unknown,true,false,false,false\n"
        "ECad,Extra,unknown,false,true,false,false\n"
    )

    project_path = temp_dir / "intra_extra_project"
    project = Project.create(
        path=str(project_path),
        project_id="INTRA_EXTRA_TEST",
        project_name="Intra Extra Test",
        standard_marker_file=str(standard_file),
        marker_alias_file=str(example_marker_aliases_path),
    )

    run = project.add_run(run_id="run_001")

    csv_path = temp_dir / "sample.csv"
    csv_path.write_text("H3,ECad\n100,300\n150,350\n")

    sample_meta = temp_dir / "sample_metadata.csv"
    sample_meta.write_text("file_name,sample_id,line_id\nsample.csv,S001,MCF7\n")

    run.ingest(
        files=[str(csv_path)],
        sample_metadata=str(sample_meta),
        copy_raw=False,
        strict_markers=True,
    )

    assert "H3" in run.markers_intracellular
    assert "ECad" in run.markers_extracellular
