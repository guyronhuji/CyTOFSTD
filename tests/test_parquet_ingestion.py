"""Tests for Parquet ingestion."""

import pandas as pd
import shutil

from cytofstandard import Project


def test_parquet_ingestion(
    temp_dir, example_standard_markers_path, example_marker_aliases_path
):
    project_path = temp_dir / "parquet_project"
    project = Project.create(
        path=str(project_path),
        project_id="PARQUET_TEST",
        project_name="Parquet Test",
        standard_marker_file=str(example_standard_markers_path),
        marker_alias_file=str(example_marker_aliases_path),
    )

    run = project.add_run(run_id="run_001")

    parquet_path = temp_dir / "sample.parquet"
    df = pd.DataFrame(
        {
            "H3": [100, 150, 120],
            "H3K27me3": [200, 250, 220],
            "ECad": [300, 350, 320],
        }
    )
    df.to_parquet(parquet_path, index=False)

    sample_meta = temp_dir / "sample_metadata.csv"
    sample_meta.write_text(
        "file_name,sample_id,line_id\n"
        f"{parquet_path.name},S001,MCF7\n"
    )

    raw_dir = project.path / "runs" / "run_001" / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy(parquet_path, raw_dir / parquet_path.name)

    run.ingest(
        files=[str(raw_dir / parquet_path.name)],
        sample_metadata=str(sample_meta),
        copy_raw=False,
        strict_markers=True,
    )

    adata = run.read_adata()
    assert run.is_ingested()
    assert adata.n_obs == 3
    assert set(adata.var_names) == {"H3", "H3K27me3", "ECad"}
