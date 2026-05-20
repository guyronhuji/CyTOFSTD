"""Provenance logging utilities for cytofstandard."""

import json
from datetime import datetime
from pathlib import Path


def get_timestamp() -> str:
    """Return current timestamp in ISO format."""
    return datetime.utcnow().isoformat()


class ProvenanceLogger:
    """Logger for provenance events."""

    def __init__(self, log_path: str):
        """Initialize logger.

        Args:
            log_path: Path to the JSONL log file
        """
        self.log_path = Path(log_path)
        self.log_path.parent.mkdir(parents=True, exist_ok=True)

    def log(self, event_type: str, data: dict):
        """Log an event.

        Args:
            event_type: Type of event (e.g., 'project_created', 'run_ingested')
            data: Event data dictionary
        """
        entry = {
            "timestamp": get_timestamp(),
            "event_type": event_type,
            **data,
        }
        with open(self.log_path, "a") as f:
            f.write(json.dumps(entry) + "\n")


def log_project_created(
    logger: ProvenanceLogger,
    project_id: str,
    project_name: str,
):
    """Log project creation event."""
    logger.log(
        event_type="project_created",
        data={
            "project_id": project_id,
            "project_name": project_name,
        },
    )


def log_project_loaded(
    logger: ProvenanceLogger,
    project_id: str,
):
    """Log project load event."""
    logger.log(
        event_type="project_loaded",
        data={
            "project_id": project_id,
        },
    )


def log_run_registered(
    logger: ProvenanceLogger,
    project_id: str,
    run_id: str,
    run_name: str,
):
    """Log run registration event."""
    logger.log(
        event_type="run_registered",
        data={
            "project_id": project_id,
            "run_id": run_id,
            "run_name": run_name,
        },
    )


def log_run_loaded(
    logger: ProvenanceLogger,
    project_id: str,
    run_id: str,
):
    """Log run load event."""
    logger.log(
        event_type="run_loaded",
        data={
            "project_id": project_id,
            "run_id": run_id,
        },
    )


def log_run_removed(
    logger: ProvenanceLogger,
    project_id: str,
    run_id: str,
    run_name: str | None,
    had_files: bool,
):
    """Log run removal event."""
    logger.log(
        event_type="run_removed",
        data={
            "project_id": project_id,
            "run_id": run_id,
            "run_name": run_name,
            "had_files": had_files,
        },
    )


def log_run_renamed(
    logger: ProvenanceLogger,
    project_id: str,
    run_id: str,
    old_run_name: str | None,
    new_run_name: str,
):
    """Log run rename event."""
    logger.log(
        event_type="run_renamed",
        data={
            "project_id": project_id,
            "run_id": run_id,
            "old_run_name": old_run_name,
            "new_run_name": new_run_name,
        },
    )


def log_ingestion_started(
    logger: ProvenanceLogger,
    project_id: str,
    run_id: str,
    n_files: int,
):
    """Log ingestion start event."""
    logger.log(
        event_type="run_ingestion_started",
        data={
            "project_id": project_id,
            "run_id": run_id,
            "n_files": n_files,
        },
    )


def log_run_ingested(
    logger: ProvenanceLogger,
    project_id: str,
    run_id: str,
    n_cells: int,
    n_markers: int,
    n_files: int,
):
    """Log successful ingestion event."""
    logger.log(
        event_type="run_ingested",
        data={
            "project_id": project_id,
            "run_id": run_id,
            "n_cells": n_cells,
            "n_markers": n_markers,
            "n_files": n_files,
        },
    )


def log_ingestion_failed(
    logger: ProvenanceLogger,
    project_id: str,
    run_id: str,
    error: str,
):
    """Log ingestion failure event."""
    logger.log(
        event_type="run_ingestion_failed",
        data={
            "project_id": project_id,
            "run_id": run_id,
            "error": error,
        },
    )


__all__ = [
    "ProvenanceLogger",
    "get_timestamp",
    "log_project_created",
    "log_project_loaded",
    "log_run_registered",
    "log_run_loaded",
    "log_run_removed",
    "log_run_renamed",
    "log_ingestion_started",
    "log_run_ingested",
    "log_ingestion_failed",
]
