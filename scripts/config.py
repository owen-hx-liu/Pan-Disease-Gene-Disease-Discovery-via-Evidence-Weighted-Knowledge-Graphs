"""Central configuration: repo-relative paths and secrets-from-environment.

This module removes hard-coded absolute Windows paths and hard-coded database
credentials from the pipeline scripts (a Tier-1 release blocker; see CLAUDE.md
and PROJECT_REPORT.md). Import it instead of pasting machine-specific paths.

Paths resolve relative to the repository root (the parent of this ``scripts/``
directory), so the code runs unchanged on any machine or OS.

Neo4j credentials are read from environment variables; nothing secret is stored
in the source tree:

    NEO4J_URI       (default bolt://localhost:7687)
    NEO4J_USER      (default neo4j)
    NEO4J_PASSWORD  (no default -- must be set to use the Neo4j scripts)
"""
from __future__ import annotations

import os
from pathlib import Path

# Repo root = parent of this file's directory (scripts/ -> repo).
REPO_ROOT = Path(__file__).resolve().parent.parent

DATA_DIR = REPO_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"

# Canonical analysis inputs.
RELATIONSHIPS_CSV = PROCESSED_DIR / "relationships.csv"
EDGES_CLEAN_CSV = PROCESSED_DIR / "edges_clean.csv"
NODES_CSV = PROCESSED_DIR / "nodes.csv"
CANONICAL_NODES_CSV = PROCESSED_DIR / "canonical_nodes.csv"


def data_path(*parts: str) -> str:
    """Return a path under data/ as a string, e.g. data_path('processed', 'x.csv')."""
    return str(DATA_DIR.joinpath(*parts))


def neo4j_credentials() -> tuple[str, tuple[str, str]]:
    """Return (uri, (user, password)) from the environment.

    Raises RuntimeError if NEO4J_PASSWORD is unset so we fail loudly instead of
    silently shipping a default password.
    """
    uri = os.environ.get("NEO4J_URI", "bolt://localhost:7687")
    user = os.environ.get("NEO4J_USER", "neo4j")
    password = os.environ.get("NEO4J_PASSWORD")
    if not password:
        raise RuntimeError(
            "NEO4J_PASSWORD is not set. Export your Neo4j password before running "
            "the Neo4j scripts, e.g.  setx NEO4J_PASSWORD <password>  (Windows) or "
            "export NEO4J_PASSWORD=<password>  (Linux/macOS)."
        )
    return uri, (user, password)
