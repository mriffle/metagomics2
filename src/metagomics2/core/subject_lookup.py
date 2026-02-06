"""Runtime lookup of subject annotations from a companion SQLite database."""

import logging
import sqlite3
from pathlib import Path

from metagomics2.core.annotation import SubjectAnnotation
from metagomics2.core.diamond import parse_uniprot_accession

logger = logging.getLogger(__name__)

# SQLite has a default variable limit of 999
_SQLITE_VAR_LIMIT = 999


def load_subject_annotations(
    db_path: Path,
    subject_ids: set[str],
) -> dict[str, SubjectAnnotation]:
    """Look up taxonomy IDs and GO terms for DIAMOND subject hits.

    Extracts bare accessions from full subject IDs (e.g.,
    ``sp|Q21HH2|RS2_SACD2`` -> ``Q21HH2``), queries the companion
    SQLite database, and returns ``SubjectAnnotation`` objects keyed
    by the **full** subject ID so they match ``protein_to_subjects``
    values in the pipeline.

    Args:
        db_path: Path to the companion ``.annotations.db`` SQLite file
        subject_ids: Set of full DIAMOND subject IDs

    Returns:
        Dictionary mapping full subject ID to SubjectAnnotation
    """
    if not subject_ids:
        return {}

    if not db_path.exists():
        logger.warning(f"Annotations database not found: {db_path}")
        return {}

    # Build mapping: bare_accession -> list of full subject IDs
    # (multiple full IDs could map to the same accession in theory)
    accession_to_full: dict[str, list[str]] = {}
    for sid in subject_ids:
        acc = parse_uniprot_accession(sid)
        accession_to_full.setdefault(acc, []).append(sid)

    all_accessions = list(accession_to_full.keys())

    logger.info(
        f"Looking up annotations for {len(all_accessions)} unique accessions "
        f"from {len(subject_ids)} subject IDs"
    )

    conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    try:
        # Fetch taxonomy
        tax_map = _batch_query_taxonomy(conn, all_accessions)

        # Fetch GO annotations
        go_map = _batch_query_go(conn, all_accessions)
    finally:
        conn.close()

    # Build SubjectAnnotation objects keyed by full subject ID
    annotations: dict[str, SubjectAnnotation] = {}
    for accession, full_ids in accession_to_full.items():
        tax_id = tax_map.get(accession)
        go_terms = go_map.get(accession, set())

        for full_id in full_ids:
            annotations[full_id] = SubjectAnnotation(
                subject_id=full_id,
                tax_id=tax_id,
                go_terms=go_terms,
            )

    n_with_tax = sum(1 for a in annotations.values() if a.tax_id is not None)
    n_with_go = sum(1 for a in annotations.values() if a.go_terms)
    logger.info(
        f"Loaded annotations: {len(annotations)} subjects, "
        f"{n_with_tax} with taxonomy, {n_with_go} with GO terms"
    )

    return annotations


def _batch_query_taxonomy(
    conn: sqlite3.Connection,
    accessions: list[str],
) -> dict[str, int]:
    """Query taxonomy table in batches.

    Returns:
        Mapping of accession -> tax_id
    """
    result: dict[str, int] = {}

    for i in range(0, len(accessions), _SQLITE_VAR_LIMIT):
        batch = accessions[i : i + _SQLITE_VAR_LIMIT]
        placeholders = ",".join("?" * len(batch))
        cursor = conn.execute(
            f"SELECT accession, tax_id FROM taxonomy WHERE accession IN ({placeholders})",
            batch,
        )
        for row in cursor:
            result[row[0]] = row[1]

    return result


def _batch_query_go(
    conn: sqlite3.Connection,
    accessions: list[str],
) -> dict[str, set[str]]:
    """Query go_annotations table in batches.

    Returns:
        Mapping of accession -> set of GO term IDs
    """
    result: dict[str, set[str]] = {}

    for i in range(0, len(accessions), _SQLITE_VAR_LIMIT):
        batch = accessions[i : i + _SQLITE_VAR_LIMIT]
        placeholders = ",".join("?" * len(batch))
        cursor = conn.execute(
            f"SELECT accession, go_id FROM go_annotations WHERE accession IN ({placeholders})",
            batch,
        )
        for row in cursor:
            result.setdefault(row[0], set()).add(row[1])

    return result
