"""Tests for runtime subject annotation lookup from companion SQLite database."""

import sqlite3

from metagomics2.core.subject_lookup import load_subject_annotations


def _create_test_db(db_path):
    """Create a small test annotations database."""
    conn = sqlite3.connect(str(db_path))
    conn.executescript("""
        CREATE TABLE taxonomy (
            accession TEXT PRIMARY KEY,
            tax_id INTEGER NOT NULL
        );
        CREATE TABLE go_annotations (
            accession TEXT NOT NULL,
            go_id TEXT NOT NULL,
            aspect TEXT NOT NULL,
            evidence_code TEXT NOT NULL
        );
        CREATE INDEX idx_go_accession ON go_annotations(accession);
    """)
    conn.executemany(
        "INSERT INTO taxonomy (accession, tax_id) VALUES (?, ?)",
        [("Q21HH2", 266940), ("P12345", 9606), ("A0B1C2", 7227)],
    )
    conn.executemany(
        "INSERT INTO go_annotations (accession, go_id, aspect, evidence_code) VALUES (?, ?, ?, ?)",
        [
            ("Q21HH2", "GO:0003735", "F", "IDA"),
            ("Q21HH2", "GO:0005840", "C", "IEA"),
            ("P12345", "GO:0006412", "P", "IEA"),
        ],
    )
    conn.commit()
    conn.close()


def test_load_annotations_basic(tmp_path):
    db_path = tmp_path / "test.annotations.db"
    _create_test_db(db_path)

    subject_ids = {"sp|Q21HH2|RS2_SACD2", "sp|P12345|MYP_HUMAN"}
    result = load_subject_annotations(db_path, subject_ids)

    assert len(result) == 2

    ann_q = result["sp|Q21HH2|RS2_SACD2"]
    assert ann_q.subject_id == "sp|Q21HH2|RS2_SACD2"
    assert ann_q.tax_id == 266940
    assert ann_q.go_terms == {"GO:0003735", "GO:0005840"}

    ann_p = result["sp|P12345|MYP_HUMAN"]
    assert ann_p.subject_id == "sp|P12345|MYP_HUMAN"
    assert ann_p.tax_id == 9606
    assert ann_p.go_terms == {"GO:0006412"}


def test_load_annotations_missing_accession(tmp_path):
    """Subject IDs not in the DB still get an entry with None tax_id and empty GO."""
    db_path = tmp_path / "test.annotations.db"
    _create_test_db(db_path)

    subject_ids = {"sp|XXXXXX|UNKNOWN"}
    result = load_subject_annotations(db_path, subject_ids)

    assert len(result) == 1
    ann = result["sp|XXXXXX|UNKNOWN"]
    assert ann.tax_id is None
    assert ann.go_terms == set()


def test_load_annotations_empty_set(tmp_path):
    db_path = tmp_path / "test.annotations.db"
    _create_test_db(db_path)

    result = load_subject_annotations(db_path, set())
    assert result == {}


def test_load_annotations_missing_db(tmp_path):
    db_path = tmp_path / "nonexistent.db"
    result = load_subject_annotations(db_path, {"sp|Q21HH2|RS2_SACD2"})
    assert result == {}


def test_load_annotations_bare_accession(tmp_path):
    """Bare accessions (no db|...|... format) should also work."""
    db_path = tmp_path / "test.annotations.db"
    _create_test_db(db_path)

    subject_ids = {"A0B1C2"}
    result = load_subject_annotations(db_path, subject_ids)

    assert len(result) == 1
    ann = result["A0B1C2"]
    assert ann.tax_id == 7227
    assert ann.go_terms == set()


def test_load_annotations_large_batch(tmp_path):
    """Test with more than 999 subject IDs to verify batching."""
    db_path = tmp_path / "test.annotations.db"
    conn = sqlite3.connect(str(db_path))
    conn.executescript("""
        CREATE TABLE taxonomy (
            accession TEXT PRIMARY KEY,
            tax_id INTEGER NOT NULL
        );
        CREATE TABLE go_annotations (
            accession TEXT NOT NULL,
            go_id TEXT NOT NULL,
            aspect TEXT NOT NULL,
            evidence_code TEXT NOT NULL
        );
        CREATE INDEX idx_go_accession ON go_annotations(accession);
    """)
    # Insert 1500 entries
    entries = [(f"ACC{i:05d}", i) for i in range(1500)]
    conn.executemany("INSERT INTO taxonomy (accession, tax_id) VALUES (?, ?)", entries)
    conn.commit()
    conn.close()

    subject_ids = {f"ACC{i:05d}" for i in range(1500)}
    result = load_subject_annotations(db_path, subject_ids)

    assert len(result) == 1500
    assert result["ACC00000"].tax_id == 0
    assert result["ACC01499"].tax_id == 1499
