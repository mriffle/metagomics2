"""Tests for the annotations database build script."""

import sqlite3
from pathlib import Path

from scripts.build_annotations_db import build_annotations_db


_FASTA_CONTENT = """\
>sp|Q21HH2|RS2_SACD2 30S ribosomal protein S2 OS=Saccharopolyspora erythraea OX=266940 GN=rpsB PE=3 SV=1
MKTLQIRNPQRNAARRRSRPFQFERGQTLVHISGEPVTLKECNLVGSTLNPRGVNALTK
>sp|P12345|MYP_HUMAN Myosin OS=Homo sapiens OX=9606 GN=MYH7 PE=1 SV=4
MSSDSEMAIFGEAAPYLRKSEKERIEAQNKPFDAKTSVFVAEPDEEVGALVKRQGVMYLFK
"""

_GAF_CONTENT = """\
!gaf-version: 2.2
!generated-by: test
UniProtKB\tQ21HH2\tRS2\tenables\tGO:0003735\tPMID:12345\tIDA\t\tF\t30S ribosomal protein S2\t\tprotein\ttaxon:266940\t20200101\tUniProt\t\t
UniProtKB\tQ21HH2\tRS2\tlocated_in\tGO:0005840\tPMID:12345\tIEA\t\tC\t30S ribosomal protein S2\t\tprotein\ttaxon:266940\t20200101\tUniProt\t\t
UniProtKB\tP12345\tMYP\tinvolved_in\tGO:0006412\tPMID:11111\tIEA\t\tP\tMyosin\t\tprotein\ttaxon:9606\t20200101\tUniProt\t\t
UniProtKB\tP12345\tMYP\tNOT|enables\tGO:0009999\tPMID:22222\tIEA\t\tF\tMyosin\t\tprotein\ttaxon:9606\t20200101\tInterPro\t\t
UniProtKB\tP12345\tMYP\tenables\tGO:0005575\tGO_REF:0000015\tND\t\tC\tMyosin\t\tprotein\ttaxon:9606\t20200101\tUniProt\t\t
"""


def test_build_annotations_db(tmp_path):
    fasta_path = tmp_path / "test.fasta"
    fasta_path.write_text(_FASTA_CONTENT)

    gaf_path = tmp_path / "test.gaf"
    gaf_path.write_text(_GAF_CONTENT)

    output_path = tmp_path / "test.annotations.db"

    build_annotations_db(fasta_path, gaf_path, output_path)

    assert output_path.exists()

    conn = sqlite3.connect(str(output_path))

    # Check taxonomy
    rows = conn.execute("SELECT accession, tax_id FROM taxonomy ORDER BY accession").fetchall()
    assert len(rows) == 2
    assert rows[0] == ("P12345", 9606)
    assert rows[1] == ("Q21HH2", 266940)

    # Check GO annotations (NOT qualifier and ND evidence lines should be filtered out)
    rows = conn.execute(
        "SELECT accession, go_id, aspect FROM go_annotations ORDER BY accession, go_id"
    ).fetchall()
    assert len(rows) == 3
    assert rows[0] == ("P12345", "GO:0006412", "P")
    assert rows[1] == ("Q21HH2", "GO:0003735", "F")
    assert rows[2] == ("Q21HH2", "GO:0005840", "C")

    # Check metadata
    meta = dict(conn.execute("SELECT key, value FROM metadata").fetchall())
    assert meta["taxonomy_count"] == "2"
    assert meta["go_annotation_count"] == "3"
    assert meta["gaf_filters"] == "exclude NOT qualifier, exclude ND evidence"

    conn.close()


def test_build_annotations_db_filters_by_fasta_accessions(tmp_path):
    """GAF rows for accessions not in the FASTA should be excluded."""
    fasta_path = tmp_path / "test.fasta"
    # Only Q21HH2 in the FASTA — P12345 is absent
    fasta_path.write_text(
        ">sp|Q21HH2|RS2_SACD2 OS=Saccharopolyspora erythraea OX=266940\n"
        "MKTLQIRNPQRNAARRRSRPFQFERGQTLVHISGEPVTLKECNLVGSTLNPRGVNALTK\n"
    )

    gaf_path = tmp_path / "test.gaf"
    gaf_path.write_text(_GAF_CONTENT)

    output_path = tmp_path / "test.annotations.db"
    build_annotations_db(fasta_path, gaf_path, output_path)

    conn = sqlite3.connect(str(output_path))

    # Only Q21HH2 GO annotations should be present
    rows = conn.execute(
        "SELECT accession, go_id FROM go_annotations ORDER BY accession, go_id"
    ).fetchall()
    assert len(rows) == 2
    assert all(r[0] == "Q21HH2" for r in rows)

    # P12345 rows should have been skipped
    p_rows = conn.execute(
        "SELECT COUNT(*) FROM go_annotations WHERE accession = 'P12345'"
    ).fetchone()
    assert p_rows[0] == 0

    conn.close()
