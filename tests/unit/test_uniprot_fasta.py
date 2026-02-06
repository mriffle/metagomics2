"""Tests for UniProt FASTA OX= annotation parser."""

import gzip
import io

from metagomics2.core.uniprot_fasta import (
    parse_uniprot_fasta_annotations,
    parse_uniprot_fasta_annotations_stream,
)


_FASTA_CONTENT = """\
>sp|Q21HH2|RS2_SACD2 30S ribosomal protein S2 OS=Saccharopolyspora erythraea OX=266940 GN=rpsB PE=3 SV=1
MKTLQIRNPQRNAARRRSRPFQFERGQTLVHISGEPVTLKECNLVGSTLNPRGVNALTK
>sp|P12345|MYP_HUMAN Myosin OS=Homo sapiens OX=9606 GN=MYH7 PE=1 SV=4
MSSDSEMAIFGEAAPYLRKSEKERIEAQNKPFDAKTSVFVAEPDEEVGALVKRQGVMYLFK
>sp|A0B1C2|NOOXTAG No OX tag here GN=FOO PE=1 SV=1
ACDEFGHIKLMNPQRSTVWY
"""


def test_parse_stream_extracts_accession_and_taxid():
    records = list(parse_uniprot_fasta_annotations_stream(io.StringIO(_FASTA_CONTENT)))

    assert len(records) == 2
    assert records[0] == ("Q21HH2", 266940)
    assert records[1] == ("P12345", 9606)


def test_parse_stream_skips_missing_ox():
    """Entries without OX= are skipped."""
    records = list(parse_uniprot_fasta_annotations_stream(io.StringIO(_FASTA_CONTENT)))
    accessions = [r[0] for r in records]
    assert "A0B1C2" not in accessions


def test_parse_file_plain(tmp_path):
    fasta_path = tmp_path / "test.fasta"
    fasta_path.write_text(_FASTA_CONTENT)

    records = list(parse_uniprot_fasta_annotations(fasta_path))
    assert len(records) == 2
    assert records[0] == ("Q21HH2", 266940)


def test_parse_file_gzipped(tmp_path):
    fasta_path = tmp_path / "test.fasta.gz"
    with gzip.open(fasta_path, "wt", encoding="utf-8") as f:
        f.write(_FASTA_CONTENT)

    records = list(parse_uniprot_fasta_annotations(fasta_path))
    assert len(records) == 2
    assert records[1] == ("P12345", 9606)


def test_parse_stream_trembl_format():
    content = (
        ">tr|A0A0A0MQG0|A0A0A0MQG0_HUMAN Some protein OS=Homo sapiens OX=9606 PE=4 SV=1\n"
        "ACDEFGHIKLMNPQRSTVWY\n"
    )
    records = list(parse_uniprot_fasta_annotations_stream(io.StringIO(content)))
    assert len(records) == 1
    assert records[0] == ("A0A0A0MQG0", 9606)
