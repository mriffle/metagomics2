"""Tests for GAF 2.2 parser."""

import io
import gzip
import tempfile
from pathlib import Path

from metagomics2.core.gaf_parser import GOARecord, parse_gaf_file, parse_gaf_stream


# Sample GAF lines (tab-delimited, 17 columns)
_HEADER = "!gaf-version: 2.2\n!generated-by: UniProt\n"

_LINE_NORMAL = (
    "UniProtKB\tQ21HH2\tRS2\tenables\tGO:0003735\tPMID:12345\tIDA\t\tF\t"
    "30S ribosomal protein S2\t\tprotein\ttaxon:266940\t20200101\tUniProt\t\t\n"
)

_LINE_NORMAL_2 = (
    "UniProtKB\tP12345\tMYP\tinvolved_in\tGO:0006412\tPMID:11111\tIEA\t\tP\t"
    "Myosin\t\tprotein\ttaxon:9606\t20200101\tUniProt\t\t\n"
)

_LINE_NOT_QUALIFIER = (
    "UniProtKB\tQ21HH2\tRS2\tNOT|enables\tGO:0005840\tPMID:99999\tIEA\t\tC\t"
    "30S ribosomal protein S2\t\tprotein\ttaxon:266940\t20200101\tInterPro\t\t\n"
)

_LINE_ND_EVIDENCE = (
    "UniProtKB\tP99999\tXYZ\tenables\tGO:0005575\tGO_REF:0000015\tND\t\tC\t"
    "Some protein\t\tprotein\ttaxon:9606\t20200101\tUniProt\t\t\n"
)


def test_parse_gaf_stream_default_filters():
    text = _HEADER + _LINE_NORMAL + _LINE_NOT_QUALIFIER + _LINE_ND_EVIDENCE + _LINE_NORMAL_2
    records = list(parse_gaf_stream(io.StringIO(text)))

    assert len(records) == 2
    assert records[0] == GOARecord(
        accession="Q21HH2", go_id="GO:0003735", aspect="F", evidence_code="IDA"
    )
    assert records[1] == GOARecord(
        accession="P12345", go_id="GO:0006412", aspect="P", evidence_code="IEA"
    )


def test_parse_gaf_stream_excludes_not_qualifier():
    text = _HEADER + _LINE_NORMAL + _LINE_NOT_QUALIFIER
    records = list(parse_gaf_stream(io.StringIO(text), exclude_not_qualifier=True))

    assert len(records) == 1
    assert records[0].accession == "Q21HH2"


def test_parse_gaf_stream_includes_not_qualifier_when_disabled():
    text = _HEADER + _LINE_NORMAL + _LINE_NOT_QUALIFIER
    records = list(parse_gaf_stream(io.StringIO(text), exclude_not_qualifier=False))

    assert len(records) == 2


def test_parse_gaf_stream_excludes_nd_evidence():
    text = _HEADER + _LINE_NORMAL + _LINE_ND_EVIDENCE
    records = list(parse_gaf_stream(io.StringIO(text), exclude_nd_evidence=True))

    assert len(records) == 1
    assert records[0].accession == "Q21HH2"


def test_parse_gaf_stream_includes_nd_evidence_when_disabled():
    text = _HEADER + _LINE_NORMAL + _LINE_ND_EVIDENCE
    records = list(parse_gaf_stream(io.StringIO(text), exclude_nd_evidence=False))

    assert len(records) == 2


def test_parse_gaf_stream_no_filters():
    text = _HEADER + _LINE_NORMAL + _LINE_NOT_QUALIFIER + _LINE_ND_EVIDENCE
    records = list(parse_gaf_stream(
        io.StringIO(text), exclude_not_qualifier=False, exclude_nd_evidence=False
    ))

    assert len(records) == 3


def test_parse_gaf_stream_skips_comments():
    text = "! this is a comment\n" + _LINE_NORMAL
    records = list(parse_gaf_stream(io.StringIO(text)))

    assert len(records) == 1


def test_parse_gaf_stream_skips_short_lines():
    text = "too\tfew\tcolumns\n" + _LINE_NORMAL
    records = list(parse_gaf_stream(io.StringIO(text)))

    assert len(records) == 1


def test_parse_gaf_file_plain(tmp_path):
    gaf_path = tmp_path / "test.gaf"
    gaf_path.write_text(_HEADER + _LINE_NORMAL + _LINE_NORMAL_2)

    records = list(parse_gaf_file(gaf_path))
    assert len(records) == 2


def test_parse_gaf_file_gzipped(tmp_path):
    gaf_path = tmp_path / "test.gaf.gz"
    content = (_HEADER + _LINE_NORMAL + _LINE_NORMAL_2).encode("utf-8")
    with gzip.open(gaf_path, "wb") as f:
        f.write(content)

    records = list(parse_gaf_file(gaf_path))
    assert len(records) == 2
    assert records[0].accession == "Q21HH2"
    assert records[1].accession == "P12345"
