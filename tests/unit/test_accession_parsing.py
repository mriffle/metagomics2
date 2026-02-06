"""Tests for UniProt accession parsing from DIAMOND subject IDs."""

from metagomics2.core.diamond import parse_uniprot_accession


def test_swissprot_format():
    assert parse_uniprot_accession("sp|Q21HH2|RS2_SACD2") == "Q21HH2"


def test_trembl_format():
    assert parse_uniprot_accession("tr|A0A0A0MQG0|A0A0A0MQG0_HUMAN") == "A0A0A0MQG0"


def test_bare_accession():
    assert parse_uniprot_accession("P12345") == "P12345"


def test_bare_accession_with_underscore():
    assert parse_uniprot_accession("A0A0A0MQG0") == "A0A0A0MQG0"


def test_swissprot_short_accession():
    assert parse_uniprot_accession("sp|P12345|MYP_HUMAN") == "P12345"


def test_accession_with_hyphen():
    assert parse_uniprot_accession("sp|A0A-123|ENTRY_ORG") == "A0A-123"
