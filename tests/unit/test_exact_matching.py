"""Unit tests for exact peptide to protein matching."""

from pathlib import Path

import pytest

from metagomics2.core.fasta import parse_fasta, build_protein_dict
from metagomics2.core.matching import (
    MatchResult,
    build_automaton,
    get_union_hit_proteins,
    match_peptides_to_proteins,
)


class TestBuildAutomaton:
    """Tests for Aho-Corasick automaton building."""

    def test_builds_automaton(self):
        peptides = {"PEPTIDE", "ABC"}
        automaton = build_automaton(peptides)

        # Automaton should be searchable
        assert automaton is not None

    def test_empty_peptides(self):
        automaton = build_automaton(set())
        # Should not raise
        assert automaton is not None


class TestMatchPeptidesToProteins:
    """Tests for peptide-to-protein matching."""

    def test_correct_matches_fixture(self, small_background_fasta: Path):
        """Test with the standard fixture data."""
        records = parse_fasta(small_background_fasta)
        proteins = build_protein_dict(records)

        peptides = {"PEPTIDE", "ABC", "NOMATCH"}
        result = match_peptides_to_proteins(peptides, proteins)

        # PEPTIDE should match B1 and B2
        assert result.peptide_to_proteins["PEPTIDE"] == {"B1", "B2"}

        # ABC should match B3
        assert result.peptide_to_proteins["ABC"] == {"B3"}

        # NOMATCH should match nothing
        assert result.peptide_to_proteins["NOMATCH"] == set()

    def test_matched_unmatched_sets(self, small_background_fasta: Path):
        """Test matched and unmatched peptide tracking."""
        records = parse_fasta(small_background_fasta)
        proteins = build_protein_dict(records)

        peptides = {"PEPTIDE", "ABC", "NOMATCH"}
        result = match_peptides_to_proteins(peptides, proteins)

        assert result.matched_peptides == {"PEPTIDE", "ABC"}
        assert result.unmatched_peptides == {"NOMATCH"}
        assert result.n_matched == 2
        assert result.n_unmatched == 1

    def test_hit_proteins(self, small_background_fasta: Path):
        """Test hit protein tracking."""
        records = parse_fasta(small_background_fasta)
        proteins = build_protein_dict(records)

        peptides = {"PEPTIDE", "ABC", "NOMATCH"}
        result = match_peptides_to_proteins(peptides, proteins)

        assert result.hit_proteins == {"B1", "B2", "B3"}

    def test_multiple_occurrences_in_one_protein(self):
        """Peptide appearing multiple times in one protein maps once."""
        proteins = {"P1": "PEPTIDEPEPTIDE"}
        peptides = {"PEPTIDE"}

        result = match_peptides_to_proteins(peptides, proteins)

        # Should map to P1 once (set, not list)
        assert result.peptide_to_proteins["PEPTIDE"] == {"P1"}

    def test_overlapping_peptides(self):
        """Multiple peptides can match overlapping regions."""
        proteins = {"P1": "MPEPTIDEKAAA"}
        peptides = {"PEP", "PEPTIDE", "TIDE"}

        result = match_peptides_to_proteins(peptides, proteins)

        assert result.peptide_to_proteins["PEP"] == {"P1"}
        assert result.peptide_to_proteins["PEPTIDE"] == {"P1"}
        assert result.peptide_to_proteins["TIDE"] == {"P1"}

    def test_case_insensitive_protein_sequence(self):
        """Matching should work with lowercase protein sequences."""
        proteins = {"P1": "mpeptidekaaa"}  # lowercase
        peptides = {"PEPTIDE"}  # uppercase (normalized)

        result = match_peptides_to_proteins(peptides, proteins)

        assert result.peptide_to_proteins["PEPTIDE"] == {"P1"}

    def test_empty_peptides(self):
        """Empty peptide set returns empty result."""
        proteins = {"P1": "MPEPTIDEKAAA"}

        result = match_peptides_to_proteins(set(), proteins)

        assert result.peptide_to_proteins == {}
        assert result.matched_peptides == set()
        assert result.unmatched_peptides == set()

    def test_empty_proteins(self):
        """Empty protein dict returns all peptides unmatched."""
        peptides = {"PEPTIDE", "ABC"}

        result = match_peptides_to_proteins(peptides, {})

        assert result.peptide_to_proteins["PEPTIDE"] == set()
        assert result.peptide_to_proteins["ABC"] == set()
        assert result.unmatched_peptides == {"PEPTIDE", "ABC"}
        assert result.matched_peptides == set()

    def test_peptide_longer_than_protein(self):
        """Peptide longer than protein sequence doesn't match."""
        proteins = {"P1": "ABC"}
        peptides = {"ABCDEFGH"}

        result = match_peptides_to_proteins(peptides, proteins)

        assert result.peptide_to_proteins["ABCDEFGH"] == set()

    def test_exact_match_full_protein(self):
        """Peptide that exactly matches full protein sequence."""
        proteins = {"P1": "PEPTIDE"}
        peptides = {"PEPTIDE"}

        result = match_peptides_to_proteins(peptides, proteins)

        assert result.peptide_to_proteins["PEPTIDE"] == {"P1"}

    def test_peptide_at_start(self):
        """Peptide at the start of protein."""
        proteins = {"P1": "PEPTIDEXXXX"}
        peptides = {"PEPTIDE"}

        result = match_peptides_to_proteins(peptides, proteins)

        assert result.peptide_to_proteins["PEPTIDE"] == {"P1"}

    def test_peptide_at_end(self):
        """Peptide at the end of protein."""
        proteins = {"P1": "XXXXPEPTIDE"}
        peptides = {"PEPTIDE"}

        result = match_peptides_to_proteins(peptides, proteins)

        assert result.peptide_to_proteins["PEPTIDE"] == {"P1"}


class TestGetUnionHitProteins:
    """Tests for getting union of hit proteins."""

    def test_union_hit_proteins(self, small_background_fasta: Path):
        records = parse_fasta(small_background_fasta)
        proteins = build_protein_dict(records)

        peptides = {"PEPTIDE", "ABC"}
        result = match_peptides_to_proteins(peptides, proteins)

        union = get_union_hit_proteins(result)
        assert union == {"B1", "B2", "B3"}

    def test_union_empty_result(self):
        result = MatchResult()
        union = get_union_hit_proteins(result)
        assert union == set()


class TestMatchResultDataclass:
    """Tests for MatchResult dataclass."""

    def test_default_values(self):
        result = MatchResult()

        assert result.peptide_to_proteins == {}
        assert result.matched_peptides == set()
        assert result.unmatched_peptides == set()
        assert result.hit_proteins == set()
        assert result.n_matched == 0
        assert result.n_unmatched == 0
