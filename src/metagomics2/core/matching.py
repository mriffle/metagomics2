"""Exact peptide to protein matching using Aho-Corasick algorithm."""

from dataclasses import dataclass, field

import ahocorasick


@dataclass
class MatchResult:
    """Result of peptide-to-protein matching."""

    peptide_to_proteins: dict[str, set[str]] = field(default_factory=dict)
    matched_peptides: set[str] = field(default_factory=set)
    unmatched_peptides: set[str] = field(default_factory=set)
    hit_proteins: set[str] = field(default_factory=set)

    @property
    def n_matched(self) -> int:
        """Number of peptides that matched at least one protein."""
        return len(self.matched_peptides)

    @property
    def n_unmatched(self) -> int:
        """Number of peptides that matched no proteins."""
        return len(self.unmatched_peptides)


def build_automaton(peptides: set[str]) -> ahocorasick.Automaton:
    """Build an Aho-Corasick automaton from a set of peptide sequences.

    Args:
        peptides: Set of peptide sequences (should be normalized/uppercase)

    Returns:
        Configured Aho-Corasick automaton
    """
    automaton = ahocorasick.Automaton()

    for peptide in peptides:
        automaton.add_word(peptide, peptide)

    automaton.make_automaton()
    return automaton


def match_peptides_to_proteins(
    peptides: set[str],
    proteins: dict[str, str],
) -> MatchResult:
    """Match peptides to proteins using exact substring matching.

    Uses Aho-Corasick algorithm for efficient multi-pattern matching.
    Each peptide maps to a set of protein IDs (no multiplicity).

    Args:
        peptides: Set of peptide sequences to search for
        proteins: Dictionary mapping protein ID to sequence

    Returns:
        MatchResult containing the mapping and statistics
    """
    if not peptides:
        return MatchResult(
            peptide_to_proteins={},
            matched_peptides=set(),
            unmatched_peptides=set(),
            hit_proteins=set(),
        )

    # Initialize result with all peptides as unmatched
    result = MatchResult(
        peptide_to_proteins={peptide: set() for peptide in peptides},
        matched_peptides=set(),
        unmatched_peptides=set(peptides),
        hit_proteins=set(),
    )

    if not proteins:
        return result

    # Build automaton
    automaton = build_automaton(peptides)

    # Search each protein sequence
    for protein_id, sequence in proteins.items():
        # Ensure sequence is uppercase for matching
        sequence_upper = sequence.upper()

        for _, peptide in automaton.iter(sequence_upper):
            result.peptide_to_proteins[peptide].add(protein_id)
            result.hit_proteins.add(protein_id)

    # Update matched/unmatched sets
    for peptide, protein_ids in result.peptide_to_proteins.items():
        if protein_ids:
            result.matched_peptides.add(peptide)
            result.unmatched_peptides.discard(peptide)

    return result


def get_union_hit_proteins(match_result: MatchResult) -> set[str]:
    """Get the union of all proteins hit by any peptide.

    Args:
        match_result: Result from match_peptides_to_proteins

    Returns:
        Set of all protein IDs that were hit
    """
    return match_result.hit_proteins
