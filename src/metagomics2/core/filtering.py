"""Homology hit filtering based on user-configurable policies."""

from dataclasses import dataclass, field
from typing import Literal


@dataclass
class HomologyHit:
    """A single homology search hit."""

    query_id: str
    subject_id: str
    evalue: float
    bitscore: float
    pident: float  # percent identity
    qcov: float  # query coverage
    alnlen: int  # alignment length


@dataclass
class FilterPolicy:
    """Policy for filtering homology hits."""

    max_evalue: float | None = None
    min_pident: float | None = None  # minimum percent identity
    min_qcov: float | None = None  # minimum query coverage
    min_alnlen: int | None = None  # minimum alignment length
    top_k: int | None = None  # keep only top K by bitscore

    def to_dict(self) -> dict:
        """Convert policy to dictionary for manifest."""
        return {
            "max_evalue": self.max_evalue,
            "min_pident": self.min_pident,
            "min_qcov": self.min_qcov,
            "min_alnlen": self.min_alnlen,
            "top_k": self.top_k,
        }


@dataclass
class FilterResult:
    """Result of filtering hits for a query."""

    query_id: str
    accepted_subjects: set[str] = field(default_factory=set)
    total_hits: int = 0
    passed_threshold_hits: int = 0


def passes_thresholds(hit: HomologyHit, policy: FilterPolicy) -> bool:
    """Check if a hit passes all threshold filters.

    Args:
        hit: The homology hit to check
        policy: The filter policy

    Returns:
        True if hit passes all thresholds
    """
    if policy.max_evalue is not None and hit.evalue > policy.max_evalue:
        return False

    if policy.min_pident is not None and hit.pident < policy.min_pident:
        return False

    if policy.min_qcov is not None and hit.qcov < policy.min_qcov:
        return False

    if policy.min_alnlen is not None and hit.alnlen < policy.min_alnlen:
        return False

    return True


def filter_hits_for_query(
    hits: list[HomologyHit],
    policy: FilterPolicy,
) -> FilterResult:
    """Filter hits for a single query protein.

    Filtering steps:
    1. Apply threshold filters (evalue, pident, qcov, alnlen)
    2. Apply tie-aware top_k ranking: keep the top K hits by bitscore,
       plus any additional hits tied with the Kth-best bitscore

    Args:
        hits: List of hits for a single query
        policy: The filter policy

    Returns:
        FilterResult with accepted subject IDs
    """
    if not hits:
        return FilterResult(query_id="", accepted_subjects=set(), total_hits=0)

    query_id = hits[0].query_id
    result = FilterResult(query_id=query_id, total_hits=len(hits))

    # Step 1: Apply threshold filters
    passing_hits = [h for h in hits if passes_thresholds(h, policy)]
    result.passed_threshold_hits = len(passing_hits)

    if not passing_hits:
        return result

    # Step 2: Sort by bitscore descending, then by subject_id for determinism
    passing_hits.sort(key=lambda h: (-h.bitscore, h.subject_id))

    # Step 3: Apply tie-aware top_k ranking filter
    if policy.top_k is not None and len(passing_hits) > policy.top_k:
        # Find the bitscore of the Kth hit (0-indexed: position top_k - 1)
        kth_bitscore = passing_hits[policy.top_k - 1].bitscore
        # Keep all hits whose bitscore >= the Kth-best bitscore
        passing_hits = [h for h in passing_hits if h.bitscore >= kth_bitscore]

    result.accepted_subjects = {h.subject_id for h in passing_hits}
    return result


def filter_all_hits(
    hits_by_query: dict[str, list[HomologyHit]],
    policy: FilterPolicy,
) -> dict[str, set[str]]:
    """Filter hits for all query proteins.

    Args:
        hits_by_query: Dictionary mapping query_id to list of hits
        policy: The filter policy

    Returns:
        Dictionary mapping query_id to set of accepted subject IDs
    """
    result: dict[str, set[str]] = {}

    for query_id, hits in hits_by_query.items():
        filter_result = filter_hits_for_query(hits, policy)
        result[query_id] = filter_result.accepted_subjects

    return result


def filter_all_hits_with_hits(
    hits_by_query: dict[str, list[HomologyHit]],
    policy: FilterPolicy,
) -> dict[str, dict[str, HomologyHit]]:
    """Filter hits for all query proteins, returning accepted HomologyHit objects.

    Like filter_all_hits but preserves the full HomologyHit (including evalue,
    pident, bitscore) for each accepted (query, subject) pair.

    Args:
        hits_by_query: Dictionary mapping query_id to list of hits
        policy: The filter policy

    Returns:
        Dictionary mapping query_id to a dict of {subject_id: HomologyHit}
        for all accepted hits.
    """
    result: dict[str, dict[str, HomologyHit]] = {}

    for query_id, hits in hits_by_query.items():
        filter_result = filter_hits_for_query(hits, policy)
        if filter_result.accepted_subjects:
            subject_to_hit: dict[str, HomologyHit] = {}
            for h in hits:
                if h.subject_id in filter_result.accepted_subjects:
                    subject_to_hit[h.subject_id] = h
            result[query_id] = subject_to_hit

    return result


def parse_blast_tabular(
    lines: list[str],
    columns: list[str] | None = None,
) -> dict[str, list[HomologyHit]]:
    """Parse BLAST/DIAMOND tabular output.

    Default expected columns (outfmt 6 style):
    qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

    Args:
        lines: Lines from the tabular output file
        columns: Column names if non-standard format

    Returns:
        Dictionary mapping query_id to list of HomologyHit objects
    """
    if columns is None:
        columns = [
            "qseqid", "sseqid", "pident", "length", "mismatch",
            "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"
        ]

    # Find column indices
    try:
        qseqid_idx = columns.index("qseqid")
        sseqid_idx = columns.index("sseqid")
        pident_idx = columns.index("pident")
        length_idx = columns.index("length")
        evalue_idx = columns.index("evalue")
        bitscore_idx = columns.index("bitscore")
    except ValueError as e:
        raise ValueError(f"Missing required column: {e}")

    # qcov might not be in standard output, compute if qlen available
    qcov_idx = columns.index("qcov") if "qcov" in columns else None
    qlen_idx = columns.index("qlen") if "qlen" in columns else None

    hits_by_query: dict[str, list[HomologyHit]] = {}

    for line in lines:
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        parts = line.split("\t")

        query_id = parts[qseqid_idx]
        subject_id = parts[sseqid_idx]
        pident = float(parts[pident_idx])
        alnlen = int(parts[length_idx])
        evalue = float(parts[evalue_idx])
        bitscore = float(parts[bitscore_idx])

        # Compute query coverage if possible
        if qcov_idx is not None:
            qcov = float(parts[qcov_idx])
        elif qlen_idx is not None:
            qlen = int(parts[qlen_idx])
            qcov = (alnlen / qlen) * 100 if qlen > 0 else 0.0
        else:
            qcov = 0.0  # Unknown, will pass any filter if min_qcov is None

        hit = HomologyHit(
            query_id=query_id,
            subject_id=subject_id,
            evalue=evalue,
            bitscore=bitscore,
            pident=pident,
            qcov=qcov,
            alnlen=alnlen,
        )

        if query_id not in hits_by_query:
            hits_by_query[query_id] = []
        hits_by_query[query_id].append(hit)

    return hits_by_query
