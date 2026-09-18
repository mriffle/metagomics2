"""Unit tests for homology hit filtering."""

import pytest

from metagomics2.core.filtering import (
    FilterPolicy,
    HomologyHit,
    filter_all_hits,
    filter_all_hits_with_hits,
    filter_hits_for_query,
    parse_blast_tabular,
    passes_thresholds,
)


def make_hit(
    query_id: str = "Q1",
    subject_id: str = "S1",
    evalue: float = 1e-10,
    bitscore: float = 100.0,
    pident: float = 90.0,
    qcov: float = 95.0,
    alnlen: int = 100,
) -> HomologyHit:
    """Helper to create a HomologyHit."""
    return HomologyHit(
        query_id=query_id,
        subject_id=subject_id,
        evalue=evalue,
        bitscore=bitscore,
        pident=pident,
        qcov=qcov,
        alnlen=alnlen,
    )


class TestPassesThresholds:
    """Tests for threshold filtering."""

    def test_passes_all_thresholds(self):
        hit = make_hit(evalue=1e-10, pident=90.0, qcov=95.0, alnlen=100)
        policy = FilterPolicy(
            max_evalue=1e-5,
            min_pident=80.0,
            min_qcov=90.0,
            min_alnlen=50,
        )
        assert passes_thresholds(hit, policy) is True

    def test_fails_evalue(self):
        hit = make_hit(evalue=1e-3)  # Too high
        policy = FilterPolicy(max_evalue=1e-5)
        assert passes_thresholds(hit, policy) is False

    def test_fails_pident(self):
        hit = make_hit(pident=70.0)  # Too low
        policy = FilterPolicy(min_pident=80.0)
        assert passes_thresholds(hit, policy) is False

    def test_fails_qcov(self):
        hit = make_hit(qcov=50.0)  # Too low
        policy = FilterPolicy(min_qcov=80.0)
        assert passes_thresholds(hit, policy) is False

    def test_fails_alnlen(self):
        hit = make_hit(alnlen=30)  # Too short
        policy = FilterPolicy(min_alnlen=50)
        assert passes_thresholds(hit, policy) is False

    def test_no_thresholds_always_passes(self):
        hit = make_hit()
        policy = FilterPolicy()
        assert passes_thresholds(hit, policy) is True

    def test_boundary_evalue_equal(self):
        hit = make_hit(evalue=1e-5)
        policy = FilterPolicy(max_evalue=1e-5)
        assert passes_thresholds(hit, policy) is True

    def test_boundary_pident_equal(self):
        hit = make_hit(pident=80.0)
        policy = FilterPolicy(min_pident=80.0)
        assert passes_thresholds(hit, policy) is True


class TestFilterHitsForQuery:
    """Tests for filtering hits for a single query."""

    def test_threshold_filters(self):
        hits = [
            make_hit(subject_id="S1", evalue=1e-10, bitscore=100),
            make_hit(subject_id="S2", evalue=1e-3, bitscore=50),  # Fails evalue
            make_hit(subject_id="S3", evalue=1e-8, bitscore=80),
        ]
        policy = FilterPolicy(max_evalue=1e-5)

        result = filter_hits_for_query(hits, policy)

        assert result.accepted_subjects == {"S1", "S3"}
        assert result.total_hits == 3
        assert result.passed_threshold_hits == 2

    def test_top_k_by_bitscore_no_ties(self):
        hits = [
            make_hit(subject_id="S1", bitscore=100),
            make_hit(subject_id="S2", bitscore=90),
            make_hit(subject_id="S3", bitscore=80),
            make_hit(subject_id="S4", bitscore=70),
        ]
        policy = FilterPolicy(top_k=2)

        result = filter_hits_for_query(hits, policy)

        assert result.accepted_subjects == {"S1", "S2"}

    def test_top_k_keeps_ties_at_boundary(self):
        hits = [
            make_hit(subject_id="S1", bitscore=100),
            make_hit(subject_id="S2", bitscore=90),
            make_hit(subject_id="S3", bitscore=90),  # Tied with S2 at boundary
            make_hit(subject_id="S4", bitscore=70),
        ]
        policy = FilterPolicy(top_k=2)

        result = filter_hits_for_query(hits, policy)

        # S2 and S3 are tied at the Kth-best bitscore, both kept
        assert result.accepted_subjects == {"S1", "S2", "S3"}

    def test_top_k_1_keeps_all_tied_best_hits(self):
        hits = [
            make_hit(subject_id="S1", bitscore=100),
            make_hit(subject_id="S2", bitscore=100),
            make_hit(subject_id="S3", bitscore=100),
            make_hit(subject_id="S4", bitscore=90),
        ]
        policy = FilterPolicy(top_k=1)

        result = filter_hits_for_query(hits, policy)

        # All three tied at the best bitscore are kept
        assert result.accepted_subjects == {"S1", "S2", "S3"}

    def test_top_k_ties_within_k_do_not_expand(self):
        hits = [
            make_hit(subject_id="S1", bitscore=100),
            make_hit(subject_id="S2", bitscore=100),  # Tied within top_k
            make_hit(subject_id="S3", bitscore=80),
        ]
        policy = FilterPolicy(top_k=2)

        result = filter_hits_for_query(hits, policy)

        # S1 and S2 are within top_k, S3 is below the Kth bitscore
        assert result.accepted_subjects == {"S1", "S2"}

    def test_top_k_fewer_hits_than_k(self):
        hits = [
            make_hit(subject_id="S1", bitscore=100),
        ]
        policy = FilterPolicy(top_k=5)

        result = filter_hits_for_query(hits, policy)

        assert result.accepted_subjects == {"S1"}

    def test_top_k_all_tied(self):
        hits = [
            make_hit(subject_id="S1", bitscore=100),
            make_hit(subject_id="S2", bitscore=100),
            make_hit(subject_id="S3", bitscore=100),
            make_hit(subject_id="S4", bitscore=100),
            make_hit(subject_id="S5", bitscore=100),
        ]
        policy = FilterPolicy(top_k=1)

        result = filter_hits_for_query(hits, policy)

        # All 5 are tied at the best bitscore
        assert result.accepted_subjects == {"S1", "S2", "S3", "S4", "S5"}

    def test_determinism_independent_of_input_order(self):
        hits_order1 = [
            make_hit(subject_id="S3", bitscore=80),
            make_hit(subject_id="S1", bitscore=100),
            make_hit(subject_id="S2", bitscore=90),
        ]
        hits_order2 = [
            make_hit(subject_id="S1", bitscore=100),
            make_hit(subject_id="S2", bitscore=90),
            make_hit(subject_id="S3", bitscore=80),
        ]
        policy = FilterPolicy(top_k=2)

        result1 = filter_hits_for_query(hits_order1, policy)
        result2 = filter_hits_for_query(hits_order2, policy)

        assert result1.accepted_subjects == result2.accepted_subjects

    def test_empty_hits(self):
        result = filter_hits_for_query([], FilterPolicy())

        assert result.accepted_subjects == set()
        assert result.total_hits == 0

    def test_all_hits_filtered_out(self):
        hits = [
            make_hit(subject_id="S1", evalue=1),
            make_hit(subject_id="S2", evalue=1),
        ]
        policy = FilterPolicy(max_evalue=1e-10)

        result = filter_hits_for_query(hits, policy)

        assert result.accepted_subjects == set()
        assert result.passed_threshold_hits == 0

    def test_combined_threshold_and_ranking(self):
        hits = [
            make_hit(subject_id="S1", evalue=1e-10, bitscore=100),
            make_hit(subject_id="S2", evalue=1, bitscore=200),  # Fails evalue
            make_hit(subject_id="S3", evalue=1e-8, bitscore=90),
            make_hit(subject_id="S4", evalue=1e-9, bitscore=80),
        ]
        policy = FilterPolicy(max_evalue=1e-5, top_k=2)

        result = filter_hits_for_query(hits, policy)

        # S2 filtered by evalue, then top 2 of remaining by bitscore
        assert result.accepted_subjects == {"S1", "S3"}

    def test_accepted_hits_match_accepted_subjects(self):
        hits = [
            make_hit(subject_id="S1", bitscore=100),
            make_hit(subject_id="S2", bitscore=90),
            make_hit(subject_id="S3", bitscore=80),
        ]
        result = filter_hits_for_query(hits, FilterPolicy(top_k=2))

        assert {h.subject_id for h in result.accepted_hits} == result.accepted_subjects
        assert result.accepted_subjects == {"S1", "S2"}

    def test_multiple_hsps_for_one_subject_take_one_rank_slot(self):
        # S1 aligns with two HSPs.  Ranking is over subjects, so with top_k=3
        # the three distinct subjects are kept; counting HSPs instead would
        # fill the third slot with S1's weaker HSP and drop S3.
        hits = [
            make_hit(subject_id="S1", bitscore=100),
            make_hit(subject_id="S1", bitscore=60),
            make_hit(subject_id="S2", bitscore=100),
            make_hit(subject_id="S3", bitscore=40),
        ]
        result = filter_hits_for_query(hits, FilterPolicy(top_k=3))

        assert result.accepted_subjects == {"S1", "S2", "S3"}
        assert len(result.accepted_hits) == 3

    def test_duplicate_hsps_do_not_push_a_subject_below_the_cut(self):
        # Without collapsing, S1's three HSPs occupy ranks 1-3 and S2 (rank 4)
        # is dropped by top_k=3 even though only two distinct subjects exist.
        hits = [
            make_hit(subject_id="S1", bitscore=100),
            make_hit(subject_id="S1", bitscore=95),
            make_hit(subject_id="S1", bitscore=90),
            make_hit(subject_id="S2", bitscore=50),
        ]
        result = filter_hits_for_query(hits, FilterPolicy(top_k=2))

        assert result.accepted_subjects == {"S1", "S2"}

    def test_representative_hit_is_best_passing_hsp(self):
        hits = [
            make_hit(subject_id="S1", evalue=1e-3, bitscore=30.0, pident=40.0),
            make_hit(subject_id="S1", evalue=1e-50, bitscore=200.0, pident=95.0),
            make_hit(subject_id="S1", evalue=1e-20, bitscore=120.0, pident=70.0),
        ]
        result = filter_hits_for_query(hits, FilterPolicy(max_evalue=1e-10))

        assert len(result.accepted_hits) == 1
        best = result.accepted_hits[0]
        assert best.evalue == 1e-50
        assert best.bitscore == 200.0
        assert best.pident == 95.0

    def test_representative_hit_tie_goes_to_lowest_evalue(self):
        hits = [
            make_hit(subject_id="S1", evalue=1e-20, bitscore=100.0, pident=80.0),
            make_hit(subject_id="S1", evalue=1e-30, bitscore=100.0, pident=85.0),
        ]
        result = filter_hits_for_query(hits, FilterPolicy())

        assert result.accepted_hits[0].evalue == 1e-30


class TestFilterAllHitsWithHits:
    """Tests for filter_all_hits_with_hits."""

    def test_returns_hit_objects_for_accepted_pairs(self):
        hits_by_query = {
            "Q1": [
                make_hit(query_id="Q1", subject_id="S1", bitscore=100, evalue=1e-40),
                make_hit(query_id="Q1", subject_id="S2", bitscore=50, evalue=1e-8),
            ],
        }
        result = filter_all_hits_with_hits(hits_by_query, FilterPolicy(top_k=1))

        assert set(result) == {"Q1"}
        assert set(result["Q1"]) == {"S1"}
        assert result["Q1"]["S1"].evalue == 1e-40

    def test_failing_hsp_never_overrides_the_accepted_one(self):
        # Regression: the map used to be rebuilt from every raw hit, so a later
        # HSP for the same subject that failed the thresholds replaced the
        # accepted one and its evalue/pident were reported.
        hits_by_query = {
            "Q1": [
                make_hit(query_id="Q1", subject_id="S1", evalue=1e-50, bitscore=200, pident=95),
                make_hit(query_id="Q1", subject_id="S1", evalue=1e-3, bitscore=30, pident=40),
            ],
        }
        result = filter_all_hits_with_hits(hits_by_query, FilterPolicy(max_evalue=1e-10))

        assert result["Q1"]["S1"].evalue == 1e-50
        assert result["Q1"]["S1"].pident == 95

    def test_queries_without_accepted_hits_are_omitted(self):
        hits_by_query = {
            "Q1": [make_hit(query_id="Q1", subject_id="S1", evalue=1)],
            "Q2": [make_hit(query_id="Q2", subject_id="S2", evalue=1e-20)],
        }
        result = filter_all_hits_with_hits(hits_by_query, FilterPolicy(max_evalue=1e-10))

        assert set(result) == {"Q2"}

    def test_agrees_with_filter_all_hits(self):
        hits_by_query = {
            "Q1": [
                make_hit(query_id="Q1", subject_id="S1", bitscore=100),
                make_hit(query_id="Q1", subject_id="S1", bitscore=90),
                make_hit(query_id="Q1", subject_id="S2", bitscore=100),
                make_hit(query_id="Q1", subject_id="S3", bitscore=10),
            ],
        }
        policy = FilterPolicy(top_k=1)
        subjects = filter_all_hits(hits_by_query, policy)
        with_hits = filter_all_hits_with_hits(hits_by_query, policy)

        assert subjects == {"Q1": {"S1", "S2"}}
        assert {q: set(d) for q, d in with_hits.items()} == subjects


class TestFilterPolicyValidate:
    """FilterPolicy.validate rejects values outside their meaningful range."""

    def test_all_none_is_valid(self):
        FilterPolicy().validate()

    def test_typical_policy_is_valid(self):
        FilterPolicy(max_evalue=1e-5, min_pident=50, min_qcov=80, min_alnlen=20, top_k=5).validate()

    def test_boundaries_are_valid(self):
        FilterPolicy(min_pident=0, min_qcov=100, min_alnlen=1, top_k=1).validate()

    @pytest.mark.parametrize(
        ("kwargs", "field_name"),
        [
            ({"top_k": 0}, "top_k"),
            ({"top_k": -1}, "top_k"),
            ({"min_alnlen": 0}, "min_alnlen"),
            ({"max_evalue": 0.0}, "max_evalue"),
            ({"max_evalue": -1.0}, "max_evalue"),
            ({"max_evalue": float("inf")}, "max_evalue"),
            ({"max_evalue": float("nan")}, "max_evalue"),
            ({"min_pident": -0.1}, "min_pident"),
            ({"min_pident": 100.1}, "min_pident"),
            ({"min_qcov": 101}, "min_qcov"),
            ({"min_qcov": float("nan")}, "min_qcov"),
        ],
    )
    def test_rejects_out_of_range(self, kwargs: dict, field_name: str):
        with pytest.raises(ValueError, match=field_name):
            FilterPolicy(**kwargs).validate()

    def test_reports_every_problem(self):
        with pytest.raises(ValueError) as exc_info:
            FilterPolicy(top_k=0, min_alnlen=0).validate()
        assert "top_k" in str(exc_info.value)
        assert "min_alnlen" in str(exc_info.value)


class TestFilterAllHits:
    """Tests for filtering hits across all queries."""

    def test_filters_multiple_queries(self):
        hits_by_query = {
            "Q1": [
                make_hit(query_id="Q1", subject_id="S1", bitscore=100),
                make_hit(query_id="Q1", subject_id="S2", bitscore=50),
            ],
            "Q2": [
                make_hit(query_id="Q2", subject_id="S3", bitscore=80),
            ],
        }
        policy = FilterPolicy(top_k=1)

        result = filter_all_hits(hits_by_query, policy)

        assert result["Q1"] == {"S1"}
        assert result["Q2"] == {"S3"}

    def test_empty_input(self):
        result = filter_all_hits({}, FilterPolicy())
        assert result == {}


class TestParseBlastTabular:
    """Tests for parsing BLAST/DIAMOND tabular output."""

    def test_parses_standard_format(self):
        lines = [
            "Q1\tS1\t95.0\t100\t5\t0\t1\t100\t1\t100\t1e-50\t200",
            "Q1\tS2\t90.0\t80\t8\t0\t1\t80\t1\t80\t1e-30\t150",
        ]

        hits_by_query = parse_blast_tabular(lines)

        assert "Q1" in hits_by_query
        assert len(hits_by_query["Q1"]) == 2
        assert hits_by_query["Q1"][0].subject_id == "S1"
        assert hits_by_query["Q1"][0].pident == 95.0
        assert hits_by_query["Q1"][0].bitscore == 200.0

    def test_skips_comments_and_empty_lines(self):
        lines = [
            "# Comment line",
            "",
            "Q1\tS1\t95.0\t100\t5\t0\t1\t100\t1\t100\t1e-50\t200",
        ]

        hits_by_query = parse_blast_tabular(lines)

        assert len(hits_by_query["Q1"]) == 1

    def test_multiple_queries(self):
        lines = [
            "Q1\tS1\t95.0\t100\t5\t0\t1\t100\t1\t100\t1e-50\t200",
            "Q2\tS2\t90.0\t80\t8\t0\t1\t80\t1\t80\t1e-30\t150",
        ]

        hits_by_query = parse_blast_tabular(lines)

        assert "Q1" in hits_by_query
        assert "Q2" in hits_by_query


class TestFilterPolicy:
    """Tests for FilterPolicy."""

    def test_to_dict(self):
        policy = FilterPolicy(
            max_evalue=1e-5,
            min_pident=80.0,
            top_k=10,
        )

        d = policy.to_dict()

        assert d["max_evalue"] == 1e-5
        assert d["min_pident"] == 80.0
        assert d["top_k"] == 10

    def test_qcovhsp_column_sets_query_coverage(self):
        columns = [
            "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovhsp",
        ]
        lines = ["Q1\tS1\t95.0\t100\t5\t0\t1\t100\t1\t100\t1e-50\t200\t87.5"]

        hits_by_query = parse_blast_tabular(lines, columns=columns)

        assert hits_by_query["Q1"][0].qcov == 87.5

    def test_qlen_column_computes_query_coverage(self):
        columns = ["qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "qlen"]
        lines = ["Q1\tS1\t95.0\t50\t1e-50\t200\t200"]

        hits_by_query = parse_blast_tabular(lines, columns=columns)

        assert hits_by_query["Q1"][0].qcov == 25.0

    def test_default_columns_have_no_coverage(self):
        """Plain outfmt 6 carries no coverage, so qcov is 0.0 (documented fallback)."""
        lines = ["Q1\tS1\t95.0\t100\t5\t0\t1\t100\t1\t100\t1e-50\t200"]

        hits_by_query = parse_blast_tabular(lines)

        assert hits_by_query["Q1"][0].qcov == 0.0

    def test_short_line_raises(self):
        columns = ["qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "qcovhsp"]
        lines = ["Q1\tS1\t95.0\t100\t1e-50\t200"]

        with pytest.raises(ValueError, match="Line 1: expected 7"):
            parse_blast_tabular(lines, columns=columns)

    def test_min_qcov_filters_on_parsed_coverage(self):
        """Regression: with a real coverage column, min_qcov keeps good hits."""
        columns = [
            "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovhsp",
        ]
        lines = [
            "Q1\tHIGH\t95.0\t180\t5\t0\t1\t180\t1\t180\t1e-50\t300\t90.0",
            "Q1\tLOW\t95.0\t60\t5\t0\t1\t60\t1\t60\t1e-20\t100\t30.0",
        ]

        hits_by_query = parse_blast_tabular(lines, columns=columns)
        accepted = filter_all_hits(hits_by_query, FilterPolicy(min_qcov=50.0))

        assert accepted == {"Q1": {"HIGH"}}
