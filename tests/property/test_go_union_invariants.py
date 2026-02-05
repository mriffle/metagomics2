"""Property-based tests for GO union invariants."""

import pytest
from hypothesis import given, settings, strategies as st

from metagomics2.core.go import GODAG, GOTerm, load_go_from_dict


# Strategy for generating GO DAGs
@st.composite
def go_dag_strategy(draw):
    """Generate a random GO DAG."""
    n_terms = draw(st.integers(min_value=2, max_value=20))

    terms = {}
    for i in range(n_terms):
        term_id = f"GO:{i:07d}"
        terms[term_id] = {
            "name": f"term_{i}",
            "namespace": "biological_process",
        }

    # Generate edges (ensure DAG property by only allowing edges to lower-numbered terms)
    edges = {"is_a": []}
    for i in range(1, n_terms):
        # Each term has 1-3 parents from earlier terms
        n_parents = draw(st.integers(min_value=1, max_value=min(3, i)))
        parents = draw(
            st.lists(
                st.integers(min_value=0, max_value=i - 1),
                min_size=n_parents,
                max_size=n_parents,
                unique=True,
            )
        )
        for parent_idx in parents:
            edges["is_a"].append([f"GO:{i:07d}", f"GO:{parent_idx:07d}"])

    return load_go_from_dict({"terms": terms, "edges": edges})


class TestGOUnionInvariants:
    """Property-based tests for GO union operations."""

    @given(go_dag_strategy(), st.sets(st.integers(min_value=0, max_value=19), min_size=1, max_size=5))
    @settings(max_examples=50, deadline=None)
    def test_closure_contains_self(self, dag: GODAG, term_indices: set[int]):
        """Closure with include_self=True should contain the term itself."""
        for idx in term_indices:
            term_id = f"GO:{idx:07d}"
            if term_id in dag.terms:
                closure = dag.get_closure(term_id, include_self=True)
                assert term_id in closure

    @given(go_dag_strategy(), st.sets(st.integers(min_value=0, max_value=19), min_size=1, max_size=5))
    @settings(max_examples=50, deadline=None)
    def test_closure_excludes_self_when_requested(self, dag: GODAG, term_indices: set[int]):
        """Closure with include_self=False should not contain the term itself."""
        for idx in term_indices:
            term_id = f"GO:{idx:07d}"
            if term_id in dag.terms:
                closure = dag.get_closure(term_id, include_self=False)
                assert term_id not in closure

    @given(go_dag_strategy(), st.sets(st.integers(min_value=0, max_value=19), min_size=1, max_size=5))
    @settings(max_examples=50, deadline=None)
    def test_union_is_superset_of_individual_closures(self, dag: GODAG, term_indices: set[int]):
        """Union of closures should be superset of each individual closure."""
        term_ids = {f"GO:{idx:07d}" for idx in term_indices if f"GO:{idx:07d}" in dag.terms}

        if not term_ids:
            return

        union = dag.get_closure_union(term_ids)

        for term_id in term_ids:
            individual = dag.get_closure(term_id)
            assert individual <= union, f"Individual closure not subset of union"

    @given(go_dag_strategy(), st.sets(st.integers(min_value=0, max_value=19), min_size=1, max_size=5))
    @settings(max_examples=50, deadline=None)
    def test_union_is_set(self, dag: GODAG, term_indices: set[int]):
        """Union should be a proper set with no duplicates."""
        term_ids = {f"GO:{idx:07d}" for idx in term_indices if f"GO:{idx:07d}" in dag.terms}

        if not term_ids:
            return

        union = dag.get_closure_union(term_ids)

        # A set by definition has no duplicates
        assert len(union) == len(set(union))

    @given(go_dag_strategy())
    @settings(max_examples=50, deadline=None)
    def test_closure_is_transitive(self, dag: GODAG):
        """If B is in closure(A) and C is in closure(B), then C is in closure(A)."""
        for term_id in dag.terms:
            closure_a = dag.get_closure(term_id)

            for term_b in closure_a:
                if term_b in dag.terms:
                    closure_b = dag.get_closure(term_b)
                    # All ancestors of B should also be ancestors of A
                    assert closure_b <= closure_a, (
                        f"Transitivity violated: closure({term_b}) not subset of closure({term_id})"
                    )

    @given(go_dag_strategy())
    @settings(max_examples=50, deadline=None)
    def test_root_has_minimal_closure(self, dag: GODAG):
        """Root term (GO:0000000) should have closure containing only itself."""
        root_id = "GO:0000000"
        if root_id in dag.terms:
            closure = dag.get_closure(root_id, include_self=True)
            assert closure == {root_id}

            closure_no_self = dag.get_closure(root_id, include_self=False)
            assert closure_no_self == set()
