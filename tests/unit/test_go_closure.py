"""Unit tests for GO DAG loading and closure computation."""

import pytest

from metagomics2.core.go import GODAG, GOTerm, get_all_parent_ids, load_go_from_dict


class TestLoadGOFromDict:
    """Tests for loading GO DAG from dictionary."""

    def test_loads_terms(self, small_go: dict):
        dag = load_go_from_dict(small_go)

        assert "GO:0000001" in dag.terms
        assert dag.terms["GO:0000001"].name == "root_BP"
        assert dag.terms["GO:0000001"].namespace == "biological_process"

    def test_loads_is_a_edges(self, small_go: dict):
        dag = load_go_from_dict(small_go)

        # A -> root_BP via is_a
        assert "is_a" in dag.terms["GO:0000002"].parents
        assert "GO:0000001" in dag.terms["GO:0000002"].parents["is_a"]

    def test_loads_part_of_edges(self, small_go: dict):
        dag = load_go_from_dict(small_go)

        # D part_of B
        assert "part_of" in dag.terms["GO:0000005"].parents
        assert "GO:0000003" in dag.terms["GO:0000005"].parents["part_of"]

    def test_multi_parent_term(self, small_go: dict):
        dag = load_go_from_dict(small_go)

        # C has two is_a parents: A and B
        c_parents = dag.terms["GO:0000004"].parents.get("is_a", set())
        assert "GO:0000002" in c_parents  # A
        assert "GO:0000003" in c_parents  # B


class TestGOClosure:
    """Tests for GO closure computation."""

    def test_closure_includes_self(self, small_go: dict):
        dag = load_go_from_dict(small_go)

        closure = dag.get_closure("GO:0000005", include_self=True)
        assert "GO:0000005" in closure

    def test_closure_excludes_self_when_requested(self, small_go: dict):
        dag = load_go_from_dict(small_go)

        closure = dag.get_closure("GO:0000005", include_self=False)
        assert "GO:0000005" not in closure

    def test_closure_is_a_only_for_D(self, small_go: dict):
        """D (GO:0000005) -> A (GO:0000002) -> root_BP (GO:0000001) via is_a only."""
        dag = load_go_from_dict(small_go)

        closure = dag.get_closure("GO:0000005", edge_types={"is_a"}, include_self=True)

        # D, A, root_BP
        assert closure == {"GO:0000005", "GO:0000002", "GO:0000001"}

    def test_closure_is_a_plus_part_of_for_D(self, small_go: dict):
        """D (GO:0000005) has is_a->A and part_of->B."""
        dag = load_go_from_dict(small_go)

        closure = dag.get_closure(
            "GO:0000005", edge_types={"is_a", "part_of"}, include_self=True
        )

        # D, A, B, root_BP (via both A and B)
        assert "GO:0000005" in closure  # D
        assert "GO:0000002" in closure  # A
        assert "GO:0000003" in closure  # B (via part_of)
        assert "GO:0000001" in closure  # root_BP

    def test_closure_multi_parent_C(self, small_go: dict):
        """C (GO:0000004) has is_a parents A and B."""
        dag = load_go_from_dict(small_go)

        closure = dag.get_closure("GO:0000004", edge_types={"is_a"}, include_self=True)

        # C, A, B, root_BP
        assert closure == {"GO:0000004", "GO:0000002", "GO:0000003", "GO:0000001"}

    def test_closure_F_through_C(self, small_go: dict):
        """F (GO:0000007) -> C (GO:0000004) -> A,B -> root_BP."""
        dag = load_go_from_dict(small_go)

        closure = dag.get_closure("GO:0000007", edge_types={"is_a"}, include_self=True)

        # F, C, A, B, root_BP
        expected = {"GO:0000007", "GO:0000004", "GO:0000002", "GO:0000003", "GO:0000001"}
        assert closure == expected

    def test_closure_root_term(self, small_go: dict):
        """Root term closure is just itself."""
        dag = load_go_from_dict(small_go)

        closure = dag.get_closure("GO:0000001", include_self=True)
        assert closure == {"GO:0000001"}

    def test_closure_unknown_term(self, small_go: dict):
        """Unknown term returns empty or just self."""
        dag = load_go_from_dict(small_go)

        closure = dag.get_closure("GO:9999999", include_self=True)
        assert closure == {"GO:9999999"}

        closure_no_self = dag.get_closure("GO:9999999", include_self=False)
        assert closure_no_self == set()

    def test_closure_default_edge_type_is_is_a(self, small_go: dict):
        """Default edge type should be is_a only."""
        dag = load_go_from_dict(small_go)

        # D with default (is_a only) should NOT include B via part_of
        closure = dag.get_closure("GO:0000005", include_self=True)

        assert "GO:0000005" in closure  # D
        assert "GO:0000002" in closure  # A
        assert "GO:0000001" in closure  # root_BP
        # B should NOT be included (only reachable via part_of)
        assert "GO:0000003" not in closure


class TestGOClosureUnion:
    """Tests for union of multiple closures."""

    def test_union_of_single_term(self, small_go: dict):
        dag = load_go_from_dict(small_go)

        union = dag.get_closure_union({"GO:0000005"}, edge_types={"is_a"})
        single = dag.get_closure("GO:0000005", edge_types={"is_a"})

        assert union == single

    def test_union_of_multiple_terms(self, small_go: dict):
        dag = load_go_from_dict(small_go)

        # Union of D and E closures
        # D (GO:0000005) -> A -> root_BP
        # E (GO:0000006) -> B -> root_BP
        union = dag.get_closure_union(
            {"GO:0000005", "GO:0000006"}, edge_types={"is_a"}
        )

        expected = {
            "GO:0000005",  # D
            "GO:0000006",  # E
            "GO:0000002",  # A
            "GO:0000003",  # B
            "GO:0000001",  # root_BP
        }
        assert union == expected

    def test_union_is_non_redundant(self, small_go: dict):
        """Union should be a set with no duplicates."""
        dag = load_go_from_dict(small_go)

        # Both D and E have root_BP in their closure
        union = dag.get_closure_union(
            {"GO:0000005", "GO:0000006"}, edge_types={"is_a"}
        )

        # Count should match set size (no duplicates)
        assert len(union) == len(set(union))

    def test_union_empty_set(self, small_go: dict):
        dag = load_go_from_dict(small_go)

        union = dag.get_closure_union(set())
        assert union == set()


class TestCycleProtection:
    """Tests for cycle handling in GO DAG."""

    def test_closure_handles_cycle(self):
        """Closure computation should terminate even with cycles."""
        # Create a DAG with a cycle (artificial, shouldn't happen in real GO)
        data = {
            "terms": {
                "GO:X": {"name": "X", "namespace": "test"},
                "GO:Y": {"name": "Y", "namespace": "test"},
            },
            "edges": {
                "is_a": [
                    ["GO:X", "GO:Y"],
                    ["GO:Y", "GO:X"],  # Creates cycle
                ]
            },
        }
        dag = load_go_from_dict(data)

        # Should terminate and include both
        closure = dag.get_closure("GO:X", include_self=True)
        assert closure == {"GO:X", "GO:Y"}


class TestGetAllParentIds:
    """Tests for getting all direct parent IDs."""

    def test_gets_all_parents(self, small_go: dict):
        dag = load_go_from_dict(small_go)

        # D has is_a->A and part_of->B
        parents = get_all_parent_ids(dag, "GO:0000005")
        assert parents == {"GO:0000002", "GO:0000003"}

    def test_unknown_term(self, small_go: dict):
        dag = load_go_from_dict(small_go)

        parents = get_all_parent_ids(dag, "GO:9999999")
        assert parents == set()
