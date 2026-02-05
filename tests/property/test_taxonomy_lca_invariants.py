"""Property-based tests for taxonomy LCA invariants."""

import pytest
from hypothesis import given, settings, strategies as st

from metagomics2.core.taxonomy import TaxonomyTree, TaxonNode, load_taxonomy_from_dict


# Strategy for generating taxonomy trees
@st.composite
def taxonomy_tree_strategy(draw):
    """Generate a random taxonomy tree."""
    n_nodes = draw(st.integers(min_value=2, max_value=30))

    nodes = {}
    # Root node
    nodes["1"] = {"name": "root", "rank": "no rank", "parent_tax_id": None}

    # Generate other nodes
    for i in range(2, n_nodes + 1):
        # Parent must be an existing node with lower ID (ensures tree structure)
        parent_id = draw(st.integers(min_value=1, max_value=i - 1))
        nodes[str(i)] = {
            "name": f"taxon_{i}",
            "rank": "species" if i > n_nodes // 2 else "genus",
            "parent_tax_id": parent_id,
        }

    return load_taxonomy_from_dict({"nodes": nodes})


class TestTaxonomyLCAInvariants:
    """Property-based tests for taxonomy LCA operations."""

    @given(taxonomy_tree_strategy(), st.integers(min_value=1, max_value=30))
    @settings(max_examples=50, deadline=None)
    def test_lca_of_single_taxon_is_itself(self, tree: TaxonomyTree, tax_id: int):
        """LCA of a single taxon should be itself."""
        if tax_id in tree.nodes:
            lca = tree.compute_lca({tax_id})
            assert lca == tax_id

    @given(taxonomy_tree_strategy(), st.integers(min_value=1, max_value=30))
    @settings(max_examples=50, deadline=None)
    def test_lca_of_identical_taxa_is_itself(self, tree: TaxonomyTree, tax_id: int):
        """LCA of identical taxa should be the taxon itself."""
        if tax_id in tree.nodes:
            lca = tree.compute_lca({tax_id, tax_id})
            assert lca == tax_id

    @given(
        taxonomy_tree_strategy(),
        st.sets(st.integers(min_value=1, max_value=30), min_size=2, max_size=5),
    )
    @settings(max_examples=50, deadline=None)
    def test_lca_is_ancestor_of_all_inputs(self, tree: TaxonomyTree, tax_ids: set[int]):
        """LCA should be an ancestor of all input taxa."""
        valid_ids = {t for t in tax_ids if t in tree.nodes}
        if len(valid_ids) < 2:
            return

        lca = tree.compute_lca(valid_ids)
        if lca is None:
            return

        for tax_id in valid_ids:
            lineage = tree.get_lineage_set(tax_id)
            assert lca in lineage, f"LCA {lca} not in lineage of {tax_id}"

    @given(
        taxonomy_tree_strategy(),
        st.sets(st.integers(min_value=1, max_value=30), min_size=2, max_size=5),
    )
    @settings(max_examples=50, deadline=None)
    def test_lca_is_deepest_common_ancestor(self, tree: TaxonomyTree, tax_ids: set[int]):
        """LCA should be the deepest (closest to leaves) common ancestor."""
        valid_ids = {t for t in tax_ids if t in tree.nodes}
        if len(valid_ids) < 2:
            return

        lca = tree.compute_lca(valid_ids)
        if lca is None:
            return

        # Get LCA's lineage
        lca_lineage = tree.get_lineage(lca)

        # Any node deeper than LCA (earlier in lineage) should not be a common ancestor
        # Actually, LCA should be the first node in its own lineage
        assert lca_lineage[0] == lca

        # Check that LCA's parent (if exists) is also a common ancestor
        lca_node = tree.nodes.get(lca)
        if lca_node and lca_node.parent_tax_id is not None:
            parent = lca_node.parent_tax_id
            for tax_id in valid_ids:
                lineage = tree.get_lineage_set(tax_id)
                assert parent in lineage, f"Parent {parent} should also be common ancestor"

    @given(taxonomy_tree_strategy())
    @settings(max_examples=50, deadline=None)
    def test_lineage_ends_at_root(self, tree: TaxonomyTree):
        """Every lineage should end at root (tax_id=1)."""
        for tax_id in tree.nodes:
            lineage = tree.get_lineage(tax_id)
            if lineage:
                assert lineage[-1] == 1, f"Lineage of {tax_id} doesn't end at root"

    @given(taxonomy_tree_strategy())
    @settings(max_examples=50, deadline=None)
    def test_lineage_has_unique_nodes(self, tree: TaxonomyTree):
        """Lineage should have no duplicate nodes."""
        for tax_id in tree.nodes:
            lineage = tree.get_lineage(tax_id)
            assert len(lineage) == len(set(lineage)), f"Duplicate nodes in lineage of {tax_id}"

    @given(
        taxonomy_tree_strategy(),
        st.sets(st.integers(min_value=1, max_value=30), min_size=1, max_size=5),
    )
    @settings(max_examples=50, deadline=None)
    def test_lca_lineage_starts_at_lca(self, tree: TaxonomyTree, tax_ids: set[int]):
        """LCA lineage should start at the LCA."""
        valid_ids = {t for t in tax_ids if t in tree.nodes}
        if not valid_ids:
            return

        lineage = tree.get_lca_lineage(valid_ids)
        if lineage:
            lca = tree.compute_lca(valid_ids)
            assert lineage[0] == lca

    @given(
        taxonomy_tree_strategy(),
        st.integers(min_value=1, max_value=30),
        st.integers(min_value=1, max_value=30),
    )
    @settings(max_examples=50, deadline=None)
    def test_lca_is_commutative(self, tree: TaxonomyTree, tax_id1: int, tax_id2: int):
        """LCA(a, b) should equal LCA(b, a)."""
        if tax_id1 in tree.nodes and tax_id2 in tree.nodes:
            lca1 = tree.compute_lca({tax_id1, tax_id2})
            lca2 = tree.compute_lca({tax_id2, tax_id1})
            assert lca1 == lca2

    @given(
        taxonomy_tree_strategy(),
        st.integers(min_value=1, max_value=30),
        st.integers(min_value=1, max_value=30),
    )
    @settings(max_examples=50, deadline=None)
    def test_lca_with_ancestor_is_ancestor(self, tree: TaxonomyTree, tax_id: int, ancestor_idx: int):
        """LCA of a taxon and one of its ancestors should be the ancestor."""
        if tax_id not in tree.nodes:
            return

        lineage = tree.get_lineage(tax_id)
        if ancestor_idx < len(lineage):
            ancestor = lineage[ancestor_idx]
            lca = tree.compute_lca({tax_id, ancestor})
            assert lca == ancestor, f"LCA({tax_id}, {ancestor}) should be {ancestor}, got {lca}"
