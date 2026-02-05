"""Unit tests for taxonomy tree loading and LCA computation."""

import pytest

from metagomics2.core.taxonomy import TaxonomyTree, load_taxonomy_from_dict


class TestLoadTaxonomyFromDict:
    """Tests for loading taxonomy from dictionary."""

    def test_loads_nodes(self, small_taxonomy: dict):
        tree = load_taxonomy_from_dict(small_taxonomy)

        assert 1 in tree.nodes
        assert tree.nodes[1].name == "root"
        assert tree.nodes[1].rank == "no rank"
        assert tree.nodes[1].parent_tax_id is None

    def test_loads_parent_relationships(self, small_taxonomy: dict):
        tree = load_taxonomy_from_dict(small_taxonomy)

        # KingdomA (10) -> root (1)
        assert tree.nodes[10].parent_tax_id == 1
        # SpeciesA (70) -> GenusA (60)
        assert tree.nodes[70].parent_tax_id == 60

    def test_loads_all_nodes(self, small_taxonomy: dict):
        tree = load_taxonomy_from_dict(small_taxonomy)

        # Should have 20 nodes total based on fixture
        assert len(tree.nodes) == 20


class TestGetLineage:
    """Tests for lineage computation."""

    def test_lineage_from_species_to_root(self, small_taxonomy: dict):
        tree = load_taxonomy_from_dict(small_taxonomy)

        # SpeciesA (70) lineage: 70 -> 60 -> 50 -> 40 -> 30 -> 20 -> 10 -> 1
        lineage = tree.get_lineage(70)

        assert lineage == [70, 60, 50, 40, 30, 20, 10, 1]

    def test_lineage_from_root(self, small_taxonomy: dict):
        tree = load_taxonomy_from_dict(small_taxonomy)

        lineage = tree.get_lineage(1)
        assert lineage == [1]

    def test_lineage_unknown_taxon(self, small_taxonomy: dict):
        tree = load_taxonomy_from_dict(small_taxonomy)

        lineage = tree.get_lineage(99999)
        assert lineage == [99999]

    def test_lineage_set(self, small_taxonomy: dict):
        tree = load_taxonomy_from_dict(small_taxonomy)

        lineage_set = tree.get_lineage_set(70)
        assert lineage_set == {70, 60, 50, 40, 30, 20, 10, 1}

    def test_lineage_has_unique_nodes(self, small_taxonomy: dict):
        tree = load_taxonomy_from_dict(small_taxonomy)

        lineage = tree.get_lineage(70)
        assert len(lineage) == len(set(lineage))


class TestComputeLCA:
    """Tests for LCA computation."""

    def test_lca_of_identical_taxa(self, small_taxonomy: dict):
        """LCA(70, 70) = 70"""
        tree = load_taxonomy_from_dict(small_taxonomy)

        lca = tree.compute_lca({70, 70})
        assert lca == 70

    def test_lca_of_single_taxon(self, small_taxonomy: dict):
        """LCA of single taxon is itself."""
        tree = load_taxonomy_from_dict(small_taxonomy)

        lca = tree.compute_lca({70})
        assert lca == 70

    def test_lca_of_sibling_species(self, small_taxonomy: dict):
        """LCA(70, 71) = 60 (GenusA) - both are under GenusA."""
        tree = load_taxonomy_from_dict(small_taxonomy)

        # SpeciesA (70) and SpeciesB (71) are both under GenusA (60)
        lca = tree.compute_lca({70, 71})
        assert lca == 60

    def test_lca_of_cousins_same_class(self, small_taxonomy: dict):
        """LCA(70, 72) = 30 (ClassA) - different orders but same class."""
        tree = load_taxonomy_from_dict(small_taxonomy)

        # SpeciesA (70): under OrderA (40)
        # SpeciesC (72): under OrderB (41)
        # Both orders are under ClassA (30)
        lca = tree.compute_lca({70, 72})
        assert lca == 30

    def test_lca_across_kingdoms(self, small_taxonomy: dict):
        """LCA(70, 73) = 1 (root) - different kingdoms."""
        tree = load_taxonomy_from_dict(small_taxonomy)

        # SpeciesA (70): KingdomA
        # SpeciesD (73): KingdomB
        lca = tree.compute_lca({70, 73})
        assert lca == 1

    def test_lca_of_empty_set(self, small_taxonomy: dict):
        """LCA of empty set is None."""
        tree = load_taxonomy_from_dict(small_taxonomy)

        lca = tree.compute_lca(set())
        assert lca is None

    def test_lca_three_taxa_same_genus(self, small_taxonomy: dict):
        """LCA of multiple taxa under same genus."""
        tree = load_taxonomy_from_dict(small_taxonomy)

        # Only 70 and 71 are under GenusA, but let's test with 70, 71
        lca = tree.compute_lca({70, 71})
        assert lca == 60

    def test_lca_three_taxa_different_branches(self, small_taxonomy: dict):
        """LCA(70, 71, 72) - two from GenusA, one from GenusB."""
        tree = load_taxonomy_from_dict(small_taxonomy)

        # 70, 71 under GenusA (60) under OrderA (40)
        # 72 under GenusB (61) under OrderB (41)
        # All under ClassA (30)
        lca = tree.compute_lca({70, 71, 72})
        assert lca == 30

    def test_lca_ancestor_descendant(self, small_taxonomy: dict):
        """LCA of a taxon and its ancestor is the ancestor."""
        tree = load_taxonomy_from_dict(small_taxonomy)

        # SpeciesA (70) and ClassA (30)
        lca = tree.compute_lca({70, 30})
        assert lca == 30


class TestGetLCALineage:
    """Tests for getting LCA lineage."""

    def test_lca_lineage_siblings(self, small_taxonomy: dict):
        tree = load_taxonomy_from_dict(small_taxonomy)

        # LCA(70, 71) = 60, lineage from 60 to root
        lineage = tree.get_lca_lineage({70, 71})

        assert lineage[0] == 60  # LCA
        assert lineage[-1] == 1  # root
        assert lineage == [60, 50, 40, 30, 20, 10, 1]

    def test_lca_lineage_empty_set(self, small_taxonomy: dict):
        tree = load_taxonomy_from_dict(small_taxonomy)

        lineage = tree.get_lca_lineage(set())
        assert lineage == []

    def test_lca_lineage_single_taxon(self, small_taxonomy: dict):
        tree = load_taxonomy_from_dict(small_taxonomy)

        lineage = tree.get_lca_lineage({70})
        assert lineage == [70, 60, 50, 40, 30, 20, 10, 1]


class TestCycleProtection:
    """Tests for cycle handling in taxonomy tree."""

    def test_lineage_handles_cycle(self):
        """Lineage computation should terminate even with cycles."""
        # Create a tree with a cycle (artificial, shouldn't happen in real taxonomy)
        data = {
            "nodes": {
                "1": {"name": "A", "rank": "test", "parent_tax_id": 2},
                "2": {"name": "B", "rank": "test", "parent_tax_id": 1},
            }
        }
        tree = load_taxonomy_from_dict(data)

        # Should terminate and include both
        lineage = tree.get_lineage(1)
        assert set(lineage) == {1, 2}
        assert len(lineage) == 2  # No infinite loop


class TestTaxonomyTreeEdgeCases:
    """Edge case tests for taxonomy tree."""

    def test_empty_tree(self):
        tree = TaxonomyTree()

        assert tree.compute_lca({1, 2}) is None
        assert tree.get_lineage(1) == [1]

    def test_single_node_tree(self):
        data = {"nodes": {"1": {"name": "root", "rank": "no rank", "parent_tax_id": None}}}
        tree = load_taxonomy_from_dict(data)

        assert tree.compute_lca({1}) == 1
        assert tree.get_lineage(1) == [1]
