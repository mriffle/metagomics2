"""Taxonomy tree loading and LCA (Lowest Common Ancestor) computation."""

import json
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class TaxonNode:
    """A node in the taxonomy tree."""

    tax_id: int
    name: str
    rank: str
    parent_tax_id: int | None


@dataclass
class TaxonomyTree:
    """NCBI-style taxonomy tree."""

    nodes: dict[int, TaxonNode] = field(default_factory=dict)

    def get_lineage(self, tax_id: int) -> list[int]:
        """Get the lineage from a taxon to the root.

        Args:
            tax_id: The taxonomy ID to get lineage for

        Returns:
            List of tax_ids from the given taxon to root (inclusive),
            ordered from the given taxon toward root.
        """
        lineage: list[int] = []
        visited: set[int] = set()
        current = tax_id

        while current is not None and current not in visited:
            visited.add(current)
            lineage.append(current)

            node = self.nodes.get(current)
            if node is None:
                break
            current = node.parent_tax_id

        return lineage

    def get_lineage_set(self, tax_id: int) -> set[int]:
        """Get the lineage as a set for fast membership testing.

        Args:
            tax_id: The taxonomy ID

        Returns:
            Set of tax_ids in the lineage
        """
        return set(self.get_lineage(tax_id))

    def compute_lca(self, tax_ids: set[int]) -> int | None:
        """Compute the Lowest Common Ancestor of a set of taxonomy IDs.

        The LCA is the deepest node that is an ancestor of all given tax_ids.

        Args:
            tax_ids: Set of taxonomy IDs

        Returns:
            The LCA tax_id, or None if no common ancestor exists
        """
        if not tax_ids:
            return None

        tax_id_list = list(tax_ids)

        if len(tax_id_list) == 1:
            return tax_id_list[0]

        # Get lineage of first taxon as the initial candidate set
        first_lineage = self.get_lineage(tax_id_list[0])
        if not first_lineage:
            return None

        # Convert to set for intersection
        common_ancestors = set(first_lineage)

        # Intersect with lineages of all other taxa
        for tax_id in tax_id_list[1:]:
            lineage_set = self.get_lineage_set(tax_id)
            common_ancestors &= lineage_set

            if not common_ancestors:
                return None

        # Find the deepest common ancestor (closest to leaves)
        # The deepest is the one that appears first in any lineage
        # (since lineages go from leaf toward root)
        for tax_id in first_lineage:
            if tax_id in common_ancestors:
                return tax_id

        return None

    def get_lca_lineage(self, tax_ids: set[int]) -> list[int]:
        """Compute LCA and return lineage from LCA to root.

        Args:
            tax_ids: Set of taxonomy IDs

        Returns:
            Lineage from LCA to root, or empty list if no LCA
        """
        lca = self.compute_lca(tax_ids)
        if lca is None:
            return []
        return self.get_lineage(lca)


def load_taxonomy_from_json(file_path: Path | str) -> TaxonomyTree:
    """Load taxonomy tree from a JSON file.

    Expected format:
    {
        "nodes": {
            "1": {"name": "root", "rank": "no rank", "parent_tax_id": null},
            "10": {"name": "KingdomA", "rank": "kingdom", "parent_tax_id": 1},
            ...
        }
    }

    Args:
        file_path: Path to the JSON file

    Returns:
        TaxonomyTree object
    """
    file_path = Path(file_path)

    with open(file_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    return load_taxonomy_from_dict(data)


def load_taxonomy_from_dict(data: dict) -> TaxonomyTree:
    """Load taxonomy tree from a dictionary.

    Args:
        data: Dictionary with 'nodes' key

    Returns:
        TaxonomyTree object
    """
    tree = TaxonomyTree()

    for tax_id_str, node_data in data.get("nodes", {}).items():
        tax_id = int(tax_id_str)
        parent_tax_id = node_data.get("parent_tax_id")
        if parent_tax_id is not None:
            parent_tax_id = int(parent_tax_id)

        tree.nodes[tax_id] = TaxonNode(
            tax_id=tax_id,
            name=node_data.get("name", ""),
            rank=node_data.get("rank", ""),
            parent_tax_id=parent_tax_id,
        )

    return tree
