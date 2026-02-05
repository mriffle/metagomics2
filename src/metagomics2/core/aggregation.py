"""Aggregation of peptide quantities into taxonomy nodes and GO terms."""

from dataclasses import dataclass, field
from typing import TypeVar

from metagomics2.core.annotation import PeptideAnnotation


@dataclass
class NodeAggregate:
    """Aggregated data for a single node (taxonomy or GO)."""

    node_id: str | int
    quantity: float = 0.0
    n_peptides: int = 0  # Distinct peptides contributing
    contributing_peptides: set[str] = field(default_factory=set)

    @property
    def ratio_total(self) -> float:
        """Ratio relative to total peptide quantity. Set externally."""
        return self._ratio_total

    @ratio_total.setter
    def ratio_total(self, value: float) -> None:
        self._ratio_total = value

    @property
    def ratio_annotated(self) -> float | None:
        """Ratio relative to annotated peptide quantity. Set externally."""
        return self._ratio_annotated

    @ratio_annotated.setter
    def ratio_annotated(self, value: float | None) -> None:
        self._ratio_annotated = value

    def __post_init__(self) -> None:
        self._ratio_total: float = 0.0
        self._ratio_annotated: float | None = None


@dataclass
class CoverageStats:
    """Coverage statistics for a peptide list."""

    total_peptide_quantity: float = 0.0
    annotated_peptide_quantity: float = 0.0
    unannotated_peptide_quantity: float = 0.0
    n_peptides_total: int = 0
    n_peptides_annotated: int = 0
    n_peptides_unannotated: int = 0

    @property
    def annotation_coverage_ratio(self) -> float:
        """Ratio of annotated to total quantity."""
        if self.total_peptide_quantity == 0:
            return 0.0
        return self.annotated_peptide_quantity / self.total_peptide_quantity


@dataclass
class AggregationResult:
    """Result of aggregating peptide annotations."""

    taxonomy_nodes: dict[int, NodeAggregate] = field(default_factory=dict)
    go_terms: dict[str, NodeAggregate] = field(default_factory=dict)
    coverage: CoverageStats = field(default_factory=CoverageStats)


def aggregate_peptide_annotations(
    annotations: list[PeptideAnnotation],
) -> AggregationResult:
    """Aggregate peptide annotations into node totals.

    For each taxonomy node t:
        quantity(t) = Σ_p q(p) * 1[t ∈ TAX(p)]

    For each GO term g:
        quantity(g) = Σ_p q(p) * 1[g ∈ GO(p)]

    Args:
        annotations: List of PeptideAnnotation objects

    Returns:
        AggregationResult with taxonomy nodes, GO terms, and coverage stats
    """
    result = AggregationResult()

    # Compute coverage stats
    for ann in annotations:
        result.coverage.total_peptide_quantity += ann.quantity
        result.coverage.n_peptides_total += 1

        if ann.is_annotated:
            result.coverage.annotated_peptide_quantity += ann.quantity
            result.coverage.n_peptides_annotated += 1
        else:
            result.coverage.unannotated_peptide_quantity += ann.quantity
            result.coverage.n_peptides_unannotated += 1

    # Aggregate taxonomy nodes
    for ann in annotations:
        if not ann.is_annotated:
            continue

        for tax_id in ann.taxonomy_nodes:
            if tax_id not in result.taxonomy_nodes:
                result.taxonomy_nodes[tax_id] = NodeAggregate(node_id=tax_id)

            node = result.taxonomy_nodes[tax_id]
            node.quantity += ann.quantity
            node.contributing_peptides.add(ann.peptide)

    # Update n_peptides for taxonomy
    for node in result.taxonomy_nodes.values():
        node.n_peptides = len(node.contributing_peptides)

    # Aggregate GO terms
    for ann in annotations:
        if not ann.is_annotated:
            continue

        for go_id in ann.go_terms:
            if go_id not in result.go_terms:
                result.go_terms[go_id] = NodeAggregate(node_id=go_id)

            node = result.go_terms[go_id]
            node.quantity += ann.quantity
            node.contributing_peptides.add(ann.peptide)

    # Update n_peptides for GO
    for node in result.go_terms.values():
        node.n_peptides = len(node.contributing_peptides)

    # Compute ratios
    total_qty = result.coverage.total_peptide_quantity
    annotated_qty = result.coverage.annotated_peptide_quantity

    for node in result.taxonomy_nodes.values():
        node.ratio_total = node.quantity / total_qty if total_qty > 0 else 0.0
        node.ratio_annotated = (
            node.quantity / annotated_qty if annotated_qty > 0 else None
        )

    for node in result.go_terms.values():
        node.ratio_total = node.quantity / total_qty if total_qty > 0 else 0.0
        node.ratio_annotated = (
            node.quantity / annotated_qty if annotated_qty > 0 else None
        )

    return result


def validate_aggregation_invariants(result: AggregationResult) -> list[str]:
    """Validate that aggregation invariants hold.

    Invariants:
    - For any node: 0 ≤ quantity(n) ≤ total_peptide_quantity
    - When annotated > 0: 0 ≤ ratio_total(n) ≤ ratio_annotated(n) ≤ 1

    Args:
        result: The aggregation result to validate

    Returns:
        List of violation messages (empty if all invariants hold)
    """
    violations: list[str] = []
    total_qty = result.coverage.total_peptide_quantity
    annotated_qty = result.coverage.annotated_peptide_quantity

    # Check taxonomy nodes
    for tax_id, node in result.taxonomy_nodes.items():
        if node.quantity < 0:
            violations.append(f"Taxonomy {tax_id}: negative quantity {node.quantity}")
        if node.quantity > total_qty:
            violations.append(
                f"Taxonomy {tax_id}: quantity {node.quantity} > total {total_qty}"
            )
        if node.ratio_total < 0 or node.ratio_total > 1:
            violations.append(
                f"Taxonomy {tax_id}: ratio_total {node.ratio_total} not in [0,1]"
            )
        if annotated_qty > 0 and node.ratio_annotated is not None:
            if node.ratio_annotated < 0 or node.ratio_annotated > 1:
                violations.append(
                    f"Taxonomy {tax_id}: ratio_annotated {node.ratio_annotated} not in [0,1]"
                )
            if node.ratio_total > node.ratio_annotated:
                violations.append(
                    f"Taxonomy {tax_id}: ratio_total {node.ratio_total} > ratio_annotated {node.ratio_annotated}"
                )

    # Check GO terms
    for go_id, node in result.go_terms.items():
        if node.quantity < 0:
            violations.append(f"GO {go_id}: negative quantity {node.quantity}")
        if node.quantity > total_qty:
            violations.append(
                f"GO {go_id}: quantity {node.quantity} > total {total_qty}"
            )
        if node.ratio_total < 0 or node.ratio_total > 1:
            violations.append(
                f"GO {go_id}: ratio_total {node.ratio_total} not in [0,1]"
            )
        if annotated_qty > 0 and node.ratio_annotated is not None:
            if node.ratio_annotated < 0 or node.ratio_annotated > 1:
                violations.append(
                    f"GO {go_id}: ratio_annotated {node.ratio_annotated} not in [0,1]"
                )
            if node.ratio_total > node.ratio_annotated:
                violations.append(
                    f"GO {go_id}: ratio_total {node.ratio_total} > ratio_annotated {node.ratio_annotated}"
                )

    return violations
