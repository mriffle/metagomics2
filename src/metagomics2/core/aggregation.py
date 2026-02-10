"""Aggregation of peptide quantities into taxonomy nodes and GO terms."""

from dataclasses import dataclass, field

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


@dataclass
class ComboAggregate:
    """Aggregated data for a (taxonomy node, GO term) pair."""

    tax_id: int
    go_id: str
    quantity: float = 0.0
    n_peptides: int = 0
    contributing_peptides: set[str] = field(default_factory=set)
    fraction_of_taxon: float = 0.0
    fraction_of_go: float = 0.0
    ratio_total_taxon: float = 0.0
    ratio_total_go: float = 0.0


def aggregate_go_taxonomy_combos(
    annotations: list[PeptideAnnotation],
    aggregation: AggregationResult,
) -> dict[tuple[int, str], ComboAggregate]:
    """Cross-tabulate taxonomy nodes and GO terms.

    For each annotated peptide p, for every (tax_id, go_id) in
    TAX(p) × GO(p), accumulate q(p).

    Then compute:
        fraction_of_taxon = combo_quantity / taxon_total_quantity
        fraction_of_go    = combo_quantity / go_total_quantity

    Args:
        annotations: List of PeptideAnnotation objects
        aggregation: The already-computed AggregationResult (for per-node totals)

    Returns:
        Dictionary mapping (tax_id, go_id) to ComboAggregate
    """
    combos: dict[tuple[int, str], ComboAggregate] = {}

    for ann in annotations:
        if not ann.is_annotated:
            continue
        if not ann.taxonomy_nodes or not ann.go_terms:
            continue

        for tax_id in ann.taxonomy_nodes:
            for go_id in ann.go_terms:
                key = (tax_id, go_id)
                if key not in combos:
                    combos[key] = ComboAggregate(tax_id=tax_id, go_id=go_id)

                combo = combos[key]
                combo.quantity += ann.quantity
                combo.contributing_peptides.add(ann.peptide)

    # Finalize n_peptides and fractions
    for combo in combos.values():
        combo.n_peptides = len(combo.contributing_peptides)

        tax_node = aggregation.taxonomy_nodes.get(combo.tax_id)
        if tax_node and tax_node.quantity > 0:
            combo.fraction_of_taxon = combo.quantity / tax_node.quantity
        if tax_node:
            combo.ratio_total_taxon = tax_node.ratio_total

        go_node = aggregation.go_terms.get(combo.go_id)
        if go_node and go_node.quantity > 0:
            combo.fraction_of_go = combo.quantity / go_node.quantity
        if go_node:
            combo.ratio_total_go = go_node.ratio_total

    return combos
