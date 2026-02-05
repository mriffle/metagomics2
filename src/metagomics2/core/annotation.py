"""Peptide annotation semantics: taxonomy LCA and GO union."""

from dataclasses import dataclass, field

from metagomics2.core.go import GODAG
from metagomics2.core.taxonomy import TaxonomyTree


@dataclass
class SubjectAnnotation:
    """Annotation data for a subject protein."""

    subject_id: str
    tax_id: int | None = None
    go_terms: set[str] = field(default_factory=set)


@dataclass
class PeptideAnnotation:
    """Computed annotation for a peptide."""

    peptide: str
    quantity: float
    is_annotated: bool = False

    # Taxonomy annotation
    lca_tax_id: int | None = None
    taxonomy_nodes: set[int] = field(default_factory=set)  # Lineage from LCA to root

    # GO annotation
    go_terms: set[str] = field(default_factory=set)  # Union of closures

    # Provenance
    implied_subjects: set[str] = field(default_factory=set)
    background_proteins: set[str] = field(default_factory=set)


def get_implied_subjects(
    peptide: str,
    peptide_to_proteins: dict[str, set[str]],
    protein_to_subjects: dict[str, set[str]],
) -> set[str]:
    """Get all implied subject proteins for a peptide.

    S(p) = ⋃_{b∈B(p)} S(b)

    Args:
        peptide: The peptide sequence
        peptide_to_proteins: Mapping from peptide to background proteins
        protein_to_subjects: Mapping from background protein to accepted subjects

    Returns:
        Set of implied subject IDs
    """
    background_proteins = peptide_to_proteins.get(peptide, set())
    implied: set[str] = set()

    for protein_id in background_proteins:
        subjects = protein_to_subjects.get(protein_id, set())
        implied |= subjects

    return implied


def annotate_peptide_taxonomy(
    implied_subjects: set[str],
    subject_annotations: dict[str, SubjectAnnotation],
    taxonomy_tree: TaxonomyTree,
) -> tuple[int | None, set[int]]:
    """Compute taxonomy annotation for a peptide.

    Uses LCA intersection: find the lowest common ancestor of all
    implied subjects' tax_ids, then return lineage from LCA to root.

    Args:
        implied_subjects: Set of subject IDs
        subject_annotations: Mapping from subject ID to annotation
        taxonomy_tree: The taxonomy tree

    Returns:
        Tuple of (LCA tax_id or None, set of taxonomy nodes in lineage)
    """
    if not implied_subjects:
        return None, set()

    # Collect tax_ids from all subjects
    tax_ids: set[int] = set()
    for subject_id in implied_subjects:
        annotation = subject_annotations.get(subject_id)
        if annotation and annotation.tax_id is not None:
            tax_ids.add(annotation.tax_id)

    if not tax_ids:
        return None, set()

    # Compute LCA
    lca = taxonomy_tree.compute_lca(tax_ids)
    if lca is None:
        return None, set()

    # Get lineage from LCA to root
    lineage = taxonomy_tree.get_lineage(lca)
    return lca, set(lineage)


def annotate_peptide_go(
    implied_subjects: set[str],
    subject_annotations: dict[str, SubjectAnnotation],
    go_dag: GODAG,
    edge_types: set[str] | None = None,
    include_self: bool = True,
) -> set[str]:
    """Compute GO annotation for a peptide.

    Uses non-redundant union: union of GO closures across all implied subjects.

    GO(p) = ⋃_{s∈S(p)} go_closure(s)

    Args:
        implied_subjects: Set of subject IDs
        subject_annotations: Mapping from subject ID to annotation
        go_dag: The GO DAG
        edge_types: Edge types to follow in closure (default: {"is_a"})
        include_self: Whether to include terms themselves in closure

    Returns:
        Set of GO term IDs (non-redundant union of closures)
    """
    if not implied_subjects:
        return set()

    # Collect all direct GO terms from subjects
    all_direct_terms: set[str] = set()
    for subject_id in implied_subjects:
        annotation = subject_annotations.get(subject_id)
        if annotation:
            all_direct_terms |= annotation.go_terms

    if not all_direct_terms:
        return set()

    # Compute union of closures
    return go_dag.get_closure_union(all_direct_terms, edge_types, include_self)


def annotate_peptide(
    peptide: str,
    quantity: float,
    peptide_to_proteins: dict[str, set[str]],
    protein_to_subjects: dict[str, set[str]],
    subject_annotations: dict[str, SubjectAnnotation],
    taxonomy_tree: TaxonomyTree,
    go_dag: GODAG,
    go_edge_types: set[str] | None = None,
    go_include_self: bool = True,
) -> PeptideAnnotation:
    """Compute full annotation for a peptide.

    Args:
        peptide: The peptide sequence
        quantity: The peptide quantity
        peptide_to_proteins: Mapping from peptide to background proteins
        protein_to_subjects: Mapping from background protein to accepted subjects
        subject_annotations: Mapping from subject ID to annotation
        taxonomy_tree: The taxonomy tree
        go_dag: The GO DAG
        go_edge_types: Edge types for GO closure
        go_include_self: Whether to include GO terms themselves

    Returns:
        PeptideAnnotation with computed taxonomy and GO annotations
    """
    background_proteins = peptide_to_proteins.get(peptide, set())
    implied_subjects = get_implied_subjects(
        peptide, peptide_to_proteins, protein_to_subjects
    )

    is_annotated = len(implied_subjects) > 0

    # Taxonomy annotation
    lca_tax_id, taxonomy_nodes = annotate_peptide_taxonomy(
        implied_subjects, subject_annotations, taxonomy_tree
    )

    # GO annotation
    go_terms = annotate_peptide_go(
        implied_subjects,
        subject_annotations,
        go_dag,
        go_edge_types,
        go_include_self,
    )

    return PeptideAnnotation(
        peptide=peptide,
        quantity=quantity,
        is_annotated=is_annotated,
        lca_tax_id=lca_tax_id,
        taxonomy_nodes=taxonomy_nodes,
        go_terms=go_terms,
        implied_subjects=implied_subjects,
        background_proteins=background_proteins,
    )


def load_subject_annotations_from_dict(data: dict) -> dict[str, SubjectAnnotation]:
    """Load subject annotations from a dictionary.

    Expected format:
    {
        "subjects": {
            "U1": {"tax_id": 70, "go_terms": ["GO:0000004"]},
            ...
        }
    }

    Args:
        data: Dictionary with 'subjects' key

    Returns:
        Dictionary mapping subject ID to SubjectAnnotation
    """
    annotations: dict[str, SubjectAnnotation] = {}

    for subject_id, subject_data in data.get("subjects", {}).items():
        annotations[subject_id] = SubjectAnnotation(
            subject_id=subject_id,
            tax_id=subject_data.get("tax_id"),
            go_terms=set(subject_data.get("go_terms", [])),
        )

    return annotations
