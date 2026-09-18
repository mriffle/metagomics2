"""Gene Ontology (GO) DAG loading and closure computation."""
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

# Relationship types that can appear as edges in a GO release: ``is_a`` plus
# every ``relationship:`` predicate used by the ontology.  Anything else in a
# user-supplied edge-type list is a typo and must be rejected, because an
# unknown type silently matches no edges and quietly shrinks every closure.
GO_EDGE_TYPES = frozenset({
    "is_a",
    "part_of",
    "regulates",
    "positively_regulates",
    "negatively_regulates",
    "has_part",
    "occurs_in",
    "happens_during",
    "ends_during",
})
DEFAULT_GO_EDGE_TYPES = frozenset({"is_a", "part_of"})


def parse_go_edge_types(raw: str) -> set[str]:
    """Parse a comma-separated list of GO edge types.

    Entries are stripped of whitespace and empty entries are ignored, so
    ``"is_a, part_of"`` and ``"is_a,part_of,"`` both parse.  Unknown types
    and an empty result raise ``ValueError`` with the allowed list.
    """
    types = {part.strip() for part in raw.split(",") if part.strip()}
    allowed = ", ".join(sorted(GO_EDGE_TYPES))
    if not types:
        raise ValueError(f"No GO edge types given; allowed types are {allowed}")
    unknown = sorted(types - GO_EDGE_TYPES)
    if unknown:
        raise ValueError(
            f"Unknown GO edge type(s) {', '.join(unknown)}; allowed types are {allowed}"
        )
    return types


@dataclass
class GOTerm:
    """A Gene Ontology term."""

    id: str
    name: str
    namespace: str
    parents: dict[str, set[str]] = field(default_factory=dict)  # edge_type -> set of parent IDs


@dataclass
class GODAG:
    """Gene Ontology Directed Acyclic Graph."""

    terms: dict[str, GOTerm] = field(default_factory=dict)
    obsolete_terms: dict[str, GOTerm] = field(default_factory=dict)
    # Secondary IDs (OBO ``alt_id``) mapped to the primary term that absorbed
    # them.  Annotation files may still use the old ID after two terms merge.
    alt_ids: dict[str, str] = field(default_factory=dict)

    def resolve_term_id(self, term_id: str) -> str | None:
        """Map a GO ID to the term it names today, or None if the DAG has no such term.

        Returns ``term_id`` itself when it is a live term, the primary ID when
        ``term_id`` is a secondary (alt) ID of a live term, ``term_id`` when
        it is a known obsolete term (so it can still be reported as such), and
        None otherwise.
        """
        if term_id in self.terms:
            return term_id
        primary = self.alt_ids.get(term_id)
        if primary is not None and primary in self.terms:
            return primary
        if term_id in self.obsolete_terms:
            return term_id
        return None

    def get_closure(
        self,
        term_id: str,
        edge_types: set[str] | None = None,
        include_self: bool = True,
    ) -> set[str]:
        """Compute the transitive closure (ancestors) of a term.

        Args:
            term_id: The GO term ID to compute closure for
            edge_types: Set of edge types to follow (default: {"is_a"})
            include_self: Whether to include the term itself in the closure

        Returns:
            Set of GO term IDs in the closure
        """
        if edge_types is None:
            edge_types = {"is_a"}

        # A secondary ID resolves to its primary term, which is what the closure
        # (including the self term) is computed for.
        term_id = self.alt_ids.get(term_id, term_id)

        closure: set[str] = set()
        if include_self:
            closure.add(term_id)

        if term_id not in self.terms:
            return closure

        # BFS/DFS to find all ancestors
        visited: set[str] = set()
        stack = [term_id]

        while stack:
            current = stack.pop()
            if current in visited:
                continue
            visited.add(current)

            if current != term_id:
                closure.add(current)

            term = self.terms.get(current)
            if term is None:
                continue

            for edge_type in edge_types:
                parents = term.parents.get(edge_type, set())
                for parent_id in parents:
                    if parent_id not in visited:
                        stack.append(parent_id)

        return closure

    def get_closure_union(
        self,
        term_ids: set[str],
        edge_types: set[str] | None = None,
        include_self: bool = True,
    ) -> set[str]:
        """Compute the union of closures for multiple terms.

        Args:
            term_ids: Set of GO term IDs
            edge_types: Set of edge types to follow
            include_self: Whether to include the terms themselves

        Returns:
            Union of all closures as a set
        """
        result: set[str] = set()
        for term_id in term_ids:
            result |= self.get_closure(term_id, edge_types, include_self)
        return result


def load_go_from_json(file_path: Path | str) -> GODAG:
    """Load GO DAG from a JSON file.

    Expected format:
    {
        "terms": {
            "GO:0000001": {"name": "...", "namespace": "..."},
            ...
        },
        "edges": {
            "is_a": [["child", "parent"], ...],
            "part_of": [["child", "parent"], ...]
        },
        "alt_ids": {"GO:0019952": "GO:0000003", ...},      # optional: secondary -> primary
        "obsolete_terms": {"GO:0000002": {"name": "...", "namespace": "..."}}  # optional
    }

    Args:
        file_path: Path to the JSON file

    Returns:
        GODAG object
    """
    file_path = Path(file_path)

    with open(file_path, encoding="utf-8") as f:
        data = json.load(f)

    return load_go_from_dict(data)


def load_go_from_dict(data: dict[str, Any]) -> GODAG:
    """Load GO DAG from a dictionary.

    Args:
        data: Dictionary with 'terms' and 'edges' keys

    Returns:
        GODAG object
    """
    dag = GODAG()

    # Load terms
    for term_id, term_data in data.get("terms", {}).items():
        dag.terms[term_id] = GOTerm(
            id=term_id,
            name=term_data.get("name", ""),
            namespace=term_data.get("namespace", ""),
            parents={},
        )

    # Load edges
    for edge_type, edges in data.get("edges", {}).items():
        for child_id, parent_id in edges:
            if child_id in dag.terms:
                if edge_type not in dag.terms[child_id].parents:
                    dag.terms[child_id].parents[edge_type] = set()
                dag.terms[child_id].parents[edge_type].add(parent_id)

    # Load secondary IDs
    for alt_id, primary_id in data.get("alt_ids", {}).items():
        dag.alt_ids[alt_id] = primary_id

    # Load obsolete terms
    for term_id, term_data in data.get("obsolete_terms", {}).items():
        dag.obsolete_terms[term_id] = GOTerm(
            id=term_id,
            name=term_data.get("name", ""),
            namespace=term_data.get("namespace", ""),
            parents={},
        )

    return dag


def get_all_parent_ids(dag: GODAG, term_id: str) -> set[str]:
    """Get all direct parent IDs for a term across all edge types.

    Args:
        dag: The GO DAG
        term_id: The term ID

    Returns:
        Set of direct parent IDs
    """
    term = dag.terms.get(term_id)
    if term is None:
        return set()

    parents: set[str] = set()
    for parent_set in term.parents.values():
        parents |= parent_set
    return parents
