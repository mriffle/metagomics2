"""Parser for OBO (Open Biomedical Ontologies) format files."""

from pathlib import Path
from typing import TextIO

from metagomics2.core.go import GODAG, GOTerm


class OBOParsingError(Exception):
    """Raised when OBO parsing fails."""

    pass


def parse_obo_file(file_path: Path | str) -> GODAG:
    """Parse a GO OBO file into a GODAG.

    Args:
        file_path: Path to the OBO file

    Returns:
        GODAG object

    Raises:
        OBOParsingError: If parsing fails
    """
    file_path = Path(file_path)

    if not file_path.exists():
        raise OBOParsingError(f"File not found: {file_path}")

    with open(file_path, "r", encoding="utf-8") as f:
        return parse_obo_from_handle(f)


def parse_obo_from_handle(handle: TextIO) -> GODAG:
    """Parse OBO format from a file handle.

    Args:
        handle: File handle to read from

    Returns:
        GODAG object

    Raises:
        OBOParsingError: If parsing fails
    """
    dag = GODAG()
    current_stanza: dict[str, list[str]] = {}
    in_term = False

    for line_num, line in enumerate(handle, 1):
        line = line.strip()

        # Skip empty lines and comments
        if not line or line.startswith("!"):
            continue

        # Start of a new stanza
        if line == "[Term]":
            # Process previous term if exists
            if current_stanza:
                _process_term_stanza(current_stanza, dag)
            current_stanza = {}
            in_term = True
            continue

        # Other stanza types (Typedef, etc.) - skip for now
        if line.startswith("[") and line.endswith("]"):
            if current_stanza:
                _process_term_stanza(current_stanza, dag)
            current_stanza = {}
            in_term = False
            continue

        # Parse tag-value pairs
        if in_term and ":" in line:
            # Split on first colon
            tag, _, value = line.partition(":")
            tag = tag.strip()
            value = value.strip()

            # Remove trailing comments and modifiers
            if "!" in value:
                value = value.split("!")[0].strip()

            # Handle escaped characters
            value = value.replace("\\n", "\n").replace("\\t", "\t")

            if tag not in current_stanza:
                current_stanza[tag] = []
            current_stanza[tag].append(value)

    # Process last term
    if current_stanza:
        _process_term_stanza(current_stanza, dag)

    return dag


def _process_term_stanza(stanza: dict[str, list[str]], dag: GODAG) -> None:
    """Process a [Term] stanza and add to DAG.

    Args:
        stanza: Dictionary of tag -> list of values
        dag: GODAG to add term to
    """
    # Store obsolete terms with metadata but not in the main DAG
    if "is_obsolete" in stanza and stanza["is_obsolete"][0].lower() == "true":
        if "id" in stanza:
            term_id = stanza["id"][0]
            dag.obsolete_terms[term_id] = GOTerm(
                id=term_id,
                name=stanza.get("name", [""])[0],
                namespace=stanza.get("namespace", [""])[0],
                parents={},
            )
        return

    # Extract required fields
    if "id" not in stanza:
        return  # Skip terms without ID

    term_id = stanza["id"][0]

    # Get name
    name = stanza.get("name", [""])[0]

    # Get namespace
    namespace = stanza.get("namespace", [""])[0]

    # Create term
    term = GOTerm(
        id=term_id,
        name=name,
        namespace=namespace,
        parents={},
    )

    # Process relationships
    # is_a relationships
    if "is_a" in stanza:
        term.parents["is_a"] = set()
        for is_a_value in stanza["is_a"]:
            # Extract GO ID (format: "GO:0000001 ! name")
            parent_id = is_a_value.split()[0]
            term.parents["is_a"].add(parent_id)

    # relationship: part_of GO:XXXXXXX
    if "relationship" in stanza:
        for rel_value in stanza["relationship"]:
            parts = rel_value.split()
            if len(parts) >= 2:
                rel_type = parts[0]
                parent_id = parts[1]

                if rel_type not in term.parents:
                    term.parents[rel_type] = set()
                term.parents[rel_type].add(parent_id)

    dag.terms[term_id] = term


def convert_obo_to_json_dict(obo_path: Path | str) -> dict:
    """Convert an OBO file to the JSON dictionary format.

    Args:
        obo_path: Path to OBO file

    Returns:
        Dictionary in the format expected by load_go_from_dict
    """
    dag = parse_obo_file(obo_path)

    # Build the JSON structure
    result = {
        "terms": {},
        "edges": {},
    }

    # Add terms
    for term_id, term in dag.terms.items():
        result["terms"][term_id] = {
            "name": term.name,
            "namespace": term.namespace,
        }

    # Add edges
    all_edge_types = set()
    for term in dag.terms.values():
        all_edge_types.update(term.parents.keys())

    for edge_type in all_edge_types:
        result["edges"][edge_type] = []

    for term_id, term in dag.terms.items():
        for edge_type, parent_ids in term.parents.items():
            for parent_id in parent_ids:
                result["edges"][edge_type].append([term_id, parent_id])

    # Add obsolete terms
    if dag.obsolete_terms:
        result["obsolete_terms"] = {}
        for term_id, term in dag.obsolete_terms.items():
            result["obsolete_terms"][term_id] = {
                "name": term.name,
                "namespace": term.namespace,
            }

    return result
