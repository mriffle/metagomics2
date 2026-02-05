"""Parser for NCBI taxonomy dump files."""

from pathlib import Path

from metagomics2.core.taxonomy import TaxonomyTree, TaxonNode


class NCBIParsingError(Exception):
    """Raised when NCBI taxonomy parsing fails."""

    pass


def parse_ncbi_taxonomy_dump(dump_dir: Path | str) -> TaxonomyTree:
    """Parse NCBI taxonomy dump files into a TaxonomyTree.

    Expected files in dump_dir:
    - nodes.dmp: taxonomy nodes with parent relationships
    - names.dmp: taxonomy names

    Args:
        dump_dir: Directory containing NCBI taxonomy dump files

    Returns:
        TaxonomyTree object

    Raises:
        NCBIParsingError: If parsing fails
    """
    dump_dir = Path(dump_dir)

    if not dump_dir.exists():
        raise NCBIParsingError(f"Directory not found: {dump_dir}")

    nodes_file = dump_dir / "nodes.dmp"
    names_file = dump_dir / "names.dmp"

    if not nodes_file.exists():
        raise NCBIParsingError(f"nodes.dmp not found in {dump_dir}")
    if not names_file.exists():
        raise NCBIParsingError(f"names.dmp not found in {dump_dir}")

    # Parse nodes.dmp for structure and ranks
    nodes_data = _parse_nodes_dmp(nodes_file)

    # Parse names.dmp for scientific names
    names_data = _parse_names_dmp(names_file)

    # Build taxonomy tree
    tree = TaxonomyTree()

    for tax_id, node_info in nodes_data.items():
        parent_tax_id = node_info["parent_tax_id"]
        rank = node_info["rank"]
        name = names_data.get(tax_id, f"taxid_{tax_id}")

        # Root node has itself as parent
        if tax_id == parent_tax_id:
            parent_tax_id = None

        tree.nodes[tax_id] = TaxonNode(
            tax_id=tax_id,
            name=name,
            rank=rank,
            parent_tax_id=parent_tax_id,
        )

    return tree


def _parse_nodes_dmp(file_path: Path) -> dict[int, dict]:
    """Parse nodes.dmp file.

    Format: tax_id | parent_tax_id | rank | ... (pipe-delimited)

    Args:
        file_path: Path to nodes.dmp

    Returns:
        Dictionary mapping tax_id to node info
    """
    nodes = {}

    with open(file_path, "r", encoding="utf-8") as f:
        for line in f:
            # Split on pipe and strip whitespace
            parts = [p.strip() for p in line.split("|")]

            if len(parts) < 3:
                continue

            tax_id = int(parts[0])
            parent_tax_id = int(parts[1])
            rank = parts[2]

            nodes[tax_id] = {
                "parent_tax_id": parent_tax_id,
                "rank": rank,
            }

    return nodes


def _parse_names_dmp(file_path: Path) -> dict[int, str]:
    """Parse names.dmp file to get scientific names.

    Format: tax_id | name | unique_name | name_class | (pipe-delimited)

    We only keep entries where name_class is "scientific name".

    Args:
        file_path: Path to names.dmp

    Returns:
        Dictionary mapping tax_id to scientific name
    """
    names = {}

    with open(file_path, "r", encoding="utf-8") as f:
        for line in f:
            # Split on pipe and strip whitespace
            parts = [p.strip() for p in line.split("|")]

            if len(parts) < 4:
                continue

            tax_id = int(parts[0])
            name = parts[1]
            name_class = parts[3]

            # Only use scientific names
            if name_class == "scientific name":
                names[tax_id] = name

    return names


def convert_ncbi_dump_to_json_dict(dump_dir: Path | str) -> dict:
    """Convert NCBI taxonomy dump to JSON dictionary format.

    Args:
        dump_dir: Directory containing NCBI taxonomy dump files

    Returns:
        Dictionary in the format expected by load_taxonomy_from_dict
    """
    tree = parse_ncbi_taxonomy_dump(dump_dir)

    # Build the JSON structure
    result = {"nodes": {}}

    for tax_id, node in tree.nodes.items():
        result["nodes"][str(tax_id)] = {
            "name": node.name,
            "rank": node.rank,
            "parent_tax_id": node.parent_tax_id,
        }

    return result
