"""Integration tests for NCBI taxonomy parsing using real data.

Downloads the same NCBI taxonomy dump used in the Dockerfile and verifies
parsing produces a structurally valid TaxonomyTree. Files are cached in
tests/.ref_cache/ to avoid re-downloading on every test run.

These tests are marked 'slow' since the initial download is ~60 MB.
"""

import shutil
import tarfile
import urllib.request
from pathlib import Path

import pytest

from metagomics2.core.ncbi_parser import parse_ncbi_taxonomy_dump
from metagomics2.core.taxonomy import TaxonomyTree

# Same URL as Dockerfile
NCBI_TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"

CACHE_DIR = Path(__file__).parent.parent / ".ref_cache"

USER_AGENT = "metagomics2-test/1.0"


def _download(url: str, dest: Path) -> None:
    """Download a URL to a local file, setting a User-Agent header."""
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req) as resp, open(dest, "wb") as out:
        shutil.copyfileobj(resp, out)


@pytest.fixture(scope="module")
def ncbi_dump_dir() -> Path:
    """Download and cache the real NCBI taxonomy dump.

    Uses the same URL as the Dockerfile:
      https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
    """
    cache_path = CACHE_DIR / "taxonomy"
    nodes_path = cache_path / "nodes.dmp"
    names_path = cache_path / "names.dmp"

    if not nodes_path.exists() or not names_path.exists():
        cache_path.mkdir(parents=True, exist_ok=True)
        tarball = cache_path / "taxdump.tar.gz"

        print(f"\nDownloading NCBI taxonomy dump from {NCBI_TAXDUMP_URL} ...")
        _download(NCBI_TAXDUMP_URL, tarball)
        print(f"  Downloaded {tarball.stat().st_size / 1e6:.1f} MB, extracting...")

        with tarfile.open(tarball, "r:gz") as tar:
            tar.extractall(path=cache_path, filter="data")

        tarball.unlink()  # Remove tarball after extraction
        print(f"  Extracted to {cache_path}")

    return cache_path


@pytest.fixture(scope="module")
def taxonomy_tree(ncbi_dump_dir: Path) -> TaxonomyTree:
    """Parse the real NCBI taxonomy dump into a TaxonomyTree (cached per module)."""
    return parse_ncbi_taxonomy_dump(ncbi_dump_dir)


@pytest.mark.slow
class TestRealNCBITaxonomyParsing:
    """Tests that parse the real NCBI taxonomy dump and check structural invariants."""

    def test_parses_without_error(self, taxonomy_tree: TaxonomyTree):
        """Parsing should complete without raising."""
        assert taxonomy_tree is not None

    def test_has_many_nodes(self, taxonomy_tree: TaxonomyTree):
        """Real NCBI taxonomy has millions of nodes."""
        assert len(taxonomy_tree.nodes) > 2_000_000, (
            f"Expected >2M taxonomy nodes, got {len(taxonomy_tree.nodes)}"
        )

    def test_root_node_exists(self, taxonomy_tree: TaxonomyTree):
        """Tax ID 1 (root) must exist."""
        assert 1 in taxonomy_tree.nodes
        root = taxonomy_tree.nodes[1]
        assert root.name == "root"
        assert root.rank == "no rank"
        assert root.parent_tax_id is None

    def test_root_is_only_parentless_node(self, taxonomy_tree: TaxonomyTree):
        """Only the root node should have parent_tax_id == None."""
        parentless = [
            tid for tid, node in taxonomy_tree.nodes.items()
            if node.parent_tax_id is None
        ]
        assert parentless == [1], (
            f"Expected only root (1) to be parentless, got {parentless[:10]}"
        )

    def test_all_parents_exist(self, taxonomy_tree: TaxonomyTree):
        """Every node's parent_tax_id should reference an existing node."""
        dangling = []
        for tid, node in taxonomy_tree.nodes.items():
            if node.parent_tax_id is not None and node.parent_tax_id not in taxonomy_tree.nodes:
                dangling.append((tid, node.parent_tax_id))

        assert len(dangling) == 0, (
            f"{len(dangling)} nodes have dangling parent refs, e.g. {dangling[:5]}"
        )

    def test_no_self_parent_except_root(self, taxonomy_tree: TaxonomyTree):
        """No node should be its own parent (root's self-ref is resolved to None)."""
        self_parents = [
            tid for tid, node in taxonomy_tree.nodes.items()
            if node.parent_tax_id == tid
        ]
        assert len(self_parents) == 0, (
            f"Self-parent nodes: {self_parents[:10]}"
        )

    def test_all_nodes_have_names(self, taxonomy_tree: TaxonomyTree):
        """Every node should have a non-empty name."""
        nameless = [
            tid for tid, node in taxonomy_tree.nodes.items()
            if not node.name
        ]
        assert len(nameless) == 0, (
            f"{len(nameless)} nodes have no name, e.g. {nameless[:10]}"
        )

    def test_all_nodes_have_ranks(self, taxonomy_tree: TaxonomyTree):
        """Every node should have a non-empty rank."""
        rankless = [
            tid for tid, node in taxonomy_tree.nodes.items()
            if not node.rank
        ]
        assert len(rankless) == 0, (
            f"{len(rankless)} nodes have no rank, e.g. {rankless[:10]}"
        )

    def test_well_known_taxa_present(self, taxonomy_tree: TaxonomyTree):
        """Some well-known taxa should be present with correct names."""
        well_known = {
            9606: ("Homo sapiens", "species"),
            9913: ("Bos taurus", "species"),
            562: ("Escherichia coli", "species"),
            2: ("Bacteria", "domain"),
            2759: ("Eukaryota", "domain"),
            2157: ("Archaea", "domain"),
            10239: ("Viruses", "acellular root"),
        }
        for tax_id, (expected_name, expected_rank) in well_known.items():
            assert tax_id in taxonomy_tree.nodes, (
                f"Well-known taxon {tax_id} ({expected_name}) missing"
            )
            node = taxonomy_tree.nodes[tax_id]
            assert node.name == expected_name, (
                f"Tax {tax_id} name mismatch: "
                f"expected '{expected_name}', got '{node.name}'"
            )
            assert node.rank == expected_rank, (
                f"Tax {tax_id} rank mismatch: "
                f"expected '{expected_rank}', got '{node.rank}'"
            )

    def test_lineage_to_root(self, taxonomy_tree: TaxonomyTree):
        """Lineage of well-known species should reach root."""
        # Homo sapiens -> ... -> Eukaryota -> root
        lineage = taxonomy_tree.get_lineage(9606)
        assert lineage[0] == 9606, "Lineage should start with the query taxon"
        assert lineage[-1] == 1, "Lineage should end at root"
        assert len(lineage) > 5, (
            f"Homo sapiens lineage suspiciously short: {len(lineage)} nodes"
        )

        # Eukaryota should be in the lineage
        lineage_set = set(lineage)
        assert 2759 in lineage_set, "Eukaryota should be in Homo sapiens lineage"

    def test_lineage_no_cycles(self, taxonomy_tree: TaxonomyTree):
        """Lineages for a sample of nodes should have no cycles (no repeated IDs)."""
        sample = list(taxonomy_tree.nodes.keys())[:500]
        for tax_id in sample:
            lineage = taxonomy_tree.get_lineage(tax_id)
            assert len(lineage) == len(set(lineage)), (
                f"Cycle detected in lineage of {tax_id}"
            )

    def test_lca_of_human_and_cow(self, taxonomy_tree: TaxonomyTree):
        """LCA of Homo sapiens (9606) and Bos taurus (9913) should be a mammalian ancestor."""
        lca = taxonomy_tree.compute_lca({9606, 9913})
        assert lca is not None

        # LCA should be in both lineages
        human_lineage = taxonomy_tree.get_lineage_set(9606)
        cow_lineage = taxonomy_tree.get_lineage_set(9913)
        assert lca in human_lineage
        assert lca in cow_lineage

        # LCA should be above species level
        lca_node = taxonomy_tree.nodes[lca]
        assert lca_node.rank != "species", (
            f"LCA of human and cow should not be a species, got {lca_node.name}"
        )

    def test_lca_of_human_and_ecoli(self, taxonomy_tree: TaxonomyTree):
        """LCA of Homo sapiens and E. coli should be root or cellular organisms."""
        lca = taxonomy_tree.compute_lca({9606, 562})
        assert lca is not None

        # Should be a very high-level ancestor
        lca_lineage = taxonomy_tree.get_lineage(lca)
        # LCA lineage to root should be very short (1-3 nodes)
        assert len(lca_lineage) <= 3, (
            f"LCA of human and E. coli is too deep: "
            f"{taxonomy_tree.nodes[lca].name} (lineage length {len(lca_lineage)})"
        )

    def test_major_ranks_present(self, taxonomy_tree: TaxonomyTree):
        """Standard taxonomic ranks should appear in the tree."""
        ranks = set(node.rank for node in taxonomy_tree.nodes.values())
        expected_ranks = {
            "domain", "phylum", "class", "order",
            "family", "genus", "species",
        }
        for rank in expected_ranks:
            assert rank in ranks, f"Expected rank '{rank}' not found"

    def test_species_count(self, taxonomy_tree: TaxonomyTree):
        """There should be a large number of species-rank nodes."""
        species_count = sum(
            1 for node in taxonomy_tree.nodes.values()
            if node.rank == "species"
        )
        assert species_count > 1_000_000, (
            f"Expected >1M species, got {species_count}"
        )

    def test_domain_count(self, taxonomy_tree: TaxonomyTree):
        """There should be a small number of domain-rank nodes."""
        domains = [
            (tid, node.name)
            for tid, node in taxonomy_tree.nodes.items()
            if node.rank == "domain"
        ]
        assert 2 <= len(domains) <= 10, (
            f"Expected 2-10 domains, got {len(domains)}: {domains}"
        )

    def test_bacteria_subtree_is_large(self, taxonomy_tree: TaxonomyTree):
        """Bacteria (tax_id=2) should have many descendants."""
        # Count direct children of Bacteria
        bacteria_children = [
            tid for tid, node in taxonomy_tree.nodes.items()
            if node.parent_tax_id == 2
        ]
        assert len(bacteria_children) > 10, (
            f"Bacteria has only {len(bacteria_children)} direct children"
        )
