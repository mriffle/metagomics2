"""Integration tests for GO OBO parsing using real data from Gene Ontology.

Downloads the same GO OBO file used in the Dockerfile and verifies parsing
produces a structurally valid GODAG. Files are cached in tests/.ref_cache/
to avoid re-downloading on every test run.

These tests are marked 'slow' since the initial download is ~35 MB.
"""

import shutil
import urllib.request
from pathlib import Path

import pytest

from metagomics2.core.go import GODAG
from metagomics2.core.obo_parser import parse_obo_file

# Same URL and version as Dockerfile
GO_VERSION = "2024-01-17"
GO_OBO_URL = f"http://purl.obolibrary.org/obo/go/releases/{GO_VERSION}/go.obo"

CACHE_DIR = Path(__file__).parent.parent / ".ref_cache"

USER_AGENT = "metagomics2-test/1.0"


def _download(url: str, dest: Path) -> None:
    """Download a URL to a local file, setting a User-Agent header."""
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req) as resp, open(dest, "wb") as out:
        shutil.copyfileobj(resp, out)


@pytest.fixture(scope="module")
def go_obo_path() -> Path:
    """Download and cache the real GO OBO file.

    Uses the same URL as the Dockerfile:
      http://purl.obolibrary.org/obo/go/releases/2024-01-17/go.obo
    """
    cache_path = CACHE_DIR / "go" / f"go-{GO_VERSION}.obo"

    if not cache_path.exists():
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        print(f"\nDownloading GO OBO from {GO_OBO_URL} ...")
        _download(GO_OBO_URL, cache_path)
        print(f"  Saved to {cache_path} ({cache_path.stat().st_size / 1e6:.1f} MB)")

    return cache_path


@pytest.fixture(scope="module")
def go_dag(go_obo_path: Path) -> GODAG:
    """Parse the real GO OBO file into a GODAG (cached per module)."""
    return parse_obo_file(go_obo_path)


@pytest.mark.slow
class TestRealGOOBOParsing:
    """Tests that parse the real GO OBO release and check structural invariants."""

    def test_parses_without_error(self, go_dag: GODAG):
        """Parsing should complete without raising."""
        assert go_dag is not None

    def test_has_many_terms(self, go_dag: GODAG):
        """Real GO has tens of thousands of terms."""
        assert len(go_dag.terms) > 40_000, (
            f"Expected >40k GO terms, got {len(go_dag.terms)}"
        )

    def test_all_terms_have_go_id_format(self, go_dag: GODAG):
        """Every term ID should match GO:NNNNNNN format."""
        for term_id in go_dag.terms:
            assert term_id.startswith("GO:"), f"Bad ID prefix: {term_id}"
            numeric_part = term_id[3:]
            assert numeric_part.isdigit(), f"Non-numeric suffix: {term_id}"
            assert len(numeric_part) == 7, f"Wrong ID length: {term_id}"

    def test_all_terms_have_names(self, go_dag: GODAG):
        """Every non-obsolete term should have a non-empty name."""
        nameless = [tid for tid, t in go_dag.terms.items() if not t.name]
        assert len(nameless) == 0, (
            f"{len(nameless)} terms have no name, e.g. {nameless[:5]}"
        )

    def test_all_terms_have_namespace(self, go_dag: GODAG):
        """Every term should belong to one of the three GO namespaces."""
        valid_namespaces = {
            "biological_process",
            "molecular_function",
            "cellular_component",
        }
        bad = [
            (tid, t.namespace)
            for tid, t in go_dag.terms.items()
            if t.namespace not in valid_namespaces
        ]
        assert len(bad) == 0, (
            f"{len(bad)} terms have invalid namespace, e.g. {bad[:5]}"
        )

    def test_three_root_terms_exist(self, go_dag: GODAG):
        """The three root terms should be present."""
        # GO:0008150 = biological_process
        # GO:0003674 = molecular_function
        # GO:0005575 = cellular_component
        roots = {
            "GO:0008150": "biological_process",
            "GO:0003674": "molecular_function",
            "GO:0005575": "cellular_component",
        }
        for root_id, expected_ns in roots.items():
            assert root_id in go_dag.terms, f"Root term {root_id} missing"
            assert go_dag.terms[root_id].namespace == expected_ns

    def test_root_terms_have_no_is_a_parents(self, go_dag: GODAG):
        """Root terms should have no is_a parents."""
        roots = ["GO:0008150", "GO:0003674", "GO:0005575"]
        for root_id in roots:
            term = go_dag.terms[root_id]
            is_a_parents = term.parents.get("is_a", set())
            assert len(is_a_parents) == 0, (
                f"Root {root_id} has is_a parents: {is_a_parents}"
            )

    def test_most_terms_have_at_least_one_parent(self, go_dag: GODAG):
        """Non-root terms should have at least one is_a parent."""
        roots = {"GO:0008150", "GO:0003674", "GO:0005575"}
        orphans = []
        for term_id, term in go_dag.terms.items():
            if term_id in roots:
                continue
            all_parents = set()
            for parent_set in term.parents.values():
                all_parents.update(parent_set)
            if len(all_parents) == 0:
                orphans.append(term_id)

        # Allow a tiny number of edge cases but not many
        assert len(orphans) < 10, (
            f"{len(orphans)} non-root orphan terms: {orphans[:10]}"
        )

    def test_is_a_parents_reference_existing_terms(self, go_dag: GODAG):
        """All is_a parent IDs should reference terms that exist in the DAG."""
        dangling = []
        for term_id, term in go_dag.terms.items():
            for parent_id in term.parents.get("is_a", set()):
                if parent_id not in go_dag.terms:
                    dangling.append((term_id, parent_id))

        assert len(dangling) == 0, (
            f"{len(dangling)} dangling is_a refs, e.g. {dangling[:5]}"
        )

    def test_no_self_loops_in_is_a(self, go_dag: GODAG):
        """No term should be its own is_a parent."""
        self_loops = [
            tid for tid, t in go_dag.terms.items()
            if tid in t.parents.get("is_a", set())
        ]
        assert len(self_loops) == 0, f"Self-loops: {self_loops[:5]}"

    def test_closure_reaches_root(self, go_dag: GODAG):
        """For a sample of terms, the transitive closure should reach a root."""
        roots = {"GO:0008150", "GO:0003674", "GO:0005575"}
        # Test a sample of terms
        sample = list(go_dag.terms.keys())[:200]
        for term_id in sample:
            closure = go_dag.get_closure(term_id, edge_types={"is_a"}, include_self=True)
            hits_root = closure & roots
            assert len(hits_root) > 0, (
                f"Term {term_id} closure doesn't reach any root. "
                f"Closure size: {len(closure)}"
            )

    def test_closure_is_transitive(self, go_dag: GODAG):
        """For a sample: if B in closure(A), then closure(B) ⊆ closure(A)."""
        sample = list(go_dag.terms.keys())[:50]
        for term_id in sample:
            closure_a = go_dag.get_closure(term_id, edge_types={"is_a"})
            for ancestor_id in list(closure_a)[:10]:
                if ancestor_id in go_dag.terms:
                    closure_b = go_dag.get_closure(ancestor_id, edge_types={"is_a"})
                    assert closure_b <= closure_a, (
                        f"Transitivity violated: closure({ancestor_id}) "
                        f"not subset of closure({term_id})"
                    )

    def test_well_known_terms_present(self, go_dag: GODAG):
        """Some well-known GO terms should be present."""
        well_known = {
            "GO:0006915": "apoptotic process",
            "GO:0007049": "cell cycle",
            "GO:0006412": "translation",
            "GO:0005634": "nucleus",
            "GO:0005524": "ATP binding",
        }
        for go_id, expected_name in well_known.items():
            assert go_id in go_dag.terms, f"Well-known term {go_id} missing"
            assert go_dag.terms[go_id].name == expected_name, (
                f"Term {go_id} name mismatch: "
                f"expected '{expected_name}', got '{go_dag.terms[go_id].name}'"
            )

    def test_no_obsolete_terms(self, go_dag: GODAG):
        """Parser should have filtered out obsolete terms."""
        # Obsolete terms typically have "obsolete" in their name
        # but the parser filters by is_obsolete tag, not name.
        # We just verify the count is reasonable (real GO has ~50k active terms).
        assert len(go_dag.terms) < 60_000, (
            f"Suspiciously many terms ({len(go_dag.terms)}), "
            "obsolete terms may not be filtered"
        )

    def test_namespace_distribution(self, go_dag: GODAG):
        """All three namespaces should have substantial representation."""
        ns_counts = {"biological_process": 0, "molecular_function": 0, "cellular_component": 0}
        for term in go_dag.terms.values():
            if term.namespace in ns_counts:
                ns_counts[term.namespace] += 1

        for ns, count in ns_counts.items():
            assert count > 1000, (
                f"Namespace '{ns}' has only {count} terms, expected >1000"
            )

    def test_edge_types_present(self, go_dag: GODAG):
        """DAG should contain is_a edges and likely part_of edges."""
        edge_types = set()
        for term in go_dag.terms.values():
            edge_types.update(term.parents.keys())

        assert "is_a" in edge_types, "No is_a edges found"
        # part_of is common but stored under 'relationship' parsing
        # Just verify is_a is the dominant edge type
        is_a_count = sum(
            len(t.parents.get("is_a", set())) for t in go_dag.terms.values()
        )
        assert is_a_count > 50_000, f"Only {is_a_count} is_a edges, expected >50k"
