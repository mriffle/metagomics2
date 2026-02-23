"""Unit tests for write_peptide_mapping_parquet."""

from pathlib import Path

import polars as pl
import pytest

from metagomics2.core.annotation import PeptideAnnotation
from metagomics2.core.reporting import write_peptide_mapping_parquet


def make_annotation(
    peptide: str,
    is_annotated: bool = True,
    lca_tax_id: int | None = None,
    go_terms: set[str] | None = None,
    background_proteins: set[str] | None = None,
) -> PeptideAnnotation:
    """Helper to create a PeptideAnnotation."""
    return PeptideAnnotation(
        peptide=peptide,
        quantity=1.0,
        is_annotated=is_annotated,
        lca_tax_id=lca_tax_id,
        go_terms=go_terms or set(),
        background_proteins=background_proteins or set(),
    )


class TestWritePeptideMappingParquet:
    """Tests for write_peptide_mapping_parquet."""

    def test_single_annotated_peptide_single_subject(self, tmp_path: Path) -> None:
        """Annotated peptide with one background protein and one subject → one row."""
        annotations = [
            make_annotation(
                "PEPTIDE",
                is_annotated=True,
                lca_tax_id=9606,
                go_terms={"GO:0000001"},
                background_proteins={"prot_A"},
            ),
        ]
        peptide_to_proteins = {"PEPTIDE": {"prot_A"}}
        protein_to_subjects = {"prot_A": {"subj_X"}}

        output = tmp_path / "peptide_mapping.parquet"
        write_peptide_mapping_parquet(annotations, peptide_to_proteins, protein_to_subjects, output)

        df = pl.read_parquet(output)
        assert len(df) == 1
        row = df.row(0, named=True)
        assert row["peptide"] == "PEPTIDE"
        assert row["peptide_tax_id"] == 9606
        assert row["peptide_go_terms"] == ["GO:0000001"]
        assert row["background_protein"] == "prot_A"
        assert row["annotated_protein"] == "subj_X"

    def test_annotated_peptide_two_background_proteins(self, tmp_path: Path) -> None:
        """Annotated peptide with two background proteins → two rows."""
        annotations = [
            make_annotation(
                "PEPTIDE",
                is_annotated=True,
                lca_tax_id=9606,
                go_terms={"GO:0000001"},
                background_proteins={"prot_A", "prot_B"},
            ),
        ]
        peptide_to_proteins = {"PEPTIDE": {"prot_A", "prot_B"}}
        protein_to_subjects = {"prot_A": {"subj_X"}, "prot_B": {"subj_Y"}}

        output = tmp_path / "peptide_mapping.parquet"
        write_peptide_mapping_parquet(annotations, peptide_to_proteins, protein_to_subjects, output)

        df = pl.read_parquet(output)
        assert len(df) == 2
        bg_proteins = set(df["background_protein"].to_list())
        assert bg_proteins == {"prot_A", "prot_B"}
        ann_proteins = set(df["annotated_protein"].to_list())
        assert ann_proteins == {"subj_X", "subj_Y"}

    def test_unannotated_peptide_excluded(self, tmp_path: Path) -> None:
        """Unannotated peptides are not written to the Parquet file."""
        annotations = [
            make_annotation("ANNOTATED", is_annotated=True, lca_tax_id=9606, background_proteins={"prot_A"}),
            make_annotation("UNANNOTATED", is_annotated=False),
        ]
        peptide_to_proteins = {"ANNOTATED": {"prot_A"}, "UNANNOTATED": set()}
        protein_to_subjects = {"prot_A": {"subj_X"}}

        output = tmp_path / "peptide_mapping.parquet"
        write_peptide_mapping_parquet(annotations, peptide_to_proteins, protein_to_subjects, output)

        df = pl.read_parquet(output)
        assert len(df) == 1
        assert df["peptide"][0] == "ANNOTATED"

    def test_empty_annotations_produces_schema_valid_file(self, tmp_path: Path) -> None:
        """Empty annotations → schema-valid empty Parquet file."""
        output = tmp_path / "peptide_mapping.parquet"
        write_peptide_mapping_parquet([], {}, {}, output)

        assert output.exists()
        df = pl.read_parquet(output)
        assert len(df) == 0

    def test_schema_column_names(self, tmp_path: Path) -> None:
        """Output file has exactly the expected column names."""
        output = tmp_path / "peptide_mapping.parquet"
        write_peptide_mapping_parquet([], {}, {}, output)

        df = pl.read_parquet(output)
        assert df.columns == [
            "peptide",
            "peptide_tax_id",
            "peptide_go_terms",
            "background_protein",
            "annotated_protein",
        ]

    def test_schema_column_dtypes(self, tmp_path: Path) -> None:
        """Output file has the expected column dtypes."""
        output = tmp_path / "peptide_mapping.parquet"
        write_peptide_mapping_parquet([], {}, {}, output)

        df = pl.read_parquet(output)
        assert df.dtypes[0] == pl.Utf8
        assert df.dtypes[1] == pl.Int64
        assert df.dtypes[2] == pl.List(pl.Utf8)
        assert df.dtypes[3] == pl.Utf8
        assert df.dtypes[4] == pl.Utf8

    def test_nullable_tax_id(self, tmp_path: Path) -> None:
        """Peptide with no LCA tax_id writes None for peptide_tax_id."""
        annotations = [
            make_annotation(
                "PEPTIDE",
                is_annotated=True,
                lca_tax_id=None,
                go_terms={"GO:0000001"},
                background_proteins={"prot_A"},
            ),
        ]
        peptide_to_proteins = {"PEPTIDE": {"prot_A"}}
        protein_to_subjects = {"prot_A": {"subj_X"}}

        output = tmp_path / "peptide_mapping.parquet"
        write_peptide_mapping_parquet(annotations, peptide_to_proteins, protein_to_subjects, output)

        df = pl.read_parquet(output)
        assert len(df) == 1
        assert df["peptide_tax_id"][0] is None

    def test_multiple_subjects_per_background_protein(self, tmp_path: Path) -> None:
        """One background protein with two subjects → two rows."""
        annotations = [
            make_annotation(
                "PEPTIDE",
                is_annotated=True,
                lca_tax_id=9606,
                background_proteins={"prot_A"},
            ),
        ]
        peptide_to_proteins = {"PEPTIDE": {"prot_A"}}
        protein_to_subjects = {"prot_A": {"subj_X", "subj_Y"}}

        output = tmp_path / "peptide_mapping.parquet"
        write_peptide_mapping_parquet(annotations, peptide_to_proteins, protein_to_subjects, output)

        df = pl.read_parquet(output)
        assert len(df) == 2
        assert set(df["annotated_protein"].to_list()) == {"subj_X", "subj_Y"}

    def test_go_terms_sorted(self, tmp_path: Path) -> None:
        """GO terms in each row are sorted."""
        annotations = [
            make_annotation(
                "PEPTIDE",
                is_annotated=True,
                go_terms={"GO:0000003", "GO:0000001", "GO:0000002"},
                background_proteins={"prot_A"},
            ),
        ]
        peptide_to_proteins = {"PEPTIDE": {"prot_A"}}
        protein_to_subjects = {"prot_A": {"subj_X"}}

        output = tmp_path / "peptide_mapping.parquet"
        write_peptide_mapping_parquet(annotations, peptide_to_proteins, protein_to_subjects, output)

        df = pl.read_parquet(output)
        go_terms = df["peptide_go_terms"][0].to_list()
        assert go_terms == ["GO:0000001", "GO:0000002", "GO:0000003"]

    def test_all_unannotated_produces_empty_file(self, tmp_path: Path) -> None:
        """All unannotated peptides → empty but schema-valid Parquet file."""
        annotations = [
            make_annotation("P1", is_annotated=False),
            make_annotation("P2", is_annotated=False),
        ]
        output = tmp_path / "peptide_mapping.parquet"
        write_peptide_mapping_parquet(annotations, {}, {}, output)

        df = pl.read_parquet(output)
        assert len(df) == 0
        assert df.columns == [
            "peptide",
            "peptide_tax_id",
            "peptide_go_terms",
            "background_protein",
            "annotated_protein",
        ]
