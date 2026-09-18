"""Unit tests for the runner's homology stage: DIAMOND hit cap and wiring."""

import logging
from pathlib import Path
from unittest.mock import patch

import pytest

from metagomics2.core.diamond import (
    DEFAULT_DIAMOND_EVALUE,
    DEFAULT_MAX_TARGET_SEQS,
    DiamondResult,
)
from metagomics2.core.filtering import FilterPolicy, HomologyHit
from metagomics2.pipeline.runner import PipelineConfig, PipelineRunner


def _hit(query: str, subject: str, bitscore: float = 100.0) -> HomologyHit:
    return HomologyHit(
        query_id=query, subject_id=subject, evalue=1e-20, bitscore=bitscore,
        pident=90.0, qcov=95.0, alnlen=100,
    )


def _runner(tmp_path: Path, **config_kwargs) -> PipelineRunner:
    subset = tmp_path / "work" / "hit_proteins.fasta"
    subset.parent.mkdir(parents=True)
    subset.write_text(">protA\nACDEFGHIK\n")
    config = PipelineConfig(
        fasta_path=tmp_path / "bg.fasta",
        peptide_list_paths=[],
        output_dir=tmp_path / "results",
        annotated_db_path=tmp_path / "db.dmnd",
        **config_kwargs,
    )
    runner = PipelineRunner(config)
    runner.subset_fasta_path = subset
    return runner


def _diamond_result(hits_by_query: dict[str, list[HomologyHit]]) -> DiamondResult:
    return DiamondResult(
        hits_by_query=hits_by_query,
        output_path=Path("diamond_results.tsv"),
        n_queries=len(hits_by_query),
        n_hits=sum(len(h) for h in hits_by_query.values()),
        command=["diamond", "blastp"],
    )


class TestEffectiveMaxTargetSeqs:
    def _config(self, **kwargs) -> PipelineConfig:
        return PipelineConfig(
            fasta_path=Path("bg.fasta"), peptide_list_paths=[], output_dir=Path("out"), **kwargs
        )

    def test_default(self):
        assert self._config().effective_max_target_seqs() == DEFAULT_MAX_TARGET_SEQS

    def test_raised_to_top_k(self):
        config = self._config(filter_policy=FilterPolicy(top_k=DEFAULT_MAX_TARGET_SEQS + 50))
        assert config.effective_max_target_seqs() == DEFAULT_MAX_TARGET_SEQS + 50

    def test_top_k_below_cap_leaves_cap(self):
        config = self._config(filter_policy=FilterPolicy(top_k=5))
        assert config.effective_max_target_seqs() == DEFAULT_MAX_TARGET_SEQS

    def test_zero_is_unlimited_regardless_of_top_k(self):
        config = self._config(diamond_max_target_seqs=0, filter_policy=FilterPolicy(top_k=9999))
        assert config.effective_max_target_seqs() == 0


class TestEffectiveDiamondEvalue:
    def _config(self, **kwargs) -> PipelineConfig:
        return PipelineConfig(
            fasta_path=Path("bg.fasta"), peptide_list_paths=[], output_dir=Path("out"), **kwargs
        )

    def test_default_when_policy_unset(self):
        assert self._config().effective_diamond_evalue() == DEFAULT_DIAMOND_EVALUE

    def test_policy_value_used(self):
        config = self._config(filter_policy=FilterPolicy(max_evalue=1e-5))
        assert config.effective_diamond_evalue() == 1e-5


class TestHomologySearchCap:
    @patch("metagomics2.pipeline.runner.run_diamond")
    def test_cap_passed_to_diamond(self, mock_run, tmp_path):
        mock_run.return_value = _diamond_result({"protA": [_hit("protA", "S1")]})
        runner = _runner(tmp_path, diamond_max_target_seqs=40, filter_policy=FilterPolicy(top_k=60))

        runner._run_homology_search()

        assert mock_run.call_args.kwargs["max_target_seqs"] == 60
        assert mock_run.call_args.kwargs["evalue"] == DEFAULT_DIAMOND_EVALUE
        assert runner.protein_to_subjects == {"protA": {"S1"}}
        assert runner.diamond_command == "diamond blastp"

    @patch("metagomics2.pipeline.runner.run_diamond")
    def test_unlimited_passes_zero(self, mock_run, tmp_path):
        mock_run.return_value = _diamond_result({})
        runner = _runner(tmp_path, diamond_max_target_seqs=0)

        runner._run_homology_search()

        assert mock_run.call_args.kwargs["max_target_seqs"] == 0

    @patch("metagomics2.pipeline.runner.run_diamond")
    def test_warns_when_a_query_hits_the_cap(self, mock_run, tmp_path, caplog):
        mock_run.return_value = _diamond_result({
            "protA": [_hit("protA", f"S{i}") for i in range(3)],
            "protB": [_hit("protB", "S9")],
        })
        runner = _runner(tmp_path, diamond_max_target_seqs=3)

        with caplog.at_level(logging.WARNING, logger="metagomics2.pipeline.runner"):
            runner._run_homology_search()

        warnings = [r.message for r in caplog.records if r.levelno == logging.WARNING]
        assert len(warnings) == 1
        assert warnings[0].startswith(
            "1 of 2 query proteins returned the DIAMOND per-query hit cap of 3"
        )
        assert "METAGOMICS_DIAMOND_MAX_TARGET_SEQS" in warnings[0]

    @patch("metagomics2.pipeline.runner.run_diamond")
    def test_no_warning_below_cap(self, mock_run, tmp_path, caplog):
        mock_run.return_value = _diamond_result(
            {"protA": [_hit("protA", "S1"), _hit("protA", "S2")]}
        )
        runner = _runner(tmp_path, diamond_max_target_seqs=3)

        with caplog.at_level(logging.WARNING, logger="metagomics2.pipeline.runner"):
            runner._run_homology_search()

        assert not [r for r in caplog.records if r.levelno == logging.WARNING]

    def test_missing_annotated_db_raises(self, tmp_path):
        runner = _runner(tmp_path)
        runner.config.annotated_db_path = None
        with pytest.raises(ValueError, match="No annotated database"):
            runner._run_homology_search()


class TestSilentAnnotationWarning:
    """Peptides with hits but no usable annotation are reported, not hidden."""

    def _runner_with_subjects(self, tmp_path, subject_annotations):
        from metagomics2.core.annotation import SubjectAnnotation
        from metagomics2.core.go import load_go_from_dict
        from metagomics2.core.peptides import Peptide
        from metagomics2.core.taxonomy import TaxonomyTree

        runner = _runner(tmp_path)
        runner.taxonomy_tree = TaxonomyTree()
        runner.go_dag = load_go_from_dict(
            {"terms": {"GO:0000001": {"name": "root", "namespace": "bp"}}, "edges": {}}
        )
        runner.protein_to_subjects = {"B1": {"S1"}}
        runner.subject_annotations = {
            sid: SubjectAnnotation(subject_id=sid, **fields)
            for sid, fields in subject_annotations.items()
        }
        peptides = [Peptide("PEP", 2.0), Peptide("NOHIT", 1.0)]
        peptide_to_proteins = {"PEP": {"B1"}, "NOHIT": set()}
        return runner, peptides, peptide_to_proteins

    def test_warns_when_hits_carry_nothing(self, tmp_path, caplog):
        runner, peptides, p2p = self._runner_with_subjects(tmp_path, {"S1": {}})

        with caplog.at_level(logging.WARNING, logger="metagomics2.pipeline.runner"):
            annotations = runner._annotate_peptides(peptides, p2p)

        assert [a.is_annotated for a in annotations] == [False, False]
        messages = [r.message for r in caplog.records if r.levelno == logging.WARNING]
        assert len(messages) == 1
        assert messages[0].startswith("1 of 2 peptides matched homology hits but received no")

    def test_no_warning_when_hits_annotate(self, tmp_path, caplog):
        runner, peptides, p2p = self._runner_with_subjects(
            tmp_path, {"S1": {"go_terms": {"GO:0000001"}}}
        )

        with caplog.at_level(logging.WARNING, logger="metagomics2.pipeline.runner"):
            annotations = runner._annotate_peptides(peptides, p2p)

        assert annotations[0].is_annotated is True
        assert not [r for r in caplog.records if r.levelno == logging.WARNING]


def test_runner_class_has_docstring():
    """The class attributes must not precede the docstring, or Python drops it."""
    assert PipelineRunner.__doc__ is not None
    assert "Orchestrates" in PipelineRunner.__doc__
