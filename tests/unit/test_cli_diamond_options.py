"""Unit tests for the CLI's DIAMOND tuning options."""

from pathlib import Path
from unittest.mock import MagicMock, patch

from metagomics2.cli import cmd_run, create_parser
from metagomics2.core.diamond import DEFAULT_MAX_TARGET_SEQS


def _base_args(tmp_path: Path) -> list[str]:
    fasta = tmp_path / "bg.fasta"
    fasta.write_text(">p\nACDEFGHIK\n")
    peptides = tmp_path / "peps.tsv"
    peptides.write_text("peptide_sequence\tquantity\nACDEF\t1\n")
    hits = tmp_path / "hits.json"
    hits.write_text("{}")
    ann = tmp_path / "ann.json"
    ann.write_text("{}")
    return [
        "run",
        "--fasta", str(fasta),
        "--peptides", str(peptides),
        "--outdir", str(tmp_path / "out"),
        "--mock-hits", str(hits),
        "--mock-annotations", str(ann),
    ]


class TestParser:
    def test_defaults_are_none(self, tmp_path: Path):
        args = create_parser().parse_args(_base_args(tmp_path))
        assert args.diamond_block_size is None
        assert args.diamond_index_chunks is None
        assert args.diamond_tmpdir is None
        assert args.diamond_max_target_seqs is None

    def test_values_parsed(self, tmp_path: Path):
        args = create_parser().parse_args(
            _base_args(tmp_path)
            + ["--diamond-block-size", "8", "--diamond-index-chunks", "1",
               "--diamond-tmpdir", "/dev/shm", "--diamond-max-target-seqs", "0"]
        )
        assert args.diamond_block_size == 8.0
        assert args.diamond_index_chunks == 1
        assert args.diamond_tmpdir == "/dev/shm"
        assert args.diamond_max_target_seqs == 0

    def test_help_mentions_memory(self, capsys):
        try:
            create_parser().parse_args(["run", "--help"])
        except SystemExit:
            pass
        out = capsys.readouterr().out
        assert "--diamond-block-size" in out
        assert "6 GB" in out


class TestCmdRun:
    def test_options_reach_pipeline_config(self, tmp_path: Path):
        args = create_parser().parse_args(
            _base_args(tmp_path)
            + ["--diamond-block-size", "8", "--diamond-index-chunks", "1",
               "--diamond-tmpdir", str(tmp_path / "shm")]
        )
        result = MagicMock(success=True)
        with patch("metagomics2.cli.run_pipeline", return_value=result) as mock_run:
            assert cmd_run(args) == 0
        config = mock_run.call_args[0][0]
        assert config.diamond_block_size == 8.0
        assert config.diamond_index_chunks == 1
        assert config.diamond_tmpdir == tmp_path / "shm"
        assert config.diamond_max_target_seqs == DEFAULT_MAX_TARGET_SEQS  # not given

    def test_max_target_seqs_reaches_pipeline_config(self, tmp_path: Path):
        args = create_parser().parse_args(
            _base_args(tmp_path) + ["--diamond-max-target-seqs", "2000"]
        )
        result = MagicMock(success=True)
        with patch("metagomics2.cli.run_pipeline", return_value=result) as mock_run:
            assert cmd_run(args) == 0
        assert mock_run.call_args[0][0].diamond_max_target_seqs == 2000

    def test_rejects_negative_max_target_seqs(self, tmp_path: Path, capsys):
        args = create_parser().parse_args(
            _base_args(tmp_path) + ["--diamond-max-target-seqs", "-1"]
        )
        with patch("metagomics2.cli.run_pipeline") as mock_run:
            assert cmd_run(args) == 1
        mock_run.assert_not_called()
        assert "--diamond-max-target-seqs" in capsys.readouterr().err

    def test_rejects_non_positive_block_size(self, tmp_path: Path, capsys):
        args = create_parser().parse_args(_base_args(tmp_path) + ["--diamond-block-size", "0"])
        with patch("metagomics2.cli.run_pipeline") as mock_run:
            assert cmd_run(args) == 1
        mock_run.assert_not_called()
        assert "--diamond-block-size" in capsys.readouterr().err

    def test_rejects_zero_index_chunks(self, tmp_path: Path, capsys):
        args = create_parser().parse_args(
            _base_args(tmp_path) + ["--diamond-index-chunks", "0"]
        )
        with patch("metagomics2.cli.run_pipeline") as mock_run:
            assert cmd_run(args) == 1
        mock_run.assert_not_called()
        assert "--diamond-index-chunks" in capsys.readouterr().err


class TestGoEdgeTypesOption:
    def test_spaces_tolerated(self, tmp_path: Path):
        args = create_parser().parse_args(
            _base_args(tmp_path) + ["--go-edge-types", "is_a, part_of"]
        )
        result = MagicMock(success=True)
        with patch("metagomics2.cli.run_pipeline", return_value=result) as mock_run:
            assert cmd_run(args) == 0
        assert mock_run.call_args[0][0].go_edge_types == {"is_a", "part_of"}

    def test_unknown_type_rejected(self, tmp_path: Path, capsys):
        args = create_parser().parse_args(
            _base_args(tmp_path) + ["--go-edge-types", "is_a,partof"]
        )
        with patch("metagomics2.cli.run_pipeline") as mock_run:
            assert cmd_run(args) == 1
        mock_run.assert_not_called()
        err = capsys.readouterr().err
        assert "--go-edge-types" in err
        assert "partof" in err
