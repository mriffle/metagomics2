"""Tests for input sanitization of JobParams fields."""

import math

import pytest
from pydantic import ValidationError

from metagomics2.models.job import JobParams


class TestMaxEvalue:
    """Validation of max_evalue."""

    def test_none_allowed(self):
        assert JobParams(max_evalue=None).max_evalue is None

    def test_valid_small(self):
        assert JobParams(max_evalue=1e-10).max_evalue == 1e-10

    def test_valid_boundary_1000(self):
        assert JobParams(max_evalue=1000).max_evalue == 1000

    def test_reject_zero(self):
        with pytest.raises(ValidationError):
            JobParams(max_evalue=0)

    def test_reject_negative(self):
        with pytest.raises(ValidationError):
            JobParams(max_evalue=-1e-5)

    def test_reject_nan(self):
        with pytest.raises(ValidationError):
            JobParams(max_evalue=float("nan"))

    def test_reject_inf(self):
        with pytest.raises(ValidationError):
            JobParams(max_evalue=float("inf"))

    def test_reject_neg_inf(self):
        with pytest.raises(ValidationError):
            JobParams(max_evalue=float("-inf"))

    def test_reject_over_1000(self):
        with pytest.raises(ValidationError):
            JobParams(max_evalue=1001)


class TestMinPident:
    """Validation of min_pident."""

    def test_none_allowed(self):
        assert JobParams(min_pident=None).min_pident is None

    def test_valid_zero(self):
        assert JobParams(min_pident=0).min_pident == 0

    def test_valid_100(self):
        assert JobParams(min_pident=100).min_pident == 100

    def test_valid_80(self):
        assert JobParams(min_pident=80).min_pident == 80

    def test_reject_negative(self):
        with pytest.raises(ValidationError):
            JobParams(min_pident=-1)

    def test_reject_over_100(self):
        with pytest.raises(ValidationError):
            JobParams(min_pident=101)

    def test_reject_nan(self):
        with pytest.raises(ValidationError):
            JobParams(min_pident=float("nan"))

    def test_reject_inf(self):
        with pytest.raises(ValidationError):
            JobParams(min_pident=float("inf"))


class TestMinQcov:
    """Validation of min_qcov."""

    def test_none_allowed(self):
        assert JobParams(min_qcov=None).min_qcov is None

    def test_valid_zero(self):
        assert JobParams(min_qcov=0).min_qcov == 0

    def test_valid_100(self):
        assert JobParams(min_qcov=100).min_qcov == 100

    def test_reject_negative(self):
        with pytest.raises(ValidationError):
            JobParams(min_qcov=-0.1)

    def test_reject_over_100(self):
        with pytest.raises(ValidationError):
            JobParams(min_qcov=100.1)

    def test_reject_nan(self):
        with pytest.raises(ValidationError):
            JobParams(min_qcov=float("nan"))


class TestMinAlnlen:
    """Validation of min_alnlen."""

    def test_none_allowed(self):
        assert JobParams(min_alnlen=None).min_alnlen is None

    def test_valid_1(self):
        assert JobParams(min_alnlen=1).min_alnlen == 1

    def test_valid_large(self):
        assert JobParams(min_alnlen=500).min_alnlen == 500

    def test_reject_zero(self):
        with pytest.raises(ValidationError):
            JobParams(min_alnlen=0)

    def test_reject_negative(self):
        with pytest.raises(ValidationError):
            JobParams(min_alnlen=-1)


class TestTopK:
    """Validation of top_k."""

    def test_none_allowed(self):
        assert JobParams(top_k=None).top_k is None

    def test_valid_1(self):
        assert JobParams(top_k=1).top_k == 1

    def test_valid_large(self):
        assert JobParams(top_k=100).top_k == 100

    def test_reject_zero(self):
        with pytest.raises(ValidationError):
            JobParams(top_k=0)

    def test_reject_negative(self):
        with pytest.raises(ValidationError):
            JobParams(top_k=-5)


class TestDbChoice:
    """Validation of db_choice (path traversal protection)."""

    def test_empty_allowed(self):
        assert JobParams(db_choice="").db_choice == ""

    def test_valid_filename(self):
        assert JobParams(db_choice="uniprot_sprot.dmnd").db_choice == "uniprot_sprot.dmnd"

    def test_reject_dotdot(self):
        with pytest.raises(ValidationError):
            JobParams(db_choice="../etc/passwd")

    def test_reject_dotdot_only(self):
        with pytest.raises(ValidationError):
            JobParams(db_choice="..")

    def test_reject_forward_slash(self):
        with pytest.raises(ValidationError):
            JobParams(db_choice="subdir/file.dmnd")

    def test_reject_backslash(self):
        with pytest.raises(ValidationError):
            JobParams(db_choice="subdir\\file.dmnd")

    def test_reject_complex_traversal(self):
        with pytest.raises(ValidationError):
            JobParams(db_choice="../../etc/shadow")


class TestSearchTool:
    """Validation of search_tool (Literal allowlist)."""

    def test_diamond_accepted(self):
        assert JobParams(search_tool="diamond").search_tool == "diamond"

    def test_reject_unknown(self):
        with pytest.raises(ValidationError):
            JobParams(search_tool="blast")

    def test_reject_empty(self):
        with pytest.raises(ValidationError):
            JobParams(search_tool="")

    def test_reject_injection(self):
        with pytest.raises(ValidationError):
            JobParams(search_tool="diamond; rm -rf /")


