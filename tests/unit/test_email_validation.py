"""Tests for email validation in JobParams."""

import pytest
from pydantic import ValidationError

from metagomics2.models.job import JobParams


class TestEmailValidation:
    """Tests for the notification_email field validator."""

    def test_valid_email_accepted(self):
        params = JobParams(notification_email="user@example.com")
        assert params.notification_email == "user@example.com"

    def test_empty_email_allowed(self):
        params = JobParams(notification_email="")
        assert params.notification_email == ""

    def test_default_email_is_empty(self):
        params = JobParams()
        assert params.notification_email == ""

    def test_whitespace_stripped(self):
        params = JobParams(notification_email="  user@example.com  ")
        assert params.notification_email == "user@example.com"

    def test_whitespace_only_becomes_empty(self):
        params = JobParams(notification_email="   ")
        assert params.notification_email == ""

    def test_invalid_email_no_at(self):
        with pytest.raises(ValidationError):
            JobParams(notification_email="userexample.com")

    def test_invalid_email_no_domain(self):
        with pytest.raises(ValidationError):
            JobParams(notification_email="user@")

    def test_invalid_email_no_tld(self):
        with pytest.raises(ValidationError):
            JobParams(notification_email="user@example")

    def test_invalid_email_spaces_in_address(self):
        with pytest.raises(ValidationError):
            JobParams(notification_email="user @example.com")

    def test_valid_email_with_subdomain(self):
        params = JobParams(notification_email="user@mail.example.com")
        assert params.notification_email == "user@mail.example.com"

    def test_valid_email_with_plus(self):
        params = JobParams(notification_email="user+tag@example.com")
        assert params.notification_email == "user+tag@example.com"
