"""Tests for logging setup."""

import logging

import pytest

from healer.logging_config import configure_logging
from healer.utils import progress as progress_mod
from healer.utils.progress import progress_bar


class _FakeBar(list):
    """Stands in for tqdm: iterable, with the methods progress_bar calls."""

    def close(self):
        pass


def test_library_logger_has_a_null_handler():
    """Importing healer must not print anything on its own."""
    import healer  # noqa: F401

    handlers = logging.getLogger("healer").handlers
    assert any(isinstance(h, logging.NullHandler) for h in handlers)


def test_configure_logging_does_not_stack_handlers():
    before = len(logging.getLogger("healer").handlers)
    configure_logging("info")
    once = len(logging.getLogger("healer").handlers)
    configure_logging("debug")
    twice = len(logging.getLogger("healer").handlers)

    assert once == twice
    assert once >= before
    assert logging.getLogger("healer").level == logging.DEBUG


def test_configure_logging_rejects_unknown_level():
    with pytest.raises(ValueError, match="Unknown level"):
        configure_logging("chatty")


def test_records_are_not_duplicated_while_a_bar_is_drawn(monkeypatch, capsys):
    """
    The redirect must target whichever logger owns the handlers. Redirecting one
    that only propagates prints every record twice.
    """
    configure_logging("info")
    monkeypatch.setattr(progress_mod, "tqdm", lambda iterable, **kw: _FakeBar(iterable))

    with progress_bar([1], desc="x", show_progress=True) as bar:
        for _ in bar:
            logging.getLogger("healer.test").info("only-once-marker")

    err = capsys.readouterr().err
    assert err.count("only-once-marker") == 1
