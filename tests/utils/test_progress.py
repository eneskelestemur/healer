"""Tests for progress bar policy."""

import pytest

from healer.utils import progress as progress_mod
from healer.utils.progress import progress_bar, progress_enabled


def test_explicit_choice_wins_over_environment(monkeypatch):
    monkeypatch.setenv("HEALER_PROGRESS", "1")
    assert progress_enabled(False) is False
    monkeypatch.setenv("HEALER_PROGRESS", "0")
    assert progress_enabled(True) is True


@pytest.mark.parametrize(
    "value,expected",
    [
        ("0", False),
        ("false", False),
        ("no", False),
        ("off", False),
        ("1", True),
        ("yes", True),
    ],
)
def test_environment_variable_is_honoured(monkeypatch, value, expected):
    monkeypatch.setenv("HEALER_PROGRESS", value)
    assert progress_enabled() is expected


def test_falls_back_to_tty_check(monkeypatch):
    monkeypatch.delenv("HEALER_PROGRESS", raising=False)
    monkeypatch.setattr("sys.stderr.isatty", lambda: False, raising=False)
    assert progress_enabled() is False


class _FakeBar(list):
    """Stands in for tqdm: iterable, with the methods progress_bar calls."""

    def close(self):
        pass


def test_inner_bar_is_suppressed_while_outer_is_active(monkeypatch):
    """
    Only the outermost caller draws a bar. A nested bar would tear down and
    redraw the one wrapping it on every iteration.
    """
    drawn = []

    def fake_tqdm(iterable, **kw):
        drawn.append(kw["desc"])
        return _FakeBar(iterable)

    monkeypatch.setattr(progress_mod, "tqdm", fake_tqdm)

    with progress_bar([1, 2], desc="outer", show_progress=True) as outer:
        for _ in outer:
            with progress_bar([1], desc="inner", show_progress=True) as inner:
                list(inner)

    assert drawn == ["outer"]


def test_bar_depth_is_released_on_exception():
    with pytest.raises(RuntimeError):
        with progress_bar([1], desc="boom", show_progress=True):
            raise RuntimeError("boom")
    assert progress_mod._active_bars == 0


def test_suppressed_bar_still_supports_tqdm_calls():
    """The stand-in accepts the calls the enumeration loops make on a bar."""
    with progress_bar([1, 2], desc="off", show_progress=False) as bar:
        bar.set_postfix(evals=1, refresh=False)
        bar.update(1)
        assert list(bar) == [1, 2]
        bar.close()
