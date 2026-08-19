"""
Progress bar helpers.
"""

import logging
import os
import sys
from contextlib import contextmanager
from typing import Any, Iterable, Iterator, Optional

from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm

logger = logging.getLogger(__name__)

# Depth of currently drawn bars. Only the outermost caller draws one, so an
# inner loop does not tear down and redraw the bar wrapping it.
_active_bars = 0


def _redirect_targets() -> list:
    """
    Pick the loggers whose handlers should write through `tqdm.write` while a
    bar is drawn. Redirecting a logger that only propagates would print each
    record twice, once from the redirect and once from the handler upstream.

    Returns:
        The loggers that own the stream handlers records reach.
    """
    healer_logger = logging.getLogger("healer")
    owns_handlers = any(
        not isinstance(h, logging.NullHandler) for h in healer_logger.handlers
    )
    if owns_handlers and not healer_logger.propagate:
        return [healer_logger]
    return [logging.getLogger()]


def progress_enabled(show_progress: Optional[bool] = None) -> bool:
    """
    Decide whether progress bars should be drawn.

    Args:
        show_progress: explicit choice, which takes priority. None consults
            the HEALER_PROGRESS environment variable, then falls back to
            whether stderr is a terminal.

    Returns:
        True if bars should be drawn.
    """
    if show_progress is not None:
        return bool(show_progress)

    env = os.environ.get("HEALER_PROGRESS")
    if env is not None:
        return env.strip().lower() not in ("0", "false", "no", "off", "")

    return sys.stderr.isatty()


class _NullBar:
    """
    Stand-in used when no bar is drawn. Iterates the wrapped iterable and
    accepts the tqdm calls made by the enumeration loops as no-ops.
    """

    def __init__(self, iterable: Iterable[Any]) -> None:
        self._iterable = iterable

    def __iter__(self) -> Iterator[Any]:
        return iter(self._iterable)

    def set_postfix(self, *args: Any, **kwargs: Any) -> None:
        pass

    def update(self, n: int = 1) -> None:
        pass

    def close(self) -> None:
        pass


@contextmanager
def progress_bar(
    iterable: Iterable[Any],
    desc: str,
    total: Optional[int] = None,
    show_progress: Optional[bool] = None,
    unit: str = "it",
) -> Iterator[Any]:
    """
    Wrap an iterable in a progress bar, unless one is already being drawn or
    progress is switched off. While a bar is drawn, log records from the
    `healer` logger are routed through `tqdm.write` so they scroll above it
    instead of breaking it apart.

    Args:
        iterable: the iterable to wrap.
        desc: bar label.
        total: item count, for iterables without a length.
        show_progress: explicit override, see `progress_enabled`.
        unit: unit label for the rate display.

    Yields:
        The iterable, wrapped in a `tqdm` bar when one is drawn.
    """
    global _active_bars

    if _active_bars > 0 or not progress_enabled(show_progress):
        yield _NullBar(iterable)
        return

    _active_bars += 1
    try:
        with logging_redirect_tqdm(loggers=_redirect_targets()):
            bar = tqdm(iterable, desc=desc, total=total, unit=unit)
            try:
                yield bar
            finally:
                bar.close()
    finally:
        _active_bars -= 1
