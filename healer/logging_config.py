'''
    Logging setup helpers for scripts and notebooks.
'''
import logging
import os
from typing import Optional, Union

_LEVELS = {
    'debug': logging.DEBUG,
    'info': logging.INFO,
    'warning': logging.WARNING,
    'error': logging.ERROR,
}

_DEFAULT_FORMAT = '%(asctime)s [%(levelname)s] %(message)s'
_DEFAULT_DATEFMT = '%H:%M:%S'


def configure_logging(
    level: Union[str, int] = 'info',
    fmt: str = _DEFAULT_FORMAT,
    datefmt: str = _DEFAULT_DATEFMT,
) -> logging.Logger:
    '''
        Attach a stderr handler to the `healer` logger, for scripts and notebooks
        that have not set up logging themselves. Applications with their own
        logging setup do not need this; HEALER logs through the standard
        hierarchy and stays silent until a handler is attached.

        Calling this repeatedly replaces the handler rather than stacking them.

        Args:
            level: 'debug', 'info', 'warning', 'error', or a logging constant.
            fmt: format string for the handler.
            datefmt: date format for the handler.

        Returns:
            The configured `healer` logger.

        Raises:
            ValueError: if `level` is not a recognized name.
    '''
    if isinstance(level, str):
        try:
            level = _LEVELS[level.lower()]
        except KeyError:
            raise ValueError(
                f"Unknown level {level!r}. Use one of {sorted(_LEVELS)}."
            ) from None

    logger = logging.getLogger('healer')
    for handler in [h for h in logger.handlers if getattr(h, '_healer_handler', False)]:
        logger.removeHandler(handler)

    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter(fmt, datefmt=datefmt))
    handler._healer_handler = True

    logger.addHandler(handler)
    logger.setLevel(level)
    logger.propagate = False
    return logger


def _default_level() -> Optional[int]:
    '''
        Read a level from the HEALER_LOG_LEVEL environment variable.

        Returns:
            The level, or None if the variable is unset or unrecognized.
    '''
    name = os.environ.get('HEALER_LOG_LEVEL')
    if not name:
        return None
    return _LEVELS.get(name.strip().lower())
