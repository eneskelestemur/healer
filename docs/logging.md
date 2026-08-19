# Logging & Progress

HEALER logs through the standard `logging` module and stays silent until you
attach a handler. Progress bars are separate from log verbosity and are
controlled independently.

## Turning logging on

In a script or notebook:

```python
import healer

healer.configure_logging('info')     # 'debug', 'info', 'warning', 'error'
```

This attaches a stderr handler to the `healer` logger. Calling it again replaces
the handler rather than stacking duplicates.

You can also set a level without touching code:

```bash
export HEALER_LOG_LEVEL=debug
```

Applications with their own logging setup need neither — HEALER logs into the
standard hierarchy and inherits whatever you have configured:

```python
import logging

logging.basicConfig(level=logging.INFO)
logging.getLogger('healer.domain').setLevel(logging.WARNING)   # quieten a subtree
```

Logger names follow the package layout, so `healer.application.healer`,
`healer.domain.bb_repository`, and so on can each be tuned separately.

## What each level means

| Level | Contents |
|-------|----------|
| `ERROR` | An operation failed and you lost work, such as one molecule in a batch |
| `WARNING` | Results differ from what you asked for: a reaction template skipped, building block pools truncated, a query that could not be fragmented, an optimizer that abandoned a composition |
| `INFO` | Milestones — library loaded, enumeration finished with counts and timing |
| `DEBUG` | Per-composition and per-stage detail: stop reasons, composition dumps, per-round counts |

Nothing is logged per product or per candidate above `DEBUG`. Where a loop needs
to report a problem it counts occurrences and emits one summary, so a long run
cannot bury you in repeated lines.

`INFO` is a handful of lines per run:

```
[INFO] Loaded 100 building blocks indexed for 82 reaction types
[INFO] Enumerated 50 product(s) from 2 composition(s) in 0.1s
```

## Progress bars

Bars are drawn when stderr is a terminal, so piped output, log files, CI, and
executed notebooks stay clean. Override per instance:

```python
healer = MoleculeHEALER(bb_source='test', show_progress=True)    # always
healer = MoleculeHEALER(bb_source='test', show_progress=False)   # never
```

or globally:

```bash
export HEALER_PROGRESS=0
```

**Only the outermost bar is drawn.** Enumerating a batch shows one bar over
molecules; the per-composition bars inside each `enumerate()` suppress
themselves. Nothing nests, and nothing is redrawn by an inner loop.

While a bar is live, log records are routed through `tqdm.write`, so they scroll
above a bar that stays pinned to the bottom line. Running at `DEBUG` with bars on
is legible.

> `verbose` is still accepted as a deprecated alias for `show_progress`, where
> `verbose >= 1` enables bars. It no longer affects log verbosity — use
> `configure_logging` for that.

## Command line

| Flag | Effect |
|------|--------|
| *(none)* | `INFO` |
| `-v` | `DEBUG` |
| `-q`, `--quiet` | `WARNING`, and no progress bars |

```bash
healer molecule inputs.csv --bb-source test -q -o out.csv    # errors only
healer molecule inputs.csv --bb-source test -v -o out.csv    # full detail
```

## Web interface

The web app logs through the same hierarchy and leaves formatting to whatever is
hosting it, so uvicorn and Celery control the output.
