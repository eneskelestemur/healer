# HEALER Documentation

HEALER generates synthetically accessible analogs of a query molecule. It
retrosynthetically fragments the molecule, matches each fragment to commercially
available building blocks, and reassembles new molecules through validated
reaction templates.

## Getting started

| Page | Contents |
|------|----------|
| [Installation](installation.md) | Installing the package and its optional extras |
| [Building Blocks](building-blocks.md) | Downloading, preprocessing, and selecting BB libraries |
| [Quick Start](quickstart.md) | A first enumeration in Python and on the command line |

## Using HEALER

| Page | Contents |
|------|----------|
| [HEALER Classes](healer-classes.md) | `MoleculeHEALER`, `FragmentHEALER`, `SiteHEALER` and their parameters |
| [Guided Enumeration](guided-enumeration.md) | Steering enumeration toward an objective with optimizers |
| [Results](results.md) | Output columns, similarity, and property profiling |
| [Logging & Progress](logging.md) | Log levels, progress bars, and quieting output |
| [Command Line](cli.md) | The `healer` CLI |
| [Web Interface](web-interface.md) | The `healer-ui` local app and server mode |

## Reference

| Page | Contents |
|------|----------|
| [Reactions](reactions.md) | The reaction library, tags, and contributing templates |
| [Architecture](architecture.md) | How the pieces fit together |

## Demo

[`notebooks/healer_demo.ipynb`](../notebooks/healer_demo.ipynb) walks through the
Python API end to end with visualizations.
