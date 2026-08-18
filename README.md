<p align="center">
  <img src="assets/healer_logo.png" alt="HEALER Logo" width="400"/>
</p>

<h1 align="center">HEALER</h1>
<h3 align="center">Hit Expansion to Advanced Leads Using Enumerated Reactions</h3>

<p align="center">
  <a href="docs/installation.md">Installation</a> •
  <a href="docs/quickstart.md">Quick Start</a> •
  <a href="docs/README.md">Documentation</a> •
  <a href="notebooks/healer_demo.ipynb">Demo</a> •
  <a href="#citation">Citation</a>
</p>

---

HEALER generates synthetically accessible molecular analogs by combining
retrosynthetic fragmentation with commercially available building blocks and
validated reaction templates. It bridges the gap between computational design and
laboratory synthesis.

## Features

- **Molecule HEALER** — retrosynthetically fragment a molecule and re-enumerate with similar building blocks
- **Fragment HEALER** — enumerate from pre-fragmented molecules
- **Site HEALER** — targeted enumeration at specific reactive sites with property filters
- **Guided enumeration** — steer the search with beam search, genetic algorithms, or Bayesian optimization
- **Synthetically accessible** — every product comes from a validated reaction template
- **Flexible** — works with any building block library in SDF format

## Installation

```bash
pip install mol-healer
```

Optional extras: `mol-healer[web]` for the browser interface, `mol-healer[opt]`
for guided enumeration. See [Installation](docs/installation.md).

## Quick Start

```python
from healer import MoleculeHEALER

healer = MoleculeHEALER(
    bb_source='test',                              # bundled 100-BB test set
    reaction_tags=['amide coupling', 'N-arylation'],
    sim_threshold=0.5,
)

healer.set_query_mol("CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O", n_compositions=10)
healer.enumerate(max_evals_per_comp=500)

results = healer.get_results(calc_similarity=True, calc_properties=True)
print(f"Generated {len(results)} analogs")
```

```bash
healer molecule "CCO" --bb-source test -o results.csv
```

Production runs need a real building block library — see
[Building Blocks](docs/building-blocks.md).

## Documentation

Full documentation is in [`docs/`](docs/README.md).

| Guide | |
|-------|---|
| [Installation](docs/installation.md) | Package, extras, console scripts |
| [Building Blocks](docs/building-blocks.md) | Downloading, preprocessing, selecting libraries |
| [Quick Start](docs/quickstart.md) | A first enumeration |
| [HEALER Classes](docs/healer-classes.md) | The three modes and their parameters |
| [Guided Enumeration](docs/guided-enumeration.md) | Optimizing toward an objective |
| [Results](docs/results.md) | Output columns and properties |
| [Command Line](docs/cli.md) | The `healer` CLI |
| [Web Interface](docs/web-interface.md) | `healer-ui` and server mode |
| [Reactions](docs/reactions.md) | Reaction library and contributing templates |
| [Architecture](docs/architecture.md) | How the pieces fit together |

[`notebooks/healer_demo.ipynb`](notebooks/healer_demo.ipynb) walks through the
Python API with visualizations.

## Contributing

Reaction templates are welcome — see [Reactions](docs/reactions.md#contributing)
for the entry format, then open an
[issue](https://github.com/eneskelestemur/healer/issues).

For code contributions:

```bash
git clone https://github.com/eneskelestemur/healer.git
cd healer
conda env create -f environment.yml && conda activate healer
pip install -e '.[dev,web,opt]'
pytest
```

## Citation

If you use HEALER in your research, please cite:

```bibtex
@article{healer2025,
  title={HEALER: Hit Expansion to Advanced Leads Using Enumerated Reactions},
  author={Kelestemur, Enes and ...},
  journal={...},
  year={2025}
}
```

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.

## Acknowledgments

- Reaction template formats mainly adapted from [datamol](https://github.com/datamol-io/datamol)
- Building block preprocessing inspired by retrosynthesis literature
- [Ketcher](https://github.com/epam/ketcher) for molecular drawing
