# Installation

## Requirements

- Python ≥ 3.11
- [Conda](https://docs.conda.io/en/latest/miniconda.html) is recommended but not required

## From PyPI

```bash
pip install mol-healer
```

## From source

```bash
git clone https://github.com/eneskelestemur/healer.git
cd healer
pip install -e .
```

RDKit installs from PyPI on all common platforms. If you prefer conda:

```bash
conda create -n healer python=3.12 -c conda-forge rdkit
conda activate healer
pip install -e .
```

## Optional extras

| Extra | Install | Adds |
|-------|---------|------|
| `web` | `pip install 'mol-healer[web]'` | FastAPI/Celery web interface (`healer-ui`) |
| `opt` | `pip install 'mol-healer[opt]'` | PyGAD and BayBE for [guided enumeration](guided-enumeration.md) |
| `dev` | `pip install 'mol-healer[dev]'` | pytest, coverage, ruff |

Extras combine: `pip install 'mol-healer[web,opt]'`.

## Verifying the installation

```python
import healer

print(healer.__version__)
```

The package ships a 100-compound building block set for testing, so this works
without downloading anything:

```python
from healer import MoleculeHEALER

h = MoleculeHEALER(bb_source='test')
```

For real work you need a full building block library — see
[Building Blocks](building-blocks.md).

## Console scripts

Installing the package puts three commands on your PATH:

| Command | Purpose |
|---------|---------|
| `healer` | [Command line enumeration](cli.md) |
| `healer-ui` | [Web interface](web-interface.md) (requires the `web` extra) |
| `preprocess-bb` | [Building block preprocessing](building-blocks.md) |
