# HEALER

**Hit Expansion to Advanced Leads Using Enumerated Reactions**

HEALER generates novel chemical analogs by fragmenting query molecules and recombining them with commercially available building blocks using validated reaction templates.

## Installation

### Prerequisites
- Python ≥ 3.11
- [Conda](https://docs.conda.io/en/latest/miniconda.html) (recommended)

### Setup

```bash
git clone https://github.com/eneskelestemur/healer.git
cd healer

# Create environment
conda env create -f environment.yml
conda activate healer

# Install package
pip install -e .
```

### Building Block Setup

HEALER requires pre-processed building blocks. Download catalogs from [Enamine](https://enamine.net/building-blocks/building-blocks-catalog) and process them:

```bash
preprocess-bb data/buildingblocks/YOUR_FILE.zip --verbose
```

This creates `*_processed.sdf` files with reaction annotations.

## Usage

### Command Line

```bash
# Enumerate analogs for a molecule
healer molecule "CCO" --bb-source US_stock -o results.csv

# With more options
healer molecule input.csv --bb-source US_stock --reactions "amide coupling,N-arylation" \
    --n-compositions 20 --max-evals 500 --similarity -o results.csv

# Site-specific enumeration (use 'healer view' to find atom indices)
healer view "c1ccccc1N"
healer site "c1ccccc1N" --reactive-sites "[5]" --bb-source US_stock

# Fragment-based enumeration (dot-separated SMILES)
healer fragment "c1ccccc1.CC(=O)O" --bb-source US_stock

# Parallel processing
healer molecule input.csv --workers 4 -o results.csv
```

### Python API

```python
from healer import MoleculeHEALER

healer = MoleculeHEALER(
    bb_source='US_stock',
    reaction_tags=['amide coupling', 'N-arylation'],
    sim_threshold=0.3,
)

healer.set_query_mol("CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O", n_compositions=10)
healer.enumerate(max_evals_per_comp=500)
results = healer.get_results(calc_similarity=True)
```

**Available HEALER Classes:**
- `MoleculeHEALER` — Full retrosynthetic enumeration
- `SiteHEALER` — Enumeration at specific atom sites with property filters
- `FragmentHEALER` — Enumeration from pre-fragmented molecules

### Web Interface

#### Local UI

```bash
pip install healer[web]
healer-ui
```

Then open http://localhost:8000

> **Note:** The frontend must be pre-built. If you cloned from source, run:
> `cd web_client && npm install && npm run build && cd ..`

#### Server Deployment (Docker)

For multi-user deployments with background job processing:

```bash
export HEALER_USE_CELERY=true
docker-compose up
```
- API: http://localhost:8000

## Project Structure

```
healer/
├── healer/              # Core Python package
│   ├── application/     # HEALER classes (MoleculeHEALER, SiteHEALER, etc.)
│   ├── domain/          # Domain models (reactions, building blocks)
│   ├── scripts/         # Preprocessing scripts
│   ├── utils/           # Utility functions
│   ├── web/             # FastAPI server and Celery workers
│   └── cli.py           # Command-line interface
├── data/                # Building blocks and reaction templates
├── tests/               # Unit tests
└── web_client/          # React/TypeScript frontend
```

## License

MIT License. See [LICENSE](LICENSE).

## Citation

If you use HEALER in your research, please cite:

```
@software{healer,
  author = {Kelestemur, Enes},
  title = {HEALER: Hit Expansion to Advanced Leads Using Enumerated Reactions},
  url = {https://github.com/eneskelestemur/healer}
}
```
