# HEALER: Hit Expansion to Advanced Leads Using Enumerated Reactions

HEALER is a powerful computational chemistry tool designed for **hit expansion** and **lead optimization** through systematic enumeration of chemical reactions. It enables researchers to explore chemical space by generating novel molecules from query structures using retrosynthetic analysis and building block databases.

## 🎯 Overview

HEALER performs **retrosynthetic fragmentation** of query molecules and systematically **recombines** the fragments with commercially available building blocks to generate chemically feasible analogs. The tool supports three main enumeration strategies:

- **MoleculeHEALER**: Full retrosynthetic analysis and enumeration
- **SiteHEALER**: Site-specific enumeration with property-based filtering  
- **FragmentHEALER**: Fragment-based enumeration

## 🏗️ Repository Structure

```
healer/
├── healer/                    # Core package
│   ├── application/           # Main HEALER classes
│   ├── domain/               # Domain models (reactions, building blocks, etc.)
│   ├── utils/                # Utility functions
│   ├── data/                 # Data files
│   │   ├── buildingblocks/   # Building block databases (ZIP files)
│   │   └── reactions/        # Reaction templates
│   └── cli.py               # Command-line interface
├── scripts/                  # Preprocessing and utility scripts
├── webserver/               # Dash web application
├── benchmark/               # Benchmarking data and results
├── external/                # External/optional dependencies
├── tests/                   # Unit tests
└── notebooks/               # Analysis and visualization notebook
```

## 🚀 Installation

### Prerequisites

- Python ≥ 3.11
- [Anaconda](https://www.anaconda.com/products/individual) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html)

### Quick Install

1. **Clone the repository:**
   ```bash
   git clone https://github.com/eneskelestemur/healer.git
   cd healer
   ```

2. **Create and activate the conda environment:**
   ```bash
   conda env create -f environment.yml
   conda activate healer
   ```

3. **Install HEALER in development mode:**
   ```bash
   pip install -e .
   ```

### Building Block Setup

HEALER needs a pre-processed building block database. You can use your own library of building blocks or use commercial building block databases from **Enamine**, which are freely available:

1. **Download building blocks** from [Enamine Building Blocks Catalog](https://enamine.net/building-blocks/building-blocks-catalog)
   - Download the SDF files for your desired catalog (US Stock, EU Stock, or Global Stock)
   - Place the downloaded ZIP files in `healer/data/buildingblocks/`

2. **Preprocess building blocks:**
   
   This step might take some time depending on the number of files and their sizes. For US Stock, it takes around 10-12 minutes.

   ```bash
   # Process a single building block file (automatically extracts ZIP files)
   # Change the zip file name as needed
   python scripts/preprocess_bb_source.py healer/data/buildingblocks/Enamine_building_blocks_stock.zip --verbose
   
   # Process all building block files
   for file in healer/data/buildingblocks/*.zip; do
       python scripts/preprocess_bb_source.py "$file" --verbose
   done
   ```

   The preprocessing script:
   - Automatically extracts ZIP files
   - Removes salt/solvent fragments
   - Annotates molecules with compatible reactions
   - Creates `*_processed.sdf` files ready for use

## 💻 Usage

### Python API

#### MoleculeHEALER - Full Retrosynthetic Enumeration

```python
from healer.application.healer import MoleculeHEALER

# Initialize with building block source and reaction filters
healer = MoleculeHEALER(
    bb_supplier='US_stock',  # Options: 'US_stock', 'EU_stock', 'Global_stock'
    reaction_tags=['amide coupling', 'N-arylation', 'alkylation'],  # or 'all'
    max_evals_per_comp=1000,
    n_compositions=10,
    sim_threshold=0.30,
    max_bbs_per_comp=10,
    verbose=1
)

# Set query molecule
query_smiles = "CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O"  # Penicillin G
healer.set_query_mol(
    query_mol=query_smiles,
    retro_tree_depth=2,     # Depth of retrosynthetic analysis
    min_frag_size=3,        # Minimum fragment size
    custom_split_sites=None # Optional: custom reaction sites
)

# Perform enumeration
healer.enumerate()

# Get results
results = healer.get_results(
    calc_similarity=True,    # Calculate Tanimoto similarity to query
    calc_stoplight=False,    # Calculate stoplight scores
    calc_cns_mpo=False      # Calculate CNS-MPO scores
)

print(f"Generated {len(results)} analogs")
print(results.head())

```

#### SiteHEALER - Site-Specific Enumeration

```python
from healer.application.healer import SiteHEALER

# Initialize with property-based filters
healer = SiteHEALER(
    bb_supplier='EU_stock',
    reaction_tags=['amide coupling', 'C-N bond formation'],
    max_evals_per_comp=500,
    rules={
        'MW': (100, 500),      # Molecular weight range
        'HBD': (0, 5),         # H-bond donors
        'HBA': (0, 10),        # H-bond acceptors  
        'TPSA': (0, 140),      # Topological polar surface area
        'RotB': (0, 10),       # Rotatable bonds
        'Rings': (1, 4),       # Number of rings
        'ArRings': (0, 3),     # Aromatic rings
        'Chiral': (0, 2)       # Chiral centers
    },
    struct_rules=[],           # Optional SMARTS patterns
    verbose=1
)

# Set query with specific reactive sites
healer.set_query_mol(
    query_mol=query_smiles,
    reactive_sites=[5, 12, 18]  # Atom indices for reaction sites
)

healer.enumerate()
results = healer.get_results(calc_similarity=True)
```

#### FragmentHEALER - Fragment-Based Enumeration

```python
from healer.application.healer import FragmentHEALER

# Initialize for fragment-based enumeration
healer = FragmentHEALER(
    bb_supplier='Global_stock',
    reaction_tags='all',
    max_evals_per_comp=2000,
    sim_threshold=0.25,
    max_bbs_per_comp=15,
    verbose=1
)

# Set query molecule containing multiple fragments
query_smiles = "CC2(C)SC1C(N)C(=O)N1C2C(=O)O.O=C(O)Cc1ccccc1"
healer.set_query_mol(query_mol=query_smiles)
healer.enumerate()
results = healer.get_results(calc_similarity=True)
```

### Command-Line Interface

The CLI supports batch processing and parallel execution:

```bash
# Single molecule enumeration
healer molecule "CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O" \
    --bb_source US_stock \
    --reaction_tags "amide coupling,N-arylation" \
    --max_evals_per_comp 1000 \
    --calculate_similarity \
    --output results.csv

# Batch processing from file
healer molecule molecules.smi \
    --header \
    --column_name smiles \
    --bb_source EU_stock \
    --workers 4 \
    --output batch_results.csv

# Site enumeration with filters
healer site "CN1C(=O)C2C(c3ccccc3)NNC2C1=O" \
    --bb_source Global_stock \
    --MW_range 150:400 \
    --HBD_range 0:3 \
    --TPSA_range 0:100 \
    --output site_results.csv

# Fragment enumeration
healer fragment "CC2(C)SC1C(N)C(=O)N1C2C(=O)O.O=C(O)Cc1ccccc1" \
    --bb_source US_stock \
    --sim_threshold 0.40 \
    --output fragment_results.csv

# Available options
healer --help
```

### Web Application

Launch the interactive Dash web interface:

```bash
# Navigate to webserver directory
cd webserver

# Local Mode (default - unlimited parameters)
python app.py

# Server Mode (enforced limits for shared deployment) 
export HEALER_SERVER_MODE=true && python app.py
```

Then navigate to `http://localhost:8053` in your browser for the modern web interface featuring:
- **Automatic Fragment Detection**: Smart switching between MoleculeHEALER and FragmentHEALER
- **Server/Local Mode**: Configurable computational limits for deployment
- **Enhanced UI**: Bootstrap Minty theme with comprehensive tooltips
- **Tab-based Interface**: Separate Molecule HEALER and Site HEALER workflows
- **Real-time Enumeration**: Background processing with progress indicators
- **Results Export**: CSV download with building block URLs
- **Parameter Validation**: Real-time validation with user feedback

## 🔧 Configuration

### Building Block Sources

- **US_stock**: US-available Enamine building blocks
- **EU_stock**: EU-available Enamine building blocks  
- **Global_stock**: Complete global Enamine catalog
- **Custom path**: Path to your own processed SDF file

### Reaction Tags

Available reaction types include:
- `amide coupling`
- `N-arylation` 
- `alkylation`
- `azole formation`
- `C-C coupling`
- `cyclization`
- And many more... (use `'all'` for complete list)

### Optimization (Still in development)

HEALER supports integration with optimization algorithms for objective-driven enumeration:

```python
from healer.application.optimizers import BaseStagewiseOptimizer

# Custom optimizer example
class MyOptimizer(BaseStagewiseOptimizer):
    def target_fn(self, mol):
        # Your objective function (e.g., docking score, QSAR prediction)
        return calculate_score(mol)
    
    def filter(self, candidates, depth):
        # Filter candidates at each stage
        return top_candidates

healer.enumerate(optimizer=MyOptimizer())
```

## 📊 Output Format

Results are returned as pandas DataFrames with the following columns:

- **ID**: Unique identifier (HEAL_XXXXXX)
- **Product**: SMILES of the enumerated molecule
- **BB1, BB2, ...**: Building blocks used
- **Reaction1_name, Reaction2_name, ...**: Reactions applied
- **URL1, URL2, ...**: Enamine catalog URLs for building blocks
- **Similarity_to_query**: Tanimoto similarity (if calculated)
- **Additional properties**: Stoplight scores, CNS-MPO, custom metrics

## 🧪 Examples

See the `benchmark/` directory for comprehensive examples and the `figures.ipynb` notebook for analysis workflows.

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📚 Citation

If you use HEALER in your research, please cite:

```bibtex
@software{healer2024,
  title={HEALER: Hit Expansion to Advanced Leads Using Enumerated Reactions},
  author={Kelestemur, Enes},
  year={2024},
  url={https://github.com/eneskelestemur/healer}
}
```

## 🆘 Support

- **Issues**: Report bugs and request features via [GitHub Issues](https://github.com/eneskelestemur/healer/issues)
- **Documentation**: Additional examples in `figures.ipynb`
- **Contact**: enesk@email.unc.edu

---

**Note**: HEALER uses freely available building block databases from Enamine. Please comply with Enamine's terms of use when downloading and using their catalogs.
