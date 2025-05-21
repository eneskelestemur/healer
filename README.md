# HEALER: Hit Expansion to Advanced Leads Using Enumerated Reactions

## Overview

- Brief description of the project.
- Objectives and goals.

## Repository Structure

- List of primary directories and files (e.g., scripts, notebooks, environment config).

-*Note: Exclude details for files that are being ignored by Git.*

## Installation

1. **Prerequisites**
   - Install [Anaconda](https://www.anaconda.com/products/individual) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) for Python package management.

2. **Clone the Repository**
   Open a terminal and run:
   ```sh
   git clone <repository-url>
   cd <repository-directory>
   ```

3. **Set Up the Conda Environment**
   Use the provided `environment.yml` file to create the Conda environment:
   ```sh
   conda env create -f environment.yml
   conda activate healer
   ```
   *Note: Replace `healer` with the actual environment name if different.*

4. **Extract Building Blocks**
   The zip files located inside the `buildingblocks` folder should be unzipped into that same folder:
   ```sh
   unzip "buildingblocks/*.zip" -d buildingblocks/
   ```


## Usage

### MoleculeEnumerator

```python
from enumerator import MoleculeEnumerator

# Define your query molecule as a SMILES string (or an RDKit Mol object)
molecule_smiles = "CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C"  # Example SMILES

# Initialize the MoleculeEnumerator with desired parameters:
enumerator = MoleculeEnumerator(
    building_blocks='US_stock',  # Options: 'US_stock', 'EU_stock', 'Global_stock', or a custom file path
    reaction_tags=['amide coupling', 'C-N bond formation'],  # List of reaction tags to filter reactions
    custom_comp_sites=[],  # List of tuples for custom composition sites (if any)
    n_compositions=10,     # Number of compositions to enumerate
    sim_threshold=0.5      # Similarity threshold for filtering
)

# Run the enumeration process
enumerator.enumerate(molecule_smiles)

# Retrieve the results as a pandas DataFrame and print them
results = enumerator.get_results()
print(results)

# Optionally, save the results to a CSV file
enumerator.save_results("enumerated_molecules.csv")
```

### SiteEnumerator

```python
from enumerator import SiteEnumerator

# Define your query molecule as a SMILES string (or an RDKit Mol object)
molecule_smiles = "CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C"  # Example SMILES

# Initialize the SiteEnumerator with desired parameters:
enumerator = SiteEnumerator(
    building_blocks='US_stock',  # Options: 'US_stock', 'EU_stock', 'Global_stock', or a custom file path
    reaction_sites=[1, 2, 3],  # List of atom indices to consider as reaction sites
    reaction_tags=['amide coupling', 'C-N bond formation'],  # List of reaction tags to filter reactions
    rules={
        'MW': (0, 500),       # Molecular weight range
        'HBD': (0, 5),        # Hydrogen bond donors
        'HBA': (0, 10),       # Hydrogen bond acceptors
        'TPSA': (0, 200),     # Topological polar surface area
        'RotB': (0, 10),      # Rotatable bonds
        'Rings': (0, 10),     # Number of rings
        'ArRings': (0, 5),    # Number of aromatic rings
        'Chiral': (0, 5)      # Number of chiral centers
    },
    struct_rules=[]  # Optionally add structure-based rules as a list of SMARTS strings
)

# Run the enumeration process
enumerator.enumerate(molecule_smiles)

# Retrieve the results as a pandas DataFrame and print them
results = enumerator.get_results()
print(results)

# Optionally, save the results to a CSV file
enumerator.save_results("site_enumerated_molecules.csv")
```

### Web Application 

Alternatively, you can host the web application locally by running the following command line. This will also allow you to visualize the enuermated molecules.

```sh
python webserver/EnumeratorApp.py
```

### CLI 


## Key Files

-**enumerator.py**: Describe the main functionality.

-**reaction.py**: Describe its purpose.

-**utils.py**: Provide a brief overview.

- Others as needed.

## Contributing

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
