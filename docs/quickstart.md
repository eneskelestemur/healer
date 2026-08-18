# Quick Start

## Python

```python
from healer import MoleculeHEALER

healer = MoleculeHEALER(
    bb_source='test',                              # bundled 100-BB set
    reaction_tags=['amide coupling', 'N-arylation'],
    sim_threshold=0.5,
)

healer.set_query_mol(
    "CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O",  # penicillin G
    n_compositions=10,
)
healer.enumerate(max_evals_per_comp=500)

results = healer.get_results(calc_similarity=True, calc_properties=True)
print(f"Generated {len(results)} analogs")
```

The three calls map to the three stages of a run:

| Call | Does | Reused across |
|------|------|---------------|
| `MoleculeHEALER(...)` | Loads building blocks and reactions | Everything |
| `set_query_mol(...)` | Fragments one query molecule | Many `enumerate()` calls |
| `enumerate(...)` | Reassembles products | — |

The split exists so the building block library is loaded once and reused across
many query molecules:

```python
for smiles in my_molecules:
    healer.set_query_mol(smiles, n_compositions=5)
    healer.enumerate(max_total_products=100)
    healer.save_results(f'{smiles}.csv')
```

## Command line

```bash
# one molecule
healer molecule "CC(=O)Nc1ccccc1" --bb-source test -o results.csv

# a file of molecules, using all cores
healer molecule inputs.csv --n-jobs -1 -o results.csv

# view atom indices, needed for site-specific runs
healer view "c1ccccc1N"

# enumerate only at atom 6
healer site "c1ccccc1N" --reactive-sites "[6]" --bb-source test
```

## Next steps

- [HEALER Classes](healer-classes.md) — choosing between the three modes
- [Guided Enumeration](guided-enumeration.md) — optimizing toward an objective
- [Results](results.md) — what comes out
