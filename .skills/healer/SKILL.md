---
name: healer
description: Generate synthetically accessible analogs of a molecule with HEALER (mol-healer). Use when asked to expand a hit, enumerate analogs, decorate a scaffold, or search chemical space under a scoring function. Covers MoleculeHEALER, FragmentHEALER, SiteHEALER, guided enumeration with optimizers, and the CLI.
---

# HEALER

HEALER builds analogs that come with a proposed synthesis: it fragments a query
molecule retrosynthetically, matches each fragment to purchasable building
blocks, and reassembles new molecules through validated reaction templates.

Every product is one or two real reactions away from real starting material.

## Pick the right class

| Goal | Class |
|------|-------|
| Expand a hit into analogs | `MoleculeHEALER` |
| Decorate a scaffold at chosen atoms, keeping it intact | `SiteHEALER` |
| Recombine fragments you already have | `FragmentHEALER` |

## The three-call workflow

```python
from healer import MoleculeHEALER

healer = MoleculeHEALER(bb_source="US_stock", max_bbs_per_frag=100)
healer.set_query_mol(smiles, n_compositions=10)
healer.enumerate(max_total_products=5000)

df = healer.get_results(calc_similarity=True, calc_properties=True)
```

The calls are separate so one loaded building block library serves many queries.
**Reuse the instance** across molecules rather than constructing it in a loop:

```python
for smiles in molecules:
    healer.set_query_mol(smiles, n_compositions=5)
    healer.enumerate(max_total_products=200)
    healer.save_results(f"{name}.csv")
```

## Building blocks come first

Nothing works well without a real library. `bb_source="test"` is a 100-compound
set for smoke tests only — it will produce few or no analogs for most queries.

```bash
preprocess-bb ~/Downloads/Enamine_BBs.zip -o ~/.healer/buildingblocks/
export HEALER_DATA_DIR=~/.healer
```

Then use `bb_source="US_stock"`, `"EU_stock"`, `"Global_stock"`, or a path.
Preprocessing annotates each block with the reactions it can serve, so it must be
rerun if reaction templates change.

## Controlling the search

`max_bbs_per_frag` and `sim_threshold` are alternatives. Prefer the first:

- `max_bbs_per_frag=100` — keep the 100 closest blocks per fragment. Predictable
  pool sizes and run time. Required if you want `BayesianSequenceOptimizer` to
  rank blocks by similarity.
- `sim_threshold=0.5` — used **only** when `max_bbs_per_frag=-1`. Pool size
  depends on the query, so run time is unpredictable.

`retro_tree_depth=1` splits the query in two; `2` gives up to four fragments at
much higher cost. Depth above 2 is rarely synthetically practical.

The three `enumerate` limits are **approximate**. They bound run time rather than
guaranteeing exact counts, because one reaction attempt can yield several
products.

## Guided enumeration

When the combinatorial space is too large to enumerate, an optimizer steers
assembly toward a scoring function. Higher is always better.

```python
from healer import MoleculeHEALER, GeneticAlgorithmOptimizer
from rdkit.Chem import QED

opt = GeneticAlgorithmOptimizer(population_size=50, target_fn=QED.qed)
healer.enumerate(optimizer=opt, max_evals_per_comp=2000)
```

| Scorer cost | Use |
|-------------|-----|
| Cheap (descriptors) | `BeamSearchOptimizer` |
| Expensive (docking, ML model) | `GeneticAlgorithmOptimizer` |
| Very expensive, small budget | `BayesianSequenceOptimizer` |

Pass `batch_target_fn` instead of `target_fn` when scoring benefits from
batching — it takes priority and is usually the difference between usable and
not. Return `None` for molecules you cannot score.

Requires `pip install "mol-healer[opt]"`. See
[references/optimizers.md](references/optimizers.md) for custom optimizers and
budget sizing.

## Reading the results

Row one is the query itself, as a baseline. Drop it with `df[df["BB1"] != ""]`.

Columns: `ID`, `Product`, `BB1..N` (the blocks used), `Reaction1..N-1_name`,
`URL1..N` (supplier links), plus `Similarity_to_query` and a property profile
when requested. `optimization_score` appears when an optimizer was used.

`healer.enumerated_molecules` holds the underlying records if you want RDKit
`Mol` objects rather than a table.

## Command line

```bash
healer molecule input.csv --bb-source US_stock --max-bbs-per-frag 100 \
    --max-total 5000 --properties --n-jobs -1 -o analogs.csv

healer view "c1ccccc1N"                                    # atom indices
healer site "c1ccccc1N" --reactive-sites "[6]" --rules "MW:0:300"
```

Optimizers are Python-only and not exposed on the CLI.

## Logging and progress

HEALER is silent until a handler is attached. Progress bars appear only when
stderr is a terminal, so piped output and notebooks stay clean.

```python
import healer
healer.configure_logging("info")     # or "debug" for per-composition detail
```

Use `show_progress=True/False` to override bars; it is independent of log level.

## Common mistakes

- Using `bb_source="test"` for real work — it is a 100-block smoke-test set.
- Constructing a HEALER inside a loop, reloading the library every iteration.
- Setting `sim_threshold` while `max_bbs_per_frag > 0`, where it is ignored.
- Expecting `max_total_products` to be exact.
- Reading row one as an analog; it is the query.
- Forgetting to rerun `preprocess-bb` after changing reaction templates.

## Full documentation

`docs/` in the repository: `healer-classes.md`, `guided-enumeration.md`,
`building-blocks.md`, `results.md`, `cli.md`, `logging.md`, `reactions.md`,
`architecture.md`.
