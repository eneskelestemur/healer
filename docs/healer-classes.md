# HEALER Classes

Three classes share one workflow — construct, set a query, enumerate — and differ
in how the query is turned into fragments.

| Class | Input | Fragments come from |
|-------|-------|--------------------|
| [`MoleculeHEALER`](#moleculehealer) | one molecule | a retrosynthesis tree |
| [`FragmentHEALER`](#fragmenthealer) | a multi-component SMILES | the components you supply |
| [`SiteHEALER`](#sitehealer) | a molecule plus atom indices | none; the molecule is coupled directly |

## Shared constructor parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `bb_source` | Named source or path to a processed SDF | `'US_stock'` |
| `reaction_tags` | Tag, list of tags, or `'all'` | see below |
| `bb_repository` | Pre-loaded `BBRepository` to share | `None` |
| `shuffle_bb_order` | Randomize building block order | `False` |
| `show_progress` | Draw progress bars; None = only when stderr is a terminal | `None` |

`reaction_tags` defaults to `['amide coupling', 'amide', 'C-N bond formation',
'C-N', 'alkylation', 'N-arylation', 'azole', 'amination']`. See
[Reactions](reactions.md) for the full tag list.

Log verbosity is separate from `show_progress` — see
[Logging & Progress](logging.md). `verbose` is accepted as a deprecated alias
for `show_progress`.

## Shared `enumerate()` parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `optimizer` | Optimizer to steer the search — see [Guided Enumeration](guided-enumeration.md) | `None` |
| `max_evals_per_comp` | Max reaction attempts per composition | `None` |
| `max_products_per_comp` | Max products per composition | `None` |
| `max_total_products` | Stop after this many products in total | `None` |
| `n_jobs` | Parallel workers for synthesis, `-1` for all cores | `1` |

> **The three limits are approximate.** They bound run time and guard against
> combinatorial explosion; they are not exact quotas. One reaction attempt can
> yield several products, so totals may overshoot.

---

## MoleculeHEALER

Fragments the query with a retrosynthesis tree, then finds building blocks
similar to each fragment.

```python
from healer import MoleculeHEALER

healer = MoleculeHEALER(bb_source='test', sim_threshold=0.5, max_bbs_per_frag=50)
healer.set_query_mol(smiles, n_compositions=10, retro_tree_depth=1)
healer.enumerate(max_total_products=1000)
```

### Constructor

| Parameter | Description | Default |
|-----------|-------------|---------|
| `sim_threshold` | Minimum fragment-to-BB similarity | `0.5` |
| `max_bbs_per_frag` | Keep only the N most similar BBs per fragment; `-1` uses the threshold | `-1` |

`max_bbs_per_frag` gives predictable pool sizes and run times, and is required if
you want [BayBE](guided-enumeration.md#bayesiansequenceoptimizer) to rank
building blocks by similarity. `sim_threshold` is used only when
`max_bbs_per_frag` is `-1`.

Matching uses a Tversky similarity weighted by a size penalty, so building blocks
much larger than the fragment are down-ranked.

### `set_query_mol()`

| Parameter | Description | Default |
|-----------|-------------|---------|
| `query_mol` | SMILES or `Chem.Mol`, a single connected component | — |
| `n_compositions` | Number of fragment compositions to keep | `10` |
| `retro_tree_depth` | Retrosynthesis tree depth | `1` |
| `min_frag_size` | Minimum fragment size in heavy atoms | `3` |
| `max_retro_nodes` | Node budget for the retrosynthesis tree; `None` removes it | `10000` |
| `randomize_compositions` | Shuffle rather than order by fragment count | `False` |
| `random_seed` | Seed, `-1` for none | `-1` |
| `custom_split_sites` | Explicit bonds to break, skipping the tree | `None` |

Depth 1 splits the molecule in two; depth 2 splits those halves again, giving up
to four fragments. Cost grows exponentially with depth and routes beyond depth 2
are rarely practical.

`max_retro_nodes` bounds that growth. In ordinary use it never fires — fragments
fall below `min_frag_size` after a couple of levels, which prunes the tree well
before the default budget. It is a guard against pathological inputs: a large
query, `reaction_tags='all'`, and a low `min_frag_size` together. When it does
fire, enumeration continues with the smaller tree and logs a warning.

Compositions are ordered by fragment count and truncated to `n_compositions`.

To break specific bonds instead:

```python
healer.set_query_mol(smiles, custom_split_sites=[[(3, 4)], [(7, 8), (11, 12)]])
```

Each inner list is one composition, given as `(begin_atom, end_atom)` bond pairs.

---

## FragmentHEALER

Takes fragments directly, skipping retrosynthesis. Use it when you already know
how the molecule should be decomposed.

```python
from healer import FragmentHEALER

healer = FragmentHEALER(bb_source='test', max_bbs_per_frag=50)
healer.set_query_mol('c1ccccc1N.CC(=O)O')       # or a tuple of SMILES/Mols
healer.enumerate(max_total_products=500)
```

Constructor parameters match `MoleculeHEALER`. `set_query_mol()` takes only the
query, which must have at least two components.

---

## SiteHEALER

Couples a single building block onto the query at chosen atoms, leaving the rest
of the molecule untouched. There is no fragmentation, so products always contain
the intact query.

```python
from healer import SiteHEALER

healer = SiteHEALER(bb_source='test', rules={'MW': (0, 300)})
healer.set_query_mol('c1ccccc1N', reactive_sites=[6])
healer.enumerate(max_total_products=500)
```

Use `healer view "<smiles>"` to see atom indices.

### Constructor

| Parameter | Description | Default |
|-----------|-------------|---------|
| `rules` | Property windows for filtering building blocks | see below |
| `struct_rules` | SMARTS patterns building blocks must match | `[]` |

`rules` defaults to `MW (0, 500)`, `HBD (0, 5)`, `HBA (0, 10)`, `TPSA (0, 200)`,
`RotB (0, 10)`, `Rings (0, 10)`, `ArRings (0, 5)`, `Chiral (0, 5)`. Supply a dict
to override, or call `set_rules(MW=(0, 300))` afterwards.

Because products keep the whole query, filtering building blocks by property is
how you keep products in a reasonable range.

### `set_query_mol()`

| Parameter | Description | Default |
|-----------|-------------|---------|
| `query_mol` | SMILES or `Chem.Mol` | — |
| `reactive_sites` | Atom indices allowed to react | `None` (all atoms) |

Atoms outside `reactive_sites` and their immediate neighbours are protected from
reacting. Leaving it `None` allows reactions anywhere and logs a warning.
