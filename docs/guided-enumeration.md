# Guided Enumeration

By default enumeration is exhaustive: every fragment is paired with every
compatible building block. When the search space is too large for that, an
optimizer steers assembly toward a scoring function so only a promising subset is
built.

Optimizers are available through the Python API only. Install them with:

```bash
pip install 'mol-healer[opt]'
```

## Scoring functions

You supply the objective. Higher is always better — every optimizer maximizes.

```python
from rdkit.Chem import Crippen, QED

healer.enumerate(optimizer=BeamSearchOptimizer(beam_width=50, target_fn=QED.qed))
```

Scoring is usually the expensive step, so a batched form is also accepted and
takes priority when both are given:

```python
def batch_score(mols):
    return my_model.predict([Chem.MolToSmiles(m) for m in mols])

opt = BeamSearchOptimizer(beam_width=50, batch_target_fn=batch_score)
```

`batch_target_fn` must return one value per input molecule, using `None` for any
it cannot score. Failures are logged and the molecule is dropped from the
optimizer's feedback rather than ending the run.

## Two families

| Family | Interface | Prunes |
|--------|-----------|--------|
| Stagewise | `select_candidates`, `filter` | at each assembly stage |
| Sequence | `init_search`, `ask`, `tell` | proposes whole building block combinations |

### Choosing one

| Situation | Use |
|-----------|-----|
| Cheap scorer, want the best intermediates kept | `BeamSearchOptimizer` |
| Expensive scorer, large search space | `GeneticAlgorithmOptimizer` |
| Very expensive scorer, small evaluation budget | `BayesianSequenceOptimizer` |

---

## BeamSearchOptimizer

Scores assembled intermediates at each stage and carries forward the best
`beam_width`.

```python
from healer import MoleculeHEALER, BeamSearchOptimizer
from rdkit.Chem import QED

healer = MoleculeHEALER(bb_source='test', max_bbs_per_frag=50)
healer.set_query_mol(smiles, n_compositions=5)
healer.enumerate(optimizer=BeamSearchOptimizer(beam_width=100, target_fn=QED.qed))
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `beam_width` | Intermediates kept per stage | `100` |

It scores every intermediate at every stage, so it calls the scorer far more than
the sequence optimizers. Keep the scorer cheap.

## GeneticAlgorithmOptimizer

Evolves building block combinations with PyGAD. One `ask`/`tell` round is one
generation.

```python
from healer import GeneticAlgorithmOptimizer

opt = GeneticAlgorithmOptimizer(
    population_size=50,
    mutation_percent_genes=10,
    keep_elitism=2,
    random_seed=42,
    target_fn=QED.qed,
)
healer.enumerate(optimizer=opt, max_evals_per_comp=2000)
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `population_size` | Combinations proposed per generation | `50` |
| `mutation_percent_genes` | Percent of genes mutated | `10` |
| `crossover_type` | PyGAD crossover operator | `'uniform'` |
| `keep_elitism` | Best solutions carried over unchanged | `2` |
| `random_seed` | Seed for reproducibility | `None` |
| `max_domain_per_frag` | Cap on building blocks per fragment | `None` |

Genes are building block indices, one per fragment, so a two-fragment query gives
a two-gene problem. With few genes, raise `population_size` rather than the
mutation rate to keep exploring — at two genes the mutation count is already
clamped to one.

If the search converges early, lower `keep_elitism` or raise `population_size`.

## BayesianSequenceOptimizer

Bayesian optimization through BayBE, with building blocks featurized as MORDRED
descriptors. Starts with diverse random picks and switches to a surrogate model
once it has measurements.

```python
from healer import BayesianSequenceOptimizer

opt = BayesianSequenceOptimizer(batch_size=10, target_fn=expensive_score)
healer.enumerate(optimizer=opt, max_evals_per_comp=200)
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `batch_size` | Combinations recommended per round | `10` |
| `encoding` | BayBE `SubstanceParameter` encoding | `'MORDRED'` |
| `max_domain_per_frag` | Cap on building blocks per fragment | `200` |

BayBE enumerates the full product of the building block pools up front, so the
pools are capped by default. **Set `max_bbs_per_frag` on the HEALER** so the cap
keeps the building blocks most similar to each fragment; otherwise it truncates
arbitrarily and logs a warning:

```python
healer = MoleculeHEALER(bb_source='US_stock', max_bbs_per_frag=200)
```

---

## Budgets

Each composition is searched independently: `init_search` is called once per
composition, so an optimizer starts cold each time and

```
total evaluations ≈ n_compositions × max_evals_per_comp
```

Size `n_compositions` and `max_evals_per_comp` together, and remember the limits
are approximate.

## Failures

An optimizer that cannot continue raises `OptimizerError`. Enumeration logs a
warning, abandons that composition, and moves to the next — one bad composition
does not end the run. Running out of candidates is not an error; the search
simply ends for that composition.

## Custom optimizers

Subclass either base class. For stagewise, `filter` is required and
`select_candidates` is optional:

```python
from healer import BaseStagewiseOptimizer

class HeavyAtomCap(BaseStagewiseOptimizer):
    def select_candidates(self, candidates, depth):
        # runs before any reaction, so this saves the synthesis work
        for rec, bb, rxn in candidates:
            if bb.num_heavy_atoms <= 20:
                yield rec, bb, rxn

    def filter(self, records, depth):
        scores = self.evaluate_batch([r.product for r in records])
        ranked = sorted(
            ((r, s) for r, s in zip(records, scores) if s is not None),
            key=lambda x: x[1], reverse=True,
        )
        return [r for r, _ in ranked[:100]]
```

Use `select_candidates` for anything judged from the building block alone, since
it prunes before synthesis and receives a lazy generator. Anything that needs the
product must go in `filter`.

For sequence optimizers, implement `init_search`, `ask`, and `tell`. `ask`
returns building block tuples, one per fragment; `tell` receives `(bb_tuple,
score)` pairs. A combination that yields several products is reported once, at
its best product's score.
