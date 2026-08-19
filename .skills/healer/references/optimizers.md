# Guided enumeration in depth

## Budgets

`init_search` runs once per composition, so each composition is searched from a
cold start:

```
total evaluations ≈ n_compositions × max_evals_per_comp
```

Size the two together. With `n_compositions=10` and `max_evals_per_comp=200`,
expect ~2000 scoring calls, not 200.

## Choosing an optimizer

**`BeamSearchOptimizer(beam_width=100, target_fn=...)`** scores every assembled
intermediate at every stage and keeps the best `beam_width`. It calls the scorer
far more than the sequence optimizers, so the scorer must be cheap.

**`GeneticAlgorithmOptimizer(population_size=50, ...)`** evolves building block
combinations; one ask/tell round is one generation. Genes are block indices, one
per fragment, so a two-fragment query is a two-gene problem — with few genes,
raise `population_size` rather than the mutation rate, since the mutation count
is already clamped to one. If it converges early, lower `keep_elitism`.

**`BayesianSequenceOptimizer(batch_size=10, ...)`** fits a surrogate over MORDRED
descriptors. It enumerates the full product of the building block pools up front,
so pools are capped at `max_domain_per_frag=200` by default. **Set
`max_bbs_per_frag` on the HEALER** so that cap keeps the blocks most similar to
each fragment; otherwise it truncates arbitrarily and warns.

## Scoring functions

```python
def batch_score(mols):
    smiles = [Chem.MolToSmiles(m) for m in mols]
    return my_model.predict(smiles)          # None for anything unscorable

opt = GeneticAlgorithmOptimizer(batch_target_fn=batch_score)
```

`batch_target_fn` takes priority over `target_fn`. Failures are logged and the
molecule is dropped from the optimizer's feedback rather than ending the run.

A combination that yields several products is reported to the optimizer once, at
its best product's score.

## Custom optimizers

Stagewise optimizers have two hooks. `select_candidates` runs **before** any
reaction and receives a lazy generator, so building-block-level heuristics prune
without paying for synthesis. `filter` runs after, on assembled products.

```python
from healer import BaseStagewiseOptimizer

class CheapAndGood(BaseStagewiseOptimizer):
    def select_candidates(self, candidates, depth):
        for rec, bb, rxn in candidates:
            if bb.num_heavy_atoms <= 20:
                yield rec, bb, rxn

    def filter(self, records, depth):
        scores = self.evaluate_batch([r.product for r in records])
        ranked = sorted(
            ((r, s) for r, s in zip(records, scores) if s is not None),
            key=lambda x: x[1],
            reverse=True,
        )
        return [r for r, _ in ranked[:100]]
```

Anything requiring the product must go in `filter` — you cannot score a molecule
you have not built.

Sequence optimizers implement `init_search`, `ask`, and `tell`. `ask` returns
building block tuples, one block per fragment; returning an empty list ends the
search for that composition. Raise `OptimizerError` for a genuine failure —
enumeration logs a warning, abandons that composition, and continues.
