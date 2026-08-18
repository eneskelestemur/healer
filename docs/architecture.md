# Architecture

## Layers

```
healer/
├── domain/         value objects and building block storage
├── application/    orchestration
├── web/            FastAPI app, Celery worker
├── scripts/        preprocessing
└── utils/          fingerprints, similarity, drawing
```

`domain` knows nothing about `application`; the interfaces depend on both.

### domain

| Module | Holds |
|--------|-------|
| `reaction_template.py` | `ReactionTemplate21` — a 2 → 1 template, forward and retro |
| `building_block.py` | `BuildingBlock` — an SDF record with parsed properties |
| `bb_repository.py` | `BBRepository` — loading, indexing, and caching a library |
| `composition.py` | `CompositionPath`, `CompositionWithBBs` — fragments and their pools |
| `retro_step.py` | `RetroStep` — one retrosynthetic split |
| `enumeration_record.py` | `EnumerationRecord` — one product and its provenance |

### application

| Module | Holds |
|--------|-------|
| `healer.py` | `_BaseHEALER` and the three concrete classes |
| `tree_builder.py` | `RetrosynthesisTree` — recursive decomposition |
| `optimizers.py` | Optimizer interfaces and implementations |

## The enumeration pipeline

```
query molecule
   │
   ├─ _process_query_mol()        → CompositionPath objects (fragments)
   │
   ├─ _process_building_blocks()  → CompositionWithBBs (fragments + BB pools)
   │
   └─ _enumerate_*()              → EnumerationRecord objects (products)
```

`_BaseHEALER` is a template method: it owns the assembly loop, while subclasses
supply the two `_process_*` steps. That is the whole difference between the three
modes — `MoleculeHEALER` builds a retrosynthesis tree, `FragmentHEALER` takes the
components as given, `SiteHEALER` produces a single "fragment" that is the intact
query.

Assembly seeds one record per building block in the first pool, then for each
subsequent pool builds `(record, building block, reaction)` candidates and
applies them. Candidates are generated lazily, so a wide stage is not
materialized before it is pruned.

## Two enumeration paths

`enumerate()` dispatches on the optimizer:

| Optimizer | Path |
|-----------|------|
| `None` or a `BaseStagewiseOptimizer` | `_enumerate_stagewise` |
| A `BaseSequenceOptimizer` | `_enumerate_sequence` |

Unoptimized runs use `PassthroughOptimizer`, whose hooks are the identity and
whose scores are all `None`, so exhaustive and stagewise-guided enumeration share
one code path rather than duplicating the loop.

The sequence path is separate because the optimizer drives the outer loop through
`ask`/`tell` rather than being called at fixed points.

## Building block indexing

`BBRepository.load()` reads the SDF once, computes a Morgan fingerprint per
record, and builds an inverted index from reaction name to compatible building
block indices — populated from the `rxn_annotations` written during
preprocessing. Filtering by reaction is then a set union over that index.

Repositories are cached per resolved path at module level, so several HEALER
instances on the same library share one copy.

`BuildingBlock` keeps its SMILES and rebuilds the RDKit `Mol` on demand;
`evict()` drops the cached `Mol`. Records whose atoms carry properties that
SMILES cannot round-trip, such as the protection flags `SiteHEALER` sets, keep
their original `Mol`.

## Fragment matching

`MoleculeHEALER._process_building_blocks` scores every fragment against every
building block in chunks, as a `float16` matrix. The Tversky similarity is
weighted by a size penalty, so building blocks much larger than the fragment rank
lower. Pools are then taken as the top `max_bbs_per_frag` per fragment, or
everything above `sim_threshold`.

Similarities are retained alongside the pools only in the capped case, where
their size is bounded; under a threshold the pools are unbounded and the matrix
is discarded.

## Parallelism

`n_jobs` parallelizes the synthesis step with joblib's loky backend, dispatching
candidates in chunks of 1000 to a top-level picklable worker.

Worker processes re-import the RDKit patch described below, since monkey patches
do not survive process boundaries.

## The RDKit patch

`BuildingBlock` is not an RDKit `Mol`, and `Mol` cannot practically be subclassed
through its Boost bindings. `utils/rdkit_monkey_patch.py` therefore wraps the
RDKit functions HEALER calls so they unwrap `BuildingBlock` arguments
transparently.

It is applied on `import healer`, and independently by the CLI, web app, Celery
worker, and joblib workers. Anything that reaches RDKit with a `BuildingBlock`
before the patch is applied will fail with a Boost converter error.
