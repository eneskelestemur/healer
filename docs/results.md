# Results

## Retrieving results

```python
df = healer.get_results(calc_similarity=True, calc_properties=True)
healer.save_results('analogs.csv')          # same arguments, writes CSV
records = healer.get_results(as_dict=True)  # list of dicts
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `as_dict` | Return a list of dicts instead of a DataFrame | `False` |
| `calc_similarity` | Add Tanimoto similarity to the query | `True` |
| `calc_properties` | Add the property profile | `True` |
| `skip_cns_mpo` | Skip CNS MPO when profiling | `True` |

Both calculations cost time on large result sets; turn them off for a quick look.

## Columns

| Column | Description |
|--------|-------------|
| `ID` | `HEAL_000000`-style identifier |
| `Product` | Canonical SMILES of the analog |
| `BB1`, `BB2`, … | Building blocks used |
| `Reaction1_name`, … | Reaction template applied at each step |
| `URL1`, `URL2`, … | Supplier URL per building block, when present in the SDF |
| `Similarity_to_query` | Tanimoto similarity, with `calc_similarity=True` |
| `optimization_score` | Objective value, when an optimizer was used |

The number of `BB`/`Reaction` columns comes from the mode: `2 ** retro_tree_depth`
for `MoleculeHEALER`, 2 for `SiteHEALER`, and the component count for
`FragmentHEALER`. Unused slots are empty strings.

With `calc_properties=True`, [prop-profiler](https://pypi.org/project/prop-profiler/)
adds `mw`, `logp`, `hba`, `hbd`, `tpsa`, `num_rotatable_bonds`, `fsp3`, `qed`,
`esol_mg/L`, `stoplight_score`, and `stoplight_color`.

## The query molecule is the first row

Row one is the query itself, with empty building block columns and a
`Similarity_to_query` of 1.001 so it sorts to the top. It gives a baseline to
compare analogs against; drop it with `df[df['BB1'] != '']`.

## Deduplication

Rows are deduplicated on everything except the `URL` and `Reaction` columns, so
the same product reached by two different routes appears once, keeping the first
route found.

## Raw records

`healer.enumerated_molecules` holds the underlying `EnumerationRecord` objects
before tabulation:

| Attribute | Description |
|-----------|-------------|
| `product` | RDKit `Mol` |
| `bbs` | `BuildingBlock` objects used |
| `reaction_names` | Templates applied, in order |
| `props` | Extra values, including `optimization_score` |
| `origin` | Index of the proposed combination, for sequence optimizers |

Useful when you want the `Mol` objects directly:

```python
best = max(
    (r for r in healer.enumerated_molecules if 'optimization_score' in r.props),
    key=lambda r: r.props['optimization_score'],
)
```
