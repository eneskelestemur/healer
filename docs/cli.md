# Command Line

```bash
healer molecule <input> [options]   # retrosynthetic enumeration
healer fragment <input> [options]   # enumeration from pre-split fragments
healer site <input> [options]       # site-specific enumeration
healer view <smiles>                # render with atom indices
```

Run `healer <command> --help` for the full list. Optimizers are Python-only and
are not exposed on the command line — see [Guided Enumeration](guided-enumeration.md).

## Input

`input` is a SMILES string or a file. `.csv` files are read from the `smiles`
column unless `--column` says otherwise; `.smi` and `.sdf` are also accepted.
Invalid molecules are logged and skipped.

## Common options

| Option | Description | Default |
|--------|-------------|---------|
| `-o`, `--output` | Output CSV | `healer_results.csv` |
| `--column` | SMILES column for CSV input | `smiles` |
| `--config` | JSON config file of the same options | — |
| `--bb-source` | Named source or path | `US_stock` |
| `--reactions` | Comma-separated tags, or `all` | `all` |
| `--shuffle` | Shuffle building block order | off |
| `--max-evals` | Max reaction attempts per composition | — |
| `--max-products` | Max products per composition | — |
| `--max-total` | Max products in total | — |
| `--similarity` | Add similarity to the query | off |
| `--properties` | Add the property profile | off |
| `--n-jobs` | Parallel workers, `-1` for all cores | `1` |
| `-v`, `-vv` | Info, then debug logging | info |

The three limits are approximate — see [HEALER Classes](healer-classes.md).

## `molecule`

| Option | Description | Default |
|--------|-------------|---------|
| `--sim-threshold` | Minimum fragment-to-BB similarity | `0.5` |
| `--max-bbs-per-frag` | Cap BBs per fragment, `-1` for unlimited | `-1` |
| `--n-compositions` | Compositions to consider | `10` |
| `--retro-depth` | Retrosynthesis tree depth | `1` |
| `--min-frag-size` | Minimum fragment size in heavy atoms | `3` |
| `--randomize` | Shuffle composition order | off |
| `--seed` | Random seed | `-1` |

```bash
healer molecule inputs.csv --bb-source US_stock --max-bbs-per-frag 100 \
    --max-total 5000 --properties --n-jobs -1 -o analogs.csv
```

## `fragment`

Takes `--sim-threshold` and `--max-bbs-per-frag`. Input must be a
dot-separated SMILES:

```bash
healer fragment "c1ccccc1N.CC(=O)O" --bb-source test
```

## `site`

| Option | Description |
|--------|-------------|
| `--reactive-sites` | Atom indices as a JSON list, e.g. `"[1,2,5]"` |
| `--rules` | Property windows, `"MW:0:500,HBD:0:5"` |
| `--struct-rules` | Comma-separated SMARTS a BB must match |

```bash
healer view "c1ccccc1N"     # find the atom index first
healer site "c1ccccc1N" --reactive-sites "[6]" --rules "MW:0:300"
```

## `view`

Renders a molecule with atom indices, opening a browser or writing SVG with
`-o out.svg`. Use it to pick indices for `site`.

## Config files

Any option can live in JSON, with command line flags taking precedence:

```json
{
  "bb_source": "US_stock",
  "reactions": "amide coupling,N-arylation",
  "max_bbs_per_frag": 100,
  "max_total": 5000
}
```

```bash
healer molecule inputs.csv --config run.json -o analogs.csv
```

## Batches

One HEALER is constructed and reused for every molecule in the file, so the
building block library is loaded once. Results are appended as each molecule
finishes, and a failure on one molecule is logged without stopping the run.
