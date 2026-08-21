# Building Blocks

HEALER assembles products from a library of purchasable building blocks. Any SDF
file works once it has been preprocessed.

## Preprocessing

Raw SDF files must be annotated before use. Preprocessing removes salts and small
disconnected fragments, then records which reaction templates each building block
can participate in and at which reactant position. Storing these annotations up
front is what keeps enumeration fast — at run time the reactant role is a
dictionary lookup rather than a substructure match.

```bash
preprocess-bb ~/Downloads/Enamine_BBs.zip -o ~/.healer/buildingblocks/ --verbose
```

| Option | Description |
|--------|-------------|
| `input_file` | SDF or ZIP to process |
| `-o`, `--output-dir` | Output directory (defaults to the input's directory) |
| `-w`, `--workers` | Worker processes, `-1` for all cores (default: 1) |
| `-v`, `--verbose` | Verbose output |

The result is written alongside the input with a `_processed` suffix. Only
processed files can be used as a `bb_source`.

## Where HEALER looks

Set `HEALER_DATA_DIR` to the directory holding your libraries:

```bash
export HEALER_DATA_DIR=~/.healer
```

HEALER then resolves the named sources below against
`$HEALER_DATA_DIR/buildingblocks/`.

## Named sources

| Name | Expected subdirectory |
|------|----------------------|
| `US_stock` | `Enamine_Rush-Delivery_Building_Blocks-US/` |
| `EU_stock` | `Enamine_Rush-Delivery_Building_Blocks-EU/` |
| `Global_stock` | `Enamine_Building_Blocks_Stock/` |
| `test` | The bundled 100-compound set |

Each directory is globbed for `*_processed.sdf`, and the most recently modified
match wins — so a dated Enamine download can be dropped in without renaming.

Catalogs are available from
[Enamine](https://enamine.net/building-blocks/building-blocks-catalog).

## Custom libraries

Pass a path instead of a name:

```python
from healer import MoleculeHEALER

h = MoleculeHEALER(bb_source='/path/to/my_library_processed.sdf')
```

Relative paths resolve against `$HEALER_DATA_DIR/buildingblocks/`, and glob
patterns are accepted.

SDF properties carried on a building block are picked up by name. An `id` field
(matched case-insensitively) is reported in the `BBID` result columns, and a
`URL` field in the `URL` columns; both are left empty when absent. Adding an
`id` to an internal library is enough to trace every analog back to a catalog
entry.

## Sharing a library across instances

Loading a large SDF takes time and memory, so repositories are cached per
resolved path. Constructing several HEALERs with the same `bb_source` reuses one
copy automatically:

```python
from healer import MoleculeHEALER, SiteHEALER, get_repository

repo = get_repository('US_stock')
repo.load()

mol_healer  = MoleculeHEALER(bb_repository=repo)
site_healer = SiteHEALER(bb_repository=repo)
```

`clear_repository_cache()` releases the memory when you are done.

## Memory

The whole library is loaded into memory with a fingerprint per building block.
Budget roughly 1–2 GB per million building blocks. If you are working on a
constrained machine, use a stock subset rather than the full catalog.
