# Reactions

HEALER ships 102 reaction templates in
[`healer/data/reactions/reactions.json`](../healer/data/reactions/reactions.json),
each defined in both directions: a retrosynthetic transform for splitting the
query, and a forward transform for reassembly.

Templates take **two reactants and give one product**. Anything else is rejected
at load time with a warning naming the template, so an invalid entry is skipped
rather than silently dropping out of the library.

## Selecting reactions

Filter by tag when constructing a HEALER:

```python
healer = MoleculeHEALER(reaction_tags=['amide coupling', 'N-arylation'])
healer = MoleculeHEALER(reaction_tags='all')
```

A template is included if it carries **any** of the requested tags. The full list
of 116 tags is in
[`reaction_tags.txt`](../healer/data/reactions/reaction_tags.txt), or:

```python
from healer.utils.utils import get_reaction_tags

print(get_reaction_tags())
```

Tags cover reaction families (`amide coupling`, `Suzuki`, `Buchwald-Hartwig`),
bond types (`C-N`, `C-C`, `C-O`), and functional groups (`azole`, `sulfonamide`).

The default selection is amide- and amine-focused, covering the most reliable and
widely used couplings. Widening to `'all'` gives more diversity at the cost of
run time and some less routine chemistry.

Reactions can be changed on an existing instance, which resets any compositions:

```python
healer.set_reactions(['Suzuki', 'amide coupling'])
```

## Entry format

```json
{
  "amide-1": {
    "description": "Amide coupling between a carboxylic acid and a primary amine.",
    "long_name": "Amide coupling",
    "syn_smarts": "[#6:101]-C(=O)-[OH].[#7;H2:102]>>[#6:101]-C(=O)-[#7:102]",
    "retro_smarts": "[#6:1]-C(=O)-[#7:2]>>[#6:1]-C(=O)-O.[#7:2]",
    "rhs_classes": ["carboxylic-acids", "amines-prim"],
    "tags": ["amide coupling", "amide", "C-N bond formation", "C-N"],
    "tier": 1
  }
}
```

| Field | Purpose |
|-------|---------|
| `syn_smarts` | Forward reaction, two reactants to one product |
| `retro_smarts` | Retrosynthetic split used to fragment the query |
| `rhs_classes` | Functional group classes of the two reactants |
| `tags` | Selection tags |
| `tier` | 1 for routine chemistry, 2 for more specialized |
| `display_smarts` | Optional annotated form used for rendering |
| `description`, `long_name` | Human-readable labels |

Atom maps matter. Templates are sanitized at load with RDKit's
`RemoveUnmappedReactantTemplates`, so a reactant whose atoms are largely unmapped
is stripped and the template is dropped as no longer 2 → 1. If a template you add
is skipped, check that both reactants carry mapped atoms.

## Building block annotations

Which templates a building block can serve, and at which reactant position, is
computed once during [preprocessing](building-blocks.md) and stored on the SDF
record. Adding or changing templates therefore requires reprocessing your
building block libraries.

## Contributing

Open a [GitHub issue](https://github.com/eneskelestemur/healer/issues) with the
proposed SMARTS, or contact enesk@email.unc.edu. Include `syn_smarts`,
`retro_smarts`, `rhs_classes`, and `tags`.
