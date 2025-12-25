'''
    Monkey patch used rdkit functions to handle BuildingBlock objects.
'''
from functools import wraps

from rdkit import Chem
from rdkit.Chem import (
    rdChemReactions, rdFingerprintGenerator, 
    rdMolDescriptors, Descriptors
)

from healer.domain.building_block import BuildingBlock


# List of tuples with target and method names to patch
PATCH_TARGETS = [
    (Chem, [
        "MolToSmiles", "MolToSmarts", "SanitizeMol",
        "GetMolFrags"
    ]),
    (rdChemReactions.ChemicalReaction, [
        "RunReactants", "IsMoleculeReactant",
        "IsMoleculeProduct", "GetReactants",
        "GetProducts"
    ]),
    (rdFingerprintGenerator.FingerprintGenerator64, [
        "GetFingerprint", "GetFingerprints"
    ]),
    (Descriptors, [
        "MolWt", "NumHDonors", "NumHAcceptors",
        "TPSA", "NumRotatableBonds", "RingCount",
        "NumAromaticRings"
    ]),
    (rdMolDescriptors, [
        "CalcNumAtomStereoCenters", "_CalcMolWt"
    ]),
]

# Guard against double-patching on reimport
_PATCHED = set()


def _unwrap(obj):
    '''
        Recursively unwrap BuildingBlock -> raw Mol; leave everything else.
        Optimized with fast paths for common types.
    '''
    # Fast path: most common case
    if type(obj) is BuildingBlock:
        return obj._mol
    if type(obj) is list:
        return [_unwrap(x) for x in obj]
    if type(obj) is tuple:
        return tuple(_unwrap(x) for x in obj)
    # Slow path: subclasses of BuildingBlock (if any)
    if isinstance(obj, BuildingBlock):
        return obj._mol
    return obj

def _wrap(fn):
    '''
        Wrap any function or method so that all args/kwargs are unwrapped first.
        Works for single-Mol, sequences of Mol, reaction tuples, etc.
    '''
    @wraps(fn)
    def inner(*args, **kwargs):
        new_args = tuple(_unwrap(a) for a in args)
        new_kwargs = {k: _unwrap(v) for k, v in kwargs.items()}
        return fn(*new_args, **new_kwargs)
    return inner


for target, names in PATCH_TARGETS:
    for name in names:
        key = (id(target), name)
        if key in _PATCHED:
            continue
        orig = getattr(target, name, None)
        if orig:
            setattr(target, name, _wrap(orig))
            _PATCHED.add(key)

