'''
    Test cases for the retrosynthesis tree builder.
'''
from pathlib import Path

import pytest
from rdkit import Chem

from healer.application.tree_builder import RetrosynthesisTree
from healer.domain.composition import CompositionPath
import healer.utils.utils as utils

HEALER_PKG = Path(__file__).parent.parent / 'healer'
REACTIONS_PATH = HEALER_PKG / 'data' / 'reactions' / 'reactions.json'
PEN_SMILES = "CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O"

@pytest.fixture(scope="module")
def reactions():
    # Load and filter valid reactions
    rxns = utils.load_reactions_from_json(REACTIONS_PATH)
    return [r for r in rxns if r.is_valid()]

@pytest.fixture(scope="module")
def penicillin():
    return Chem.MolFromSmiles(PEN_SMILES)

def test_zero_depth(penicillin, reactions):
    # At depth 0, no splitting; don't return the root molecule as a path
    # since it has no steps.
    tree = RetrosynthesisTree(penicillin, reactions, max_depth=0)
    tree.build()
    paths = tree.get_composition_paths()
    assert len(paths) == 0, "Expected no paths at depth 0"

def test_depth1_splits(penicillin, reactions):
    # Depth 1: expect at least one split into two fragments
    tree = RetrosynthesisTree(penicillin, reactions, max_depth=1)
    tree.build()
    paths = tree.get_composition_paths()
    # There should be at least one path splitting into 2 fragments
    assert any(len(p.fragments) == 2 for p in paths)
    # All paths should have at most one step and up to two fragments
    for path in paths:
        assert len(path.steps) <= 1
        assert 1 <= len(path.fragments) <= 2

def test_min_heavy_atoms_filter(penicillin, reactions):
    # Using a very high min_heavy_atoms should prevent any splits
    tree = RetrosynthesisTree(penicillin, reactions, max_depth=2, min_heavy_atoms=100)
    tree.build()
    paths = tree.get_composition_paths()
    # No splits
    assert len(paths) == 0

def test_custom_fragments_constructor():
    # You can create a CompositionPath directly from fragments
    dummy = Chem.MolFromSmiles('CCO')
    cp = CompositionPath.from_fragments((dummy,))
    assert cp.steps is None
    assert len(cp.fragments) == 1
    assert Chem.MolToSmiles(cp.fragments[0]) == 'CCO'
