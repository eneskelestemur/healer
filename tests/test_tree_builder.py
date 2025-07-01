'''
    Test cases for the retrosynthesis tree builder.
'''
import pytest
from rdkit import Chem

from healer.application.tree_builder import RetrosynthesisTree
from healer.domain.composition import CompositionPath
import healer.utils.utils as utils


PEN_SMILES = "CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O"

@pytest.fixture(scope="module")
def reactions():
    # Load and filter valid reactions
    rxns = utils.load_reactions_from_json('reactions/reactions.json')
    return [r for r in rxns if r.is_valid()]

@pytest.fixture(scope="module")
def penicillin():
    return Chem.MolFromSmiles(PEN_SMILES)

def test_zero_depth(penicillin, reactions):
    # At depth 0, no splitting; expect one path with the original molecule
    tree = RetrosynthesisTree(penicillin, reactions, max_depth=0)
    tree.build()
    paths = tree.get_composition_paths()
    assert len(paths) == 1
    path = paths[0]
    # No steps, one fragment equal to original
    assert path.steps == ()
    assert len(path.fragments) == 1
    assert Chem.MolToSmiles(path.fragments[0]) == PEN_SMILES

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
    # No splits => only the root molecule
    assert len(paths) == 1
    assert paths[0].fragments == (penicillin,)

def test_custom_fragments_constructor():
    # You can create a CompositionPath directly from fragments
    dummy = Chem.MolFromSmiles('CCO')
    cp = CompositionPath.from_fragments((dummy,))
    assert cp.steps is None
    assert len(cp.fragments) == 1
    assert Chem.MolToSmiles(cp.fragments[0]) == 'CCO'
