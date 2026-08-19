"""
Test cases for the retrosynthesis tree builder.
"""

import logging
from pathlib import Path

import pytest
from rdkit import Chem

import healer.utils.utils as utils
from healer.application.tree_builder import RetrosynthesisTree
from healer.domain.composition import CompositionPath

HEALER_PKG = Path(__file__).parent.parent / "healer"
REACTIONS_PATH = HEALER_PKG / "data" / "reactions" / "reactions.json"
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


@pytest.fixture(scope="module")
def unbounded_tree(penicillin, reactions):
    """The full depth-2 tree, used to size budgets that are certain to bite."""
    tree = RetrosynthesisTree(
        penicillin, reactions, max_depth=2, min_heavy_atoms=1, max_nodes=None
    )
    tree.build()
    return tree


def test_node_budget_truncates_the_tree(penicillin, reactions, unbounded_tree, caplog):
    """A budget below the full size stops expansion early and says so."""
    budget = max(3, unbounded_tree.n_nodes // 2)

    with caplog.at_level(logging.WARNING, logger="healer.application.tree_builder"):
        tree = RetrosynthesisTree(
            penicillin, reactions, max_depth=2, min_heavy_atoms=1, max_nodes=budget
        )
        tree.build()

    assert tree.budget_exhausted
    assert tree.n_nodes < unbounded_tree.n_nodes
    assert "node budget" in caplog.text


def test_truncated_tree_still_yields_usable_paths(
    penicillin, reactions, unbounded_tree
):
    """Truncation must leave a smaller valid tree, not a broken one."""
    budget = max(3, unbounded_tree.n_nodes // 2)
    tree = RetrosynthesisTree(
        penicillin, reactions, max_depth=2, min_heavy_atoms=1, max_nodes=budget
    )
    tree.build()
    paths = tree.get_composition_paths()

    assert paths
    assert len(paths) <= len(unbounded_tree.get_composition_paths())
    for path in paths:
        assert len(path.fragments) >= 2
        assert all(frag is not None for frag in path.fragments)


def test_budget_none_does_not_truncate(penicillin, reactions):
    tree = RetrosynthesisTree(penicillin, reactions, max_depth=2, max_nodes=None)
    tree.build()
    assert tree.budget_exhausted is False


def test_generous_budget_matches_unbounded_tree(penicillin, reactions):
    """The default budget must not change results for an ordinary query."""
    bounded = RetrosynthesisTree(penicillin, reactions, max_depth=2, max_nodes=10000)
    bounded.build()
    unbounded = RetrosynthesisTree(penicillin, reactions, max_depth=2, max_nodes=None)
    unbounded.build()

    assert bounded.budget_exhausted is False
    assert len(bounded.get_composition_paths()) == len(
        unbounded.get_composition_paths()
    )


def test_custom_fragments_constructor():
    # You can create a CompositionPath directly from fragments
    dummy = Chem.MolFromSmiles("CCO")
    cp = CompositionPath.from_fragments((dummy,))
    assert cp.steps is None
    assert len(cp.fragments) == 1
    assert Chem.MolToSmiles(cp.fragments[0]) == "CCO"
