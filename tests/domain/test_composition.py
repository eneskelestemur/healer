"""Tests for composition data structures."""

import pytest
from rdkit import Chem

from healer.domain.building_block import BuildingBlock
from healer.domain.composition import CompositionPath, CompositionWithBBs


def mols(*smiles: str) -> tuple[Chem.Mol, ...]:
    return tuple(Chem.MolFromSmiles(s) for s in smiles)


class TestConstruction:
    def test_fragments_can_be_given_directly(self):
        path = CompositionPath(fragments=mols("CCO", "c1ccccc1"))
        assert len(path) == 2

    def test_fragments_are_flattened_from_steps(self, amide_reaction):
        from healer.domain.retro_step import RetroStep

        left, right = mols("CC(=O)O", "NCc1ccccc1")
        step = RetroStep(
            product=Chem.MolFromSmiles("CC(=O)NCc1ccccc1"),
            reaction=amide_reaction,
            reactants=(left, right),
        )
        path = CompositionPath(steps=(step,))
        assert len(path.fragments) == 2

    def test_neither_steps_nor_fragments_is_an_error(self):
        with pytest.raises(ValueError, match="requires at least steps or fragments"):
            CompositionPath()

    def test_from_fragments_carries_no_reaction_history(self):
        path = CompositionPath.from_fragments(mols("CCO"))
        assert path.steps is None
        assert len(path) == 1


class TestIdentity:
    def test_same_fragments_compare_equal(self):
        assert CompositionPath(fragments=mols("CCO", "c1ccccc1")) == CompositionPath(
            fragments=mols("OCC", "c1ccccc1")
        )

    def test_fragment_order_matters(self):
        assert CompositionPath(fragments=mols("CCO", "c1ccccc1")) != CompositionPath(
            fragments=mols("c1ccccc1", "CCO")
        )

    def test_equal_paths_hash_alike(self):
        a = CompositionPath(fragments=mols("CCO", "c1ccccc1"))
        b = CompositionPath(fragments=mols("OCC", "c1ccccc1"))
        assert len({a, b}) == 1

    def test_comparison_with_other_types_is_not_implemented(self):
        assert CompositionPath(fragments=mols("CCO")).__eq__("nope") is NotImplemented


class TestCompositionWithBBs:
    def test_similarities_default_to_absent(self):
        comp = CompositionWithBBs(
            comp=CompositionPath(fragments=mols("CCO")),
            fragment_bbs=([BuildingBlock(Chem.MolFromSmiles("CCO"))],),
        )
        assert comp.fragment_sims is None

    def test_similarities_align_with_the_pools(self):
        pool = [BuildingBlock(Chem.MolFromSmiles(s)) for s in ("CCO", "CCC")]
        comp = CompositionWithBBs(
            comp=CompositionPath(fragments=mols("CCO")),
            fragment_bbs=(pool,),
            fragment_sims=([0.9, 0.4],),
        )
        assert len(comp.fragment_sims[0]) == len(comp.fragment_bbs[0])
