"""Tests for the RDKit patch that unwraps BuildingBlock arguments."""

from rdkit import Chem
from rdkit.Chem import Descriptors

import healer.utils.rdkit_monkey_patch as patch
from healer.domain.building_block import BuildingBlock
from healer.utils.fingerprints import get_fingerprint_generator


class TestUnwrapping:
    def test_building_block_becomes_its_mol(self):
        bb = BuildingBlock(Chem.MolFromSmiles("CCO"))
        assert patch._unwrap(bb) is bb.mol

    def test_lists_and_tuples_are_unwrapped_recursively(self):
        bb = BuildingBlock(Chem.MolFromSmiles("CCO"))
        assert patch._unwrap([bb])[0] is bb.mol
        assert patch._unwrap((bb,))[0] is bb.mol

    def test_other_values_pass_through(self):
        assert patch._unwrap("CCO") == "CCO"
        assert patch._unwrap(7) == 7


class TestPatchedFunctions:
    def test_moltosmiles_accepts_a_building_block(self):
        assert Chem.MolToSmiles(BuildingBlock(Chem.MolFromSmiles("OCC"))) == "CCO"

    def test_descriptors_accept_a_building_block(self):
        bb = BuildingBlock(Chem.MolFromSmiles("CCO"))
        assert Descriptors.MolWt(bb) == Descriptors.MolWt(bb.mol)

    def test_fingerprints_accept_a_building_block(self):
        bb = BuildingBlock(Chem.MolFromSmiles("CCO"))
        gen = get_fingerprint_generator()
        assert list(gen.GetFingerprint(bb).GetOnBits()) == list(
            gen.GetFingerprint(bb.mol).GetOnBits()
        )

    def test_sanitize_accepts_a_building_block(self):
        bb = BuildingBlock(Chem.MolFromSmiles("CCO"))
        assert Chem.SanitizeMol(bb) == Chem.SanitizeFlags.SANITIZE_NONE


class TestIdempotence:
    def test_reimporting_does_not_stack_wrappers(self):
        """Workers re-import the patch; double wrapping would be silent overhead."""
        import importlib

        before = Chem.MolToSmiles
        importlib.reload(patch)
        assert Chem.MolToSmiles is before
