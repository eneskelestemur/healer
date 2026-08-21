"""Tests for the building block wrapper."""

import pickle

import pytest
from rdkit import Chem

from healer.domain.building_block import BuildingBlock


@pytest.fixture
def ethanol() -> BuildingBlock:
    return BuildingBlock(Chem.MolFromSmiles("CCO"))


class TestConstruction:
    def test_canonical_smiles_is_stored(self):
        assert BuildingBlock(Chem.MolFromSmiles("OCC")).get_smiles() == "CCO"

    def test_heavy_atom_count_is_captured(self, ethanol):
        assert ethanol.num_heavy_atoms == 3

    def test_sdf_properties_are_parsed(self):
        mol = Chem.MolFromSmiles("CCO")
        mol.SetProp("rxn_annotations", '{"amide-1": [0]}')
        mol.SetProp("URL", "https://example.com/bb")
        bb = BuildingBlock(mol)

        assert bb.get_parsed_prop("rxn_annotations") == {"amide-1": [0]}
        assert bb.get_parsed_prop("URL") == "https://example.com/bb"

    def test_missing_property_returns_empty_string(self, ethanol):
        assert ethanol.get_parsed_prop("nope") == ""


class TestIdentifier:
    @pytest.mark.parametrize("prop_name", ["id", "ID", "Id"])
    def test_the_identifier_is_found_whatever_its_case(self, prop_name):
        mol = Chem.MolFromSmiles("CCO")
        mol.SetProp(prop_name, "EN300-12345")

        assert BuildingBlock(mol).get_id() == "EN300-12345"

    def test_a_source_without_one_gives_an_empty_string(self, ethanol):
        assert ethanol.get_id() == ""

    def test_other_properties_are_not_mistaken_for_it(self):
        mol = Chem.MolFromSmiles("CCO")
        mol.SetProp("MDLNUMBER", "MFCD00003399")

        assert BuildingBlock(mol).get_id() == ""


class TestLazyMol:
    def test_mol_is_rebuilt_from_smiles(self, ethanol):
        assert Chem.MolToSmiles(ethanol.mol) == "CCO"

    def test_evict_frees_the_mol_and_it_comes_back(self, ethanol):
        _ = ethanol.mol
        assert ethanol._mol is not None

        ethanol.evict()
        assert ethanol._mol is None
        assert Chem.MolToSmiles(ethanol.mol) == "CCO"

    def test_atom_properties_survive_because_the_mol_is_kept(self):
        """SMILES cannot round-trip atom properties, so the original is retained."""
        mol = Chem.MolFromSmiles("c1ccccc1N")
        mol.GetAtomWithIdx(6).SetProp("_protected", "1")
        bb = BuildingBlock(mol)

        bb.evict()
        assert bb.mol.GetAtomWithIdx(6).HasProp("_protected")

    def test_plain_molecules_do_not_retain_the_original(self, ethanol):
        assert ethanol._mol_with_atom_props is None


class TestProperties:
    def test_set_prop_updates_both_views(self, ethanol):
        ethanol.SetProp("price", 42)
        assert ethanol.get_parsed_prop("price") == 42
        assert ethanol.mol.HasProp("price")

    def test_set_prop_accepts_strings_unchanged(self, ethanol):
        ethanol.SetProp("vendor", "Enamine")
        assert ethanol.get_parsed_prop("vendor") == "Enamine"

    def test_clear_prop_removes_from_both_views(self, ethanol):
        ethanol.SetProp("price", 42)
        ethanol.ClearProp("price")
        assert ethanol.get_parsed_prop("price") == ""
        assert not ethanol.mol.HasProp("price")


class TestDelegation:
    def test_rdkit_methods_reach_the_underlying_mol(self, ethanol):
        assert ethanol.GetNumAtoms() == 3
        assert ethanol.HasSubstructMatch(Chem.MolFromSmarts("[OH]"))

    def test_unknown_attributes_raise(self, ethanol):
        with pytest.raises(AttributeError):
            ethanol.definitely_not_a_method()


class TestPickling:
    def test_round_trips_through_pickle(self, ethanol):
        """Parallel synthesis ships building blocks to worker processes."""
        ethanol.SetProp("price", 42)
        restored = pickle.loads(pickle.dumps(ethanol))

        assert restored.get_smiles() == ethanol.get_smiles()
        assert restored.get_parsed_prop("price") == 42
        assert restored.num_heavy_atoms == ethanol.num_heavy_atoms

    def test_delegation_still_works_after_unpickling(self, ethanol):
        assert pickle.loads(pickle.dumps(ethanol)).GetNumAtoms() == 3


class TestHashing:
    def test_same_structure_hashes_alike(self):
        assert hash(BuildingBlock(Chem.MolFromSmiles("CCO"))) == hash(
            BuildingBlock(Chem.MolFromSmiles("OCC"))
        )

    def test_different_structures_differ(self):
        assert hash(BuildingBlock(Chem.MolFromSmiles("CCO"))) != hash(
            BuildingBlock(Chem.MolFromSmiles("CCC"))
        )
