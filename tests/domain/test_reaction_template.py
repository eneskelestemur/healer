"""Tests for the reaction template wrapper."""

import logging

import pytest
from rdkit import Chem

from healer.domain.building_block import BuildingBlock
from healer.domain.reaction_template import ReactionTemplate21

AMIDE_SYN = "[#6:101]-C(=O)-[OH].[#7;H2:102]>>[#6:101]-C(=O)-[#7:102]"
AMIDE_RETRO = "[#6:1]-C(=O)-[#7:2]>>[#6:1]-C(=O)-O.[#7:2]"


def make_amide() -> ReactionTemplate21:
    return ReactionTemplate21(
        name="amide-test",
        reaction_smarts=AMIDE_SYN,
        retro_smarts=AMIDE_RETRO,
        rhs_classes=["carboxylic-acids", "amines-prim"],
        tags=["amide coupling", "C-N"],
    )


class TestValidation:
    def test_two_to_one_template_is_valid(self):
        assert make_amide().is_valid()

    def test_one_to_one_template_is_rejected(self, caplog):
        with caplog.at_level(logging.WARNING, logger="healer.domain.reaction_template"):
            rxn = ReactionTemplate21(
                name="not-two-to-one",
                reaction_smarts="[c:1]Br>>[c:1]C",
                retro_smarts="[c:1]C>>[c:1]Br",
            )
        assert not rxn.is_valid()
        assert "not-two-to-one" in caplog.text

    def test_malformed_smarts_is_reported_not_raised(self, caplog):
        with caplog.at_level(logging.WARNING, logger="healer.domain.reaction_template"):
            rxn = ReactionTemplate21(
                name="broken", reaction_smarts="this is not smarts", retro_smarts=""
            )
        assert not rxn.is_valid()
        assert "broken" in caplog.text

    def test_invalid_templates_do_not_stop_a_library_load(self, tmp_path):
        """One bad entry must not prevent the rest of the file from loading."""
        import json

        import healer.utils.utils as utils

        path = tmp_path / "rxns.json"
        path.write_text(
            json.dumps(
                {
                    "good": {"syn_smarts": AMIDE_SYN, "retro_smarts": AMIDE_RETRO},
                    "bad": {"syn_smarts": "nonsense", "retro_smarts": ""},
                }
            )
        )
        loaded = utils.load_reactions_from_json(str(path))
        assert len(loaded) == 2
        assert [r.name for r in loaded if r.is_valid()] == ["good"]


class TestSmartsAccessors:
    def test_get_reaction_smarts_returns_the_constructor_value(self):
        assert make_amide().get_reaction_smarts() == AMIDE_SYN

    def test_set_reaction_smarts_updates_the_accessor(self):
        rxn = make_amide()
        rxn.set_reaction_smarts(AMIDE_SYN)
        assert rxn.get_reaction_smarts() == AMIDE_SYN
        assert rxn.is_valid()

    def test_set_reaction_smarts_invalidates_a_bad_template(self):
        rxn = make_amide()
        rxn.set_reaction_smarts("[c:1]Br>>[c:1]C")
        assert not rxn.is_valid()

    def test_reactants_and_products_are_reported(self):
        rxn = make_amide()
        assert len(rxn.get_reactants()) == 2
        assert len(rxn.get_products()) == 1
        assert len(rxn.get_reactants_smarts()) == 2
        assert len(rxn.get_products_smiles()) == 1

    def test_reactants_sort_by_weight_when_asked(self):
        rxn = make_amide()
        from rdkit.Chem import Descriptors

        weights = [Descriptors.MolWt(r) for r in rxn.get_reactants(sort_by_mw=True)]
        assert weights == sorted(weights, reverse=True)


class TestFromReactionJson:
    def test_syn_smarts_key_is_accepted(self):
        rxn = ReactionTemplate21.from_reaction_json(
            "json-amide", {"syn_smarts": AMIDE_SYN, "retro_smarts": AMIDE_RETRO}
        )
        assert rxn.name == "json-amide"
        assert rxn.is_valid()

    def test_unknown_keys_are_ignored(self):
        rxn = ReactionTemplate21.from_reaction_json(
            "extra",
            {
                "syn_smarts": AMIDE_SYN,
                "retro_smarts": AMIDE_RETRO,
                "tier": 1,
                "not_a_parameter": "ignored",
            },
        )
        assert rxn.tier == 1
        assert not hasattr(rxn, "not_a_parameter")


class TestRunRetro:
    def test_splits_an_amide_into_acid_and_amine(self):
        rxn = make_amide()
        product = Chem.MolFromSmiles("CC(=O)NCc1ccccc1")
        pairs = rxn.run_retro(product)

        assert pairs, "the amide bond should be found"
        smiles = {tuple(sorted(Chem.MolToSmiles(m) for m in pair)) for pair in pairs}
        assert any("CC(=O)O" in pair for pair in smiles)

    def test_returns_nothing_when_the_motif_is_absent(self):
        assert make_amide().run_retro(Chem.MolFromSmiles("c1ccccc1")) == []


class TestRunSyn:
    def test_couples_an_acid_and_an_amine(self):
        rxn = make_amide()
        acid = self._annotated("CC(=O)O", rxn.name, [0])
        amine = self._annotated("NCc1ccccc1", rxn.name, [1])

        products = rxn.run_syn(acid, amine)
        assert products
        for product in products:
            Chem.SanitizeMol(product)
        assert any("C(=O)N" in Chem.MolToSmiles(p) for p in products)

    def test_requires_exactly_two_reactants(self):
        rxn = make_amide()
        with pytest.raises(AssertionError):
            rxn.run_syn(Chem.MolFromSmiles("CC(=O)O"))

    def test_annotations_decide_the_reactant_order(self):
        """
        Both building blocks claiming the same slot leaves no valid ordering, so
        no products are produced.
        """
        rxn = make_amide()
        same_slot_a = self._annotated("CC(=O)O", rxn.name, [0])
        same_slot_b = self._annotated("NCc1ccccc1", rxn.name, [0])
        assert rxn.run_syn(same_slot_a, same_slot_b) == []

    def test_plain_mols_fall_back_to_substructure_matching(self):
        rxn = make_amide()
        products = rxn.run_syn(
            Chem.MolFromSmiles("CC(=O)O"), Chem.MolFromSmiles("NCc1ccccc1")
        )
        assert products

    @staticmethod
    def _annotated(smiles: str, rxn_name: str, positions: list[int]) -> BuildingBlock:
        mol = Chem.MolFromSmiles(smiles)
        bb = BuildingBlock(mol)
        bb.SetProp("rxn_annotations", {rxn_name: positions})
        return bb


class TestMembership:
    def test_is_reactant_and_is_product(self):
        rxn = make_amide()
        assert rxn.is_reactant(Chem.MolFromSmiles("CC(=O)O"))
        assert not rxn.is_reactant(Chem.MolFromSmiles("c1ccccc1"))
        assert rxn.is_product(Chem.MolFromSmiles("CC(=O)NC"))

    def test_get_reactant_index_locates_the_slot(self):
        rxn = make_amide()
        acid_slots = rxn.get_reactant_index(Chem.MolFromSmiles("CC(=O)O"))
        amine_slots = rxn.get_reactant_index(Chem.MolFromSmiles("NCc1ccccc1"))
        assert acid_slots and amine_slots
        assert set(acid_slots) != set(amine_slots)


class TestDunders:
    def test_str_and_repr_are_the_name(self):
        rxn = make_amide()
        assert str(rxn) == "amide-test"
        assert repr(rxn) == "amide-test"

    def test_equal_smarts_hash_alike(self):
        assert hash(make_amide()) == hash(make_amide())
