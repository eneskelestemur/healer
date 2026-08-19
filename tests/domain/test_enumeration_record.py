"""Tests for the enumeration record."""

from rdkit import Chem

from healer.domain.building_block import BuildingBlock
from healer.domain.enumeration_record import EnumerationRecord


class TestDefaults:
    def test_only_the_product_is_required(self):
        record = EnumerationRecord(product=Chem.MolFromSmiles("CCO"))
        assert record.bbs == []
        assert record.reaction_names == []
        assert record.props == {}
        assert record.origin is None

    def test_defaults_are_not_shared_between_records(self):
        first = EnumerationRecord(product=Chem.MolFromSmiles("CCO"))
        second = EnumerationRecord(product=Chem.MolFromSmiles("CCC"))
        first.bbs.append(BuildingBlock(Chem.MolFromSmiles("CCO")))
        first.props["score"] = 1.0

        assert second.bbs == []
        assert second.props == {}


class TestProvenance:
    def test_origin_tracks_the_proposed_combination(self):
        record = EnumerationRecord(product=Chem.MolFromSmiles("CCO"), origin=3)
        assert record.origin == 3

    def test_building_blocks_and_reactions_are_recorded_in_order(self):
        bbs = [BuildingBlock(Chem.MolFromSmiles(s)) for s in ("CC(=O)O", "NCc1ccccc1")]
        record = EnumerationRecord(
            product=Chem.MolFromSmiles("CC(=O)NCc1ccccc1"),
            bbs=bbs,
            reaction_names=["amide-1"],
        )
        assert [bb.get_smiles() for bb in record.bbs] == ["CC(=O)O", "NCc1ccccc1"]
        assert record.reaction_names == ["amide-1"]
