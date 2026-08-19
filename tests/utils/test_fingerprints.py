"""Tests for the fingerprint generator factory."""

import pytest
from rdkit import Chem

from healer.utils.fingerprints import get_fingerprint_generator


class TestFactory:
    def test_default_generator_produces_the_configured_size(self):
        fp = get_fingerprint_generator().GetFingerprint(Chem.MolFromSmiles("CCO"))
        assert fp.GetNumBits() == 2048

    def test_parameters_can_be_overridden(self):
        fp = get_fingerprint_generator(fpSize=512).GetFingerprint(
            Chem.MolFromSmiles("CCO")
        )
        assert fp.GetNumBits() == 512

    def test_unsupported_type_is_rejected(self):
        with pytest.raises(ValueError, match="Unsupported fingerprint type"):
            get_fingerprint_generator("not-a-fingerprint")


class TestDeterminism:
    def test_same_molecule_gives_the_same_fingerprint(self):
        gen = get_fingerprint_generator()
        first = gen.GetFingerprint(Chem.MolFromSmiles("CCO"))
        second = gen.GetFingerprint(Chem.MolFromSmiles("OCC"))
        assert list(first.GetOnBits()) == list(second.GetOnBits())

    def test_separate_generators_agree(self):
        mol = Chem.MolFromSmiles("c1ccccc1N")
        a = get_fingerprint_generator().GetFingerprint(mol)
        b = get_fingerprint_generator().GetFingerprint(mol)
        assert list(a.GetOnBits()) == list(b.GetOnBits())

    def test_chirality_is_distinguished(self):
        """The default generator includes chirality, so enantiomers differ."""
        gen = get_fingerprint_generator()
        left = gen.GetFingerprint(Chem.MolFromSmiles("C[C@H](N)C(=O)O"))
        right = gen.GetFingerprint(Chem.MolFromSmiles("C[C@@H](N)C(=O)O"))
        assert list(left.GetOnBits()) != list(right.GetOnBits())

    def test_batch_matches_single(self):
        gen = get_fingerprint_generator()
        mols = [Chem.MolFromSmiles(s) for s in ("CCO", "CCC")]
        batch = list(gen.GetFingerprints(mols))
        assert list(batch[0].GetOnBits()) == list(
            gen.GetFingerprint(mols[0]).GetOnBits()
        )
