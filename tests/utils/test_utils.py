"""Tests for the shared helper functions."""

import json

import numpy as np
import pytest
from rdkit import Chem

import healer.utils.utils as utils
from healer.utils.fingerprints import get_fingerprint_generator

AMIDE_SYN = "[#6:101]-C(=O)-[OH].[#7;H2:102]>>[#6:101]-C(=O)-[#7:102]"
AMIDE_RETRO = "[#6:1]-C(=O)-[#7:2]>>[#6:1]-C(=O)-O.[#7:2]"


@pytest.fixture(scope="module")
def fps():
    gen = get_fingerprint_generator()
    mols = [Chem.MolFromSmiles(s) for s in ("CCO", "CCC", "c1ccccc1", "CCOCC")]
    return [gen.GetFingerprint(m) for m in mols]


class TestReactionLoading:
    def test_every_shipped_template_is_returned(self, all_reactions):
        assert len(all_reactions) > 50
        assert all(r.is_valid() for r in all_reactions)

    def test_results_are_cached_per_path(self, tmp_path):
        path = tmp_path / "rxns.json"
        path.write_text(
            json.dumps({"a": {"syn_smarts": AMIDE_SYN, "retro_smarts": AMIDE_RETRO}})
        )

        first = utils.load_reactions_from_json(str(path))
        second = utils.load_reactions_from_json(str(path))
        assert first is second

    def test_cache_can_be_cleared(self, tmp_path):
        path = tmp_path / "rxns.json"
        path.write_text(
            json.dumps({"a": {"syn_smarts": AMIDE_SYN, "retro_smarts": AMIDE_RETRO}})
        )

        first = utils.load_reactions_from_json(str(path))
        utils.load_reactions_from_json.cache_clear()
        assert utils.load_reactions_from_json(str(path)) is not first


class TestReactionTags:
    def test_tags_are_read_from_the_packaged_file(self):
        tags = utils.get_reaction_tags()
        assert "amide coupling" in tags
        assert all(tag.strip() == tag for tag in tags)


class TestSanitizeMol:
    def test_reports_success_for_a_valid_molecule(self):
        ok, flags = utils.sanitize_mol(Chem.MolFromSmiles("CCO"))
        assert ok
        assert flags == Chem.SanitizeFlags.SANITIZE_NONE

    def test_reports_failure_for_a_bad_valence(self):
        mol = Chem.MolFromSmiles("CCO")
        mol.GetAtomWithIdx(0).SetNumExplicitHs(10)
        ok, flags = utils.sanitize_mol(mol)
        assert not ok
        assert flags != Chem.SanitizeFlags.SANITIZE_NONE


class TestSimilarity:
    def test_tanimoto_matrix_has_the_expected_shape(self, fps):
        sims = utils.get_batch_tani_sims(fps[:2], fps)
        assert sims.shape == (2, 4)

    def test_a_molecule_is_identical_to_itself(self, fps):
        assert utils.get_batch_tani_sims([fps[0]], [fps[0]])[0][0] == pytest.approx(1.0)

    def test_tanimoto_stays_in_range(self, fps):
        sims = utils.get_batch_tani_sims(fps, fps)
        assert sims.min() >= 0.0 and sims.max() <= 1.0

    def test_tanimoto_is_symmetric(self, fps):
        sims = utils.get_batch_tani_sims(fps, fps)
        assert np.allclose(sims, sims.T)

    def test_similar_molecules_score_above_dissimilar_ones(self, fps):
        ethanol, propane, benzene = fps[0], fps[1], fps[2]
        sims = utils.get_batch_tani_sims([ethanol], [propane, benzene])[0]
        assert sims[0] > sims[1]

    def test_tversky_weights_the_query_by_default(self, fps):
        """
        The default weighting favours the query, so a substructure of a larger
        molecule scores higher than the reverse comparison.
        """
        small, large = fps[0], fps[3]
        forward = utils.get_batch_tversky_sims([small], [large])[0][0]
        backward = utils.get_batch_tversky_sims([large], [small])[0][0]
        assert forward > backward

    def test_tversky_matrix_shape_and_range(self, fps):
        sims = utils.get_batch_tversky_sims(fps[:2], fps)
        assert sims.shape == (2, 4)
        assert sims.min() >= 0.0
        assert sims.max() == pytest.approx(1.0, abs=1e-9) or sims.max() < 1.0

    def test_pairwise_matches_the_batch_result(self, fps):
        pairwise = utils.get_tani_sim_fp(fps[0], fps[1])
        batch = utils.get_batch_tani_sims([fps[0]], [fps[1]])[0][0]
        assert pairwise == pytest.approx(batch)


class TestColours:
    def test_full_alpha_returns_the_foreground(self):
        assert utils.make_rgb_transparent((1.0, 0.0, 0.0), (1.0, 1.0, 1.0), 1.0) == (
            1.0,
            0.0,
            0.0,
        )

    def test_zero_alpha_returns_the_background(self):
        assert utils.make_rgb_transparent((1.0, 0.0, 0.0), (1.0, 1.0, 1.0), 0.0) == (
            1.0,
            1.0,
            1.0,
        )


class TestRendering:
    def test_svg_is_returned_as_a_data_uri(self):
        uri = utils.get_svg_mol(Chem.MolFromSmiles("CCO"))
        assert uri.startswith("data:image/svg+xml;base64,")

    def test_atom_indices_can_be_shown(self):
        assert utils.get_svg_mol(Chem.MolFromSmiles("CCO"), show_idx=True)

    def test_smiles_input_is_accepted(self):
        assert utils.get_svg_mol("CCO").startswith("data:image/svg+xml;base64,")

    def test_building_blocks_can_be_highlighted(self):
        uri = utils.get_svg_mol_with_bbs("CC(=O)NCc1ccccc1", ["CC(=O)O", "NCc1ccccc1"])
        assert uri.startswith("data:image/svg+xml;base64,")

    def test_highlighting_tolerates_a_non_matching_block(self):
        """A block that is not part of the product must not break rendering."""
        assert utils.get_svg_mol_with_bbs("CC(=O)NCc1ccccc1", ["c1ccncc1"])


class TestReadCxsmiles:
    def test_header_is_skipped_by_default(self, tmp_path):
        path = tmp_path / "mols.smi"
        path.write_text("smiles\nCCO\nCCC\n")
        assert utils.read_cxsmiles_file(str(path)) == ["CCO", "CCC"]

    def test_header_can_be_disabled(self, tmp_path):
        path = tmp_path / "mols.smi"
        path.write_text("CCO\nCCC\n")
        assert utils.read_cxsmiles_file(str(path), header=False) == ["CCO", "CCC"]

    def test_only_the_first_field_is_taken(self, tmp_path):
        path = tmp_path / "mols.smi"
        path.write_text("CCO ethanol extra\n")
        assert utils.read_cxsmiles_file(str(path), header=False) == ["CCO"]

    def test_missing_file_raises(self):
        with pytest.raises(FileNotFoundError):
            utils.read_cxsmiles_file("/nonexistent/path.smi")
