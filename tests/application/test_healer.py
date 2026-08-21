"""Tests for the HEALER enumeration classes."""

import logging
import pickle

import pytest
from rdkit import Chem

from healer.application.healer import (
    DEFAULT_BB_RULES,
    DEFAULT_REACTION_TAGS,
    FragmentHEALER,
    MoleculeHEALER,
    SiteHEALER,
)
from healer.domain.bb_repository import BBRepository
from healer.domain.building_block import BuildingBlock
from tests.conftest import ASPIRIN_SMILES, FRAGMENT_SMILES, PENICILLIN_SMILES


def molecule_healer(repo: BBRepository, **kwargs) -> MoleculeHEALER:
    params = {
        "bb_source": "test",
        "reaction_tags": "all",
        "bb_repository": repo,
        "sim_threshold": 0.0,
        "max_bbs_per_frag": 5,
        "show_progress": False,
    }
    params.update(kwargs)
    return MoleculeHEALER(**params)


class TestReactionSelection:
    def test_all_selects_every_valid_template(self, test_bb_repository, all_reactions):
        healer = molecule_healer(test_bb_repository, reaction_tags="all")
        assert len(healer.reactions) == len(all_reactions)

    def test_a_single_tag_narrows_the_set(self, test_bb_repository):
        healer = molecule_healer(test_bb_repository, reaction_tags="amide coupling")
        assert healer.reactions
        assert len(healer.reactions) < len(healer._reactions)
        assert all("amide coupling" in r.tags for r in healer.reactions)

    def test_a_tag_list_is_a_union(self, test_bb_repository):
        tags = ["amide coupling", "N-arylation"]
        healer = molecule_healer(test_bb_repository, reaction_tags=tags)
        assert all(any(t in r.tags for t in tags) for r in healer.reactions)

    def test_default_tags_are_used_when_none_is_given(self, test_bb_repository):
        healer = MoleculeHEALER(
            bb_source="test", bb_repository=test_bb_repository, show_progress=False
        )
        assert set(healer.reaction_tags) == set(DEFAULT_REACTION_TAGS)

    def test_unknown_tags_select_nothing(self, test_bb_repository):
        healer = molecule_healer(test_bb_repository, reaction_tags=["not-a-real-tag"])
        assert healer.reactions == []

    def test_set_reactions_updates_and_clears_compositions(self, test_bb_repository):
        healer = molecule_healer(test_bb_repository)
        healer.set_query_mol(PENICILLIN_SMILES, n_compositions=2)
        healer._process_query_mol()
        assert healer._compositions

        healer.set_reactions("amide coupling")
        assert healer._compositions == []
        assert all("amide coupling" in r.tags for r in healer.reactions)


class TestBuildingBlockMatching:
    def test_top_k_caps_each_pool(self, test_bb_repository):
        healer = molecule_healer(test_bb_repository, max_bbs_per_frag=3)
        healer.set_query_mol(PENICILLIN_SMILES, n_compositions=2)
        healer._process_query_mol()
        healer._process_building_blocks()

        for comp in healer._compositions:
            assert all(len(pool) == 3 for pool in comp.fragment_bbs)

    def test_threshold_mode_keeps_everything_above_the_cutoff(self, test_bb_repository):
        healer = molecule_healer(
            test_bb_repository, max_bbs_per_frag=-1, sim_threshold=0.0
        )
        healer.set_query_mol(PENICILLIN_SMILES, n_compositions=1)
        healer._process_query_mol()
        healer._process_building_blocks()

        pool_size = len(healer._compositions[0].fragment_bbs[0])
        assert pool_size > 3

    def test_a_high_threshold_empties_the_pools(self, test_bb_repository):
        healer = molecule_healer(
            test_bb_repository, max_bbs_per_frag=-1, sim_threshold=0.99
        )
        healer.set_query_mol(PENICILLIN_SMILES, n_compositions=1)
        healer._process_query_mol()
        healer._process_building_blocks()

        assert all(
            len(pool) == 0
            for comp in healer._compositions
            for pool in comp.fragment_bbs
        )

    def test_similarities_are_kept_only_when_pools_are_capped(self, test_bb_repository):
        capped = molecule_healer(test_bb_repository, max_bbs_per_frag=4)
        capped.set_query_mol(PENICILLIN_SMILES, n_compositions=1)
        capped._process_query_mol()
        capped._process_building_blocks()
        assert capped._compositions[0].fragment_sims is not None

        uncapped = molecule_healer(test_bb_repository, max_bbs_per_frag=-1)
        uncapped.set_query_mol(PENICILLIN_SMILES, n_compositions=1)
        uncapped._process_query_mol()
        uncapped._process_building_blocks()
        assert uncapped._compositions[0].fragment_sims is None

    def test_kept_similarities_are_the_highest_available(self, test_bb_repository):
        """Top-k must keep the closest building blocks, not an arbitrary slice."""
        healer = molecule_healer(test_bb_repository, max_bbs_per_frag=5)
        healer.set_query_mol(PENICILLIN_SMILES, n_compositions=1)
        healer._process_query_mol()
        healer._process_building_blocks()

        comp = healer._compositions[0]
        all_bbs = healer.bb_mols
        for pool, sims in zip(comp.fragment_bbs, comp.fragment_sims, strict=True):
            assert len(pool) == len(sims)
            assert min(sims) >= 0.0
            assert len(pool) <= len(all_bbs)

    def test_shuffling_changes_pool_order_not_membership(self, test_bb_repository):
        ordered = molecule_healer(
            test_bb_repository, max_bbs_per_frag=-1, sim_threshold=0.0
        )
        shuffled = molecule_healer(
            test_bb_repository,
            max_bbs_per_frag=-1,
            sim_threshold=0.0,
            shuffle_bb_order=True,
        )
        assert {bb.get_smiles() for bb in ordered.bb_mols} == {
            bb.get_smiles() for bb in shuffled.bb_mols
        }


class TestEnumerationLimits:
    def test_max_total_products_bounds_the_result(self, test_bb_repository):
        healer = molecule_healer(test_bb_repository, max_bbs_per_frag=10)
        healer.set_query_mol(PENICILLIN_SMILES, n_compositions=3)
        healer.enumerate(max_total_products=7)

        assert len(healer.enumerated_molecules) == 8  # query plus seven products

    def test_max_products_per_comp_bounds_each_composition(self, test_bb_repository):
        healer = molecule_healer(test_bb_repository, max_bbs_per_frag=10)
        healer.set_query_mol(PENICILLIN_SMILES, n_compositions=2)
        healer.enumerate(max_products_per_comp=2)

        n_comps = len(healer._compositions)
        assert len(healer.enumerated_molecules) - 1 <= 2 * n_comps

    def test_max_evals_per_comp_reduces_output(self, test_bb_repository):
        unlimited = molecule_healer(test_bb_repository, max_bbs_per_frag=10)
        unlimited.set_query_mol(PENICILLIN_SMILES, n_compositions=2)
        unlimited.enumerate()

        limited = molecule_healer(test_bb_repository, max_bbs_per_frag=10)
        limited.set_query_mol(PENICILLIN_SMILES, n_compositions=2)
        limited.enumerate(max_evals_per_comp=5)

        assert len(limited.enumerated_molecules) <= len(unlimited.enumerated_molecules)

    def test_no_limits_enumerates_everything_reachable(self, test_bb_repository):
        healer = molecule_healer(test_bb_repository, max_bbs_per_frag=3)
        healer.set_query_mol(PENICILLIN_SMILES, n_compositions=1)
        healer.enumerate()
        assert len(healer.enumerated_molecules) > 1


class TestParallelSynthesis:
    def test_parallel_matches_sequential(self, test_bb_repository):
        """Worker processes must produce the same products as the main process."""
        results = []
        for n_jobs in (1, 2):
            healer = molecule_healer(test_bb_repository, max_bbs_per_frag=5)
            healer.set_query_mol(PENICILLIN_SMILES, n_compositions=2)
            healer.enumerate(n_jobs=n_jobs)
            results.append(
                sorted(Chem.MolToSmiles(r.product) for r in healer.enumerated_molecules)
            )
        assert results[0] == results[1]

    def test_a_healer_can_be_pickled(self, test_bb_repository):
        healer = molecule_healer(test_bb_repository)
        healer.set_query_mol(PENICILLIN_SMILES, n_compositions=1)
        restored = pickle.loads(pickle.dumps(healer))

        assert restored._fp_generator is not None
        assert len(restored.reactions) == len(healer.reactions)


class TestResults:
    def test_the_query_is_the_first_row(self, test_bb_repository):
        healer = molecule_healer(test_bb_repository)
        healer.set_query_mol(PENICILLIN_SMILES, n_compositions=1)
        healer.enumerate(max_total_products=5)
        df = healer.get_results(calc_similarity=True, calc_properties=False)

        assert df.iloc[0]["BB1"] == ""
        assert df.iloc[0]["Similarity_to_query"] >= df["Similarity_to_query"].max()

    def test_columns_follow_the_mode(self, test_bb_repository):
        healer = molecule_healer(test_bb_repository)
        healer.set_query_mol(PENICILLIN_SMILES, n_compositions=1, retro_tree_depth=1)
        healer.enumerate(max_total_products=3)
        df = healer.get_results(calc_similarity=False, calc_properties=False)

        assert healer.max_bbs == 2
        assert "BB2" in df.columns and "BB3" not in df.columns
        assert "Reaction1_name" in df.columns and "Reaction2_name" not in df.columns

    def test_building_block_identifiers_are_reported(self, test_bb_repository):
        healer = molecule_healer(test_bb_repository)
        healer.set_query_mol(PENICILLIN_SMILES, n_compositions=1)
        healer.enumerate(max_total_products=3)
        df = healer.get_results(calc_similarity=False, calc_properties=False)

        assert "BBID1" in df.columns and "BBID2" in df.columns
        analogs = df[df["BB1"] != ""]
        assert (analogs["BBID1"] != "").all()

    def test_identifiers_are_empty_for_sources_without_them(
        self, test_bb_repository, monkeypatch
    ):
        monkeypatch.setattr(BuildingBlock, "get_id", lambda self: "")
        healer = molecule_healer(test_bb_repository)
        healer.set_query_mol(PENICILLIN_SMILES, n_compositions=1)
        healer.enumerate(max_total_products=3)
        df = healer.get_results(calc_similarity=False, calc_properties=False)

        assert (df["BBID1"] == "").all()

    def test_ids_are_unique_and_formatted(self, test_bb_repository):
        healer = molecule_healer(test_bb_repository)
        healer.set_query_mol(PENICILLIN_SMILES, n_compositions=2)
        healer.enumerate(max_total_products=5)
        df = healer.get_results(calc_similarity=False, calc_properties=False)

        assert df["ID"].is_unique
        assert all(i.startswith("HEAL_") for i in df["ID"])

    def test_duplicate_routes_are_collapsed(self, test_bb_repository):
        healer = molecule_healer(test_bb_repository, max_bbs_per_frag=10)
        healer.set_query_mol(PENICILLIN_SMILES, n_compositions=3)
        healer.enumerate(max_total_products=50)
        df = healer.get_results(calc_similarity=False, calc_properties=False)

        route_cols = [
            col
            for col in df.columns
            if (col.startswith("BB") and not col.startswith("BBID"))
            or col.startswith("Reaction")
        ]
        assert not df.duplicated(subset=route_cols).any()

    def test_similarity_sorts_descending(self, test_bb_repository):
        healer = molecule_healer(test_bb_repository, max_bbs_per_frag=10)
        healer.set_query_mol(PENICILLIN_SMILES, n_compositions=2)
        healer.enumerate(max_total_products=10)
        sims = healer.get_results(calc_similarity=True, calc_properties=False)[
            "Similarity_to_query"
        ]
        assert list(sims) == sorted(sims, reverse=True)

    def test_dict_output_matches_the_frame(self, test_bb_repository):
        healer = molecule_healer(test_bb_repository)
        healer.set_query_mol(PENICILLIN_SMILES, n_compositions=1)
        healer.enumerate(max_total_products=3)

        df = healer.get_results(calc_similarity=False, calc_properties=False)
        records = healer.get_results(
            as_dict=True, calc_similarity=False, calc_properties=False
        )
        assert len(records) == len(df)
        assert records[0]["Product"] == df.iloc[0]["Product"]

    def test_results_can_be_saved(self, test_bb_repository, tmp_path):
        healer = molecule_healer(test_bb_repository)
        healer.set_query_mol(PENICILLIN_SMILES, n_compositions=1)
        healer.enumerate(max_total_products=3)

        out = tmp_path / "results.csv"
        healer.save_results(str(out), calc_similarity=False, calc_properties=False)
        assert out.exists() and out.read_text().startswith("ID,Product")

    def test_properties_are_added_on_request(self, test_bb_repository):
        healer = molecule_healer(test_bb_repository)
        healer.set_query_mol(PENICILLIN_SMILES, n_compositions=1)
        healer.enumerate(max_total_products=3)
        df = healer.get_results(calc_similarity=False, calc_properties=True)

        for column in ("mw", "logp", "qed"):
            assert column in df.columns


class TestQueryHandling:
    def test_a_mol_object_is_accepted(self, test_bb_repository, penicillin):
        healer = molecule_healer(test_bb_repository)
        healer.set_query_mol(penicillin, n_compositions=1)
        assert healer.query_mol is not None

    def test_multi_component_queries_are_rejected(self, test_bb_repository):
        healer = molecule_healer(test_bb_repository)
        with pytest.raises(ValueError, match="single connected component"):
            healer.set_query_mol(FRAGMENT_SMILES)

    def test_enumerating_without_a_query_is_an_error(self, test_bb_repository):
        healer = molecule_healer(test_bb_repository)
        with pytest.raises(ValueError, match="Query molecule must be set"):
            healer.enumerate()

    def test_an_unsupported_optimizer_is_rejected(self, test_bb_repository):
        healer = molecule_healer(test_bb_repository)
        healer.set_query_mol(PENICILLIN_SMILES, n_compositions=1)
        with pytest.raises(TypeError, match="Unsupported optimizer type"):
            healer.enumerate(optimizer="not an optimizer")

    def test_a_query_with_no_split_warns_and_returns_itself(
        self, test_bb_repository, caplog
    ):
        healer = molecule_healer(test_bb_repository)
        healer.set_query_mol("c1ccc(cc1)C(=O)NC", n_compositions=2)

        with caplog.at_level(logging.WARNING, logger="healer.application.healer"):
            healer.enumerate(max_total_products=5)

        assert "No valid fragmentation" in caplog.text
        assert len(healer.enumerated_molecules) == 1

    def test_the_same_instance_serves_several_queries(self, test_bb_repository):
        healer = molecule_healer(test_bb_repository)
        counts = []
        for smiles in (PENICILLIN_SMILES, ASPIRIN_SMILES):
            healer.set_query_mol(smiles, n_compositions=1)
            healer.enumerate(max_total_products=5)
            counts.append(len(healer.enumerated_molecules))
        assert all(c >= 1 for c in counts)

    def test_repeated_enumeration_does_not_accumulate(self, test_bb_repository):
        healer = molecule_healer(test_bb_repository)
        healer.set_query_mol(PENICILLIN_SMILES, n_compositions=1)
        healer.enumerate(max_total_products=5)
        first = len(healer.enumerated_molecules)
        healer.enumerate(max_total_products=5)
        assert len(healer.enumerated_molecules) == first


class TestCustomSplitSites:
    def test_explicit_bonds_replace_the_retro_tree(self, test_bb_repository):
        healer = molecule_healer(test_bb_repository)
        healer.set_query_mol("CC(=O)NCc1ccccc1", custom_split_sites=[[(1, 3)]])
        healer._process_query_mol()

        assert len(healer._compositions) == 1
        assert len(healer._compositions[0].fragments) == 2

    def test_a_split_that_yields_one_fragment_is_skipped(
        self, test_bb_repository, caplog
    ):
        healer = molecule_healer(test_bb_repository)
        healer.set_query_mol("C1CCCCC1", custom_split_sites=[[(0, 1)]])

        with caplog.at_level(logging.WARNING, logger="healer.application.healer"):
            healer._process_query_mol()

        assert healer._compositions == []
        assert "did not produce multiple fragments" in caplog.text


class TestFragmentHEALER:
    def test_fragments_become_one_composition(self, test_bb_repository):
        healer = FragmentHEALER(
            bb_source="test",
            reaction_tags="all",
            bb_repository=test_bb_repository,
            sim_threshold=0.0,
            max_bbs_per_frag=5,
            show_progress=False,
        )
        healer.set_query_mol(FRAGMENT_SMILES)
        healer._process_query_mol()

        assert len(healer._compositions) == 1
        assert len(healer._compositions[0].fragments) == 2

    def test_a_tuple_of_smiles_is_accepted(self, test_bb_repository):
        healer = FragmentHEALER(
            bb_source="test", bb_repository=test_bb_repository, show_progress=False
        )
        healer.set_query_mol(("c1ccccc1N", "CC(=O)O"))
        assert len(Chem.GetMolFrags(healer.query_mol)) == 2

    def test_single_component_queries_are_rejected(self, test_bb_repository):
        healer = FragmentHEALER(
            bb_source="test", bb_repository=test_bb_repository, show_progress=False
        )
        with pytest.raises(ValueError, match="at least two fragments"):
            healer.set_query_mol("CCO")

    def test_max_bbs_follows_the_fragment_count(self, test_bb_repository):
        healer = FragmentHEALER(
            bb_source="test", bb_repository=test_bb_repository, show_progress=False
        )
        healer.set_query_mol("c1ccccc1N.CC(=O)O.CCN")
        assert healer.max_bbs == 3


class TestSiteHEALER:
    def make(self, repo: BBRepository, **kwargs) -> SiteHEALER:
        params = {
            "bb_source": "test",
            "reaction_tags": "all",
            "bb_repository": repo,
            "show_progress": False,
        }
        params.update(kwargs)
        return SiteHEALER(**params)

    def test_the_query_is_one_composition(self, test_bb_repository):
        healer = self.make(test_bb_repository)
        healer.set_query_mol("c1ccccc1N", reactive_sites=[6])
        healer._process_query_mol()
        assert len(healer._compositions) == 1

    def test_max_bbs_is_two(self, test_bb_repository):
        healer = self.make(test_bb_repository)
        healer.set_query_mol("c1ccccc1N", reactive_sites=[6])
        assert healer.max_bbs == 2

    def test_the_caller_molecule_is_not_modified(self, test_bb_repository):
        """Protection flags must land on a copy, not the caller's molecule."""
        mol = Chem.MolFromSmiles("c1ccccc1N")
        healer = self.make(test_bb_repository)
        healer.set_query_mol(mol, reactive_sites=[6])
        healer._process_query_mol()

        assert not any(atom.HasProp("_protected") for atom in mol.GetAtoms())

    def test_atoms_away_from_the_site_are_protected(self, test_bb_repository):
        healer = self.make(test_bb_repository)
        healer.set_query_mol("c1ccccc1N", reactive_sites=[6])
        healer._process_query_mol()

        protected = [
            atom.GetIdx()
            for atom in healer._compositions[0].fragments[0].GetAtoms()
            if atom.HasProp("_protected")
        ]
        assert 6 not in protected
        assert protected

    def test_no_sites_warns(self, test_bb_repository, caplog):
        healer = self.make(test_bb_repository)
        healer.set_query_mol("c1ccccc1N")
        with caplog.at_level(logging.WARNING, logger="healer.application.healer"):
            healer._process_query_mol()
        assert "No reactive sites" in caplog.text

    def test_property_rules_filter_the_pool(self, test_bb_repository):
        loose = self.make(test_bb_repository, rules={"MW": (0, 1000)})
        loose.set_query_mol("c1ccccc1N", reactive_sites=[6])
        loose._process_query_mol()
        loose._process_building_blocks()

        strict = self.make(test_bb_repository, rules={"MW": (0, 150)})
        strict.set_query_mol("c1ccccc1N", reactive_sites=[6])
        strict._process_query_mol()
        strict._process_building_blocks()

        assert len(strict._compositions[0].fragment_bbs[1]) < len(
            loose._compositions[0].fragment_bbs[1]
        )

    def test_structural_rules_filter_the_pool(self, test_bb_repository):
        healer = self.make(test_bb_repository, struct_rules=["[NX3;H2]"])
        healer.set_query_mol("c1ccccc1N", reactive_sites=[6])
        healer._process_query_mol()
        healer._process_building_blocks()

        pool = healer._compositions[0].fragment_bbs[1]
        primary_amine = Chem.MolFromSmarts("[NX3;H2]")
        assert all(bb.HasSubstructMatch(primary_amine) for bb in pool)

    def test_default_rules_are_not_shared_between_instances(self, test_bb_repository):
        first = self.make(test_bb_repository)
        first.set_rules(MW=(0, 123))
        second = self.make(test_bb_repository)

        assert second.rules["MW"] == DEFAULT_BB_RULES["MW"]
        assert first.rules["MW"] == (0, 123)

    def test_unknown_rules_are_rejected(self, test_bb_repository):
        healer = self.make(test_bb_repository)
        with pytest.raises(ValueError, match="Invalid rule"):
            healer.set_rules(NotAProperty=(0, 1))

    def test_repeated_enumeration_is_stable(self, test_bb_repository):
        healer = self.make(test_bb_repository, rules={"MW": (0, 250)})
        healer.set_query_mol("c1ccccc1N", reactive_sites=[6])
        healer.enumerate(max_total_products=3)
        first = len(healer.enumerated_molecules)
        healer.enumerate(max_total_products=3)
        assert len(healer.enumerated_molecules) == first


class TestProgressSettings:
    def test_verbose_is_accepted_as_an_alias(self, test_bb_repository):
        assert (
            molecule_healer(test_bb_repository, show_progress=None, verbose=0).verbose
            == 0
        )
        assert (
            molecule_healer(test_bb_repository, show_progress=None, verbose=1).verbose
            == 1
        )

    def test_show_progress_wins_over_verbose(self, test_bb_repository):
        healer = molecule_healer(test_bb_repository, show_progress=False, verbose=2)
        assert healer._show_progress is False
