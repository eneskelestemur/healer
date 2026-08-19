"""Tests for building block path resolution and repository loading."""

import pickle
from pathlib import Path

import pytest

from healer.domain.bb_repository import (
    BBRepository,
    _build_bb_paths,
    clear_repository_cache,
    get_repository,
    resolve_bb_path,
)


class TestNamedPaths:
    def test_the_expected_sources_are_defined(self):
        paths = _build_bb_paths()
        for key in ("US_stock", "EU_stock", "Global_stock", "test"):
            assert key in paths

    def test_stock_sources_are_globs(self):
        paths = _build_bb_paths()
        for key in ("US_stock", "EU_stock", "Global_stock"):
            assert "*" in paths[key]

    def test_the_test_source_is_a_concrete_file(self):
        assert "*" not in _build_bb_paths()["test"]


class TestResolution:
    def test_a_named_source_resolves_to_an_existing_file(self):
        resolved = Path(resolve_bb_path("test"))
        assert resolved.is_absolute() and resolved.exists()

    def test_an_absolute_path_is_returned_as_given(self, test_bb_path):
        assert resolve_bb_path(test_bb_path) == test_bb_path

    def test_a_missing_file_raises(self):
        with pytest.raises(FileNotFoundError):
            resolve_bb_path("/nonexistent/blocks.sdf")

    def test_an_unknown_name_raises(self):
        with pytest.raises(FileNotFoundError):
            resolve_bb_path("Mars_stock")

    def test_a_glob_picks_the_newest_match(self, tmp_path):
        older, newer = tmp_path / "a_processed.sdf", tmp_path / "b_processed.sdf"
        older.write_text("")
        newer.write_text("")
        import os

        os.utime(older, (1, 1))
        assert resolve_bb_path(str(tmp_path / "*_processed.sdf")) == str(newer)

    def test_a_glob_matching_nothing_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="No file matches"):
            resolve_bb_path(str(tmp_path / "*_processed.sdf"))

    def test_a_relative_path_resolves_against_a_base(self, tmp_path):
        target = tmp_path / "blocks.sdf"
        target.write_text("")
        assert resolve_bb_path("blocks.sdf", base_dir=tmp_path) == str(target)

    def test_a_relative_path_without_a_base_raises(self):
        with pytest.raises(FileNotFoundError):
            resolve_bb_path("blocks.sdf")


class TestLoading:
    def test_blocks_and_the_reaction_index_are_populated(self, test_bb_repository):
        assert test_bb_repository.is_loaded
        assert len(test_bb_repository) > 0
        assert test_bb_repository._reaction_bb_indices

    def test_counts_agree(self, test_bb_repository):
        assert test_bb_repository.loaded_count == len(test_bb_repository)
        assert test_bb_repository.total_count >= test_bb_repository.loaded_count

    def test_every_block_carries_a_fingerprint(self, test_bb_repository):
        assert all(bb.fingerprint is not None for bb in test_bb_repository)

    def test_loading_twice_is_a_no_op(self, test_bb_path):
        repo = BBRepository(source_path=test_bb_path)
        repo.load(show_progress=False)
        first = repo.get_all_bbs()
        repo.load(show_progress=False)
        assert repo.get_all_bbs() is first

    def test_using_a_repository_before_loading_raises(self, test_bb_path):
        repo = BBRepository(source_path=test_bb_path)
        with pytest.raises(RuntimeError, match="not loaded"):
            repo.get_all_bbs()
        with pytest.raises(RuntimeError, match="not loaded"):
            list(repo)


class TestReactionFiltering:
    def test_all_reactions_select_at_least_one_block(
        self, test_bb_repository, all_reactions
    ):
        assert test_bb_repository.get_bbs_for_reactions(all_reactions)

    def test_an_unused_reaction_selects_nothing(
        self, test_bb_repository, all_reactions
    ):
        class Unused:
            name = "definitely-not-a-real-reaction"

        assert test_bb_repository.get_bbs_for_reactions([Unused()]) == []

    def test_filtering_is_a_union_over_reactions(
        self, test_bb_repository, all_reactions
    ):
        pair = all_reactions[:2]
        union = {
            bb.get_smiles() for bb in test_bb_repository.get_bbs_for_reactions(pair)
        }
        singles = set()
        for rxn in pair:
            singles |= {
                bb.get_smiles()
                for bb in test_bb_repository.get_bbs_for_reactions([rxn])
            }
        assert union == singles

    def test_the_iterator_matches_the_list(self, test_bb_repository, all_reactions):
        subset = all_reactions[:3]
        assert [
            bb.get_smiles() for bb in test_bb_repository.iter_bbs_for_reactions(subset)
        ] == [
            bb.get_smiles() for bb in test_bb_repository.get_bbs_for_reactions(subset)
        ]

    def test_blocks_are_shared_not_copied(self, test_bb_repository, all_reactions):
        """Pools reference the same objects, which is what keeps memory flat."""
        first = test_bb_repository.get_bbs_for_reactions(all_reactions)
        second = test_bb_repository.get_bbs_for_reactions(all_reactions)
        assert first[0] is second[0]


class TestAccess:
    def test_indexing_returns_the_same_object(self, test_bb_repository):
        assert (
            test_bb_repository.get_bb_by_index(0) is test_bb_repository.get_all_bbs()[0]
        )

    def test_membership_is_supported(self, test_bb_repository):
        assert test_bb_repository.get_bb_by_index(0) in test_bb_repository


class TestCaching:
    def test_the_same_source_returns_the_same_repository(self, test_bb_path):
        assert get_repository(test_bb_path) is get_repository(test_bb_path)

    def test_a_name_and_its_path_share_one_repository(self, test_bb_path):
        assert get_repository("test") is get_repository(test_bb_path)

    def test_clearing_the_cache_creates_a_new_repository(self, test_bb_path):
        first = get_repository(test_bb_path)
        clear_repository_cache()
        assert get_repository(test_bb_path) is not first


class TestPickling:
    def test_a_repository_survives_a_round_trip(self, test_bb_path):
        """Repositories travel to worker processes during parallel synthesis."""
        repo = BBRepository(source_path=test_bb_path)
        repo.load(show_progress=False)

        restored = pickle.loads(pickle.dumps(repo))
        assert len(restored) == len(repo)
        assert restored.is_loaded
        assert restored._supplier is not None
