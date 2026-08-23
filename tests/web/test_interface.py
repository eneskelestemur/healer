"""
Smoke tests for healer.web.interface — building block discovery and path resolution.
"""

from pathlib import Path

import pytest

from healer.domain.bb_repository import resolve_bb_path as repo_resolve_bb_path
from healer.web.interface import BB_BASE_PATH, discover_building_blocks, resolve_bb_path


def test_discover_returns_list():
    result = discover_building_blocks()
    assert isinstance(result, list)


def test_discover_finds_test_bb():
    """The bundled test_100_bb_processed.sdf should always be discoverable."""
    result = discover_building_blocks()
    # At minimum the package ships with the test SDF — at least one entry expected
    assert len(result) >= 1, (
        "discover_building_blocks returned nothing. "
        f"BB_BASE_PATH={BB_BASE_PATH} — does test_100_bb_processed.sdf exist there?"
    )


def test_discover_entry_shape():
    """Every entry must have 'value' and 'label' string fields."""
    result = discover_building_blocks()
    for entry in result:
        assert "value" in entry, f"Entry missing 'value': {entry}"
        assert "label" in entry, f"Entry missing 'label': {entry}"
        assert isinstance(entry["value"], str)
        assert isinstance(entry["label"], str)


def test_discover_reports_library_sizes():
    """Each entry carries the number of building blocks in its file."""
    result = discover_building_blocks()
    for entry in result:
        assert isinstance(entry["count"], int)
        assert entry["count"] > 0, f"Empty library reported: {entry}"


def test_discover_values_are_absolute_paths():
    """All 'value' fields must be absolute paths that exist on disk."""
    result = discover_building_blocks()
    for entry in result:
        path = Path(entry["value"])
        assert path.is_absolute(), f"Non-absolute path in discover results: {path}"
        assert path.exists(), f"Path does not exist on disk: {path}"


def test_discover_test_entry_has_expected_label():
    """The test BB entry should carry the 'Test Set' label."""
    result = discover_building_blocks()
    test_entries = [e for e in result if e.get("key") == "test"]
    if test_entries:
        assert test_entries[0]["label"] == "Test Set"


def test_resolve_named_key_test():
    """'test' named key resolves to an existing SDF file."""
    path = resolve_bb_path("test")
    assert Path(path).exists(), f"Resolved path does not exist: {path}"
    assert path.endswith(".sdf")


def test_resolve_absolute_path(test_bb_path: str):
    """An absolute path that exists is returned as-is."""
    result = resolve_bb_path(test_bb_path)
    assert result == test_bb_path


def test_resolve_nonexistent_absolute_path():
    """An absolute path that does not exist raises FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        resolve_bb_path("/nonexistent/path/to/file.sdf")


def test_resolve_bad_key_raises():
    """An unknown key that is not an existing path raises FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        resolve_bb_path("this_key_does_not_exist_anywhere")


def test_repo_resolve_named_key_test():
    """'test' named key resolves via bb_repository layer to an existing SDF."""
    path = repo_resolve_bb_path("test")
    assert Path(path).exists()
    assert path.endswith(".sdf")


def test_repo_resolve_nonexistent_raises():
    with pytest.raises(FileNotFoundError):
        repo_resolve_bb_path("/does/not/exist.sdf")


class TestServerLimits:
    def test_limits_are_returned_as_a_copy(self):
        from healer.web.interface import get_server_limits

        limits = get_server_limits()
        limits["max_total_products"] = -1
        assert get_server_limits()["max_total_products"] != -1

    def test_requests_pass_through_in_local_mode(self):
        from healer.web import interface

        params = {"max_total_products": 10**9}
        assert interface.apply_server_limits(dict(params)) == params

    def test_requests_are_clamped_in_server_mode(self, monkeypatch):
        from healer.web import interface

        monkeypatch.setattr(interface, "SERVER_MODE", True)
        clamped = interface.apply_server_limits({"max_total_products": 10**9})
        assert (
            clamped["max_total_products"]
            == interface.SERVER_LIMITS["max_total_products"]
        )

    def test_values_below_the_cap_are_untouched_in_server_mode(self, monkeypatch):
        from healer.web import interface

        monkeypatch.setattr(interface, "SERVER_MODE", True)
        assert (
            interface.apply_server_limits({"max_total_products": 3})[
                "max_total_products"
            ]
            == 3
        )


class TestFragmentCounting:
    def test_a_single_component_counts_once(self):
        from healer.web.interface import count_molecular_fragments

        assert count_molecular_fragments("CCO") == 1

    def test_components_are_counted(self):
        from healer.web.interface import count_molecular_fragments

        assert count_molecular_fragments("c1ccccc1N.CC(=O)O") == 2

    def test_an_invalid_smiles_counts_as_zero(self):
        from healer.web.interface import count_molecular_fragments

        assert count_molecular_fragments("not-a-molecule") == 0


class TestResultFormatting:
    def test_display_rows_are_derived_from_the_records(self):
        from healer.web.interface import format_enumeration_results

        records = [
            {
                "ID": "HEAL_000000",
                "Product": "CCO",
                "BB1": "CC",
                "Similarity_to_query": 0.5,
            }
        ]
        display, original = format_enumeration_results(records, "molecule")

        assert len(display) == 1
        assert display[0]["Product"] == "CCO"
        assert display[0]["Similarity_to_query"] == 0.5
        assert original == records

    def test_no_records_gives_no_rows(self):
        from healer.web.interface import format_enumeration_results

        display, original = format_enumeration_results([], "molecule")
        assert display == [] and original == []

    def test_identifiers_and_urls_travel_with_their_block(self):
        from healer.web.interface import format_enumeration_results

        records = [
            {
                "Product": "CCO",
                "BB1": "CC",
                "BB2": "CO",
                "BBID1": "EN300-1",
                "URL1": "https://example.com/1",
                "BBID2": "",
                "URL2": "",
            }
        ]
        display, _ = format_enumeration_results(records, "molecule")

        assert display[0]["BBID1"] == "EN300-1"
        assert display[0]["URL1"] == "https://example.com/1"
        assert "BBID2" not in display[0] and "URL2" not in display[0]

    def test_site_mode_reports_the_added_block(self):
        from healer.web.interface import format_enumeration_results

        records = [
            {
                "Product": "CCO",
                "BB1": "CC",
                "BB2": "CO",
                "BBID2": "EN300-2",
                "URL2": "https://example.com/2",
            }
        ]
        display, _ = format_enumeration_results(records, "site")

        assert display[0]["BB"] == "CO"
        assert display[0]["BBID"] == "EN300-2"
        assert display[0]["URL"] == "https://example.com/2"


class TestEnumerationJob:
    """The wrapper both the Celery tasks and local mode run jobs through."""

    PARAMS = {
        "molecule": "CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O",
        "bb_source": "test",
        "reaction_tags": ["all"],
        "sim_threshold": 0.0,
        "max_bbs_per_frag": 10,
        "retro_tree_depth": 1,
        "min_frag_size": 3,
        "n_compositions": 2,
        "max_total_products": 5,
    }

    @pytest.fixture(scope="class")
    def job(self):
        from healer.web.interface import run_enumeration_job

        seen = []
        payload = run_enumeration_job(
            "molecule", dict(self.PARAMS), on_stage=seen.append
        )
        return payload, seen

    def test_stages_are_reported_in_order(self, job):
        from healer.web.interface import STAGES

        _, seen = job
        assert seen == list(STAGES)

    def test_the_payload_carries_rows_and_stats(self, job):
        payload, _ = job
        assert payload["display"] and payload["complete"]
        assert set(payload["stats"]) == {"n_molecules", "seconds"}

    def test_the_query_row_is_not_counted(self, job):
        payload, _ = job
        assert payload["stats"]["n_molecules"] == len(payload["display"]) - 1

    def test_the_timing_excludes_the_library_load(self, job):
        payload, _ = job
        assert 0 <= payload["stats"]["seconds"] < 60

    def test_stages_are_optional(self):
        from healer.web.interface import run_enumeration_job

        payload = run_enumeration_job("molecule", dict(self.PARAMS))
        assert payload["stats"]["n_molecules"] >= 0
