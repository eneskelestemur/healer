"""
Smoke tests for healer.web.interface — building block discovery and path resolution.
"""

from pathlib import Path

import pytest

from healer.domain.bb_repository import resolve_bb_path as repo_resolve_bb_path
from healer.web.interface import BB_BASE_PATH, discover_building_blocks, resolve_bb_path

# ---------------------------------------------------------------------------
# discover_building_blocks
# ---------------------------------------------------------------------------


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


def test_discover_values_are_absolute_paths():
    """All 'value' fields must be absolute paths that exist on disk."""
    result = discover_building_blocks()
    for entry in result:
        path = Path(entry["value"])
        assert path.is_absolute(), f"Non-absolute path in discover results: {path}"
        assert path.exists(), f"Path does not exist on disk: {path}"


def test_discover_test_entry_has_expected_label():
    """The test BB entry should carry the 'Test Set (100 BBs)' label."""
    result = discover_building_blocks()
    test_entries = [e for e in result if e.get("key") == "test"]
    if test_entries:
        assert test_entries[0]["label"] == "Test Set (100 BBs)"


# ---------------------------------------------------------------------------
# resolve_bb_path (interface layer)
# ---------------------------------------------------------------------------


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


# ---------------------------------------------------------------------------
# resolve_bb_path (domain / bb_repository layer)
# ---------------------------------------------------------------------------


def test_repo_resolve_named_key_test():
    """'test' named key resolves via bb_repository layer to an existing SDF."""
    path = repo_resolve_bb_path("test")
    assert Path(path).exists()
    assert path.endswith(".sdf")


def test_repo_resolve_nonexistent_raises():
    with pytest.raises(FileNotFoundError):
        repo_resolve_bb_path("/does/not/exist.sdf")
