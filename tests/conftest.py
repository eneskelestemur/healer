"""
Shared pytest fixtures for HEALER tests.
"""

import pytest

# Must be imported before any RDKit chemistry runs so BuildingBlock objects
# are transparently unwrapped by RunReactants, MolToSmiles, etc.
import healer.utils.rdkit_monkey_patch  # noqa: F401
from healer.domain.bb_repository import BBRepository, get_repository, resolve_bb_path

# Penicillin G — has amide bonds that the tree builder reliably splits at depth 1.
# Used across multiple test modules.
PENICILLIN_SMILES = "CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O"

# Aspirin — simpler; useful for routes/API tests that need a quick valid molecule.
ASPIRIN_SMILES = "CC(=O)Oc1ccccc1C(=O)O"


@pytest.fixture(scope="session")
def test_bb_path() -> str:
    """Return the absolute path to the bundled 100-BB test SDF."""
    return resolve_bb_path("test")


@pytest.fixture(scope="session")
def test_bb_repository(test_bb_path: str) -> BBRepository:
    """
    A session-scoped BBRepository loaded from the bundled test SDF.
    Loaded once and reused across all tests — mimics the module-level cache
    that avoids reloading large BB libraries in production.
    """
    repo = get_repository(test_bb_path)
    if not repo.is_loaded:
        repo.load(show_progress=False)
    return repo
