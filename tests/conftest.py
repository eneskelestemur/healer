"""Shared fixtures for the HEALER test suite."""

from pathlib import Path

import pytest
from rdkit import Chem

# Must be imported before any RDKit chemistry runs so BuildingBlock objects
# are transparently unwrapped by RunReactants, MolToSmiles, etc.
import healer.utils.rdkit_monkey_patch  # noqa: F401
import healer.utils.utils as utils
from healer.domain.bb_repository import BBRepository, get_repository, resolve_bb_path
from healer.domain.reaction_template import ReactionTemplate21

HEALER_PKG = Path(__file__).parent.parent / "healer"
REACTIONS_PATH = HEALER_PKG / "data" / "reactions" / "reactions.json"

# Penicillin G — amide bonds that the tree builder reliably splits at depth 1.
PENICILLIN_SMILES = "CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O"

# Aspirin — simpler; useful where any valid molecule will do.
ASPIRIN_SMILES = "CC(=O)Oc1ccccc1C(=O)O"

# A two-fragment query for FragmentHEALER.
FRAGMENT_SMILES = "c1ccccc1N.CC(=O)O"


@pytest.fixture(scope="session")
def test_bb_path() -> str:
    """Absolute path to the bundled 100-BB test SDF."""
    return resolve_bb_path("test")


@pytest.fixture(scope="session")
def test_bb_repository(test_bb_path: str) -> BBRepository:
    """
    Session-scoped BBRepository loaded from the bundled test SDF, mirroring the
    module-level cache that avoids reloading large libraries in production.
    """
    repo = get_repository(test_bb_path)
    if not repo.is_loaded:
        repo.load(show_progress=False)
    return repo


@pytest.fixture(scope="session")
def all_reactions() -> list[ReactionTemplate21]:
    """Every valid reaction template shipped with the package."""
    return [
        r for r in utils.load_reactions_from_json(str(REACTIONS_PATH)) if r.is_valid()
    ]


@pytest.fixture(scope="session")
def amide_reaction(all_reactions: list[ReactionTemplate21]) -> ReactionTemplate21:
    """A representative amide coupling template."""
    for rxn in all_reactions:
        if "amide coupling" in rxn.tags:
            return rxn
    pytest.skip("no amide coupling template available")


@pytest.fixture
def penicillin() -> Chem.Mol:
    return Chem.MolFromSmiles(PENICILLIN_SMILES)


@pytest.fixture
def aspirin() -> Chem.Mol:
    return Chem.MolFromSmiles(ASPIRIN_SMILES)


@pytest.fixture
def bb_pool(test_bb_repository: BBRepository):
    """A handful of real building blocks from the bundled set."""
    return test_bb_repository.get_all_bbs()[:10]
