"""HEALER: Hit Expansion to Advanced Leads Using Enumerated Reactions"""

from importlib.metadata import PackageNotFoundError, version as _pkg_version

# Must be applied before any BuildingBlock is passed to an RDKit function.
import healer.utils.rdkit_monkey_patch  # noqa: F401

from healer.application.healer import MoleculeHEALER, SiteHEALER, FragmentHEALER
from healer.domain.bb_repository import get_repository, clear_repository_cache, _build_bb_paths

try:
    __version__ = _pkg_version("mol-healer")
except PackageNotFoundError:
    __version__ = "0.0.0.dev0"

__all__ = [
    "MoleculeHEALER",
    "SiteHEALER",
    "FragmentHEALER",
    "get_repository",
    "clear_repository_cache",
    "_build_bb_paths",
    "__version__",
]
