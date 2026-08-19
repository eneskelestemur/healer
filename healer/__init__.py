"""HEALER: Hit Expansion to Advanced Leads Using Enumerated Reactions"""

import logging
from importlib.metadata import PackageNotFoundError, version as _pkg_version

# Must be applied before any BuildingBlock is passed to an RDKit function.
import healer.utils.rdkit_monkey_patch  # noqa: F401

from healer.logging_config import configure_logging, _default_level

# Stay silent until the application attaches a handler.
logging.getLogger("healer").addHandler(logging.NullHandler())
if _default_level() is not None:
    configure_logging(_default_level())

from healer.application.healer import MoleculeHEALER, SiteHEALER, FragmentHEALER
from healer.application.optimizers import (
    BaseOptimizer,
    BaseStagewiseOptimizer,
    BaseSequenceOptimizer,
    BeamSearchOptimizer,
    GeneticAlgorithmOptimizer,
    BayesianSequenceOptimizer,
    OptimizerError,
)
from healer.domain.bb_repository import get_repository, clear_repository_cache, _build_bb_paths
from healer.utils.progress import progress_enabled

try:
    __version__ = _pkg_version("mol-healer")
except PackageNotFoundError:
    __version__ = "0.0.0.dev0"

__all__ = [
    "MoleculeHEALER",
    "SiteHEALER",
    "FragmentHEALER",
    "BaseOptimizer",
    "BaseStagewiseOptimizer",
    "BaseSequenceOptimizer",
    "BeamSearchOptimizer",
    "GeneticAlgorithmOptimizer",
    "BayesianSequenceOptimizer",
    "OptimizerError",
    "configure_logging",
    "progress_enabled",
    "get_repository",
    "clear_repository_cache",
    "_build_bb_paths",
    "__version__",
]
