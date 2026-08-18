"""
Smoke tests for optimizer classes.
"""
import pytest
from rdkit import Chem
from rdkit.Chem import QED

from healer.application.healer import MoleculeHEALER, FragmentHEALER
from healer.application.optimizers import (
    BaseOptimizer,
    BeamSearchOptimizer,
    GeneticAlgorithmOptimizer,
    BayesianSequenceOptimizer,
    qed_batch_scorer,
    stoplight_batch_scorer,
)
from healer.domain.bb_repository import BBRepository
from healer.domain.building_block import BuildingBlock
from tests.conftest import PENICILLIN_SMILES


def qed_fn(mol: Chem.Mol) -> float:
    return QED.qed(mol)


def failing_fn(mol: Chem.Mol) -> float:
    raise RuntimeError("always fails")


def batch_qed(mols):
    return [QED.qed(m) for m in mols]


# ---------------------------------------------------------------------------
# BaseOptimizer
# ---------------------------------------------------------------------------

def test_base_optimizer_requires_fn():
    with pytest.raises(ValueError):
        class _Dummy(BaseOptimizer):
            pass
        _Dummy()


def test_evaluate_batch_with_target_fn():
    opt = BeamSearchOptimizer(beam_width=5, target_fn=qed_fn)
    mols = [Chem.MolFromSmiles("C"), Chem.MolFromSmiles("CC"), Chem.MolFromSmiles("CCC")]
    scores = opt._evaluate_batch(mols)
    assert len(scores) == 3
    assert all(isinstance(s, float) for s in scores)


def test_evaluate_batch_with_batch_fn():
    opt = BeamSearchOptimizer(beam_width=5, batch_target_fn=batch_qed)
    mols = [Chem.MolFromSmiles("C"), Chem.MolFromSmiles("CC")]
    scores = opt._evaluate_batch(mols)
    assert len(scores) == 2
    assert all(isinstance(s, float) for s in scores)


def test_evaluate_batch_handles_failure():
    opt = BeamSearchOptimizer(beam_width=5, target_fn=failing_fn)
    mols = [Chem.MolFromSmiles("C")]
    scores = opt._evaluate_batch(mols)
    assert scores == [None]


# ---------------------------------------------------------------------------
# qed_batch_scorer / stoplight_batch_scorer
# ---------------------------------------------------------------------------

def test_qed_batch_scorer():
    mols = [Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O"), Chem.MolFromSmiles("c1ccccc1")]
    scores = qed_batch_scorer(mols)
    assert len(scores) == 2
    assert all(s is None or isinstance(s, float) for s in scores)


def test_stoplight_batch_scorer():
    mols = [Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")]
    scores = stoplight_batch_scorer(mols)
    assert len(scores) == 1
    assert scores[0] is None or isinstance(scores[0], float)


# ---------------------------------------------------------------------------
# BeamSearchOptimizer — full pipeline
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def beam_healer(test_bb_repository: BBRepository) -> MoleculeHEALER:
    return MoleculeHEALER(
        bb_source="test",
        reaction_tags="all",
        bb_repository=test_bb_repository,
        sim_threshold=0.0,
        max_bbs_per_frag=10,
        verbose=0,
    )


def test_beam_search_optimizer_pipeline(beam_healer: MoleculeHEALER):
    opt = BeamSearchOptimizer(beam_width=5, target_fn=qed_fn)
    beam_healer.set_query_mol(PENICILLIN_SMILES, n_compositions=3, retro_tree_depth=1, min_frag_size=3)
    beam_healer.enumerate(optimizer=opt, max_total_products=10)
    results = beam_healer.get_results(as_dict=True, calc_similarity=False, calc_properties=False)
    assert len(results) >= 1
    assert all(Chem.MolFromSmiles(r["Product"]) is not None for r in results)
    # At least some products should have an optimization_score
    scored = [r for r in results if "optimization_score" in r and r["optimization_score"] is not None]
    assert len(scored) > 0


def test_beam_search_with_batch_fn(beam_healer: MoleculeHEALER):
    opt = BeamSearchOptimizer(beam_width=5, batch_target_fn=qed_batch_scorer)
    beam_healer.set_query_mol(PENICILLIN_SMILES, n_compositions=2, retro_tree_depth=1, min_frag_size=3)
    beam_healer.enumerate(optimizer=opt, max_total_products=5)
    assert len(beam_healer.enumerated_molecules) >= 1


# ---------------------------------------------------------------------------
# GeneticAlgorithmOptimizer — full pipeline
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def ga_healer(test_bb_repository: BBRepository) -> MoleculeHEALER:
    return MoleculeHEALER(
        bb_source="test",
        reaction_tags="all",
        bb_repository=test_bb_repository,
        sim_threshold=0.0,
        max_bbs_per_frag=10,
        verbose=0,
    )


def test_ga_optimizer_pipeline(ga_healer: MoleculeHEALER):
    opt = GeneticAlgorithmOptimizer(
        population_size=10,
        mutation_percent_genes=20,
        random_seed=42,
        target_fn=qed_fn,
    )
    ga_healer.set_query_mol(PENICILLIN_SMILES, n_compositions=2, retro_tree_depth=1, min_frag_size=3)
    ga_healer.enumerate(optimizer=opt, max_evals_per_comp=30, max_total_products=20)
    results = ga_healer.get_results(as_dict=True, calc_similarity=False, calc_properties=False)
    assert len(results) >= 1
    assert all(Chem.MolFromSmiles(r["Product"]) is not None for r in results)


def test_ga_optimizer_evolves(ga_healer: MoleculeHEALER):
    """GA should run multiple generations when budget allows."""
    opt = GeneticAlgorithmOptimizer(population_size=8, random_seed=0, target_fn=qed_fn)
    ga_healer.set_query_mol(PENICILLIN_SMILES, n_compositions=1, retro_tree_depth=1, min_frag_size=3)
    ga_healer.enumerate(optimizer=opt, max_evals_per_comp=50, max_total_products=40)
    assert opt._ga.generations_completed >= 1


# ---------------------------------------------------------------------------
# BayesianSequenceOptimizer — full pipeline
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def bayes_healer(test_bb_repository: BBRepository) -> MoleculeHEALER:
    return MoleculeHEALER(
        bb_source="test",
        reaction_tags="all",
        bb_repository=test_bb_repository,
        sim_threshold=0.0,
        max_bbs_per_frag=10,
        verbose=0,
    )


def test_bayesian_optimizer_pipeline(bayes_healer: MoleculeHEALER):
    opt = BayesianSequenceOptimizer(batch_size=3, target_fn=qed_fn)
    bayes_healer.set_query_mol(PENICILLIN_SMILES, n_compositions=1, retro_tree_depth=1, min_frag_size=3)
    bayes_healer.enumerate(optimizer=opt, max_evals_per_comp=15, max_total_products=10)
    results = bayes_healer.get_results(as_dict=True, calc_similarity=False, calc_properties=False)
    assert len(results) >= 1
    assert all(Chem.MolFromSmiles(r["Product"]) is not None for r in results)


def test_bayesian_optimizer_with_batch_fn(bayes_healer: MoleculeHEALER):
    opt = BayesianSequenceOptimizer(batch_size=3, batch_target_fn=batch_qed)
    bayes_healer.set_query_mol(PENICILLIN_SMILES, n_compositions=1, retro_tree_depth=1, min_frag_size=3)
    bayes_healer.enumerate(optimizer=opt, max_evals_per_comp=12, max_total_products=8)
    assert len(bayes_healer.enumerated_molecules) >= 1


# ---------------------------------------------------------------------------
# Failure handling: scoring failures excluded from tell() feedback
# ---------------------------------------------------------------------------

def test_failing_scorer_does_not_crash_ga(ga_healer: MoleculeHEALER):
    """A target_fn that always fails should not crash GA enumeration."""
    opt = GeneticAlgorithmOptimizer(population_size=6, random_seed=1, target_fn=failing_fn)
    ga_healer.set_query_mol(PENICILLIN_SMILES, n_compositions=1, retro_tree_depth=1, min_frag_size=3)
    ga_healer.enumerate(optimizer=opt, max_total_products=5)
    assert len(ga_healer.enumerated_molecules) >= 1


# ---------------------------------------------------------------------------
# Domain capping
# ---------------------------------------------------------------------------

def _dummy_pool(n: int):
    return [BuildingBlock(Chem.MolFromSmiles("C" * (i + 1))) for i in range(n)]


def test_prepare_domain_without_cap_returns_pools_unchanged():
    pool = _dummy_pool(10)
    opt = GeneticAlgorithmOptimizer(target_fn=qed_fn, max_domain_per_frag=None)
    assert opt.prepare_domain([pool]) == [pool]


def test_prepare_domain_keeps_most_similar_bbs():
    pool = _dummy_pool(6)
    sims = ([0.1, 0.9, 0.2, 0.8, 0.3, 0.95],)
    opt = BayesianSequenceOptimizer(target_fn=qed_fn, max_domain_per_frag=3)
    kept = opt.prepare_domain([pool], sims)[0]
    assert [bb.get_smiles() for bb in kept] == ["CCCCCC", "CC", "CCCC"]


def test_prepare_domain_falls_back_to_order_without_sims():
    pool = _dummy_pool(6)
    opt = BayesianSequenceOptimizer(target_fn=qed_fn, max_domain_per_frag=3)
    kept = opt.prepare_domain([pool], None)[0]
    assert kept == pool[:3]


def test_capped_healer_retains_sims_uncapped_does_not(test_bb_repository: BBRepository):
    for max_bbs, expect_sims in ((8, True), (-1, False)):
        healer = MoleculeHEALER(
            bb_source="test",
            reaction_tags="all",
            bb_repository=test_bb_repository,
            sim_threshold=0.0,
            max_bbs_per_frag=max_bbs,
            verbose=0,
        )
        healer.set_query_mol(PENICILLIN_SMILES, n_compositions=2)
        healer._process_query_mol()
        healer._process_building_blocks()
        comp = healer._compositions[0]
        if expect_sims:
            assert comp.fragment_sims is not None
            assert [len(s) for s in comp.fragment_sims] == [len(p) for p in comp.fragment_bbs]
        else:
            assert comp.fragment_sims is None


# ---------------------------------------------------------------------------
# Sequence optimizers must assemble exactly what they proposed
# ---------------------------------------------------------------------------

def test_sequence_optimizer_assembles_proposed_tuples(test_bb_repository: BBRepository):
    """
    With three or more fragments, records must be paired with building blocks
    through their origin tag. Pairing by list position instead lets products
    drift onto combinations the optimizer never asked for, because applying a
    reaction can yield several products per candidate or none at all.
    """
    class SpyGA(GeneticAlgorithmOptimizer):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.proposed = set()

        def ask(self):
            bb_tuples = super().ask()
            self.proposed.update(
                tuple(bb.get_smiles() for bb in bb_tuple) for bb_tuple in bb_tuples
            )
            return bb_tuples

    healer = FragmentHEALER(
        bb_source="test",
        reaction_tags="all",
        bb_repository=test_bb_repository,
        sim_threshold=0.0,
        max_bbs_per_frag=6,
        verbose=0,
    )
    healer.set_query_mol("c1ccccc1N.CC(=O)O.CCN")
    opt = SpyGA(population_size=6, random_seed=0, target_fn=qed_fn)
    healer.enumerate(optimizer=opt, max_evals_per_comp=18, max_total_products=25)

    assembled = [r for r in healer.enumerated_molecules if r.origin is not None]
    assert assembled, "no products were assembled"
    assert all(len(r.bbs) == 3 for r in assembled)
    for rec in assembled:
        assert tuple(bb.get_smiles() for bb in rec.bbs) in opt.proposed
