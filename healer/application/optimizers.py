'''
    Optimizer interfaces and implementations for guided enumeration.
'''
import logging
from abc import ABC, abstractmethod
from typing import List, Tuple, Callable, Iterable, Optional

from rdkit import Chem

from healer.domain.enumeration_record import EnumerationRecord
from healer.domain.building_block import BuildingBlock
from healer.domain.reaction_template import ReactionTemplate21

logger = logging.getLogger(__name__)

Candidate = Tuple[EnumerationRecord, BuildingBlock, ReactionTemplate21]


class BaseOptimizer(ABC):
    '''Base class holding the scoring function(s) for all optimizers.'''

    def __init__(
        self,
        target_fn: Optional[Callable[[Chem.Mol], float]] = None,
        batch_target_fn: Optional[Callable[[List[Chem.Mol]], List[Optional[float]]]] = None,
    ) -> None:
        '''
            Args:
                target_fn: function that scores a single Mol; called per-molecule when
                    batch_target_fn is not provided.
                batch_target_fn: function that scores a list of Mols in one call.
                    Takes priority over target_fn when both are given. Should return
                    a list of the same length as the input, with None for any molecule
                    that could not be scored.
        '''
        if target_fn is None and batch_target_fn is None:
            raise ValueError("Provide target_fn or batch_target_fn.")
        self.target_fn = target_fn
        self.batch_target_fn = batch_target_fn

    def _evaluate(self, mol: Chem.Mol) -> Optional[float]:
        return self._evaluate_batch([mol])[0]

    def _evaluate_batch(self, mols: List[Chem.Mol]) -> List[Optional[float]]:
        '''Score a batch of molecules; returns None for any that fail.'''
        if self.batch_target_fn is not None:
            try:
                return [float(s) if s is not None else None for s in self.batch_target_fn(mols)]
            except Exception as e:
                logger.warning("batch_target_fn failed: %s", e)
                return [None] * len(mols)
        results = []
        for mol in mols:
            try:
                s = self.target_fn(mol)
                results.append(float(s) if s is not None else None)
            except Exception as e:
                logger.warning("target_fn failed on %s: %s", Chem.MolToSmiles(mol), e)
                results.append(None)
        return results


class BaseStagewiseOptimizer(BaseOptimizer, ABC):
    '''
        Interface for stagewise optimizers that prune or reorder assembled products
        at each fragment-assembly stage.
    '''
    @abstractmethod
    def filter(self, records: List[EnumerationRecord], depth: int) -> List[EnumerationRecord]:
        '''
            Prune or reorder assembled reaction products at this stage depth.

            Called AFTER reactions have been applied; return the subset of records
            to carry forward as seeds for the next reaction stage.
        '''
        ...


class BaseSequenceOptimizer(BaseOptimizer, ABC):
    '''Interface for full-sequence optimizers using ask/tell over BB tuples.'''

    def __init__(
        self,
        target_fn: Optional[Callable[[Chem.Mol], float]] = None,
        batch_target_fn: Optional[Callable[[List[Chem.Mol]], List[Optional[float]]]] = None,
        max_domain_per_frag: Optional[int] = None,
    ) -> None:
        '''
            Args:
                target_fn: function that scores a single Mol.
                batch_target_fn: function that scores a list of Mols in one call.
                max_domain_per_frag: maximum building blocks per fragment to search
                    over. None searches the full pools.
        '''
        super().__init__(target_fn=target_fn, batch_target_fn=batch_target_fn)
        self.max_domain_per_frag = max_domain_per_frag

    def prepare_domain(
        self,
        domain: List[List[BuildingBlock]],
        sims: Optional[Tuple[List[float], ...]] = None,
    ) -> List[List[BuildingBlock]]:
        '''
            Truncate each fragment's building block pool to `max_domain_per_frag`.
            The search space is the product of the pool sizes, so uncapped pools can
            grow past what a surrogate model can enumerate.

            Args:
                domain: building block pool for each fragment.
                sims: fragment-to-BB similarity for each pool, used to keep the
                    closest building blocks. Pools are truncated in their existing
                    order when not available.

            Returns:
                The pools, each no longer than `max_domain_per_frag`.
        '''
        cap = self.max_domain_per_frag
        if cap is None or all(len(pool) <= cap for pool in domain):
            return domain

        prepared: List[List[BuildingBlock]] = []
        for i, pool in enumerate(domain):
            if len(pool) <= cap:
                prepared.append(pool)
                continue
            if sims is not None:
                order = sorted(range(len(pool)), key=lambda j: sims[i][j], reverse=True)
                prepared.append([pool[j] for j in order[:cap]])
            else:
                prepared.append(pool[:cap])

        logger.info(
            "%s: truncated building block pools from %s to %s (max_domain_per_frag=%d)%s.",
            type(self).__name__,
            [len(p) for p in domain], [len(p) for p in prepared], cap,
            "" if sims is not None else "; set max_bbs_per_frag to truncate by similarity",
        )
        return prepared

    @abstractmethod
    def init_search(self, domain: List[List[BuildingBlock]], budget: int) -> None:
        '''Initialize search with the BB domain and evaluation budget.'''
        ...

    @abstractmethod
    def ask(self) -> List[Tuple[BuildingBlock, ...]]:
        '''Propose a list of BB tuples (one BB per fragment) to evaluate next.'''
        ...

    @abstractmethod
    def tell(self, results: List[Tuple[Tuple[BuildingBlock, ...], float]]) -> None:
        '''Receive (bb_tuple, score) pairs for the last asked sequences.'''
        ...



##### Example batch scorers using prop_profiler #####
# TODO: Remove these example batch scorers

def _profile_batch(mols: List[Chem.Mol], column: str, skip_cns_mpo: bool) -> List[Optional[float]]:
    from prop_profiler import profile_molecules
    smiles_in = [Chem.MolToSmiles(m) for m in mols]
    try:
        df = profile_molecules(smiles_in, skip_cns_mpo=skip_cns_mpo, verbose=False)
        score_map = dict(zip(df['smiles'], df[column]))
        return [score_map.get(smi) for smi in smiles_in]
    except Exception as e:
        logger.warning("profile_molecules failed: %s", e)
        return [None] * len(mols)


def qed_batch_scorer(mols: List[Chem.Mol]) -> List[Optional[float]]:
    '''Batch scorer that returns QED scores using prop_profiler.'''
    return _profile_batch(mols, 'qed', skip_cns_mpo=True)


def stoplight_batch_scorer(mols: List[Chem.Mol]) -> List[Optional[float]]:
    '''Batch scorer that returns CNS STOPLIGHT scores using prop_profiler.'''
    return _profile_batch(mols, 'stoplight_score', skip_cns_mpo=True)



class BeamSearchOptimizer(BaseStagewiseOptimizer):
    '''
        Stagewise beam-search optimizer. At each fragment-assembly stage it scores
        the current intermediate molecules and keeps only the top `beam_width`
        candidates, discarding the rest before applying the next reaction.

        Note: scores the assembled intermediate at every stage, so it multiplies
        calls to target_fn. Not suitable for expensive scoring functions — use
        GeneticAlgorithmOptimizer or BayesianSequenceOptimizer instead.
    '''

    def __init__(
        self,
        beam_width: int = 100,
        target_fn: Optional[Callable[[Chem.Mol], float]] = None,
        batch_target_fn: Optional[Callable[[List[Chem.Mol]], List[Optional[float]]]] = None,
    ):
        super().__init__(target_fn=target_fn, batch_target_fn=batch_target_fn)
        self.beam_width = beam_width

    def filter(self, records: List[EnumerationRecord], depth: int) -> List[EnumerationRecord]:
        if not records:
            return records
        mols = [rec.product for rec in records]
        scores = self._evaluate_batch(mols)
        scored = [(rec, s) for rec, s in zip(records, scores) if s is not None]
        if not scored:
            logger.warning("BeamSearchOptimizer: all scores None at depth %d; keeping all %d records", depth, len(records))
            return records[:self.beam_width]
        scored.sort(key=lambda x: x[1], reverse=True)
        return [rec for rec, _ in scored[:self.beam_width]]


class GeneticAlgorithmOptimizer(BaseSequenceOptimizer):
    '''
        Sequence optimizer using a genetic algorithm (PyGAD) over building-block
        index tuples. Each generation corresponds to one ask/tell round.
    '''

    def __init__(
        self,
        population_size: int = 50,
        mutation_percent_genes: float = 10,
        crossover_type: str = 'uniform',
        keep_elitism: int = 2,
        random_seed: Optional[int] = None,
        target_fn: Optional[Callable[[Chem.Mol], float]] = None,
        batch_target_fn: Optional[Callable[[List[Chem.Mol]], List[Optional[float]]]] = None,
        max_domain_per_frag: Optional[int] = None,
    ):
        super().__init__(
            target_fn=target_fn,
            batch_target_fn=batch_target_fn,
            max_domain_per_frag=max_domain_per_frag,
        )
        self.population_size = population_size
        self.mutation_percent_genes = mutation_percent_genes
        self.crossover_type = crossover_type
        self.keep_elitism = keep_elitism
        self.random_seed = random_seed
        self._ga = None
        self._domain: Optional[List[List[BuildingBlock]]] = None
        self._bb_to_idx: Optional[List[dict]] = None
        self._score_cache: dict = {}

    def init_search(self, domain: List[List[BuildingBlock]], budget: int) -> None:
        try:
            import pygad
        except ImportError:
            raise ImportError("Install pygad: pip install 'mol-healer[opt]'")

        self._domain = domain
        self._bb_to_idx = [{bb.get_smiles(): i for i, bb in enumerate(bbs)} for bbs in domain]
        self._score_cache = {}
        gene_space = [list(range(len(bbs))) for bbs in domain]

        def fitness_func(ga, solution, solution_idx):
            return self._score_cache.get(tuple(int(x) for x in solution), -1e9)

        n_genes = len(domain)
        # Compute mutation_num_genes explicitly to avoid PyGAD's rounding-to-zero warning
        # when mutation_percent_genes * n_genes < 1 (e.g. 15% of 2 genes = 0.3 → 0).
        mutation_num_genes = max(1, round(self.mutation_percent_genes / 100 * n_genes))

        self._ga = pygad.GA(
            num_generations=1,
            num_parents_mating=max(2, self.population_size // 2),
            sol_per_pop=self.population_size,
            num_genes=n_genes,
            gene_space=gene_space,
            gene_type=int,
            fitness_func=fitness_func,
            crossover_type=self.crossover_type,
            mutation_num_genes=mutation_num_genes,
            keep_elitism=self.keep_elitism,
            random_seed=self.random_seed,
        )

    def ask(self) -> List[Tuple[BuildingBlock, ...]]:
        return [
            tuple(self._domain[j][int(idx)] for j, idx in enumerate(sol))
            for sol in self._ga.population
        ]

    def tell(self, results: List[Tuple[Tuple[BuildingBlock, ...], float]]) -> None:
        self._score_cache = {}
        for bb_tuple, score in results:
            try:
                key = tuple(self._bb_to_idx[j][bb.get_smiles()] for j, bb in enumerate(bb_tuple))
                self._score_cache[key] = score
            except KeyError:
                logger.warning("BB not found in domain during tell(); skipping.")
        self._ga.run()


class BayesianSequenceOptimizer(BaseSequenceOptimizer):
    '''
        Sequence optimizer using Bayesian optimization (BayBE) with a chemistry-aware
        surrogate model. Uses MORDRED molecular descriptors (via BayBE SubstanceParameter)
        to featurize building blocks; no separate featurization library needed.

        Starts with diversity-based random recommendations (FPS) and switches to
        Bayesian recommendations (BotorchRecommender) after initial measurements.

        The search space is the product of the building block pools, which BayBE
        enumerates up front, so pools are capped at `max_domain_per_frag` by default.
        Setting `max_bbs_per_frag` on the HEALER lets that cap keep the building
        blocks most similar to each fragment instead of an arbitrary subset.
    '''

    def __init__(
        self,
        batch_size: int = 10,
        encoding: str = 'MORDRED',
        target_fn: Optional[Callable[[Chem.Mol], float]] = None,
        batch_target_fn: Optional[Callable[[List[Chem.Mol]], List[Optional[float]]]] = None,
        max_domain_per_frag: Optional[int] = 200,
    ):
        super().__init__(
            target_fn=target_fn,
            batch_target_fn=batch_target_fn,
            max_domain_per_frag=max_domain_per_frag,
        )
        self.batch_size = batch_size
        self.encoding = encoding
        self._campaign = None
        self._domain: Optional[List[List[BuildingBlock]]] = None
        self._bb_to_idx: Optional[List[dict]] = None

    def init_search(self, domain: List[List[BuildingBlock]], budget: int) -> None:
        try:
            from baybe import Campaign
            from baybe.targets import NumericalTarget
            from baybe.objectives import SingleTargetObjective
            from baybe.parameters import SubstanceParameter
            from baybe.searchspace import SearchSpace
            from baybe.recommenders import TwoPhaseMetaRecommender, FPSRecommender, BotorchRecommender
        except ImportError:
            raise ImportError("Install baybe[chem]: pip install 'mol-healer[opt]'")

        self._domain = domain
        self._bb_to_idx = [{bb.get_smiles(): i for i, bb in enumerate(bbs)} for bbs in domain]

        params = [
            SubstanceParameter(
                name=f'BB{j}',
                data={str(i): bb.get_smiles() for i, bb in enumerate(bbs)},
                encoding=self.encoding,
            )
            for j, bbs in enumerate(domain)
        ]
        target = NumericalTarget(name='Score')
        objective = SingleTargetObjective(target=target)
        searchspace = SearchSpace.from_product(params)
        recommender = TwoPhaseMetaRecommender(
            initial_recommender=FPSRecommender(),
            recommender=BotorchRecommender(),
        )
        self._campaign = Campaign(searchspace, objective, recommender)

    def ask(self) -> List[Tuple[BuildingBlock, ...]]:
        try:
            df = self._campaign.recommend(batch_size=self.batch_size)
        except Exception as e:
            logger.debug("BayesianSequenceOptimizer: recommendation failed (%s); stopping search for this composition.", e)
            return []
        return [
            tuple(self._domain[j][int(df.at[idx, f'BB{j}'])] for j in range(len(self._domain)))
            for idx in df.index
        ]

    def tell(self, results: List[Tuple[Tuple[BuildingBlock, ...], float]]) -> None:
        if not results:
            return
        import pandas as pd
        rows = []
        for bb_tuple, score in results:
            try:
                row = {f'BB{j}': str(self._bb_to_idx[j][bb.get_smiles()]) for j, bb in enumerate(bb_tuple)}
                row['Score'] = score
                rows.append(row)
            except KeyError:
                logger.warning("BB not found in domain during tell(); skipping.")
        if rows:
            self._campaign.add_measurements(pd.DataFrame(rows))

