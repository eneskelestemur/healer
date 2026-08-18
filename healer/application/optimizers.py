'''
    Optimizer interfaces and implementations for guided enumeration.
'''
import logging
from abc import ABC, abstractmethod
from typing import List, Tuple, Callable, Iterable, Optional

import pandas as pd
from rdkit import Chem

from healer.domain.enumeration_record import EnumerationRecord
from healer.domain.building_block import BuildingBlock
from healer.domain.reaction_template import ReactionTemplate21

logger = logging.getLogger(__name__)

Candidate = Tuple[EnumerationRecord, BuildingBlock, ReactionTemplate21]


class OptimizerError(RuntimeError):
    '''
        Raised when an optimizer cannot continue for reasons other than having
        exhausted its search space. Enumeration stops for the current composition
        and moves on to the next one.
    '''


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

    def evaluate(self, mol: Chem.Mol) -> Optional[float]:
        '''
            Score a single molecule.

            Args:
                mol: molecule to score.

            Returns:
                float score, or None if the molecule could not be scored.
        '''
        return self.evaluate_batch([mol])[0]

    def evaluate_batch(self, mols: List[Chem.Mol]) -> List[Optional[float]]:
        '''
            Score a batch of molecules, using batch_target_fn when it is available
            and falling back to per-molecule calls to target_fn. Scoring failures
            are logged rather than raised, so one bad molecule does not end a run.

            Args:
                mols: molecules to score.

            Returns:
                list of scores, one per input molecule, with None for any that
                could not be scored.
        '''
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
        Interface for stagewise optimizers that prune or reorder candidates at each
        fragment-assembly stage. Pruning happens at two points: `select_candidates`
        drops reactions before they are run, and `filter` drops assembled products
        after. Scoring a molecule requires building it, so anything score-driven
        belongs in `filter`, while cheap building block heuristics belong in
        `select_candidates`, where they save the synthesis work entirely.
    '''

    def select_candidates(self, candidates: Iterable[Candidate], depth: int) -> Iterable[Candidate]:
        '''
            Prune candidate reactions before they are applied. The default keeps
            every candidate.

            Args:
                candidates: (record, building block, reaction) triples, streamed
                    lazily. Return a generator to keep them from being materialized.
                depth: index of the current fragment-assembly stage, starting at 1.

            Returns:
                the candidates to apply at this stage.
        '''
        return candidates

    @abstractmethod
    def filter(self, records: List[EnumerationRecord], depth: int) -> List[EnumerationRecord]:
        '''
            Prune or reorder assembled reaction products at this stage depth.
            Called after reactions have been applied.

            Args:
                records: products assembled at this stage.
                depth: index of the current fragment-assembly stage, starting at 1.

            Returns:
                the records to carry forward as seeds for the next stage.
        '''
        ...


class PassthroughOptimizer(BaseStagewiseOptimizer):
    '''
        Stagewise optimizer that prunes nothing and scores nothing. Used when
        `enumerate` is called without an optimizer, so exhaustive enumeration and
        guided enumeration share one code path.
    '''

    def __init__(self) -> None:
        pass

    def evaluate_batch(self, mols: List[Chem.Mol]) -> List[Optional[float]]:
        '''
            Args:
                mols: molecules to score.

            Returns:
                None for every molecule, leaving records without a score.
        '''
        return [None] * len(mols)

    def filter(self, records: List[EnumerationRecord], depth: int) -> List[EnumerationRecord]:
        '''
            Args:
                records: products assembled at this stage.
                depth: index of the current fragment-assembly stage.

            Returns:
                every record, unchanged.
        '''
        return records


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
        self._domain: Optional[List[List[BuildingBlock]]] = None
        self._bb_to_idx: Optional[List[dict]] = None

    def _set_domain(self, domain: List[List[BuildingBlock]]) -> None:
        '''
            Store the domain and index each pool by SMILES so proposals can be
            mapped back to building block positions.

            Args:
                domain: building block pool for each fragment.
        '''
        self._domain = domain
        self._bb_to_idx = [{bb.get_smiles(): i for i, bb in enumerate(bbs)} for bbs in domain]

    def _to_indices(self, bb_tuple: Tuple[BuildingBlock, ...]) -> Optional[Tuple[int, ...]]:
        '''
            Map a building block tuple to its per-fragment indices in the domain.

            Args:
                bb_tuple: one building block per fragment.

            Returns:
                the indices, or None if any building block is not in the domain.
        '''
        try:
            return tuple(self._bb_to_idx[j][bb.get_smiles()] for j, bb in enumerate(bb_tuple))
        except KeyError:
            logger.warning("BB not found in domain during tell(); skipping.")
            return None

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
        '''
            Start a new search. Called once per composition, so each composition is
            searched independently.

            Args:
                domain: building block pool for each fragment, already truncated by
                    `prepare_domain`.
                budget: evaluation budget for this composition, 0 when unlimited.
        '''
        ...

    @abstractmethod
    def ask(self) -> List[Tuple[BuildingBlock, ...]]:
        '''
            Propose the next building block combinations to assemble and score.

            Returns:
                list of tuples holding one building block per fragment. An empty
                list ends the search for the current composition.
        '''
        ...

    @abstractmethod
    def tell(self, results: List[Tuple[Tuple[BuildingBlock, ...], float]]) -> None:
        '''
            Report the scores of the last proposed combinations.

            Args:
                results: (bb_tuple, score) pairs. Combinations whose product could
                    not be scored are left out, so this may be shorter than the
                    matching `ask` result.
        '''
        ...



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
        '''
            Args:
                beam_width: number of highest scoring intermediates to keep at each
                    stage.
                target_fn: function that scores a single Mol.
                batch_target_fn: function that scores a list of Mols in one call.
        '''
        super().__init__(target_fn=target_fn, batch_target_fn=batch_target_fn)
        self.beam_width = beam_width

    def filter(self, records: List[EnumerationRecord], depth: int) -> List[EnumerationRecord]:
        '''
            Keep the `beam_width` highest scoring records. All records are kept when
            none of them could be scored.

            Args:
                records: products assembled at this stage.
                depth: index of the current fragment-assembly stage.

            Returns:
                the highest scoring records, at most `beam_width` of them.
        '''
        if not records:
            return records
        mols = [rec.product for rec in records]
        scores = self.evaluate_batch(mols)
        scored = [(rec, s) for rec, s in zip(records, scores, strict=True) if s is not None]
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
        '''
            Args:
                population_size: number of combinations proposed per generation.
                mutation_percent_genes: percentage of genes mutated per solution.
                crossover_type: PyGAD crossover operator name.
                keep_elitism: number of best solutions carried into the next
                    generation unchanged.
                random_seed: seed for reproducible runs.
                target_fn: function that scores a single Mol.
                batch_target_fn: function that scores a list of Mols in one call.
                max_domain_per_frag: maximum building blocks per fragment to search
                    over. None searches the full pools.
        '''
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
        self._score_cache: dict = {}

    def init_search(self, domain: List[List[BuildingBlock]], budget: int) -> None:
        '''
            Build a PyGAD population over building block indices, one gene per
            fragment. Fitness is served from the scores passed to `tell`, so the
            GA never calls the scoring function itself.

            Args:
                domain: building block pool for each fragment.
                budget: evaluation budget, unused since generations are driven by
                    the enumeration loop.

            Raises:
                ImportError: if pygad is not installed.
        '''
        try:
            import pygad
        except ImportError:
            raise ImportError("Install pygad: pip install 'mol-healer[opt]'") from None

        self._set_domain(domain)
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
        '''
            Map the current population onto building block combinations.

            Returns:
                one tuple per solution in the population.
        '''
        return [
            tuple(self._domain[j][int(idx)] for j, idx in enumerate(sol))
            for sol in self._ga.population
        ]

    def tell(self, results: List[Tuple[Tuple[BuildingBlock, ...], float]]) -> None:
        '''
            Load the reported scores into the fitness cache and advance the GA by
            one generation. The cache accumulates across generations so elites
            carried over from earlier ones keep their fitness. Combinations that
            have never been scored fall back to -1e9 and are selected against.

            Args:
                results: (bb_tuple, score) pairs for the last proposed combinations.
        '''
        for bb_tuple, score in results:
            key = self._to_indices(bb_tuple)
            if key is not None:
                self._score_cache[key] = score
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
        '''
            Args:
                batch_size: number of combinations recommended per round.
                encoding: BayBE SubstanceParameter encoding used to featurize
                    building blocks.
                target_fn: function that scores a single Mol.
                batch_target_fn: function that scores a list of Mols in one call.
                max_domain_per_frag: maximum building blocks per fragment to search
                    over. None searches the full pools, which BayBE may not be able
                    to enumerate.
        '''
        super().__init__(
            target_fn=target_fn,
            batch_target_fn=batch_target_fn,
            max_domain_per_frag=max_domain_per_frag,
        )
        self.batch_size = batch_size
        self.encoding = encoding
        self._campaign = None
        self._exhausted_errors: Tuple[type, ...] = ()

    def init_search(self, domain: List[List[BuildingBlock]], budget: int) -> None:
        '''
            Build a BayBE campaign whose search space is the product of the building
            block pools, with one SubstanceParameter per fragment.

            Args:
                domain: building block pool for each fragment.
                budget: evaluation budget, unused since rounds are driven by the
                    enumeration loop.

            Raises:
                ImportError: if baybe is not installed.
        '''
        try:
            from baybe import Campaign
            from baybe.targets import NumericalTarget
            from baybe.objectives import SingleTargetObjective
            from baybe.parameters import SubstanceParameter
            from baybe.searchspace import SearchSpace
            from baybe.recommenders import TwoPhaseMetaRecommender, FPSRecommender, BotorchRecommender
            from baybe.exceptions import (
                EmptySearchSpaceError, NoRecommendersLeftError, NotEnoughPointsLeftError
            )
        except ImportError:
            raise ImportError("Install baybe[chem]: pip install 'mol-healer[opt]'") from None

        self._exhausted_errors = (
            EmptySearchSpaceError, NoRecommendersLeftError, NotEnoughPointsLeftError
        )
        self._set_domain(domain)

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
        '''
            Recommend the next combinations from the campaign.

            Returns:
                `batch_size` tuples holding one building block per fragment, or an
                empty list once the search space has been exhausted.

            Raises:
                OptimizerError: if the recommendation fails for any other reason.
        '''
        try:
            df = self._campaign.recommend(batch_size=self.batch_size)
        except self._exhausted_errors as e:
            logger.debug("BayesianSequenceOptimizer: search space exhausted (%s).", e)
            return []
        except Exception as e:
            raise OptimizerError(f"BayBE recommendation failed: {e}") from e
        return [
            tuple(self._domain[j][int(df.at[idx, f'BB{j}'])] for j in range(len(self._domain)))
            for idx in df.index
        ]

    def tell(self, results: List[Tuple[Tuple[BuildingBlock, ...], float]]) -> None:
        '''
            Add the reported scores to the campaign as measurements, which the
            surrogate model is refit on before the next recommendation.

            Args:
                results: (bb_tuple, score) pairs for the last proposed combinations.
        '''
        if not results:
            return
        rows = []
        for bb_tuple, score in results:
            indices = self._to_indices(bb_tuple)
            if indices is None:
                continue
            row = {f'BB{j}': str(idx) for j, idx in enumerate(indices)}
            row['Score'] = score
            rows.append(row)
        if rows:
            self._campaign.add_measurements(pd.DataFrame(rows))

