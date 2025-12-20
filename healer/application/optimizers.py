'''
    Optimizer interfaces to interact with the enumerator.
'''
from abc import ABC, abstractmethod
from typing import List, Tuple, Callable, Iterable

from rdkit import Chem

from healer.domain.enumeration_record import EnumerationRecord
from healer.domain.building_block import BuildingBlock
from healer.domain.reaction_template import ReactionTemplate21

# Type alias for candidate tuples
Candidate = Tuple[EnumerationRecord, BuildingBlock, ReactionTemplate21]


class BaseOptimizer(ABC):
    '''Holds the expensive scoring function for all optimizers.'''
    def __init__(self, target_fn: Callable[[Chem.Mol], float]) -> None:
        """
            target_fn: a function that takes an RDKit Mol and returns a float score.
        """
        
        self.target_fn = target_fn


class BaseStagewiseOptimizer(BaseOptimizer, ABC):
    '''
        Interface for stagewise optimizers that prune or reorder candidates
        at each fragment-assembly stage.
    '''
    @abstractmethod
    def filter(
        self,
        candidates: Iterable[Candidate],
        depth: int
    ) -> Iterable[Candidate]:
        '''
            Return a (possibly pruned or reordered) iterable of candidates for this stage.
            
            The input `candidates` may be a generator (for memory efficiency) or a list.
            The output can also be either:
              - A generator: for streaming/memory-efficient filtering
              - A list: when the filter needs to sort, rank, or sample candidates
            
            The downstream consumer (_apply_candidates) accepts any iterable.
            
            Example streaming filter:
                def filter(self, candidates, depth):
                    for c in candidates:
                        if self.should_keep(c):
                            yield c
            
            Example list-based filter (for sorting):
                def filter(self, candidates, depth):
                    cands = list(candidates)  # materialize once
                    cands.sort(key=lambda c: self.score(c))
                    return cands[:self.top_k]
        '''
        ...


class BaseSequenceOptimizer(BaseOptimizer, ABC):
    '''
        Interface for full-sequence optimizers using ask/tell over BB-index tuples.
    '''
    @abstractmethod
    def init_search(
        self,
        domain: List[List[BuildingBlock]],
        budget: int
    ) -> None:
        '''
            Initialize search with given domain and evaluation budget.
        '''
        ...

    @abstractmethod
    def ask(self) -> List[Tuple[BuildingBlock, ...]]:
        '''
            Propose a list of BB-index tuples to evaluate next.
        '''
        ...

    @abstractmethod
    def tell(
        self,
        results: List[Tuple[Tuple[int, ...], float]]  # (seq, score)
    ) -> None:
        '''
            Provide the optimizer with scores for the last asked sequences.
        '''
        ...

