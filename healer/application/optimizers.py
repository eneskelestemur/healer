'''
    Optimizer interfaces to interact with the enumerator.
'''
from abc import ABC, abstractmethod
from typing import List, Tuple, Callable

from rdkit import Chem

from healer.domain.enumeration_record import EnumerationRecord
from healer.domain.building_block import BuildingBlock
from healer.domain.reaction_template import ReactionTemplate21


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
        candidates: List[Tuple[EnumerationRecord, BuildingBlock, ReactionTemplate21]],
        depth: int
    ) -> List[Tuple[EnumerationRecord, BuildingBlock, ReactionTemplate21]]:
        '''
            Return a (possibly pruned or reordered) list of candidates for this stage.
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

