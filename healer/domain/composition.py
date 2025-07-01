'''
    Data structures for retrosynthetic compositions.
'''
from dataclasses import dataclass
from typing import List, Optional, Tuple, Any

from rdkit import Chem

from healer.domain.building_block import BuildingBlock
from healer.domain.retro_step import RetroStep


@dataclass
class CompositionPath:
    '''
        Represents a retrosynthetic fragment composition:
        Either provided as a custom tuple of fragments, or as
        a sequence of RetroSteps from which fragments are flattened.
    '''
    steps: Optional[Tuple[RetroStep, ...]] = None
    fragments: Optional[Tuple[Chem.Mol, ...]] = None

    def __hash__(self) -> int:
        # hash solely on the final fragments SMILES
        sigs = tuple(Chem.MolToSmiles(m, canonical=True) for m in self.fragments)
        return hash(sigs)

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, (CompositionPath)):
            return NotImplemented
        return tuple(Chem.MolToSmiles(m, canonical=True) for m in self.fragments) == \
               tuple(Chem.MolToSmiles(m, canonical=True) for m in other.fragments)
    
    def __len__(self) -> int:
        return len(self.fragments) if self.fragments is not None else 0

    def __post_init__(self):
        if self.fragments is not None:
            return
        
        if self.steps is not None:
            frags: List[Chem.Mol] = []
            for step in self.steps:
                frags.extend(step.reactants)
            object.__setattr__(self, 'fragments', tuple(frags))
            return

        raise ValueError('CompositionPath requires at least steps or fragments')

    @classmethod
    def from_fragments(cls, fragments: Tuple[Chem.Mol, ...]) -> 'CompositionPath':
        '''
        Construct a CompositionPath directly from a tuple of fragments,
        without any reaction information.
        '''
        return cls(steps=None, fragments=fragments)


@dataclass
class CompositionWithBBs:
    '''
        Pairs a CompositionPath with a list of BuildingBlocks for each fragment.
    '''
    comp: CompositionPath
    fragment_bbs: Tuple[List[BuildingBlock], ...] 

