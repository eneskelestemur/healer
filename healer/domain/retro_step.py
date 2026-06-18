'''
    Data class for a retrosynthetic step.
'''
from dataclasses import dataclass
from typing import Tuple

from rdkit import Chem

from healer.domain.reaction_template import ReactionTemplate21


@dataclass
class RetroStep:
    '''
        Represents one retrosynthetic split:
        - product: the molecule being split
        - reaction: the ReactionTemplate21 used
        - reactants: tuple of resulting reactant fragments
    '''
    product: Chem.Mol
    reaction: ReactionTemplate21
    reactants: Tuple[Chem.Mol, Chem.Mol]

