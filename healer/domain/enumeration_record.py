'''
    Dataclass to represent an Enumeration records.
'''
from dataclasses import dataclass, field
from typing import List, Dict, Any
from rdkit import Chem
from healer.domain.building_block import BuildingBlock


@dataclass
class EnumerationRecord:
    '''
        Holds one enumerated molecule record:
            - product: SMILES of the molecule
            - bbs: list of BuildingBlock objects (length ≥ 0)
            - reaction_names: list of reaction names  (length ≥ 0)
            - props: additional properties (e.g., optimization scores, properties, etc.)
    '''
    product: Chem.Mol
    bbs: List[BuildingBlock] = field(default_factory=list)
    reaction_names: List[str] = field(default_factory=list)
    props: Dict[str, Any] = field(default_factory=dict)

