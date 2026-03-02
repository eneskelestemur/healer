'''
    Wrapper for buildingblock molecules to parse the properties automatically.
'''
import json
from typing import Any, Dict, Optional
from rdkit import Chem
from rdkit.DataStructs.cDataStructs import ExplicitBitVect


class BuildingBlock:
    def __init__(self, molecule: Chem.Mol) -> None:
        '''
            Initialize the BuildingBlock with a molecule.
        '''
        self._smiles: str = Chem.MolToSmiles(molecule)
        self._mol: Optional[Chem.Mol] = None      # lazy, reconstructed on demand
        self.num_heavy_atoms: int = molecule.GetNumHeavyAtoms()
        self.fingerprint: Optional[ExplicitBitVect] = None
        self.props: Dict[str, Any] = {
            k: self._parse_value(v)
            for k, v in molecule.GetPropsAsDict().items()
        }

    def __hash__(self) -> int:
        '''
            Hash the building block based on its SMILES representation.
        '''
        return hash(self._smiles)

    def __getattr__(self, attr: str) -> Any:
        '''
            Delegate attribute access to the underlying RDKit molecule.
            This allows us to access properties like GetNumAtoms, GetNumBonds, etc.
        '''
        # Prevent infinite recursion during pickle reconstruction
        if '_smiles' not in self.__dict__:
            raise AttributeError(f"'{type(self).__name__}' object has no attribute '{attr}'")
        return getattr(self.mol, attr)

    def get_parsed_prop(self, name: str) -> Any:
        '''
            Fetch the parsed Python object for this property.
        '''
        return self.props.get(name, '')

    @property
    def mol(self) -> Chem.Mol:
        '''
            RDKit molecule object. Lazily reconstructed from SMILES
            if it has been evicted or not yet created.
        '''
        if self._mol is None:
            self._mol = Chem.MolFromSmiles(self._smiles)
        return self._mol
    
    def evict(self) -> None:
        '''
            Drop the cached Mol to free memory.
            It will be lazily reconstructed on next access of ``mol``.
        '''
        self._mol = None

    def get_smiles(self) -> str:
        '''
            Get the canonical SMILES representation of the building block.
        '''
        return self._smiles
    
    def SetProp(self, name: str, value: Any) -> None:
        '''
            Set a property on the underlying Mol *and* update our parsed props.
        '''
        if not isinstance(value, str):
            raw = json.dumps(value)
        else:
            raw = value
        self.mol.SetProp(name, raw)
        self.props[name] = self._parse_value(raw)

    def ClearProp(self, name: str) -> None:
        '''
            Remove a property from the Mol and from parsed props.
        '''
        self.mol.ClearProp(name)
        self.props.pop(name, None)

    def _parse_value(self, val: str) -> Any:
        '''
            Parse a string value into a Python object.
        '''
        try:
            return json.loads(val)
        except (json.JSONDecodeError, TypeError):
            return val
