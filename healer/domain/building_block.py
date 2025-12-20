'''
    Wrapper for buildingblock molecules to parse the properties automatically.
'''
import json
from typing import Any, Dict
from rdkit import Chem


class BuildingBlock:
    def __init__(self, molecule: Chem.Mol) -> None:
        '''
            Initialize the BuildingBlock with a molecule.
        '''
        self._mol: Chem.Mol = molecule
        self.props: Dict[str, Any] = {
            k: self._parse_value(v)
            for k, v in molecule.GetPropsAsDict().items()
        }

    def __hash__(self) -> int:
        '''
            Hash the building block based on its SMILES representation.
        '''
        return hash(self.get_smiles())

    def __getattr__(self, attr: str) -> Any:
        '''
            Delegate attribute access to the underlying RDKit molecule.
            This allows us to access properties like GetNumAtoms, GetNumBonds, etc.
        '''
        # Prevent infinite recursion during pickle reconstruction
        # Use __dict__ to avoid triggering __getattr__ recursion
        if '_mol' not in self.__dict__ or self.__dict__['_mol'] is None:
            raise AttributeError(f"'{type(self).__name__}' object has no attribute '{attr}'")
        return getattr(self._mol, attr)

    def get_parsed_prop(self, name: str) -> Any:
        '''
            Fetch the parsed Python object for this property.
        '''
        return self.props.get(name, '')
    
    def get_mol(self) -> Chem.Mol:
        '''
            Get the RDKit molecule object.
        '''
        return self._mol
    
    def get_smiles(self) -> str:
        '''
            Get the SMILES representation of the building block.
            Cached on first access to avoid repeated computation.
        '''
        if '_smiles' not in self.__dict__:
            self.__dict__['_smiles'] = Chem.MolToSmiles(self._mol)
        return self.__dict__['_smiles']
    
    def SetProp(self, name: str, value: Any) -> None:
        '''
            Set a property on the underlying Mol *and* update our parsed props.
        '''
        if not isinstance(value, str):
            raw = json.dumps(value)
        else:
            raw = value
        self._mol.SetProp(name, raw)
        self.props[name] = self._parse_value(raw)

    def ClearProp(self, name: str) -> None:
        '''
            Remove a property from the Mol and from parsed props.
        '''
        self._mol.ClearProp(name)
        self.props.pop(name, None)

    def _parse_value(self, val: str) -> Any:
        '''
            Parse a string value into a Python object.
        '''
        try:
            return json.loads(val)
        except (json.JSONDecodeError, TypeError):
            return val
