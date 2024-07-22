'''
    This file contains helper functions for the project.
'''

import json

from rdkit import RDLogger, Chem
from rdkit.Chem import AllChem, DataStructs
from reaction import ReactionTemplate21

RDLogger.DisableLog('rdApp.*')


def load_reactions_from_json(file_path):
    '''
        Loads reactions from a json file.

        Args:
            file_path: str, path to the json file.

        Returns:
            list of ReactionTemplate21 objects.

        link to reaction data source:
            https://github.com/datamol-io/datamol/blob/9e94d026534b2a534250dfbfab924ab6f089e477/datamol/data/reactions.json
    '''
    with open(file_path, 'r') as file:
        data = json.load(file)

    reactions = []
    for key, values in data.items():
        reaction = ReactionTemplate21.from_reaction_json(name=key, reaction_json=values)
        reactions.append(reaction)

    return reactions

def get_tani_sim(mol1, mol2):
    '''
        Calculates the Tanimoto similarity between two molecules.

        Args:
            mol1: rdkit.Chem.rdchem.Mol, molecule 1.
            mol2: rdkit.Chem.rdchem.Mol, molecule 2.

        Returns:
            float, Tanimoto similarity.
    '''
    # calculate similarity
    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 3, nBits=2048)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 3, nBits=2048)
    return DataStructs.TanimotoSimilarity(fp1, fp2)

def split_molecule(mol: Chem.rdchem.Mol | str, split_site: tuple[int, int]):
    '''
        Splits the molecule at the given bond.

        Args:
            mol: rdkit mol object or smiles string.
            split_site: tuple of two integers, bond to split.

        Returns:
            tuple of rdkit mol objects, split molecules.
    '''
    if isinstance(mol, str):
        mol = Chem.MolFromSmiles(mol)
    
    # split the molecule
    with Chem.RWMol(mol) as m:
        m.RemoveBond(split_site[0], split_site[1])
    m = m.GetMol()
    frags = Chem.GetMolFrags(m, asMols=True, sanitizeFrags=True)
    return frags

def _dummy2dummy(mol: Chem.rdchem.Mol):
    '''
        Possibly a useful function, but not used in the project. 
        Helper function to replace [*] with [*H5] in the molecule.

        NOTE: [*] causes problems in the reaction SMARTS since they
            are considered as wildcard atoms.

        Args:
            mol: rdkit mol object.

        Returns:
            mol: rdkit mol object with [*] replaced by [*H5].
    '''
    if '[*]' in Chem.MolToSmiles(mol):
        return Chem.MolFromSmiles(Chem.MolToSmiles(mol).replace('[*]', '[*H5]'))
    elif '*' in Chem.MolToSmiles(mol):
        return Chem.MolFromSmiles(Chem.MolToSmiles(mol).replace('*', '[*H5]'))
    else:
        return ValueError('No dummy atom found.')

