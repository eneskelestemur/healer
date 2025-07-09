'''
    This script implements the preprocessing of the Building Block datasets.
'''

import os
import json
import healer.utils.utils as utils

from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem.FastSDMolSupplier import FastSDMolSupplier
from rdkit.Chem.rdmolfiles import SDWriter

REACTIONS = utils.load_reactions_from_json('healer/data/reactions/reactions.json')
REACTIONS = [rxn for rxn in REACTIONS if rxn.is_valid()]


def add_rxn_annotations(mol: Chem.Mol) -> Chem.Mol:
    """
        Add a new property called 'rxn_annotations' to the molecule.
        The propery is a dictionary where the keys are the reaction
        names in which the molecule can be found and the values are the
        positions of the molecule in the reaction. Example:
        ```
            {
                'amide coupling-1': [0, 1],
                'sulfoxide': [1],
            }
        ```

        Args:
            mol (Chem.Mol): The molecule to which the property will be added.

        Returns:
            Chem.Mol: The molecule with the new property added.
    """
    rxn_annotations = {}
    for rxn in REACTIONS:
        idx = rxn.get_reactant_index(mol)
        if idx is not None:
            rxn_annotations[rxn.name] = idx
    mol.SetProp('rxn_annotations', json.dumps(rxn_annotations))
    
    return mol

def remove_smaller_fragments(mol: Chem.Mol) -> Chem.Mol:
    """
        Remove the smaller fragments from the molecule.

        Args:
            mol (Chem.Mol): The molecule to process.

        Returns:
            Chem.Mol: The processed molecule.
    """
    frags = Chem.GetMolFrags(mol, asMols=True)
    largest_frag = max(frags, key=lambda x: x.GetNumAtoms())
    return largest_frag

def process_buildingblock_file(input_file: str, verbose: bool=True) -> None:
    """
        Process the building block file and add the 'rxn_annotations' property
        to each molecule.

        Args:
            input_file (str): Path to the input file.
    """
    output_file = os.path.splitext(input_file)[0] + '_processed.sdf'
    suppl = FastSDMolSupplier(input_file)
    
    count = 0
    with SDWriter(output_file) as writer:
        for mol in tqdm(suppl, desc="Processing BBs", unit="molecule", total=len(suppl), disable=not verbose):
            if mol is None:
                continue
            mol = remove_smaller_fragments(mol)
            mol = add_rxn_annotations(mol)
            annotations = mol.GetProp('rxn_annotations')
            if annotations:
                count += 1
                writer.write(mol)

    print(f"Processed {len(suppl)} molecules, annotated {count} with reactions.")
    print(f"Output written to {output_file}")

            
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Process building block files.")
    parser.add_argument("input_file", type=str, help="Path to the input file.")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output.")
    args = parser.parse_args()
    process_buildingblock_file(args.input_file, verbose=args.verbose)

