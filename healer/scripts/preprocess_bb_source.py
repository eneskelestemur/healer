'''
    This script implements the preprocessing of the Building Block datasets.
'''

import os
import json
import zipfile
from pathlib import Path
import healer.utils.utils as utils

from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem.FastSDMolSupplier import FastSDMolSupplier
from rdkit.Chem.rdmolfiles import SDWriter

# Get reactions path from package data
_HEALER_PKG = Path(__file__).parent.parent
_REACTIONS_FILE = _HEALER_PKG / 'data' / 'reactions' / 'reactions.json'

REACTIONS = utils.load_reactions_from_json(str(_REACTIONS_FILE))
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

def extract_zip_if_needed(input_file: str, verbose: bool = True) -> str:
    """
        Extract ZIP file if the input is a ZIP file, otherwise return the input file path.
        
        Args:
            input_file (str): Path to the input file (could be ZIP or SDF).
            verbose (bool): Whether to print extraction progress.
            
        Returns:
            str: Path to the extracted SDF file or original file if not a ZIP.
    """
    input_path = Path(input_file)
    
    if input_path.suffix.lower() == '.zip':
        extract_dir = input_path.parent / input_path.stem
        extract_dir.mkdir(exist_ok=True)
        
        if verbose:
            print(f"Extracting {input_file} to {extract_dir}")
            
        with zipfile.ZipFile(input_file, 'r') as zip_ref:
            zip_ref.extractall(extract_dir)
            
        # Find the SDF file in the extracted directory
        sdf_files = list(extract_dir.glob('*.sdf'))
        if not sdf_files:
            raise FileNotFoundError(f"No SDF file found in extracted ZIP: {input_file}")
        
        if len(sdf_files) > 1:
            print(f"Warning: Multiple SDF files found. Using: {sdf_files[0]}")
            
        return str(sdf_files[0])
    
    return input_file

def main(input_file: str, output_dir: str = None, verbose: bool=True) -> None:
    """
        Process the building block file and add the 'rxn_annotations' property
        to each molecule. Automatically extracts ZIP files if needed.

        Args:
            input_file (str): Path to the input file (SDF or ZIP containing SDF).
            output_dir (str): Directory to save the processed file. If None, saves
                              in the same directory as the input file.
            verbose (bool): Whether to print progress information.
    """
    # Extract ZIP file if needed
    sdf_file = extract_zip_if_needed(input_file, verbose)
    
    # Determine output path
    sdf_path = Path(sdf_file)
    if output_dir:
        out_dir = Path(output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        output_file = out_dir / (sdf_path.stem + '_processed.sdf')
    else:
        output_file = sdf_path.parent / (sdf_path.stem + '_processed.sdf')
    
    suppl = FastSDMolSupplier(sdf_file)
    
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


def cli():
    """CLI entry point for preprocess-bb command."""
    import argparse
    parser = argparse.ArgumentParser(
        prog="preprocess-bb",
        description="Preprocess building block files for HEALER."
    )
    parser.add_argument("input_file", type=str, help="Path to the input SDF or ZIP file.")
    parser.add_argument("-o", "--output-dir", type=str, default=None,
                        help="Output directory for processed file. Defaults to same directory as input.")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose output.")
    args = parser.parse_args()
    main(args.input_file, output_dir=args.output_dir, verbose=args.verbose)


if __name__ == "__main__":
    cli()
