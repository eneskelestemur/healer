'''
    Main to run enumerator on a given molecule.
'''

import os
import utils
import argparse

from enumerator import MoleculeEnumerator, SiteEnumerator
from rdkit import Chem
from rdkit.Chem.FastSDMolSupplier import FastSDMolSupplier

# TODO: add all the necessary arguments

def main():
    parser = argparse.ArgumentParser(description='Enumerate a given molecule.')
    parser.add_argument('input', type=str,
                        help='Either an .sdf file or a SMILES string.')
    parser.add_argument('--output', type=str, required=False,
                        default='enumerated_molecules.csv', 
                        help='Path to output file, csv format.')
    parser.add_argument('--bb_source', type=str, required=False, 
                        default='US_stock', 
                        help='One of the following: "US_stock", "EU_stock",'
                             '"Global_stock", "NoRush_stock", or a path to '
                             'a custom BB file.')
    parser.add_argument('--reaction_tags', type=str, required=False, nargs='+', 
                        default='"amide coupling" "amide" "amination"', 
                        help='List of reactions to use. If "all" is passed, '
                             'all reactions will be used.')
    parser.add_argument('--enumerator', type=str, required=False,
                        default='Molecule', 
                        help='One of the following: "Molecule" or "Site".')
    parser.add_argument('--max_enumerations', type=int, required=False,
                        default=100, 
                        help='Maximum number of splits to enumerate.'
                             'Used only for Molecule Enumerator.')
    parser.add_argument('--sim_threshold', type=float, required=False,
                        default=0.5, 
                        help='Similarity threshold.'
                             'Used only for Molecule Enumerator.')
    # parser.add_argument('--mw_range', type=str, required=False,)
    parser.add_argument()
    args = parser.parse_args()

    if args.input.endswith('.sdf'):
            mol = FastSDMolSupplier(args.input)[0]
    else:
        mol = Chem.MolFromSmiles(args.input)

    if args.reaction_tags == 'all':
        reactions = utils.load_reactions_from_json('reactions/reactions.json')
        reactions = [r for r in reactions if r.is_valid()]
        reaction_tags = [r.tags for r in reactions if r.is_valid()]
        reaction_tags = list(set([tag for tags in reaction_tags for tag in tags]))
    else:
        reaction_tags = args.reaction_tags

    if args.enumerator.lower() == 'molecule':
        enumerator = MoleculeEnumerator(
             mol, args.bb_source, 
             reaction_tags, [], 
             args.max_enumerations,
             args.sim_threshold
        )
        enumerator.enumerate()
        enumerator.save_results(args.output)

    elif args.enumerator.lower() == 'site':
        enumerator = SiteEnumerator(
             mol, args.bb_source, 
             [], reaction_tags,
            #  rules=rules,
            #  struct_rules=struct_rules,
        )
    

if __name__ == '__main__':
    main()
            
