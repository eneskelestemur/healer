'''
    This file contains helper functions for the project.
'''

import os
import sys
import json
import base64
from typing import List, Union, Tuple
from itertools import chain

import numpy as np
from rdkit import RDLogger, Chem
from rdkit.Chem import AllChem, DataStructs, Draw, rdFMCS

from healer.domain.reaction_template import ReactionTemplate21


RDLogger.DisableLog('rdApp.*')

rdColors = {
    'blue': (0.19, 0.51, 0.70),
    'purple': (0.68, 0.45, 0.8),
    'pink': (0.94, 0.32, 0.65),
    'green': (0.48, 0.68, 0.35),
    'yellow': (0.81, 0.82, 0.0),
    'red': (0.95, 0.42, 0.19),
    'orange': (0.93, 0.69, 0.17)
}

def make_rgb_transparent(
    rgb: Tuple[int, int, int], 
    bg_rgb: Tuple[int, int, int], 
    alpha: float
) -> Tuple[int, int, int]:
    '''
        Makes an RGB color transparent over a background color.
        Source: https://stackoverflow.com/questions/33371939/

        Args:
            rgb: Tuple[int, int, int], RGB color to make transparent.
            bg_rgb: Tuple[int, int, int], background RGB color.
            alpha: float, transparency level (0.0 to 1.0).

        Returns:
            Tuple[int, int, int], resulting RGB color.
    '''
    return tuple(alpha * c1 + (1 - alpha) * c2 for (c1, c2) in zip(rgb, bg_rgb))


def get_sascore(mol: Chem.rdchem.Mol | str) -> float:
    '''
        Calculates the synthetic accessibility score of a molecule.

        Args:
            mol: rdkit.Chem.rdchem.Mol or str, molecule.

        Returns:
            float, synthetic accessibility score, ranging from 1 to 10,
            where 1 is easy to synthesize and 10 is hard to synthesize.
    '''
    sys.path.append(os.path.join(os.environ['CONDA_PREFIX'],'share','RDKit','Contrib'))
    from SA_Score import sascorer
    
    if isinstance(mol, str):
        mol = Chem.MolFromSmiles(mol)
    if mol is None:
        raise ValueError('Invalid molecule')
    return sascorer.calculateScore(mol)

def read_cxsmiles_file(file_path: str, header: bool=True) -> list[str]:
    '''
        Reads a CXSMILES file and returns a list of SMILES strings.
        
        Args:
            file_path: str, path to the file.
            header: bool, whether the file has a header or not.

        Returns:
            list of SMILES strings.
    '''
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")
    
    with open(file_path, 'r') as file:
        smiles_list = []
        if header:
            next(file)  
        for line in file:
            smi = line.split(maxsplit=1)[0]
            smiles_list.append(smi)

    return smiles_list

def load_reactions_from_json(file_path: str) -> list[ReactionTemplate21]:
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

def sanitize_mol(mol: Chem.rdchem.Mol):
    '''
        Sanitizes the mol inplace and returns whether the sanitization was successful 
        along with the sanitization flags.

        Args:
            mol: rdkit.Chem.rdchem.Mol, molecule to sanitize.

        Returns:
            returns a tuple (bool, error_message).
    '''
    flags = Chem.SanitizeMol(mol, catchErrors=True)
    return flags == Chem.SanitizeFlags.SANITIZE_NONE, flags

def get_batch_tversky_sims(
    query_fps: List[DataStructs.ExplicitBitVect], 
    stock_fps: List[DataStructs.ExplicitBitVect], 
    query_factor: float = 0.95, 
    stock_factor: float = 0.05
) -> np.ndarray:
    '''
        Calculates the Tversky similarity between query fingerprints and
        stock fingerprints in batches.

        Args:
            query_fps: list of rdkit DataStructs.ExplicitBitVect, query fingerprints.
            stock_fps: list of rdkit DataStructs.ExplicitBitVect, stock fingerprints.
            query_factor: float, factor for the contribution of the query
                fingerprints to the similarity. Larger values will give more
                weight to the query fingerprints.
            stock_factor: float, factor for the contribution of the stock
                fingerprints to the similarity. Larger values will give more
                weight to the stock fingerprints.

        NOTE: The sum of query_factor and stock_factor should be 1.

        Returns:
            np.ndarray, Tversky similarities, shape=(n_query, n_stock).
    '''
    sims = np.zeros((len(query_fps), len(stock_fps)))
    for i, query_fp in enumerate(query_fps):
        sims[i] = DataStructs.BulkTverskySimilarity(query_fp, stock_fps,
                                                    a=query_factor,
                                                    b=stock_factor)
    return sims

def get_batch_tani_sims(
    query_fps: List[DataStructs.ExplicitBitVect], 
    stock_fps: List[DataStructs.ExplicitBitVect]
) -> np.ndarray:
    '''
        Calculates the Tanimoto similarity between query fingerprints and
        stock fingerprints in batches.

        Args:
            query_fps: list of rdkit DataStructs.ExplicitBitVect, query fingerprints.
            stock_fps: list of rdkit DataStructs.ExplicitBitVect, stock fingerprints.

        Returns:
            np.ndarray, Tanimoto similarities, shape=(n_query, n_stock).
    '''
    sims = np.zeros((len(query_fps), len(stock_fps)))
    for i, query_fp in enumerate(query_fps):
        sims[i] = DataStructs.BulkTanimotoSimilarity(query_fp, stock_fps)
    return sims

def get_tani_sim_fp(fp1, fp2):
    '''
        Calculates the Tanimoto similarity between two fingerprints.

        Args:
            fp1: rdkit.Chem.rdchem.Mol, molecule 1.
            fp2: rdkit.Chem.rdchem.Mol, molecule 2.

        Returns:
            float, Tanimoto similarity.
    '''
    return DataStructs.TanimotoSimilarity(fp1, fp2)

def get_svg_mol(
    mol: Union[Chem.Mol, str],
    sub_mol: Union[Chem.Mol, str] = None, 
    sub_mol_color: str = 'green', 
    legend: str = '', 
    show_idx: bool = False, 
    return_drawing: bool = False, 
    width: int = 350, 
    height: int = 150
) -> str:
    '''
        Get svg image of a molecule with a substructure highlighted.

        Args:
            mol: rdkit.Chem.rdchem.Mol or str, molecule to draw.
            sub_mol: rdkit.Chem.rdchem.Mol or str, substructure to highlight.
            sub_mol_color: str, color to use for highlighting.
            legend: str, legend to display.
            show_idx: bool, whether to show atom indices.
            return_drawing: bool, whether to return the raw SVG drawing text.
            width: int, width of the drawing.
            height: int, height of the drawing.

        Returns:
            str, SVG representation of the molecule.
    '''
    sub_mol_color = rdColors[sub_mol_color]
    if isinstance(mol, str):
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            raise ValueError('Invalid molecule')
    AllChem.Compute2DCoords(mol)
    if sub_mol is not None:
        if isinstance(sub_mol, str):
            sub_struct = Chem.MolFromSmiles(sub_mol)
        else:
            sub_struct = sub_mol
        assert sub_struct is not None, 'Invalid substructure'
        assert mol.HasSubstructMatch(sub_struct), 'Substructure not found'
        hit_atoms = list(mol.GetSubstructMatch(sub_struct))
        hit_bonds = []
        for bond in sub_struct.GetBonds():
            a1 = hit_atoms[bond.GetBeginAtomIdx()]
            a2 = hit_atoms[bond.GetEndAtomIdx()]
            hit_bonds.append(mol.GetBondBetweenAtoms(a1, a2).GetIdx())
    else:
        hit_atoms = []
        hit_bonds = []
    if show_idx:
        for i, atom in enumerate(mol.GetAtoms()):
            atom.SetProp('atomNote', str(i))

    draw_opts = Draw.MolDrawOptions()
    draw_opts.clearBackground = None
    draw_opts.legendFontSize = 10
    draw_opts.legendFraction = 0.15
    drawing = Draw.MolDraw2DSVG(width, height)
    drawing.SetDrawOptions(draw_opts)
    drawing.DrawMolecule(mol, highlightAtoms=hit_atoms, highlightBonds=hit_bonds,
                          highlightAtomColors={i: sub_mol_color for i in hit_atoms},
                          highlightBondColors={i: sub_mol_color for i in hit_bonds},
                          legend=legend)
    drawing.FinishDrawing()
    if return_drawing:
        return drawing.GetDrawingText()
    svg = drawing.GetDrawingText()
    svg = base64.b64encode(svg.encode('utf-8')).decode('utf-8')
    return f"data:image/svg+xml;base64,{svg}"

def get_svg_mol_with_bbs(
    mol: Union[Chem.Mol, str],
    bbs: List[Union[Chem.Mol, str]],
    bb_colors: List[Union[Tuple[int, int, int], str]] = None,
    legend: str = '',
    width: int = 350,
    height: int = 150,
    alpha: float = 1.0,
    bg_color_for_transparency: Tuple[float, float, float] = (1.0, 1.0, 1.0)
) -> str:
    '''
        Draws a molecule with highlighted building blocks.

        Args:
            mol: The molecule to draw.
            bbs: The building blocks to highlight.
            bb_colors: The colors to use for the building blocks. 
                Preset colors: 'blue', 'purple', 'pink', 'green', 'yellow', 'red', 'orange'.
            legend: The legend to display.
            width: The width of the drawing.
            height: The height of the drawing.
            alpha: The transparency level for the building block colors (0.0 to 1.0).
            bg_color_for_transparency: The background color to use when applying transparency.

        Returns:
            The SVG representation of the drawing.
    '''
    num_bbs = len(bbs)
    
    if bb_colors is None:
        rd_color_names = sorted(rdColors.keys(), reverse=True)
        bb_colors = [rdColors[rd_color_names[i % len(rd_color_names)]] for i in range(num_bbs)]
    elif isinstance(bb_colors[0], tuple):
        bb_colors = bb_colors
    elif isinstance(bb_colors[0], str):
        bb_colors = [rdColors[c] for c in bb_colors]

    if alpha < 1.0:
        bb_colors = [make_rgb_transparent(c, bg_color_for_transparency, alpha) for c in bb_colors]

    if num_bbs != len(bb_colors):
        raise ValueError(f'Number of building blocks ({num_bbs}) does not match number of colors ({len(bb_colors)})')

    if isinstance(mol, str):
       mol = Chem.MolFromSmiles(mol)
    AllChem.Compute2DCoords(mol)

    bbs = [Chem.MolFromSmiles(bb) if isinstance(bb, str) else bb for bb in bbs]
    for i, bb in enumerate(bbs):
        if bb is None:
            raise ValueError(f'Invalid building block: bb{i}')

    def _get_mcs_smarts_mol(bb):
        mcs = rdFMCS.FindMCS([mol, bb])
        if mcs.numAtoms == 0:
            raise ValueError(f'No common substructure found between mol and bb{i}')
        return Chem.MolFromSmarts(mcs.smartsString)
    new_bbs = [_get_mcs_smarts_mol(bb) for bb in bbs]

    hit_atoms_list = [list(mol.GetSubstructMatch(bb)) for bb in new_bbs]
    hit_bonds_list = []
    for i, bb in enumerate(new_bbs):
        hit_bonds = []
        for bond in bb.GetBonds():
            a1 = hit_atoms_list[i][bond.GetBeginAtomIdx()]
            a2 = hit_atoms_list[i][bond.GetEndAtomIdx()]
            try:
                hit_bonds.append(mol.GetBondBetweenAtoms(a1, a2).GetIdx())
            except:
                pass
        hit_bonds_list.append(hit_bonds)

    # Draw the molecule with highlighted building blocks
    draw_opts = Draw.MolDrawOptions()
    draw_opts.clearBackground = None
    drawing = Draw.MolDraw2DSVG(width, height)
    drawing.SetDrawOptions(draw_opts)
    highlight_atom_colors = {}
    for i, hit_atoms in enumerate(hit_atoms_list):
        for atom in hit_atoms:
            highlight_atom_colors[atom] = bb_colors[i]
    highlight_bond_colors = {}
    for i, hit_bonds in enumerate(hit_bonds_list):
        for bond in hit_bonds:
            highlight_bond_colors[bond] = bb_colors[i]
    drawing.DrawMolecule(mol, highlightAtoms=list(chain.from_iterable(hit_atoms_list)), highlightBonds=list(chain.from_iterable(hit_bonds_list)),
                          highlightAtomColors=highlight_atom_colors, highlightBondColors=highlight_bond_colors,
                          legend=legend)
    drawing.FinishDrawing()
    svg = drawing.GetDrawingText()
    svg = base64.b64encode(svg.encode('utf-8')).decode('utf-8')
    return f"data:image/svg+xml;base64,{svg}"

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

