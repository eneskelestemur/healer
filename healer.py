import os
import abc
import json
import logging
from typing import List, Union, Dict, Tuple, Any

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, rdFingerprintGenerator
from rdkit.Chem.FastSDMolSupplier import FastSDMolSupplier
from collections import namedtuple
from itertools import product as iter_product, chain, compress
from tqdm import tqdm
import utils

# Configure logging
logging.basicConfig(
    # level=logging.INFO,
    format="%(Y-%m-%d %H:%M:%S) [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger(__name__)

# Constants: mapping building-block sources to file paths
BB_PATHS: Dict[str, str] = {
    "US_stock": "buildingblocks/Enamine_Rush-Delivery_Building_Blocks-US_195312cmpd_20240610_processed.sdf",
    "EU_stock": "buildingblocks/Enamine_Rush-Delivery_Building_Blocks-EU_153230cmpd_20240806_processed.sdf",
    "Global_stock": "buildingblocks/Enamine_Building_Blocks_Stock_290951cmpd_20240806_processed.sdf",
    "test": "buildingblocks/test_100_bb_processed.sdf",
}

class _BaseHEALER(abc.ABC):
    '''
        Base HEALER.
    '''
    def __init__(
        self, 
        bb_supplier: str, 
        reaction_tags: Union[List[str], str],
        verbose: int=1,
    ) -> None:
        '''
            Initialize BaseHEALER.

            Args:
                bb_supplier: one of "US_stock", "EU_stock" or "Global_stock"; or path to an SDF file.
                reaction_tags: list of tags or 'all'.
                verbose: verbosity level.
                    - 0: only errors
                    - 1: warnings
                    - 2: info
        '''
        # Set verbosity level
        logger.setLevel(logging.ERROR if verbose == 0 else logging.WARNING if verbose == 1 else logging.INFO)

        # resolve supplier path
        self._supplier_path: str = BB_PATHS.get(bb_supplier, bb_supplier)
        self.bb_supplier = FastSDMolSupplier(self._supplier_path, sanitize=True)
        logger.info("Loaded building blocks from %s", self._supplier_path)

        # load and filter reactions
        self.set_reactions(reaction_tags)
        logger.info("Using %d reactions", len(self.reactions))

        # filter building blocks by reaction compatibility
        self._load_bbs()
        logger.info("Filtered to %d building blocks matching reactions", len(self.mols))

        # fingerprint generator
        self._fp_generator = rdFingerprintGenerator.GetMorganGenerator(
            radius=3, fpSize=2048, includeChirality=True
        )

    @abc.abstractmethod
    def enumerate(self, molecule: Union[str, Chem.Mol]) -> None:
        ...

    @abc.abstractmethod
    def get_results(self, as_dict: bool=False) -> Union[pd.DataFrame, List[Dict[str, Any]]]:
        ...

    @abc.abstractmethod
    def _process_building_blocks(self) -> None:
        ...

    @abc.abstractmethod
    def _prepare_molecule(self) -> None:
        ...
    
    def save_results(self, path: str='results.csv') -> None:
        '''
            Save results to a CSV file.
            
            Args:
                path: path to save the results.
        '''
        df = self.get_results()
        df.to_csv(path, index=False)
    
    def set_reactions(self, reaction_tags: Union[List[str], str]) -> None:
        '''
            Set reactions to use for enumeration.

            Args:
                reaction_tags: a rection tag or a list of tags to filter by.
                    'all' to use all reactions.
        '''
        # reaction data
        all_rxns = utils.load_reactions_from_json('reactions/reactions.json')
        self._reactions = [r for r in all_rxns if r.is_valid()]
        all_tags = list(set(chain(*[r.tags for r in self._reactions])))
        if isinstance(reaction_tags, str):
            if reaction_tags == 'all':
                self.reaction_tags = all_tags
                self.reactions = self._reactions
            else:
                self.reaction_tags = [reaction_tags]
                self.reactions = [r for r in self._reactions if reaction_tags in r.tags]
        else:
            self.reaction_tags = [tag for tag in reaction_tags if tag in all_tags]
            self.reactions = [r for r in self._reactions if any(tag in r.tags for tag in self.reaction_tags)]

    def _get_fingerprints(self, mols: List[Chem.Mol]) -> List[Any]:
        return [self._fp_generator.GetFingerprint(m) for m in mols]
        
    def _get_fingerprint(self, mol: Chem.Mol) -> Any:
        return self._fp_generator.GetFingerprint(mol)

    def _set_molecule(self, molecule: Union[str, Chem.Mol]) -> None:
        if isinstance(molecule, str):
            self.molecule = Chem.MolFromSmiles(molecule)
        else:
            self.molecule = molecule
        flag = Chem.SanitizeMol(self.molecule, catchErrors=True)
        assert flag == Chem.rdmolops.SanitizeFlags.SANITIZE_NONE, f"SanitizeMol failed: {flag}"

    def _load_bbs(self) -> List[Chem.Mol]:
        self.mols = [mol for mol in self.bb_supplier 
                     if set(json.loads(mol.GetProp('rxn_annotations'))).intersection(set([rxn.name for rxn in self.reactions]))]


class SiteHEALER(_BaseHEALER):
    '''
        Site HEALER.
    '''
    def __init__(
            self, 
            bb_supplier: str='US_stock',
            reaction_tags: list[str] | str=['amide coupling', 'amide', 'C-N bond formation', 'C-N',
                                            'alkylation', 'N-arylation', 'azole', 'amination'],
            reactive_sites: list[int]=[],
            rules: dict[str, tuple[int, int]]={
                'MW': (0, 500), # molecular weight
                'HBD': (0, 5), # hydrogen bond donors
                'HBA': (0, 10), # hydrogen bond acceptors
                'TPSA': (0, 200), # topological polar surface area
                'RotB': (0, 10), # rotatable bonds
                'Rings': (0, 10), # number of rings
                'ArRings': (0, 5), # number of aromatic rings
                'Chiral': (0, 5), # number of chiral centers
            },
            struct_rules: list[str]=[],
            verbose: int=1,
    ):
        '''
            Initialize SiteHEALER.

            Args:
                bb_supplier: one of "US_stock", "EU_stock" or "Global_stock"; or path to an SDF file.
                reaction_tags: list of tags or 'all'.
                reactive_sites: list of indices of reactive sites in the molecule.
                rules: dictionary of rules for filtering molecules.
                struct_rules: list of structural rules for filtering molecules.
                verbose: verbosity level, 0 for errors, 1 for warnings, 2 for info.
        '''
        super().__init__(bb_supplier, reaction_tags, verbose)
        self.reactive_sites = reactive_sites
        self.rules = rules
        self.struct_rules = struct_rules

    def enumerate(self, molecule: Union[str, Chem.Mol]) -> None:
        '''
            Enumerate the molecule.

            Args:
                molecule: SMILES string or RDKit molecule object.
        '''
        self._set_molecule(molecule)
        self._prepare_molecule()
        self._process_building_blocks()

        # enumeration loop
        query_fp = self._get_fingerprint(self.molecule)
        self.enumerated_molecules = []
        for bb in tqdm(self._filtered_bb, desc='Enumerating building blocks', total=len(self._filtered_bb), disable=not self.verbose):
            url = bb.GetProp('URL') if bb.HasProp('URL') else ''
            for reaction in self.reactions:
                products = reaction.run_syn(self._prepared_mol, bb)
                if products:
                    for p in products:
                        try:
                            p = Chem.MolFromSmiles(Chem.MolToSmiles(p))
                            product_fp = self._get_fingerprint(p)
                            tani_sim = utils.get_tani_sim_fp(query_fp, product_fp)
                        except:
                            print('Fingerprint calculation for a product failed! Skipping...', 
                                    flush=True)
                            continue
                        Enumeration = namedtuple(
                            'Enumeration', 
                            ['Product', 'Similarity_to_query', 'BB', 'Reaction_name', 'URL']
                        )
                        enum = Enumeration(Chem.MolToSmiles(p),
                                            round(tani_sim, 2),
                                            Chem.MolToSmiles(bb),
                                            reaction.name,
                                            url)
                        self.enumerated_molecules.append(enum)

    def get_results(self, as_dict: bool=False):
        '''
            Get the results as a pandas DataFrame.

            Args:
                as_dict (bool): return the results as a dictionary.

            Returns:
                pd.DataFrame: results.
        '''
        query_molecule_row = {'Product': Chem.MolToSmiles(self.molecule),
                              'Similarity_to_query': 1.0,
                              'BB': '',
                              'Reaction_name': '',
                              'URL': ''}

        if not self.enumerated_molecules:
            print('No enumerated molecules found!')
            df = pd.DataFrame([query_molecule_row])
            if as_dict:
                return df.to_dict(orient='records')
            return df
        
        column_names = ['Product', 'Similarity_to_query', 'BB', 'Reaction_name', 'URL']
        df = pd.DataFrame(self.enumerated_molecules, columns=column_names)
        df = df.drop_duplicates(subset=['Product'], keep='first', ignore_index=True)
        df = pd.concat([pd.DataFrame([query_molecule_row]), df], ignore_index=True)
        df = df.sort_values(by='Similarity_to_query', ascending=False, ignore_index=True)
        df = df.reset_index(drop=True)
        df['ID'] = [f'HEAL_{i:06d}' for i in df.index]
        if as_dict:
            return df.to_dict(orient='records')
        return df
    
    def set_rules(self, **kwargs):
        '''
            Set the rules for the building blocks.

            Args:
                kwargs: dictionary containing the rules.
                    MW: tuple (min, max) -- molecular weight
                    HBD: tuple (min, max) -- hydrogen bond donors
                    HBA: tuple (min, max) -- hydrogen bond acceptors
                    TPSA: tuple (min, max) -- topological polar surface area
                    RotB: tuple (min, max) -- rotatable bonds
                    Rings: tuple (min, max) -- number of rings
                    ArRings: tuple (min, max) -- number of aromatic rings
                    Chiral: tuple (min, max) -- number of chiral centers
        '''
        for key, value in kwargs.items():
            if key in self.rules:
                self.rules[key] = value
            else:
                raise ValueError(f'Invalid rule: {key}')
            
    def _process_building_blocks(self):
        '''
            Process the building blocks to filter out the ones that do not
            satisfy the rules.
        '''
        self._filtered_bb = []
        for bb in tqdm(self.mols, desc='Processing building blocks', total=len(self.mols), disable=not self.verbose):
            if self._check_rules(bb) and self._check_struct_rules(bb):
                for rxn in self.reactions:
                    if rxn.is_reactant(bb):
                        self._filtered_bb.append(bb)
                        break
    
    def _prepare_molecule(self, protect_neighbors: bool=False):
        '''
            Prepare the molecule by adding protection to the atoms that are not 
            part of the reaction site. If protect_neighbors is set to True, the
            neighbors of the reaction site will be protected.

            Args:
                protect_neighbors (bool): protect the neighbors of the reaction site
        '''
        if not hasattr(self, 'molecule'):
            raise ValueError('Molecule not set! Please set the molecule first.')
        
        self._prepared_mol = Chem.MolFromSmiles(Chem.MolToSmiles(self.molecule))
        if self.reaction_sites:
            dont_protect = set()
            for atom in self._prepared_mol.GetAtoms():
                if atom.GetIdx() in self.reaction_sites:
                    dont_protect.add(atom.GetIdx())
                    if not protect_neighbors:
                        for neighbor in atom.GetNeighbors():
                            dont_protect.add(neighbor.GetIdx())
            for atom in self._prepared_mol.GetAtoms():
                if atom.GetIdx() not in dont_protect:
                    atom.SetProp('_protected', '1')
        else:
            print('No reaction sites provided! All atoms will be considered reactive.', flush=True)
            pass

    def _check_struct_rules(self, building_block: Chem.rdchem.Mol):
        '''
            Check if the building block satisfies the structure-based rules.

            Args:
                building_block: rdkit mol object or SMILES string.
        '''
        if not self.struct_rules:
            return True
        for rule in self.struct_rules:
            if not building_block.HasSubstructMatch(Chem.MolFromSmarts(rule)):
                return False
        return True

    def _check_rules(self, building_block: Chem.rdchem.Mol):
        '''
            Check if the building block satisfies the rules.

            Args:
                building_block: rdkit mol or SMILES string.
        '''
        for key, value in self.rules.items():
            if key == 'MW':
                if not value[0] <= Descriptors.MolWt(building_block) <= value[1]:
                    return False
            elif key == 'HBD':
                if not value[0] <= Descriptors.NumHDonors(building_block) <= value[1]:
                    return False
            elif key == 'HBA':
                if not value[0] <= Descriptors.NumHAcceptors(building_block) <= value[1]:
                    return False
            elif key == 'TPSA':
                if not value[0] <= Descriptors.TPSA(building_block) <= value[1]:
                    return False
            elif key == 'RotB':
                if not value[0] <= Descriptors.NumRotatableBonds(building_block) <= value[1]:
                    return False
            elif key == 'Rings':
                if not value[0] <= Descriptors.RingCount(building_block) <= value[1]:
                    return False
            elif key == 'ArRings':
                if not value[0] <= Descriptors.NumAromaticRings(building_block) <= value[1]:
                    return False
            elif key == 'Chiral':
                if not value[0] <= rdMolDescriptors.CalcNumAtomStereoCenters(building_block) <= value[1]:
                    return False
            else:
                raise ValueError(f'Invalid rule: {key}')
            
        return True


class MoleculeHEALER(_BaseHEALER):
    '''
        Molecule HEALER.
    '''
    def __init__(
            self, 
            bb_supplier: str='US_stock', 
            reaction_tags: list[str]=['amide coupling', 'amide', 'C-N bond formation', 'C-N',
                                      'alkylation', 'N-arylation', 'azole', 'amination'],
            custom_comp_sites: list[tuple]=[],
            n_compositions: int=10,
            sim_threshold: float=0.5,
            max_bbs_per_comp: int=-1,
            verbose: int=1,
    ):
        '''
            Initialize MoleculeHEALER.

            Args:
                bb_supplier: one of "US_stock", "EU_stock" or "Global_stock"; or path to an SDF file.
                reaction_tags: list of tags or 'all'.
                custom_comp_sites: list of tuples containing the atom indices for splitting site of 
                    the molecule.
                n_compositions: number of compositions to consider for enumeration.
                sim_threshold: similarity threshold for filtering building blocks.
                max_bbs_per_comp: maximum number of building blocks per fragment.
                    If <= 0, all building blocks will be considered. Otherwise, the similarity
                    threshold will be asjusted to the number of building blocks.
                verbose: verbosity level, 0 for errors, 1 for warnings, 2 for info.
        '''
        super().__init__(bb_supplier, reaction_tags, verbose)
        self.custom_comp_sites = custom_comp_sites
        self.n_compositions = n_compositions
        self.sim_threshold = sim_threshold
        self.max_bbs_per_comp = max_bbs_per_comp

    def enumerate(self, molecule: str | Chem.rdchem.Mol):
        '''
            Enumerate the molecule with building blocks.
            
            Args:
                molecule (str): SMILES string or rdkit mol object. This molecule will be
                                enumerated with building blocks at the given reaction site.
        '''
        self._set_molecule(molecule)
        self._prepare_molecule()
        if not self._compositions: # check if preparation added compositions
            self.enumerated_molecules = []
            return
        self._process_building_blocks()

        # enumerate with building blocks
        query_fp = self._get_fingerprint(self.molecule)
        self.enumerated_molecules = []
        for composition in tqdm(self._filtered_bb, desc='Enumerating building blocks', total=len(self._filtered_bb), disable=not self.verbose):
            for b1, b2 in iter_product(*composition):
                url1 = b1.GetProp('URL') if b1.HasProp('URL') else ''
                url2 = b2.GetProp('URL') if b2.HasProp('URL') else ''
                for reaction in self.reactions:
                    products = reaction.run_syn(b1, b2)
                    if products:
                        for p in products:
                            try:
                                p = Chem.MolFromSmiles(Chem.MolToSmiles(p))
                                product_fp = self._get_fingerprint(p)
                                tani_sim = utils.get_tani_sim_fp(query_fp, product_fp)
                            except:
                                print('Fingerprint calculation of a product failed! Skipping...', 
                                        flush=True)
                                continue
                            # add the product to the list as a named tuple
                            Enumeration = namedtuple(
                                'Enumeration', 
                                ['Product', 'Similarity_to_query', 'BB1', 'BB2', 
                                    'Reaction_name', 'URL1', 'URL2']
                            )
                            enum = Enumeration(Chem.MolToSmiles(p),
                                                round(tani_sim, 2),
                                                Chem.MolToSmiles(b1),
                                                Chem.MolToSmiles(b2),
                                                reaction.name,
                                                url1, url2)
                            self.enumerated_molecules.append(enum)

