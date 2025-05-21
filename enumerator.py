'''
    This script contains the Enumerator class which is used to enumerate
    a given molecule with Enamine building blocks or any other building block
    source. 
'''

# libraries
import os
import abc
import json
import utils
import pandas as pd
import numpy as np

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, rdFingerprintGenerator
from rdkit.Chem.FastSDMolSupplier import FastSDMolSupplier

from collections import namedtuple
from itertools import product as iter_product, chain, compress
from tqdm import tqdm



class _BaseEnumerator(abc.ABC):
    '''
        Base Enumerator.
    '''
    def __init__(self, bb_supplier: str, reaction_tags: list[str] | str):
        '''
            Initialize the BaseEnumerator object.

            Args:
                bb_supplier (str): "US_stock", "EU_stock" or "Global_stock". A custom path to a file
                                    containing building blocks can also be provided.
                reaction_tags (list): list of reaction tags to consider for the enumeration.
                                      "all" will consider all reaction tags, but it will slow down
                                      the enumeration.
        '''
        # building blocks
        if bb_supplier == 'US_stock':
            self._supplier_path = 'buildingblocks/Enamine_Rush-Delivery_Building_Blocks-US_195312cmpd_20240610_processed.sdf'
        elif bb_supplier == 'EU_stock':
            self._supplier_path = 'buildingblocks/Enamine_Rush-Delivery_Building_Blocks-EU_153230cmpd_20240806_processed.sdf'
        elif bb_supplier == 'Global_stock':
            self._supplier_path = 'buildingblocks/Enamine_Building_Blocks_Stock_290951cmpd_20240806_processed.sdf'
        elif bb_supplier == 'test':
            self._supplier_path = 'buildingblocks/test_100_bb_processed.sdf'
        else:
            self._supplier_path = bb_supplier

        self.bb_supplier = FastSDMolSupplier(self._supplier_path, sanitize=True)

        # reaction data
        self._reactions = utils.load_reactions_from_json('reactions/reactions.json')
        self._reactions = [reaction for reaction in self._reactions if reaction.is_valid()]
        if isinstance(reaction_tags, str) and reaction_tags == 'all':
            self.reaction_tags = list(set(chain(*[reaction.tags for reaction in self._reactions])))
            self.reactions = self._reactions
        else:
            self.reaction_tags = reaction_tags
            self.reactions = [reaction for reaction in self._reactions 
                              if any(tag in reaction.tags for tag in self.reaction_tags)]
            
        # load the building blocks
        if self.reaction_tags == 'all':
            self.mols = [mol for mol in self.bb_supplier if mol is not None]
        else:
            def rxn_intersection(mol: Chem.Mol):
                return set(json.loads(mol.GetProp('rxn_annotations'))).intersection(set([rxn.name for rxn in self.reactions]))
            self.mols = [mol for mol in self.bb_supplier if rxn_intersection(mol)]

        # fingerprint generator
        self._fp_generator = rdFingerprintGenerator.GetMorganGenerator(
            radius=3, fpSize=2048, includeChirality=True
        )

    @abc.abstractmethod
    def enumerate(self):
        raise NotImplementedError
    
    @abc.abstractmethod
    def get_results(self): 
        raise NotImplementedError

    @abc.abstractmethod
    def save_results(self):
        raise NotImplementedError

    @abc.abstractmethod
    def _process_building_blocks(self):
        raise NotImplementedError

    @abc.abstractmethod
    def _prepare_molecule(self):
        raise NotImplementedError
    
    def _get_fingerprints(self, mols: list[Chem.Mol]):
        return [self._fp_generator.GetFingerprint(mol) for mol in mols]
        
    def _get_fingerprint(self, mol: Chem.Mol):
        return self._fp_generator.GetFingerprint(mol)
    
    def _set_molecule(self, molecule: str | Chem.rdchem.Mol):
        if isinstance(molecule, str):
            self.molecule = Chem.MolFromSmiles(molecule)
        else:
            self.molecule = molecule
        flag = Chem.SanitizeMol(self.molecule, catchErrors=True)
        assert flag == Chem.rdmolops.SanitizeFlags.SANITIZE_NONE, f'Molecule sanitization failed with flags: {flag}'


class SiteEnumerator(_BaseEnumerator):
    '''
        Site Enumerator class is used to enumerate a given molecule with building blocks.
    '''
    def __init__(
            self, 
            building_blocks: str='US_stock',
            reaction_sites: list[int]=[],
            reaction_tags: list[str] | str=['amide coupling', 'amide', 'C-N bond formation', 'C-N',
                                            'alkylation', 'N-arylation', 'azole', 'amination'],
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
            verbose: bool=True
    ):
        '''
            Initialize the Enumerator object.

            Args:
                building_blocks: "US_stock", "EU_stock", or "Global_stock". 
                                A custom path to a file containing building blocks can 
                                also be provided.
                reaction_sites: list of atom indices to consider for the enumeration. 
                                If no reaction site is provided, any possible reaction
                                site will be considered.
                reaction_tags: list of reaction tags to consider for the enumeration.
                rules: dictionary containing the rules if the method is rules.
                        - MW: tuple (min, max) -- molecular weight
                        - HBD: tuple (min, max) -- hydrogen bond donors
                        - HBA: tuple (min, max) -- hydrogen bond acceptors
                        - TPSA: tuple (min, max) -- topological polar surface area
                        - RotB: tuple (min, max) -- rotatable bonds
                        - Rings: tuple (min, max) -- number of rings
                        - ArRings: tuple (min, max) -- number of aromatic rings
                        - Chiral: tuple (min, max) -- number of chiral centers
                struct_rules: list of structure-based rules. List of SMILES/SMARTS to 
                                     include in the building blocks as substructure.
                verbose: print the compositions of the molecule.
        '''
        super().__init__(building_blocks, reaction_tags)
        self.reaction_sites = reaction_sites
        self.rules = rules
        self.struct_rules = struct_rules
        self.verbose = verbose

    def enumerate(self, molecule: str | Chem.Mol):
        '''
            Enumerate the molecule with building blocks.
        
            Args:
                molecule (str): SMILES string or rdkit mol object. This molecule will be
                                enumerated with building blocks at the given reaction site.
        '''
        self._set_molecule(molecule)
        self._prepare_molecule()
        self._process_building_blocks()
        
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

    def save_results(self, path: str='enumerated_molecules.csv'):
        '''
            Save the results to a file in CSV format.
            
            Columns: Product, Similarity_to_query, BB, Reaction_name

            Args:
                path (str): path to the file.
        '''
        df = self.get_results()
        df.to_csv(path, index=False)

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

    def _check_struct_rules(self, building_block: Chem.rdchem.Mol | str):
        '''
            Check if the building block satisfies the structure-based rules.

            Args:
                building_block: rdkit mol object or SMILES string.
        '''
        if isinstance(building_block, str):
            building_block = Chem.MolFromSmiles(building_block)

        if not self.struct_rules:
            return True
        for rule in self.struct_rules:
            if not building_block.HasSubstructMatch(Chem.MolFromSmarts(rule)):
                return False
        return True

    def _check_rules(self, building_block: Chem.rdchem.Mol | str):
        '''
            Check if the building block satisfies the rules.

            Args:
                building_block: rdkit mol or SMILES string.
        '''
        if isinstance(building_block, str):
            building_block = Chem.MolFromSmiles(building_block)

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
        

class MoleculeEnumerator(_BaseEnumerator):
    '''
        Molecule Enumerator class is used to enumerate a given molecule with building blocks.
    '''
    def __init__(
            self, 
            building_blocks: str='US_stock', 
            reaction_tags: list[str]=['amide coupling', 'amide', 'C-N bond formation', 'C-N',
                                      'alkylation', 'N-arylation', 'azole', 'amination'],
            custom_comp_sites: list[tuple]=[],
            n_compositions: int=10,
            sim_threshold: float=0.5,
            max_bbs_per_comp: int=-1,
            verbose: bool=True
    ):
        '''
            Initialize the MoleculeEnumerator object.

            Args:
                building_blocks (str): "US_stock", "EU_stock", or "Global_stock".
                    A custom path to a file containing building blocks can also be provided.
                reaction_tags (list): list of reaction tags to consider for the enumeration.
                    "all" will consider all reaction tags, but it will slow down the enumeration.
                custom_comp_sites (list(tuple)): list of tuples containing the atom indices for
                    splitting the molecule. Each tuple represents a composition site.
                n_compositions (int): number of compositions of the molecule to enumerate.
                sim_threshold (float): similarity threshold.
                max_bbs_per_comp (int): maximum number of building blocks per composition.
                    If <= 0, all building blocks will be considered. Otherwise, the similarity
                    threshold will be asjusted to the number of building blocks.
                verbose (bool): print the compositions of the molecule.
        '''
        super().__init__(building_blocks, reaction_tags)
        self.custom_comp_sites = custom_comp_sites
        self.n_compositions = n_compositions
        self.sim_threshold = sim_threshold
        self.max_bbs_per_comp = max_bbs_per_comp
        self.verbose = verbose

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
                              'BB1': '',
                              'BB2': '',
                              'Reaction_name': '',
                              'URL1': '',
                              'URL2': ''}
        if not self.enumerated_molecules:
            print('No enumerated molecules found! ')
            df = pd.DataFrame([query_molecule_row])
            df['ID'] = [f'HEAL_{i:06d}' for i in df.index]
            if as_dict:
                return df.to_dict(orient='records')
            return df
        
        column_names = ['Product', 'Similarity_to_query', 'BB1', 'BB2', 'Reaction_name', 'URL1', 'URL2']
        df = pd.DataFrame(self.enumerated_molecules, columns=column_names)
        df = df.drop_duplicates(subset=['Product'], keep='first', ignore_index=True)
        df = df.sort_values(by='Similarity_to_query', ascending=False, ignore_index=True)
        df = pd.concat([pd.DataFrame([query_molecule_row]), df], ignore_index=True)
        df = df.reset_index(drop=True)
        df['ID'] = [f'HEAL_{i:06d}' for i in df.index]
        if as_dict:
            return df.to_dict(orient='records')
        return df

    def save_results(self, path: str='enumerated_molecules.csv'):
        '''
            Save the results to a file in CSV format.
            
            Columns: Product, Similarity_to_query, BB1, BB2, 
                Reaction_name, URL1, URL2

            Args:
                path (str): path to the file.
        '''
        df = self.get_results()
        df.to_csv(path, index=False)

    def _process_building_blocks(self, batch_size: int=10000):
        '''
            Process the building blocks to filter out the ones that do not
            satisfy the rules. This will create a list of building blocks
            for each composition of the molecule. The final list will be
            in the following format:
            ```
                [  
                    [[bb1, bb3, ...], [bb5, bb6, ...]],     # composition 1  
                    [[bb1, bb2, ...], [bb4, bb9, ...]],     # composition 2  
                    ...
                ]
            ```
        '''
        compositions = list(chain(*self._compositions[:self.n_compositions]))
        comp_sizes = np.array([comp.GetNumHeavyAtoms() for comp in compositions])[:, None]
        composition_fps = self._get_fingerprints(compositions)
        
        n_mols = len(self.mols)
        sims = np.zeros((len(composition_fps), n_mols))
        for i in tqdm(range(0, n_mols, batch_size), desc='Processing building blocks', total=n_mols//batch_size, disable=not self.verbose):
            batch_mols = self.mols[i:i+batch_size]
            batch_sizes = np.array([mol.GetNumHeavyAtoms() for mol in batch_mols])
            batch_sim_weights = 1 - (np.clip(batch_sizes - comp_sizes, 0, None) / batch_sizes)
            batch_stock_fps = self._get_fingerprints(batch_mols)
            batch_sims = utils.get_batch_tversky_sims_rdkit(composition_fps, batch_stock_fps)
            sims[:, i:i+batch_size] = batch_sims * batch_sim_weights

        if self.max_bbs_per_comp > 0:
            quantiles = np.quantile(sims, 1-(self.max_bbs_per_comp / n_mols), axis=1)
            quantiles = np.clip(quantiles, self.sim_threshold, None)[:, None]
            sims = sims >= quantiles
        else:
            sims = sims >= self.sim_threshold
        
        self._filtered_bb = []
        for i in range(0, len(sims), 2):
            mask_row1 = sims[i]
            mask_row2 = sims[i+1]
            if (not any(mask_row1)) or (not any(mask_row2)):
                self._filtered_bb.append([[], []])
            else:
                self._filtered_bb.append([list(compress(self.mols, mask_row1)), list(compress(self.mols, mask_row2))])
    
    def _prepare_molecule(self):
        '''
            Prepare the molecule by finding possible substructure compositions with
            respect to reaction template data or custom composition sites. If fragments
            are provided, they will be used as the compositions.
        '''
        if not hasattr(self, 'molecule'):
            raise ValueError('Molecule not set! Please set the molecule first.')
        
        frags = Chem.GetMolFrags(self.molecule, asMols=True, sanitizeFrags=True)
        if len(frags) == 2:
            self._compositions = [frags]
            print(f'The molecule has 2 fragments! Using them as the compositions. If this is not intended,'
                  f'remove mixture fragments from the molecule.', flush=True)
        elif len(frags) == 1:
            self._compositions = [] 
            if self.custom_comp_sites:
                for site in self.custom_comp_sites:
                    product = self._split_molecule(site)
                    for p in product:
                        try:
                            flag = Chem.SanitizeMol(p, catchErrors=True)
                            assert flag == Chem.rdmolops.SanitizeFlags.SANITIZE_NONE
                        except AssertionError:
                            print('Sanitization failed!')
                            continue
                    self._compositions.append(product)
            else:
                for reaction in self._reactions:
                    products = reaction.run_retro(self.molecule)
                    if products:
                        for product in products:
                            try:
                                for p in product:
                                    flag = Chem.SanitizeMol(p, catchErrors=True)
                                    assert flag == Chem.rdmolops.SanitizeFlags.SANITIZE_NONE
                            except AssertionError:
                                print('Sanitization failed for a composition!')
                                continue
                            self._compositions.append(product)
                
                self._remove_duplicate_compositions()
        else:
            raise ValueError('The molecule has more than 2 fragments! Please provide a valid molecule.')

        if self.verbose:
            self.print_compositions()

    def _split_molecule(self, split_site: tuple[int, int]):
        '''
            Splits the molecule at the given bond.

            Args:
                split_site: tuple of two integers, bond to split.

            Returns:
                tuple of rdkit mol objects, split molecules.
        '''
        # split the molecule
        with Chem.RWMol(self.molecule) as m:
            m.RemoveBond(split_site[0], split_site[1])
        m = m.GetMol()
        frags = Chem.GetMolFrags(m, asMols=True, sanitizeFrags=True)
        return frags
        
    def _remove_duplicate_compositions(self):
        '''
            Remove duplicate compositions.
        '''
        composition_idx_to_remove = []
        for i, comp_i in enumerate(self._compositions):
            for j, comp_j in enumerate(self._compositions):
                if i < j:
                    all_smiles = [Chem.MolToSmiles(substruct) for substruct in comp_i]
                    all_smiles.extend([Chem.MolToSmiles(substruct) for substruct in comp_j])
                    if len(set(all_smiles)) <= 2:
                        composition_idx_to_remove.append(i)
                        break
        
        for idx in sorted(composition_idx_to_remove)[::-1]:
            self._compositions.pop(idx)

    def print_compositions(self):
        if self._compositions:
            for i, composition in enumerate(self._compositions):
                smiles = [Chem.MolToSmiles(substruct) for substruct in composition]
                print(f'Composition {i}: {smiles}')
        else:
            print('No compositions found!')


