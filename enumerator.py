'''
    This script contains the Enumerator class which is used to enumerate
    a given molecule with Enamine building blocks or any other building block
    source. 
'''

# libraries
import abc
import utils
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem.FastSDMolSupplier import FastSDMolSupplier

from itertools import product as iter_product
from tqdm import tqdm


# TODO: Speed up building block processing by using kd-trees or parallel processing

class _BaseEnumerator(abc.ABC):
    '''
        Base Enumerator.
    '''
    def __init__(self, molecule: str | Chem.rdchem.Mol, bb_supplier: str, load_reactions: bool=True):
        '''
            Initialize the BaseEnumerator object.

            Args:
                molecule (str): SMILES string or rdkit mol object.
                bb_supplier (str): "US_stocks", "EU_stocks", "Global_stocks", "No_Rush_BBs",
                                    or the path to a file containing building blocks.
        '''
        # molecule
        if isinstance(molecule, str):
            self.molecule = Chem.MolFromSmiles(molecule)
        else:
            self.molecule = molecule

        # building blocks
        if bb_supplier == 'US_stocks':
            self.bb_supplier = FastSDMolSupplier('buildingblocks/Enamine_Rush-Delivery_Building_Blocks-US_195312cmpd_20240610.sdf')
        elif bb_supplier == 'EU_stocks':
            self.bb_supplier = FastSDMolSupplier('buildingblocks/Enamine_Rush-Delivery_Building_Blocks-EU_153230cmpd_20240806.sdf')
        elif bb_supplier == 'Global_stocks':
            self.bb_supplier = FastSDMolSupplier('buildingblocks/Enamine_Building_Blocks_Stock_290951cmpd_20240806.sdf')
        elif bb_supplier == 'No_Rush_BBs':
            self.bb_supplier = FastSDMolSupplier('buildingblocks/Enamine_Building_Blocks_1333124cmpd_20240806.sdf')
        elif bb_supplier == 'test':
            self.bb_supplier = FastSDMolSupplier('buildingblocks/test_100_bb.sdf')
        else:
            self.bb_supplier = FastSDMolSupplier(bb_supplier)

        # reaction data
        if load_reactions:
            self._reactions = utils.load_reactions_from_json('reactions/reactions.json')
            self._reactions = [reaction for reaction in self._reactions if reaction.is_valid()]

    @abc.abstractmethod
    def enumerate(self):
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
    

class CustomEnumerator(_BaseEnumerator):
    '''
        Custom Enumerator class is used to enumerate a given molecule with building blocks.
    '''
    def __init__(
            self, 
            molecule: str | Chem.rdchem.Mol,
            building_blocks: str='US_stocks',
            reaction_sites: list[int]=[],
            reaction_tags: list[str]=['amide coupling', 'amide', 'C-N bond formation', 'C-N',
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
            struct_rules: list[str]=[]
    ):
        '''
            Initialize the Enumerator object.

            Args:
                molecule (str): SMILES string or rdkit mol object. This molecule will be
                                enumerated with building blocks at the given reaction site.
                reaction_sites (list): list of atom indices to consider for the enumeration.
                                      If no reaction site is provided, any possible reaction
                                      site will be considered.
                building_blocks (str): "US_stocks", "EU_stocks", "Global_stocks", "No_Rush_BBs",
                                        or the path to a file containing building blocks.
                reaction_tags (list): list of reaction tags to consider for the enumeration.
                rules (dict): dictionary containing the rules if the method is rules.
                                - MW: tuple (min, max) -- molecular weight
                                - HBD: tuple (min, max) -- hydrogen bond donors
                                - HBA: tuple (min, max) -- hydrogen bond acceptors
                                - TPSA: tuple (min, max) -- topological polar surface area
                                - RotB: tuple (min, max) -- rotatable bonds
                                - Rings: tuple (min, max) -- number of rings
                                - ArRings: tuple (min, max) -- number of aromatic rings
                                - Chiral: tuple (min, max) -- number of chiral centers
                struct_rules (list): list of structure-based rules. List of SMILES to include 
                                     in the building blocks.
        '''
        super().__init__(molecule, building_blocks, True)
        self.reaction_sites = reaction_sites
        self.rules = rules
        self.struct_rules = struct_rules
        self.reaction_tags = reaction_tags

    def enumerate(self):
        '''
            Enumerate the molecule with building blocks.
        '''
        self._prepare_molecule()
        self._process_building_blocks()
        self._reactions = [reaction for reaction in self._reactions if any(tag in reaction.tags for tag in self.reaction_tags)]
        
        # enumerate the molecule with building blocks at the reaction site(s)
        query_fp = AllChem.GetMorganFingerprintAsBitVect(self.molecule, 3, nBits=2048)
        self.enumerated_molecules = []
        for bb in tqdm(self._filtered_bb, desc='Enumerating building blocks', total=len(self._filtered_bb)):
            for reaction in self._reactions:
                products = reaction.run_syn(self._prepared_mol, bb)
                if products:
                    for product in products:
                        for p in product:
                            try:
                                p = Chem.MolFromSmiles(Chem.MolToSmiles(p))
                                product_fp = AllChem.GetMorganFingerprintAsBitVect(p, 3, nBits=2048)
                                tani_sim = utils.get_tani_sim_fp(query_fp, product_fp)
                            except RuntimeError:
                                print('Fingerprint calculation of a product failed!', flush=True)
                                tani_sim = 0.0
                            self.enumerated_molecules.append((Chem.MolToSmiles(p),
                                                                str(round(tani_sim)),
                                                                Chem.MolToSmiles(bb),
                                                                reaction.name,
                                                                str(hash(reaction))))

    def save_results(self, path: str='enumerated_molecules.csv'):
        '''
            Save the results to a file in CSV format.
            
            Columns: Product, BB1, Reaction_name, Reaction_ID

            Args:
                path (str): path to the file.
        '''
        if not self.enumerated_molecules:
            print('No enumerated molecules found! Writing an empty file...')
        
        column_names = ['Product', 'Similarity_to_query', 'BB', 'Reaction_name', 'Reaction_ID']
        df = pd.DataFrame(self.enumerated_molecules, columns=column_names)
        df = df.drop_duplicates(subset=['Product'])
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
        for bb in tqdm(self.bb_supplier, desc='Processing building blocks', total=len(self.bb_supplier)):
            if bb is None:
                continue
            try:
                flag = Chem.SanitizeMol(bb)
                assert flag == Chem.rdmolops.SanitizeFlags.SANITIZE_NONE
            except AssertionError:
                print('Sanitization failed!', flush=True)
                continue
            if self._check_rules(bb) and self._check_struct_rules(bb):
                self._filtered_bb.append(bb)
    
    def _prepare_molecule(self):
        '''
            Prepare the molecule by adding protection to the atoms that are not 
            part of the reaction site.
        '''
        self._prepared_mol = Chem.MolFromSmiles(Chem.MolToSmiles(self.molecule))
        if self.reaction_sites:
            # Only given reaction sites should be considered reactive
            dont_protect = set()
            for atom in self._prepared_mol.GetAtoms():
                if atom.GetIdx() in self.reaction_sites:
                    dont_protect.add(atom.GetIdx())
                    for neighbor in atom.GetNeighbors():
                        dont_protect.add(neighbor.GetIdx())
            for atom in self._prepared_mol.GetAtoms():
                if atom.GetIdx() not in dont_protect:
                    atom.SetProp('_protected', '1')
        else:
            pass

    def _check_struct_rules(self, building_block: Chem.rdchem.Mol):
        '''
            Check if the building block satisfies the structure-based rules.

            Args:
                building_block: rdkit mol object.
        '''
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
        

class AutomatedEnumerator(_BaseEnumerator):
    '''
        Automated Enumerator class is used to enumerate a given molecule with building blocks.
    '''
    def __init__(
            self, 
            molecule, 
            building_blocks: str='US_stocks', 
            reaction_tags: list[str]=['amide coupling', 'amide', 'C-N bond formation', 'C-N',
                                      'alkylation', 'N-arylation', 'azole', 'amination'],
            custom_comp_sites: list[tuple]=[],
            n_compositions: int=5,
            sim_threshold: float=0.5
    ):
        '''
            Initialize the AutomatedEnumerator object.

            Args:
                molecule (str): SMILES string or rdkit mol object.
                building_blocks (str): "US_stocks", "EU_stocks", "Global_stocks", "No_Rush_BBs",
                                        or the path to a file containing building blocks.
                reaction_tags (list): list of reaction tags to consider for the enumeration.
                custom_comp_sites (list(tuple)): list of tuples containing the atom indices for
                                                 splitting the molecule. Each tuple represents a
                                                 a composition site.
                n_compositions (int): number of compositions of the molecule to enumerate.
                sim_threshold (float): similarity threshold.
        '''
        super().__init__(molecule, building_blocks, True)
        self.reaction_tags = reaction_tags
        self.custom_comp_sites = custom_comp_sites
        self.n_compositions = n_compositions
        self.sim_threshold = sim_threshold

        self._compositions = [] # list of rxn based compositions of the molecule

    def enumerate(self):
        '''
            Enumerate the molecule with building blocks.
        '''
        self._prepare_molecule()
        self._process_building_blocks()
        reactions = [reaction for reaction in self._reactions if any(tag in reaction.tags for tag in self.reaction_tags)]

        # enumerate the molecule with building blocks
        query_fp = AllChem.GetMorganFingerprintAsBitVect(self.molecule, 3, nBits=2048)
        self.enumerated_molecules = []
        counter = 0
        for composition in tqdm(self._filtered_bb, desc='Enumerating building blocks', total=len(self._filtered_bb)):
            if counter == self.n_compositions:
                break
            counter += 1
            for b1, b2 in iter_product(*composition):
                for reaction in reactions:
                    products = reaction.run_syn(b1, b2)
                    if products:
                        for product in products:
                            for p in product:
                                try:
                                    p = Chem.MolFromSmiles(Chem.MolToSmiles(p))
                                    product_fp = AllChem.GetMorganFingerprintAsBitVect(p, 3, nBits=2048)
                                    tani_sim = utils.get_tani_sim_fp(query_fp, product_fp)
                                except RuntimeError:
                                    print('Fingerprint calculation of a product failed!', flush=True)
                                    tani_sim = 0.0
                                self.enumerated_molecules.append((Chem.MolToSmiles(p),
                                                                  str(round(tani_sim, 2)),
                                                                  Chem.MolToSmiles(b1),
                                                                  Chem.MolToSmiles(b2),
                                                                  reaction.name,
                                                                  str(hash(reaction))))

    def save_results(self, path: str='enumerated_molecules.csv'):
        '''
            Save the results to a file in CSV format.
            
            Columns: Product, BB1, BB2, Reaction_name, Reaction_ID

            Args:
                path (str): path to the file.
        '''
        if not self.enumerated_molecules:
            print('No enumerated molecules found! Writing an empty file...')

        column_names = ['Product', 'Similarity_to_query', 'BB1', 'BB2', 'Reaction_name', 'Reaction_ID']
        df = pd.DataFrame(self.enumerated_molecules, columns=column_names)
        df = df.drop_duplicates(subset=['Product'])
        df.to_csv(path, index=False)

    def _process_building_blocks(self):
        '''
            Process the building blocks to filter out the ones that do not
            satisfy the rules.
        '''
        self._filtered_bb = [] # [[[bb1, bb2, ...], [bb1, bb2, ...]], ...]
        for composition in tqdm(self._compositions, desc='Processing building blocks', total=len(self._compositions)):
            composition_bbs = []
            for substruct in composition:
                bb_filtered = []
                for bb in self.bb_supplier:
                    if bb is None:
                        continue
                    try:
                        flag = Chem.SanitizeMol(bb)
                        assert flag == Chem.rdmolops.SanitizeFlags.SANITIZE_NONE
                    except AssertionError:
                        print('Sanitization failed!')
                        continue
                    if utils.get_tani_sim(bb, substruct) >= self.sim_threshold:
                        bb_filtered.append(bb)
                composition_bbs.append(bb_filtered)
            if all(composition_bbs):
                self._filtered_bb.append(composition_bbs)
            else:
                comp = [Chem.MolToSmiles(substruct) for substruct in composition]
                print(f'No building blocks found for the composition: {comp}')
                
    def _prepare_molecule(self):
        '''
            Prepare the molecule by finding possible substructure compositions with
            respect to reaction template data or custom composition sites.
        '''
        if self.custom_comp_sites:
            for site in self.custom_comp_sites:
                product = utils.split_molecule(self.molecule, site)
                for p in product:
                    try:
                        flag = Chem.SanitizeMol(p)
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
                        for p in product:
                            try:
                                flag = Chem.SanitizeMol(p)
                                assert flag == Chem.rdmolops.SanitizeFlags.SANITIZE_NONE
                            except AssertionError:
                                print('Sanitization failed!')
                                continue
                        self._compositions.append(product)
            
            self._remove_duplicate_compositions()

        self.print_compositions()
        
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


