'''
    This script contains the Enumerator class which is used to enumerate
    a given molecule with Enamine building blocks or any other building block
    source. 
'''

# libraries
import utils

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors, DataStructs
from rdkit.Chem.FastSDMolSupplier import FastSDMolSupplier

from itertools import product as iter_product
from tqdm import tqdm


## Reactions
# generic reaction template for molecules with dummy atoms, R-* . R'-* >> R-R'
generic_rxn = AllChem.ReactionFromSmarts('[*:1]-[*H5:2].[*:3]-[*H5:4]>>[*:1]-[*:3]')
# amide coupling reaction
amide_coupling_rxn = AllChem.ReactionFromSmarts('[C:1](=[O:2])O.[NH2:3]>>[C:1](=[O:2])[NH:3]')

# SMART site dictionary -- key: functional group pattern, value: functional group with dummy atom [*H5]
rxn_sites = {
    'carboxylic_acid': '[C:1](=[O:2])[OH]>>[C:1](=[O:2])[*]', # carboxylic acid
    'primary_amine': '[C:1][NH2:2]>>[C:1][NH:2][*]', # amine
}


class CustomEnumerator:
    '''
        Enumerator class is used to enumerate a given molecule with building blocks.
    '''

    def __init__(self, molecule, building_blocks: str='US_stocks', method: str='similarity'):
        '''
            Initialize the Enumerator object.

            Args:
                molecule (str): SMILES string or rdkit mol object.
                building_blocks (str): "US_stocks", "EU_stocks", "Global_stocks",
                                        or the path to a file containing building blocks.
                method (str): rules, similarity, or similarity_and_rules.
        '''
        # molecule
        if isinstance(molecule, str):
            self.molecule = Chem.MolFromSmiles(molecule)
        else:
            self.molecule = molecule

        # building blocks
        if building_blocks == 'US_stocks':
            self.bb_supplier = FastSDMolSupplier('buildingblocks/Enamine_Rush-Delivery_Building_Blocks-US_195312cmpd_20240610.sdf')
        elif building_blocks == 'EU_stocks':
            return NotImplementedError
        elif building_blocks == 'Global_stocks':
            return NotImplementedError
        elif building_blocks == 'test': # remove this later
            self.bb_supplier = FastSDMolSupplier('buildingblocks/test_100_bb.sdf')
        else:
            self.bb_supplier = FastSDMolSupplier(building_blocks)

        # method
        self.method = method

        # substructs to replace
        self.substruct_to_enumerate = None

        # rules
        self.rules = {
            'MW': (0, 500), # molecular weight
            'HBD': (0, 5), # hydrogen bond donors
            'HBA': (0, 10), # hydrogen bond acceptors
            'TPSA': (0, 200), # topological polar surface area
            'RotB': (0, 10), # rotatable bonds
            'Rings': (0, 10), # number of rings
            'ArRings': (0, 5), # number of aromatic rings
            'Chiral': (0, 5), # number of chiral centers
        }

    def enumerate(self, include: list[str]=[], exclude: list[str]=[], sim_threshold: float=0.5):
        '''
            Enumerate the molecule with building blocks.
        '''
        if self.substruct_to_enumerate is None:
            raise ValueError('No substructure is added.')
        
        # prepare the molecule
        self._prepare_molecule()

        # process building blocks
        self.process_building_blocks(include, exclude, sim_threshold)
        
        # enumerate the molecule with building blocks for only one substruct
        self.enumerated_molecules = []
        for bb in tqdm(self.filtered_bb, desc='Enumerating building blocks', total=len(self.filtered_bb)):
            # enumerate the molecule with the building block
            products = generic_rxn.RunReactants((self._prepared_mol, bb))
            for product in products:
                self.enumerated_molecules.append(product[0])

    def save_enumerated_molecules(self, path: str='enumerated_molecules.smi', format: str='smi'):
        '''
            Save the enumerated molecules to a file.

            Args:
                path (str): path to the file.
        '''
        # save the enumerated molecules
        if format == 'smi':
            with open(path, 'w') as f:
                for mol in self.enumerated_molecules:
                    f.write(Chem.MolToSmiles(mol) + '\n')
        elif format == 'sdf':
            w = Chem.SDWriter(path)
            for mol in self.enumerated_molecules:
                w.write(mol)
            w.close()
        else:
            raise ValueError('Invalid format.')

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

    def add_substruct(self, substruct: Chem.rdchem.Mol | str | list[int]):
        '''
            Add a substruct to the molecule.

            Args:
                substruct: rdkit mol, SMILES or list of atom indices.
        '''
        # convert substruct to rdkit mol object
        if isinstance(substruct, Chem.rdchem.Mol):
            self.substruct_to_enumerate = substruct
        elif isinstance(substruct, str):
            self.substruct_to_enumerate = Chem.MolFromSmiles(substruct)
        elif isinstance(substruct, list):
            self.substruct_to_enumerate = Chem.PathToSubmol(self.molecule, substruct)
        else:
            raise ValueError('Invalid substruct type.')

    def process_building_blocks(self, include: list[str]=[], exclude: list[str]=[],
                                sim_threshold: float=0.7):
        '''
            Process the building blocks to filter out the ones that do not
            satisfy the rules.

            Args:
                include: list of functional groups to include.
                exclude: list of functional groups to exclude.
                sim_threshold: similarity threshold.
        '''
        # filter building blocks
        self.filtered_bb = []
        for bb in tqdm(self.bb_supplier, desc='Processing building blocks', total=len(self.bb_supplier)):
            # check rules
            if self.method == 'similarity':
                if self._check_similarity(bb, self.molecule, threshold=sim_threshold):
                    # prepare the reaction site
                    sites = self._prepare_site(bb, include, exclude)
                    for site in sites:
                        site = self._dummy2dummy(site)
                        self.filtered_bb.append(site)
            elif self.method == 'similarity_with_rules':
                if self._check_similarity(bb, self.molecule, threshold=sim_threshold) and self._check_rules(bb):
                    # prepare the reaction site
                    sites = self._prepare_site(bb, include, exclude)
                    for site in sites:
                        site = self._dummy2dummy(site)
                        self.filtered_bb.append(site)
            elif self.method == 'rules':
                if self._check_rules(bb):
                    # prepare the reaction site
                    sites = self._prepare_site(bb, include, exclude)
                    for site in sites:
                        site = self._dummy2dummy(site)
                        self.filtered_bb.append(site)
            else:
                raise ValueError('Invalid method.')

    @staticmethod  
    def _dummy2dummy(mol: Chem.rdchem.Mol):
        '''
            Helper function to replace [*] with [*H5] in the molecule.

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
    
    def _prepare_molecule(self):
        '''
            Prepare molecule by replacing the substruct with a dummy atom.
        '''
        # prepare the molecule
        if self.substruct_to_enumerate is not None:
            match = list(self.molecule.GetSubstructMatch(self.substruct_to_enumerate))
            # find the atom that is connected to the core
            atom1 = None
            for i in match:
                for neighbor in self.molecule.GetAtomWithIdx(i).GetNeighbors():
                    if neighbor.GetIdx() not in match:
                        atom1 = i

            # remove the substructure and add a dummy atom
            with Chem.RWMol(self.molecule) as new_mol:
                for i in match:
                    if i != atom1:
                        new_mol.RemoveAtom(i)
                new_mol.ReplaceAtom(atom1, Chem.Atom(0))
            prepared_mol = new_mol.GetMol()
            try:
                Chem.SanitizeMol(prepared_mol)
            except:
                return ValueError('Substructure could not be removed.')
            self._prepared_mol = self._dummy2dummy(prepared_mol)
            self._core_mol = Chem.MolFromSmiles(Chem.MolToSmiles(self._prepared_mol).replace('[*H5]', ''))
        else:
            return ValueError('No substructure is added.')

    @staticmethod
    def _prepare_site(mol: Chem.rdchem.Mol | str, include: list[str]=[], exclude: list[str]=[]):
        '''
            Prepare the reaction site for the building block.

            NOTE: If both include and exclude are empty, all functional groups
                will be included. If include is not empty, exclude will be ignored.

            Args:
                mol: rdkit mol or SMILES string.
                include: list of functional groups to include.
                exclude: list of functional groups to exclude.

            Returns:
                mols: list of rdkit mol objects with a dummy atom 
                    at the reaction site.
        '''
        # convert mol to rdkit mol object
        if isinstance(mol, str):
            mol = Chem.MolFromSmiles(mol)

        # prepare the reaction site
        mols = []
        for key, value in rxn_sites.items():
            if include:
                if key in include:
                    rxn = AllChem.ReactionFromSmarts(value)
                    products = rxn.RunReactants((mol,))
                    for product in products:
                        mols.append(product[0])
            else:
                if key not in exclude:
                    rxn = AllChem.ReactionFromSmarts(value)
                    products = rxn.RunReactants((mol,))
                    for product in products:
                        mols.append(product[0])

        return mols
    
    @staticmethod
    def _check_similarity(mol1: Chem.rdchem.Mol, mol2: Chem.rdchem.Mol, threshold: float=0.7):
        '''
            Check the similarity between two molecules.

            Args:
                mol1: rdkit mol object.
                mol2: rdkit mol object.
                threshold: similarity threshold.

            Returns:
                bool: True if the similarity is greater than the threshold.
        '''
        # calculate similarity
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 3, nBits=2048)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 3, nBits=2048)
        similarity = DataStructs.TanimotoSimilarity(fp1, fp2)

        return similarity >= threshold

    def _check_rules(self, building_block: Chem.rdchem.Mol | str):
        '''
            Check if the building block satisfies the rules.

            Args:
                building_block: rdkit mol or SMILES string.
        '''
        # convert building block to rdkit mol object
        if isinstance(building_block, str):
            building_block = Chem.MolFromSmiles(building_block)

        # check rules
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
            

class AutomatedEnumerator:
    '''
        Automated Enumerator class is used to enumerate a given molecule with building blocks.
    '''
    def __init__(self, molecule, building_blocks: str='US_stocks', sim_threshold: float=0.5):
        '''
            Initialize the AutomatedEnumerator object.

            Args:
                molecule (str): SMILES string or rdkit mol object.
                building_blocks (str): "US_stocks", "EU_stocks", "Global_stocks",
                                        or the path to a file containing building blocks.
                sim_threshold (float): similarity threshold.
        '''
        # molecule
        if isinstance(molecule, str):
            self.molecule = Chem.MolFromSmiles(molecule)
        else:
            self.molecule = molecule

        # reaction data
        self.reactions = utils.load_reactions_from_json('reactions/reactions.json')
        self.reactions = [reaction for reaction in self.reactions if reaction.is_valid()]

        # building blocks
        if building_blocks == 'US_stocks':
            self.bb_supplier = FastSDMolSupplier('buildingblocks/Enamine_Rush-Delivery_Building_Blocks-US_195312cmpd_20240610.sdf')
        elif building_blocks == 'EU_stocks':
            return NotImplementedError
        elif building_blocks == 'Global_stocks':
            return NotImplementedError
        elif building_blocks == 'test':
            self.bb_supplier = FastSDMolSupplier('buildingblocks/test_100_bb.sdf')
        else:
            self.bb_supplier = FastSDMolSupplier(building_blocks)

        # similarity threshold
        self.sim_threshold = sim_threshold

        # possible substructure compositions
        self.compositions = []

    def enumerate(self, n_compositions: int=5, 
                  reaction_tags: list[str]=['amide coupling', 'amide', 'C-N bond formation', 'C-N',
                                            'alkylation', 'N-arylation', 'azole', 'amination']):
        '''
            Enumerate the molecule with building blocks.

            Args:
                n_compositions (int): number of compositions to consider.
        '''
        # find possible compositions
        self._find_compositions()

        # process building blocks
        self._process_building_blocks()

        # filter reactions by tags
        reactions = [reaction for reaction in self.reactions if any(tag in reaction.tags for tag in reaction_tags)]

        # enumerate the molecule with building blocks
        self.enumerated_molecules = []
        counter = 0
        for composition in tqdm(self.filtered_bb, desc='Enumerating building blocks', total=len(self.filtered_bb)):
            if counter == n_compositions:
                break
            counter += 1
            for b1, b2 in iter_product(*composition):
                for reaction in reactions:
                    products = reaction.run(b1, b2)
                    if products:
                        for product in products:
                            for p in product:
                                self.enumerated_molecules.append((Chem.MolToSmiles(p),
                                                                  Chem.MolToSmiles(b1),
                                                                  Chem.MolToSmiles(b2),
                                                                  reaction.name))

    def save_enumerated_molecules(self, path: str='enumerated_molecules.csv'):
        '''
            Save the enumerated molecules to a file.

            Args:
                path (str): path to the file.
        '''
        # save the enumerated molecules
        with open(path, 'w') as f:
            f.write('Product,BB1,BB2,Reaction\n')
            for mol in self.enumerated_molecules:
                f.write(','.join(mol) + '\n')                    

    def _process_building_blocks(self):
        '''
            Process the building blocks to filter out the ones that do not
            satisfy the rules.
        '''
        # filter building blocks
        self.filtered_bb = []
        for composition in tqdm(self.compositions, desc='Processing building blocks', total=len(self.compositions)):
            bb_list = []
            for substruct in composition:
                bb_list.append([bb for bb in self.bb_supplier if self._check_similarity(bb, substruct, threshold=self.sim_threshold)])
            self.filtered_bb.append(bb_list)

    @staticmethod
    def _check_similarity(mol1: Chem.rdchem.Mol, mol2: Chem.rdchem.Mol, threshold: float=0.7):
        '''
            Check the similarity between two molecules.

            Args:
                mol1: rdkit mol object.
                mol2: rdkit mol object.
                threshold: similarity threshold.

            Returns:
                bool: True if the similarity is greater than the threshold.
        '''
        # calculate similarity
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 3, nBits=2048)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 3, nBits=2048)
        similarity = DataStructs.TanimotoSimilarity(fp1, fp2)

        return similarity >= threshold

    def _find_compositions(self):
        '''
            Find the best possible substructure compositions with
            respect to reaction template data.
        '''
        # iterate through the reactions and find possible splits
        for reaction in self.reactions:
            # run retro reaction
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
                    self.compositions.append(product)

