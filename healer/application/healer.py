import abc
import logging
from pathlib import Path
from typing import List, Union, Dict, Tuple, Any, Optional
from itertools import chain

from tqdm import tqdm
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, rdFingerprintGenerator
from rdkit.Chem.FastSDMolSupplier import FastSDMolSupplier

from healer.application.tree_builder import CompositionPath, RetrosynthesisTree
from healer.domain.composition import CompositionWithBBs
from healer.application.optimizers import BaseStagewiseOptimizer, BaseSequenceOptimizer
from healer.domain.building_block import BuildingBlock
from healer.domain.reaction_template import ReactionTemplate21
from healer.domain.enumeration_record import EnumerationRecord
import healer.utils.utils as utils

# Configure logging
logging.basicConfig(
    format="%(asctime)s [%(levelname)s]: %(message)s",
    datefmt="%H:%M:%S",
    force=True,
    level=logging.WARNING
)
logger = logging.getLogger(__name__)

# Constants: mapping building-block sources to file paths
BB_PATHS: Dict[str, str] = {
    "US_stock": "healer/data/buildingblocks/Enamine_Rush-Delivery_Building_Blocks-US_*_processed.sdf",
    "EU_stock": "healer/data/buildingblocks/Enamine_Rush-Delivery_Building_Blocks-EU_*_processed.sdf",
    "Global_stock": "healer/data/buildingblocks/Enamine_Building_Blocks_Stock_*_processed.sdf",
    "test": "healer/data/buildingblocks/test_100_bb_processed.sdf",
}


class _BaseHEALER(abc.ABC):
    '''
        Base HEALER.
    '''
    def __init__(
        self, 
        bb_supplier: str, 
        reaction_tags: Union[List[str], str],
        max_evals_per_comp: Optional[int] = None,
        verbose: int=1,
    ) -> None:
        '''
            Initialize BaseHEALER.

            Args:
                bb_supplier: one of "US_stock", "EU_stock" or "Global_stock"; or path to an SDF file.
                reaction_tags: list of tags or 'all'.
                max_evals_per_comp: maximum number of evaluations for each composition.
                verbose: verbosity level.
                    - 0: only errors
                    - 1: warnings
                    - 2: info
        '''
        # attributes
        self.verbose: int = verbose

        self.bb_supplier: FastSDMolSupplier = None
        self._supplier_path: str = ''
        self.bb_mols: List[BuildingBlock] = []
        
        self.query_mol: Chem.Mol = None
        self._compositions: List[Union[CompositionPath, CompositionWithBBs]] = []
        self.max_evals_per_comp: Optional[int] = max_evals_per_comp

        self.reactions: List[ReactionTemplate21] = []
        self.reaction_tags: List[str] = []
        self._reactions: List[ReactionTemplate21] = []  # all reactions loaded from JSON
        
        self.enumerated_molecules: List[EnumerationRecord] = []  
        
        # Set verbosity level
        logger.setLevel(logging.ERROR if verbose == 0 else logging.WARNING if verbose == 1 else logging.INFO)
        logger.info("Initializing HEALER with building block supplier: %s", bb_supplier)

        # resolve supplier path
        self._supplier_path = self._get_supplier_path(bb_supplier)
        self.bb_supplier = FastSDMolSupplier(self._supplier_path, sanitize=True)
        logger.info('Building block source has %d building blocks', len(self.bb_supplier))  

        # load and filter reactions along with compatible building blocks
        self.set_reactions(reaction_tags)

        # fingerprint generator
        self._fp_generator = rdFingerprintGenerator.GetMorganGenerator(
            radius=3, fpSize=2048, includeChirality=True
        )

    def _get_supplier_path(self, bb_supplier: str) -> str:
        '''Get the path to the building block supplier.'''
        pattern = BB_PATHS.get(bb_supplier, bb_supplier)
        p = Path(pattern)

        if any(ch in pattern for ch in ("*", "?", "[")): 
            search_dir = p.parent if p.parent != Path() else Path.cwd()
            matches = list(search_dir.glob(p.name))
            if not matches:
                raise FileNotFoundError(f"No file matches {pattern!r}")
            matches.sort(key=lambda f: f.stat().st_mtime, reverse=True)
            chosen = matches[0]
        else:
            chosen = p

        return str(chosen)

    def set_reactions(self, reaction_tags: Union[List[str], str]) -> None:
        '''
            Set reactions to use for enumeration.

            Args:
                reaction_tags: a rection tag or a list of tags to filter by.
                    'all' to use all reactions.
        '''
        # reaction data
        all_rxns = utils.load_reactions_from_json('healer/data/reactions/reactions.json')
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
        
        logger.info("Using %d reactions", len(self.reactions))

        # filter building blocks based on reactions
        self._update_bb_mols()

    def _update_bb_mols(self) -> List[Chem.Mol]:
        for mol in tqdm(self.bb_supplier, desc='Loading building blocks', total=len(self.bb_supplier), disable=self.verbose < 1):
            bb = BuildingBlock(mol)
            if any(rxn.name in bb.get_parsed_prop('rxn_annotations') for rxn in self.reactions):
                self.bb_mols.append(bb)
        logger.info("Filtered building blocks to %d based on reactions", len(self.bb_mols))

    def enumerate(
        self, 
        optimizer: Optional[Union[BaseStagewiseOptimizer, BaseSequenceOptimizer]] = None,
        max_evals_per_comp: Optional[int] = None,
    ) -> None:
        '''
            Enumerate the molecule with building blocks based on the reactions. An optimizer
            can be provided to optimize an objective function during enumeration.

            Args:
                optimizer: an optimizer object to use for enumeration.
                max_evals_per_comp: maximum number of evaluations for each composition.

            Raises:
                TypeError: if the optimizer is not of a supported type.
        '''
        if max_evals_per_comp is not None:
            self.max_evals_per_comp = max_evals_per_comp
        if not isinstance(self.query_mol, Chem.Mol):
            raise ValueError("Query molecule must be set before enumeration. Use set_query_mol() method.")
        self._process_query_mol()
        self._process_building_blocks()

        self.enumerated_molecules = [
            EnumerationRecord(
                product=self.query_mol,
                bbs=[],
                reaction_names=[],
                props={'optimization_score': optimizer.target_fn(self.query_mol)} if optimizer else {}
            )
        ]

        if optimizer is None:
            logger.info("No optimizer provided, using default enumeration.")
            self.enumerated_molecules += self._enumerate_base(self.max_evals_per_comp)
        elif isinstance(optimizer, BaseStagewiseOptimizer):
            logger.info("Using stagewise optimization for enumeration.")
            self.enumerated_molecules += self._enumerate_stagewise(optimizer, self.max_evals_per_comp)
        elif isinstance(optimizer, BaseSequenceOptimizer):
            logger.info("Using sequence optimization for enumeration.")
            self.enumerated_molecules += self._enumerate_sequence(optimizer, self.max_evals_per_comp)
        else:
            raise TypeError(f"Unsupported optimizer type: {type(optimizer)}. ")
    
    @abc.abstractmethod
    def set_query_mol(self, molecule: Union[str, Chem.Mol]) -> None:
        '''
            Set the query molecule for enumeration and set the attributes that will be used
            to process the query molecule. The purpose of this separation is to allow enumeration 
            of multiple molecules without reinitializing the HEALER instance.
        '''
        ...
        # if isinstance(molecule, str):
        #     self.query_mol = Chem.MolFromSmiles(molecule)
        # else:
        #     self.query_mol = molecule
        # flag = Chem.SanitizeMol(self.query_mol, catchErrors=True)
        # assert flag == Chem.rdmolops.SanitizeFlags.SANITIZE_NONE, f"SanitizeMol failed: {flag}"

    @abc.abstractmethod
    def _process_query_mol(self) -> None:
        '''
            Should update self._compositions with a list of `CompositionPath` objects.
        '''
        ...

    @abc.abstractmethod
    def _process_building_blocks(self) -> None:
        '''
            Should convert the `CompositionPath` to `CompositionWithBBs` objects 
            inside `self._compositions` by pairing each fragment with a list of 
            compatible building blocks.
        '''
        ...
    
    def _enumerate_base(self, max_evals: Optional[int] = None) -> List[EnumerationRecord]:
        '''
            Exhaustive enumeration without optimization.
        '''
        results: List[EnumerationRecord] = []
        for comp_bb in tqdm(self._compositions, desc="Enumerating compositions", disable=self.verbose < 1):
            eval_count = 0
            bb_lists = comp_bb.fragment_bbs 
            stage_records = self._make_seed_records(bb_lists[0])
            for bb_pool in tqdm(bb_lists[1:], desc="Enumerating stages", disable=self.verbose < 2):
                cands = self._generate_candidates(stage_records, bb_pool)
                stage_records = self._apply_candidates(cands)
            for rec in stage_records:
                results.append(rec)
                eval_count += 1
                if max_evals and eval_count >= max_evals:
                    logger.info("Reached max evaluations limit: %d", max_evals)
                    return results
        logger.info("Exhaustive enumeration completed with %d results.", len(results))
        return results
    
    def _enumerate_stagewise(
        self, 
        optimizer: BaseStagewiseOptimizer, 
        max_evals: Optional[int] = None
    ) -> List[EnumerationRecord]:
        '''
            Stagewise enumeration with optimizer.filter() hook.
        '''
        results: List[EnumerationRecord] = []
        for comp_bb in self._compositions:
            eval_count = 0
            bb_lists = comp_bb.fragment_bbs
            stage_records = self._make_seed_records(bb_lists[0])
            for depth, bb_pool in enumerate(tqdm(bb_lists[1:], desc="Enumerating stages", disable=self.verbose < 2), start=1):
                cands = self._generate_candidates(stage_records, bb_pool)
                cands = optimizer.filter(cands, depth)
                stage_records = self._apply_candidates(cands)
                if max_evals and eval_count >= max_evals:
                    logger.info("Reached max evaluations limit: %d", max_evals)
                    break
            for rec in stage_records:
                rec.props.update({'optimization_score': optimizer.target_fn(rec.product)})
                results.append(rec)
                eval_count += 1
                if max_evals and eval_count >= max_evals:
                    logger.info("Reached max evaluations limit: %d", max_evals)
                    return results
        logger.info("Stagewise enumeration completed with %d results", len(results))
        return results
    
    def _enumerate_sequence(
        self,
        optimizer: BaseSequenceOptimizer,
        max_evals: Optional[int] = None
    ) -> List[EnumerationRecord]:
        '''
            Sequence-based enumeration using optimizer.ask() and optimizer.tell().
        '''     
        results: List[EnumerationRecord] = []
        for comp_bb in self._compositions:
            eval_count = 0
            optimizer.init_search(
                domain=comp_bb.fragment_bbs, 
                budget=max_evals or 0
            )
            while True:
                bb_tuples = optimizer.ask()
                if not bb_tuples:
                    logger.info("No more building block tuples to evaluate.")
                    break

                bb_pools = [list(bb_tuple) for bb_tuple in zip(*bb_tuples)]    # transpose tuples
                stage_records = self._make_seed_records(bb_pools[0])  
                for bb_pool in tqdm(bb_pools[1:], desc="Enumerating stages", disable=self.verbose < 2):
                    cands = self._generate_candidates_positionwise(stage_records, bb_pool)
                    stage_records = self._apply_candidates(cands)

                feedback: List[Tuple[Tuple[BuildingBlock, ...], float]] = []
                for rec in stage_records:
                    score = optimizer.target_fn(rec.product)
                    rec.props.update({'optimization_score': score})
                    feedback.append((tuple(rec.bbs), score))
                    results.append(rec)
                    eval_count += 1
                    if max_evals and eval_count >= max_evals:
                        logger.info("Reached max evaluations limit: %d", max_evals)
                        break
                optimizer.tell(feedback)
                if max_evals and eval_count >= max_evals:
                    logger.info("Reached max evaluations limit: %d", max_evals)
                    break
        logger.info("Sequence-based enumeration completed with %d results", len(results))
        return results

    def _make_seed_records(self, bb0_pool: List[BuildingBlock]) -> List[EnumerationRecord]:
        '''
            Create initial enumeration records from the seed building blocks.
        '''
        return [
            EnumerationRecord(
                product=bb.get_mol(), bbs=[bb], reaction_names=[], props={}
            ) for bb in bb0_pool
        ]
    
    def _generate_candidates_positionwise(
        self,
        stage_records: List[EnumerationRecord],
        bb_pool: List[BuildingBlock]
    ) -> List[Tuple[EnumerationRecord, BuildingBlock, ReactionTemplate21]]:
        '''
            One-to-one mapping of enumeration records to building blocks. Pair 
            each current record with the BB at the same index for each reaction.
        '''
        candidates: List[Tuple[EnumerationRecord, BuildingBlock, ReactionTemplate21]] = []
        for rxn in self.reactions:
            for rec, bb in zip(stage_records, bb_pool):
                candidates.append((rec, bb, rxn))
        return candidates

    def _generate_candidates(
        self, 
        stage_records: List[EnumerationRecord], 
        bb_pool: List[BuildingBlock]
    
    ) -> List[Tuple[EnumerationRecord, BuildingBlock, ReactionTemplate21]]:
        '''
            Build all (rec, bb, rxn) triples for the next coupling.
        '''
        candidates: List[Tuple[EnumerationRecord, BuildingBlock, ReactionTemplate21]] = []
        for rxn in self.reactions:
            for rec in stage_records:
                for bb in bb_pool:
                    candidates.append((rec, bb, rxn))
        return candidates
    
    def _apply_candidates(
        self,
        candidates: List[Tuple[EnumerationRecord, BuildingBlock, ReactionTemplate21]]
    ) -> List[EnumerationRecord]:
        '''
            Batch-apply all candidates via _apply_candidate.
        '''
        next_stage_records: List[EnumerationRecord] = []
        for rec, bb, rxn in candidates:
            next_stage_records += self._apply_candidate(rec, bb, rxn)
        return next_stage_records    
    
    def _apply_candidate(
        self,
        rec: EnumerationRecord,
        bb: BuildingBlock,
        rxn: ReactionTemplate21
    ) -> List[EnumerationRecord]:
        '''
            Apply the reaction to the record and building block, returning new records.
        '''
        results: List[EnumerationRecord] = []
        mol0 = rec.product
        for product in rxn.run_syn(mol0, bb):
            flags = Chem.SanitizeMol(product, catchErrors=True)
            if flags != Chem.rdmolops.SanitizeFlags.SANITIZE_NONE:
                continue
            new_props = dict(rec.props)
            results.append(
                EnumerationRecord(
                    product=product,
                    bbs=rec.bbs + [bb],
                    reaction_names=rec.reaction_names + [rxn.name],
                    props=new_props
                )
            )
        return results

    def get_results(
        self, 
        as_dict: bool=False,
        calc_similarity: bool=False,
        calc_stoplight: bool=False,
        calc_cns_mpo: bool=False,
    ) -> Union[pd.DataFrame, List[Dict[str, Any]]]:
        '''
            Get the results of enumeration as a pandas DataFrame or a list of dicts.

            Args:
                as_dict (bool): if True, return the results as a list of dicts.
                calc_similarity (bool): if True, calculate similarity to the query molecule.
                calc_stoplight (bool): if True, calculate stoplight scores/colors.
                calc_cns_mpo (bool): if True, calculate CNS MPO scores.
        '''
        max_bb = max(len(r.bbs) for r in self.enumerated_molecules)
        max_rxn = max(len(r.reaction_names) for r in self.enumerated_molecules)

        # Build each row as a dict
        rows = []
        for r in self.enumerated_molecules:
            row = {"Product": Chem.MolToSmiles(r.product),}
            for i in range(max_bb):
                row[f"BB{i+1}"] = r.bbs[i].get_smiles() if i < len(r.bbs) else ""
            for i in range(max_rxn):
                row[f"Reaction{i+1}_name"] = (
                    r.reaction_names[i] if i < len(r.reaction_names) else ""
                )
            for i in range(max_bb):
                row[f"URL{i+1}"] = (
                    r.bbs[i].get_parsed_prop('URL') if i < len(r.bbs) else ""
                )
            row = {**row, **r.props}  # Add any additional properties
            rows.append(row)

        df = pd.DataFrame(rows)
        cols_to_consider = [col for col in df.columns if not col.startswith(("URL", "Reaction"))]
        df = df.drop_duplicates(subset=cols_to_consider, keep="first", ignore_index=True)
        
        if calc_similarity:
            enum_fps = self._get_fingerprints(df['Product'].apply(Chem.MolFromSmiles).tolist())
            query_fp = enum_fps.pop(0)  
            tani_sims = utils.get_batch_tani_sims([query_fp], enum_fps)[0]
            tani_sims = np.concat((np.array([1.001]), tani_sims))
            df['Similarity_to_query'] = tani_sims
            df = df.sort_values(by='Similarity_to_query', ascending=False, ignore_index=True)
            df = df.round({'Similarity_to_query': 2})
        
        if calc_stoplight:
            # TODO: Include stoplight calculation
            raise NotImplementedError("Stoplight calculation is not supported yet.")
        
        if calc_cns_mpo:
            # TODO: Include CNS MPO calculation
            raise NotImplementedError("CNS MPO calculation is not supported yet.")

        df.insert(0, "ID", [f"HEAL_{i:06d}" for i in df.index])
        if as_dict:
            return df.to_dict(orient="records")
        return df

    def save_results(self, path: str='results.csv') -> None:
        '''
            Save results to a CSV file.
            
            Args:
                path: path to save the results.
        '''
        df = self.get_results()
        df.to_csv(path, index=False)

    def _get_fingerprints(self, mols: List[Chem.Mol]) -> List[Any]:
        return list(self._fp_generator.GetFingerprints(mols))
        
    def _get_fingerprint(self, mol: Chem.Mol) -> Any:
        return self._fp_generator.GetFingerprint(mol)


class SiteHEALER(_BaseHEALER):
    '''
        Site HEALER: enumerates molecules by applying reactions to a query molecule
        at specified reactive sites.
    '''
    def __init__(
            self, 
            bb_supplier: str='US_stock',
            reaction_tags: list[str] | str=['amide coupling', 'amide', 'C-N bond formation', 'C-N',
                                            'alkylation', 'N-arylation', 'azole', 'amination'],
            max_evals_per_comp: Optional[int] = None,
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
                max_evals_per_comp: maximum number of evaluations for each composition.
                rules: dictionary of rules for filtering molecules.
                struct_rules: list of structural rules for filtering molecules.
                verbose: verbosity level, 0 for errors, 1 for warnings, 2 for info.
        '''
        super().__init__(bb_supplier, reaction_tags, max_evals_per_comp, verbose)
        self.rules = rules
        self.struct_rules = struct_rules
    
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
    
    def set_query_mol(
        self, query_mol: Union[str, Chem.Mol], 
        reactive_sites: Optional[List[int]] = None
    ) -> None:
        '''
            Set the query molecule for enumeration and reactive sites.

            Args:
                query_mol: a SMILES string or an RDKit Mol object.
                reactive_sites: list of indices of reactive sites in the molecule.
        '''
        if isinstance(query_mol, str):
            self.query_mol = Chem.MolFromSmiles(query_mol)
        else:
            self.query_mol = query_mol
        flag = Chem.SanitizeMol(self.query_mol, catchErrors=True)
        assert flag == Chem.rdmolops.SanitizeFlags.SANITIZE_NONE, f"SanitizeMol failed: {flag}"

        self.reactive_sites = reactive_sites if reactive_sites else []

        self._compositions = []     # reset compositions

    def _process_query_mol(self, protect_neighbors: bool=False) -> None:
        if not isinstance(self.query_mol, Chem.Mol):
            raise ValueError('Query molecule must be an RDKit Mol object.')
        
        query_mol = self.query_mol
        if self.reactive_sites:
            dont_protect = set()
            for atom in query_mol.GetAtoms():
                if atom.GetIdx() in self.reactive_sites:
                    dont_protect.add(atom.GetIdx())
                    if not protect_neighbors:
                        for neighbor in atom.GetNeighbors():
                            dont_protect.add(neighbor.GetIdx())
            for atom in query_mol.GetAtoms():
                if atom.GetIdx() not in dont_protect:
                    atom.SetProp('_protected', '1')
        else:
            logger.warning('No reactive sites provided! All atoms will be considered reactive.')

        comp = CompositionPath(fragments=(query_mol,))
        self._compositions.append(comp)

    def _process_building_blocks(self):
        filtered_bbs = []
        for bb in tqdm(self.bb_mols, desc='Processing building blocks', total=len(self.bb_mols), disable=self.verbose < 2):
            if self._check_rules(bb) and self._check_struct_rules(bb):
                filtered_bbs.append(bb)
        self._compositions = [
            CompositionWithBBs(
                comp=comp,
                fragment_bbs=([BuildingBlock(comp.fragments[0])], filtered_bbs)
            ) for comp in self._compositions
        ]


    def _check_struct_rules(self, building_block: Union[Chem.Mol, BuildingBlock]) -> bool:
        '''
            Check if the building block satisfies the structure-based rules.
        '''
        if not self.struct_rules:
            return True
        for rule in self.struct_rules:
            if not building_block.HasSubstructMatch(Chem.MolFromSmarts(rule)):
                return False
        return True

    def _check_rules(self, building_block: Union[Chem.Mol, BuildingBlock]) -> bool:
        '''
            Check if the building block satisfies the rules.
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
        Molecule HEALER: enumerates molecules by splitting them into fragments and
        applying reactions to each fragment. It can also use custom split sites to
        generate compositions from the query molecule.
    '''
    def __init__(
            self, 
            bb_supplier: str='US_stock', 
            reaction_tags: list[str]=['amide coupling', 'amide', 'C-N bond formation', 'C-N',
                                      'alkylation', 'N-arylation', 'azole', 'amination'],
            max_evals_per_comp: Optional[int] = None,
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
                max_evals_per_comp: maximum number of evaluations for each composition.
                n_compositions: number of compositions to consider for enumeration.
                sim_threshold: similarity threshold for filtering building blocks.
                max_bbs_per_comp: maximum number of building blocks per fragment.
                    If <= 0, all building blocks will be considered. Otherwise, the similarity
                    threshold will be asjusted to the number of building blocks.
                verbose: verbosity level, 0 for errors, 1 for warnings, 2 for info.
        '''
        super().__init__(bb_supplier, reaction_tags, max_evals_per_comp, verbose)
        self.n_compositions = n_compositions
        self.sim_threshold = sim_threshold
        self.max_bbs_per_comp = max_bbs_per_comp

    def set_query_mol(
        self, query_mol: Union[str, Chem.Mol], 
        custom_split_sites: Optional[List[List[Tuple[int, int]]]] = None,
        retro_tree_depth: int = 1,
        min_frag_size: int = 3,
    ) -> None:
        '''
            Set the query molecule for enumeration and custom split sites.

            Args:
                query_mol: a SMILES string or an RDKit Mol object.
                custom_split_sites: Custom split sites for the molecule. If provided, 
                    the molecule will be split into fragments based on these sites.
                    Each site is a tuple of atom indices (start, end) to split the molecule.
                    A molecule can have multiple split sites to generate multiple fragments.
                    For example, if `custom_split_sites = [[(0, 2), (3, 5)], [(1, 4)]]`, then 
                    there will be two seaparate compositions generated from the molecule. 
                retro_tree_depth: depth of retrosynthesis tree to generate compositions.
                min_frag_size: minimum number of heavy atoms in a fragment to consider it valid.
        '''
        if isinstance(query_mol, str):
            self.query_mol = Chem.MolFromSmiles(query_mol, sanitize=False)
        else:
            self.query_mol = query_mol
        flag = Chem.SanitizeMol(self.query_mol, catchErrors=True)
        assert flag == Chem.rdmolops.SanitizeFlags.SANITIZE_NONE, f"SanitizeMol failed: {flag}"

        if len(Chem.GetMolFrags(self.query_mol, sanitizeFrags=False)) > 1:
            raise ValueError('Query molecule must be a single connected component. '
                             'Use FragmentHEALER for multi-component molecules.')

        self.custom_split_sites = custom_split_sites if custom_split_sites else []
        self.retro_tree_depth = retro_tree_depth
        self.min_frag_size = min_frag_size
        
        self._compositions = []     # reset compositions

    def _process_query_mol(self) -> None:
        '''
            Process the query molecule to generate compositions based on the 
            custom split sites if provided, or create retrosynthesis tree 
            with all possible compositions.
        '''
        if not isinstance(self.query_mol, Chem.Mol):
            raise ValueError('Query molecule must be an RDKit Mol object. Set it using set_query_mol() method.')
        
        if self.custom_split_sites:
            logger.info(f"Using custom split sites: {self.custom_split_sites}")
            for split_sites in self.custom_split_sites:
                fragments = self._split_molecule(split_sites)
                if len(fragments) < 2:
                    logger.warning(f'Custom split sites {split_sites} did not produce multiple fragments. '
                                   'Skipping this composition.')
                else:
                    comp = CompositionPath(fragments=fragments)
                    self._compositions.append(comp)
        else:
            logger.info("No custom split sites provided. Generating retrosynthesis tree.")
            retro_tree = RetrosynthesisTree(
                self.query_mol, 
                self.reactions, 
                max_depth=self.retro_tree_depth, 
                min_heavy_atoms=self.min_frag_size
            )
            retro_tree.build()
            self._compositions = retro_tree.get_composition_paths()

        # print compositions
        logger.info("Compositions generated from the query molecule: "
                    f"\n{self._composition_prints()}")

        if len(self._compositions) > self.n_compositions:
            self._compositions = self._compositions[:self.n_compositions]
            logger.info(f"Truncated compositions to the first {self.n_compositions}.")
        
    def _split_molecule(self, split_sites: List[Tuple[int, int]]) -> List[Chem.Mol]:
        '''
            Split the molecule into fragments based on the provided split sites.
            Returns a list of fragments as RDKit Mol objects.
        '''
        with Chem.RWMol(self.query_mol) as rw_mol:
            for start, end in split_sites:
                rw_mol.RemoveBond(start, end)
        mol = rw_mol.GetMol()
        return Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)

    def _composition_prints(self) -> str:
        '''
            Return a loggable string representation of the compositions.
        '''
        return '\n'.join(
            f'Composition {i+1} fragments: {[Chem.MolToSmiles(frag) for frag in comp.fragments]}'
            for i, comp in enumerate(self._compositions)
        ) if self._compositions else 'No compositions found.'

    def _process_building_blocks(self, batch_size: int=10000) -> None:
        '''
            Process building blocks to filter them based on the similarity 
            to the query molecule and the number of building blocks per composition 
            if given.
        '''
        bb_sizes = np.array([bb.GetNumHeavyAtoms() for bb in self.bb_mols])
        bb_fps = self._get_fingerprints(self.bb_mols)

        frag_lists = [path.fragments for path in self._compositions]
        offsets = np.concatenate(([0], np.cumsum([len(frag_list) for frag_list in frag_lists])))
        frags_flatten = [frag for frag_list in frag_lists for frag in frag_list]
        frag_sizes = np.array([frag.GetNumHeavyAtoms() for frag in frags_flatten])[:, None]
        frag_fps = self._get_fingerprints(frags_flatten)
        
        sims = np.zeros((len(frags_flatten), len(self.bb_mols)), dtype=float)
        for start in tqdm(range(0, len(self.bb_mols), batch_size), 
                          desc='Processing building blocks', 
                          total=len(self.bb_mols)//batch_size, 
                          disable=self.verbose < 2):
            end = start + batch_size
            batch_fps = bb_fps[start:end]
            batch_sizes = bb_sizes[start:end]
            
            delta = batch_sizes[None, :] - frag_sizes
            weights = 1 - np.clip(delta, 0, None) / batch_sizes[None, :]

            sims[:, start:end] = weights * utils.get_batch_tversky_sims(frag_fps, batch_fps)

        if self.max_bbs_per_comp > 0:
            kth = np.argpartition(-sims, self.max_bbs_per_comp-1, axis=1)[:, :self.max_bbs_per_comp]
            mask = np.zeros_like(sims, dtype=bool)
            rows = np.arange(sims.shape[0])[:, None]
            mask[rows, kth] = True
        else:
            mask = sims >= self.sim_threshold

        masks_per_comp = [
            mask[offsets[i]:offsets[i+1], :] for i in range(len(self._compositions))
        ]

        orig_comps = self._compositions
        self._compositions = [
            CompositionWithBBs(
                comp=comp,
                fragment_bbs=tuple(
                    [bb for bb, keep in zip(self.bb_mols, row) if keep] for row in comp_mask
                ))
            for comp, comp_mask in zip(orig_comps, masks_per_comp)
        ]


class FragmentHEALER(MoleculeHEALER):
    '''
        Fragment HEALER: a specialized version of MoleculeHEALER that allows 
        fragment inputs that will be used as compositions for enumeration.
    '''
    def __init__(
            self, 
            bb_supplier: str='US_stock', 
            reaction_tags: list[str]=['amide coupling', 'amide', 'C-N bond formation', 'C-N',
                                      'alkylation', 'N-arylation', 'azole', 'amination'],
            max_evals_per_comp: Optional[int] = None,
            n_compositions: int=10,
            sim_threshold: float=0.5,
            max_bbs_per_comp: int=-1,
            verbose: int=1,
    ):
        '''
            Initialize FragmentHEALER.

            Args:
                bb_supplier: one of "US_stock", "EU_stock" or "Global_stock"; or path to an SDF file.
                reaction_tags: list of tags or 'all'.
                max_evals_per_comp: maximum number of evaluations for each composition.
                n_compositions: number of compositions to consider for enumeration.
                sim_threshold: similarity threshold for filtering building blocks.
                max_bbs_per_comp: maximum number of building blocks per fragment.
                    If <= 0, all building blocks will be considered. Otherwise, the similarity
                    threshold will be asjusted to the number of building blocks.
                verbose: verbosity level, 0 for errors, 1 for warnings, 2 for info.
        '''
        super().__init__(
            bb_supplier, reaction_tags, max_evals_per_comp, 
            n_compositions, sim_threshold, max_bbs_per_comp, verbose
        )

    def set_query_mol(
        self, query_mol: Union[str, Chem.Mol, tuple[Chem.Mol, ...], tuple[str, ...]], 
    ) -> None:
        '''
            Set the query molecule for enumeration. The query can be a molecule with 
            multiple fragments or a tuple of molecules. 
        '''
        if isinstance(query_mol, str):
            query_mol = Chem.MolFromSmiles(query_mol, sanitize=False)
        elif isinstance(query_mol, tuple):
            if isinstance(query_mol[0], str):
                query_mol = tuple(Chem.MolFromSmiles(smi, sanitize=False) for smi in query_mol)
            mol_out = query_mol[0]
            for m in query_mol[1:]:
                mol_out = Chem.CombineMols(mol_out, m)
            query_mol = mol_out
        self.query_mol = query_mol
        
        frags = Chem.GetMolFrags(query_mol, sanitizeFrags=False)
        if len(frags) < 2:
            raise ValueError('Query molecule must have at least two fragments. '
                             'Use MoleculeHEALER for single-component molecules.')
        
        flag = Chem.SanitizeMol(self.query_mol, catchErrors=True)
        assert flag == Chem.rdmolops.SanitizeFlags.SANITIZE_NONE, f"SanitizeMol failed: {flag}"

        self._compositions = []     # reset compositions
        
    def _process_query_mol(self) -> None:
        '''
            Process the query molecule to generate compositions based on the 
            custom split sites if provided, or create retrosynthesis tree 
            with all possible compositions.
        '''
        if not isinstance(self.query_mol, Chem.Mol):
            raise ValueError('Query molecule must be an RDKit Mol object. '
                             'Set it using set_query_mol() method.')
        
        frags = Chem.GetMolFrags(self.query_mol, asMols=True, sanitizeFrags=True)
        self._compositions.append(
            CompositionPath(fragments=frags)
        )

        # print compositions
        logger.info("Compositions generated from the query molecule: "
                    f"\n{self._composition_prints()}")

