'''
    CLI for the HEALER molecule/site healer with optional parallel execution.
'''
import healer.utils.rdkit_monkey_patch as rdkit_monkey_patch

import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import json
import logging
import argparse
import textwrap
from pathlib import Path
from multiprocessing import Pool
from itertools import chain
from functools import partial

import pandas as pd
from tqdm import tqdm

import healer.utils.utils as utils
from healer.application.healer import MoleculeHEALER, SiteHEALER, FragmentHEALER
from rdkit.Chem.FastSDMolSupplier import FastSDMolSupplier
from rdkit import Chem

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
logger = logging.getLogger(__name__)


def _load_reactions(file_path):
    reactions = utils.load_reactions_from_json(file_path)
    return [r for r in reactions if r.is_valid()]

_root = Path(__file__).parent
_rxn_data_path = Path(_root, 'data', 'reactions', 'reactions.json')
REACTIONS = _load_reactions(_rxn_data_path)
healer = None  # global for worker processes

def _init_worker(healer_type, bb_source, reaction_tags, max_evals_per_comp, shuffle_bbs,
                 sim_threshold, max_bb, rules, struct_rules, verbose):
    """Initialize the global healer in each worker."""
    global healer
    if healer_type == 'molecule':
        healer = MoleculeHEALER(
            bb_supplier=bb_source,
            reaction_tags=reaction_tags,
            max_evals_per_comp=max_evals_per_comp,
            shuffle_bbs=shuffle_bbs,
            sim_threshold=sim_threshold,
            max_bbs_per_comp=max_bb,
            verbose=verbose,
        )
    elif healer_type == 'site':
        healer = SiteHEALER(
            bb_supplier=bb_source,
            reaction_tags=reaction_tags,
            max_evals_per_comp=max_evals_per_comp,
            shuffle_bbs=shuffle_bbs,
            rules=rules,
            struct_rules=struct_rules,
            verbose=verbose,
        )
    elif healer_type == 'fragment':
        healer = FragmentHEALER(
            bb_supplier=bb_source,
            reaction_tags=reaction_tags,
            max_evals_per_comp=max_evals_per_comp,
            shuffle_bbs=shuffle_bbs,
            sim_threshold=sim_threshold,
            max_bbs_per_comp=max_bb,
            verbose=verbose,
        )
    else:
        raise ValueError(f"Unknown HEALER type: {healer_type}")

def _process_batch(smi_batch, query_kwargs: dict, results_kwargs: dict):
    '''Enumerate a batch of SMILES using the worker's healer.'''
    dfs = []
    for smi in smi_batch:
        mol = Chem.MolFromSmiles(smi) if isinstance(smi, str) else smi
        if not mol:
            logger.warning("Skipping invalid SMILES: %s", smi)
            continue
        healer.set_query_mol(query_mol=smi, **query_kwargs)
        healer.enumerate()
        dfs.append(healer.get_results(**results_kwargs))
    return pd.concat(dfs, ignore_index=True)


class HEALERCLI:
    '''Command-line interface for HEALER.'''
    def __init__(self):
        self.args = self.parse_args()
        self.reaction_tags = self.parse_reaction_tags()
        self.rules = self.parse_rules()
        self.struct_rules = self.parse_struct_rules()
        self.smiles_list = self.load_smiles()
        self.verbose = 2 if self.args.verbose else 1

        self.results_kwargs = {
            'calc_similarity': self.args.calculate_similarity,
            'calc_stoplight': self.args.calculate_stoplight,
            'calc_cns_mpo': self.args.calculate_cnsmpo,
        }

        logger.setLevel(logging.WARNING if self.verbose == 1 else logging.INFO)

    def parse_args(self):
        parser = argparse.ArgumentParser(
            description=textwrap.dedent('''\
                Enumerate molecules, sites or fragments using HEALER.
                Supports single-threaded or parallel execution.
            '''),
            formatter_class=argparse.RawTextHelpFormatter,
        )
        parser.add_argument('healer_type', choices=['molecule','site','fragment'],
                            help='Type of enumeration: molecule or site.')
        # input args
        parser.add_argument('input_smiles', help='SMILES string or input file (.smi, .csv, .sdf).')
        parser.add_argument('--header', action='store_true', help='Input CSV has a header row.')
        parser.add_argument('--column_name', default='smiles',
                            help='SMILES column in CSV if --header is set. (default: smiles)')
        
        # general enumeration args
        parser.add_argument('--output', default='healer_results.csv',
                            help='Path to output CSV. (default: healer_results.csv)')
        parser.add_argument('--bb_source', choices=['US_stock','EU_stock','Global_stock'], default='US_stock',
                            help='Building block source. (default: US_stock)')
        parser.add_argument('--reaction_tags', default='amide coupling,amide,alkylation,N-arylation,azole,amination',
                            choices=['all'] + list(set(chain(*(r.tags for r in REACTIONS)))),
                            help=f'Comma-separated reaction tags, or "all" for all valid tags. '
                                 f'(Default: amide coupling,amide,alkylation,N-arylation,azole,amination)')
        parser.add_argument('--max_evals_per_comp', type=int, default=1000,
                            help='Maximum number of evaluations for each composition. (default: 1000)')
        parser.add_argument('--shuffle_bbs', action='store_true',
                            help='Shuffle building blocks loaded from the source.')
        parser.add_argument('--calculate_similarity', action='store_true',
                            help='Calculate similarity between query and enumerated molecules.')
        parser.add_argument('--calculate_stoplight', action='store_true',
                            help='Calculate stoplight scores for enumerated molecules.')
        parser.add_argument('--calculate_cnsmpo', action='store_true',
                            help='Calculate CNS MPO scores for enumerated molecules.')
        
        # site enumeration args
        parser.add_argument('--MW_range', default='0:500',
                            help='Molecular weight range (min:max). (default: 0:500)')
        parser.add_argument('--HBD_range', default='0:5',
                            help='HBD count range (min:max). (default: 0:5)')
        parser.add_argument('--HBA_range', default='0:10',
                            help='HBA count range (min:max). (default: 0:10)')
        parser.add_argument('--TPSA_range', default='0:200',
                            help='TPSA range (min:max). (default: 0:200)')
        parser.add_argument('--RotB_range', default='0:10',
                            help='Rotatable bonds range (min:max). (default: 0:10)')
        parser.add_argument('--Rings_range', default='0:10',
                            help='Rings range (min:max). (default: 0:10)')
        parser.add_argument('--ArRings_range', default='0:5',
                            help='Aromatic rings range (min:max). (default: 0:5)')
        parser.add_argument('--Chiral_range', default='0:5',
                            help='Chiral centers range (min:max). (default: 0:5)')
        parser.add_argument('--structure_based_rules', default='',
                            help='Dot-separated SMARTS rules for building blocks. (default: None)')
        parser.add_argument('--reactive_sites', type=json.loads, default=None,
                            help='Reactive site indices for site healer (list of ints). (default: None)')
        
        # molecule enumeration args, sim_threshold and max_bb are used for both molecule and fragment healers
        parser.add_argument('--sim_threshold', type=float, default=0.40,
                            help='Similarity threshold for building block matching. (default: 0.40)')
        parser.add_argument('--max_bb', type=int, default=10,
                            help='Max building blocks per split. (default: 10)')
        parser.add_argument('--n_compositions', type=int, default=1000,
                            help='Number of compositions to consider for an enumeration. (default: 1000)')
        parser.add_argument('--randomize_compositions', action='store_true',
                            help='Randomize the order of compositions.')
        parser.add_argument('--random_seed', type=int, default=-1,
                            help='Random seed for randomizing compositions (-1 for no seed). (default: -1)')
        parser.add_argument('--custom_split_sites', type=json.loads, default=None,
                            help='Custom split sites for molecule healer (list of [ [i,j], … ]). (default: None)')
        parser.add_argument('--retro_tree_depth', type=int, default=2,
                            help='Depth of the retro synthesis tree for molecule healer. (Default is 2)')
        parser.add_argument('--min_frag_size', type=int, default=3,
                            help='Minimum fragment size for molecule healer. (Default is 3)')
        
        # other args
        parser.add_argument('--workers', type=int, default=1,
                            help='Number of parallel workers (1 for sequential). (default: 1)')
        parser.add_argument('--verbose', action='store_true',
                            help='Enable verbose output.')
        return parser.parse_args()

    def parse_reaction_tags(self):
        tags = self.args.reaction_tags.split(',')
        available = set(chain(*(r.tags for r in REACTIONS)))
        if 'all' in tags:
            return list(available)
        invalid = [t for t in tags if t not in available]
        if invalid:
            logger.warning("Invalid tags %s; they will be ignored.", invalid)
        return [t for t in tags if t in available]

    def parse_rules(self):
        rules = {}
        for prop in ['MW','HBD','HBA','TPSA','RotB','Rings','ArRings','Chiral']:
            val = getattr(self.args, f'{prop}_range')
            lo, hi = map(int, val.split(':'))
            rules[prop] = (lo, hi)
        return rules

    def parse_struct_rules(self):
        return (self.args.structure_based_rules.split('.')
                if self.args.structure_based_rules else [])

    def load_smiles(self):
        input_path = Path(self.args.input_smiles)
        suffix = input_path.suffix.lower()
        if suffix == '.sdf':
            return FastSDMolSupplier(str(input_path))
        if suffix in ('.csv','.smi'):
            df = pd.read_csv(str(input_path),
                             header=0 if self.args.header else None,
                             usecols=[self.args.column_name if self.args.header else 0])
            col = self.args.column_name if self.args.header else df.columns[0]
            return df[col].tolist()
        if Chem.MolFromSmiles(self.args.input_smiles):
            return [self.args.input_smiles]
        raise ValueError(f"Unsupported input format: {suffix}. "
                         "Please provide a valid SMILES string, SDF file, or CSV file.")
    
    def get_output_columns(self):
        sample_healer = MoleculeHEALER('test', verbose=0)
        smi = "CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O"
        sample_healer.set_query_mol(query_mol=smi)
        sample_healer.enumerate()
        sample_df = sample_healer.get_results(**self.results_kwargs)

        fixed_cols = ['ID', 'Product']
        bb_cols, rxn_cols, url_cols = [], [], []
        if self.args.healer_type == 'molecule':
            max_depth = self.args.retro_tree_depth
        elif self.args.healer_type == 'fragment':
            max_depth = 0
            for mol in self.smiles_list:
                if isinstance(mol, str):
                    mol = Chem.MolFromSmiles(mol)
                max_depth = max(max_depth, len(Chem.GetMolFrags(mol, sanitizeFrags=False)))
        elif self.args.healer_type == 'site':
            max_depth = 1
        
        bb_cols = [f'BB{i}' for i in range(1, max_depth + 2)]
        rxn_cols = [f'Reaction{i}_name' for i in range(1, max_depth + 1)]
        url_cols = [f'URL{i}' for i in range(1, max_depth + 2)]
        prop_cols = [c for c in sample_df.columns if c not in fixed_cols + bb_cols + rxn_cols + url_cols]

        return fixed_cols + bb_cols + rxn_cols + url_cols + prop_cols

    def run_sequential(self):
        logger.info("Starting sequential enumeration for %d molecule(s).", len(self.smiles_list))
        out = Path(self.args.output)
        if out.exists():
            logger.warning("Overwriting %s", out)
            out.unlink()

        if self.args.healer_type == 'molecule':
            healer = MoleculeHEALER(
                bb_supplier=self.args.bb_source,
                reaction_tags=self.reaction_tags,
                max_evals_per_comp=self.args.max_evals_per_comp,
                shuffle_bbs=self.args.shuffle_bbs,
                sim_threshold=self.args.sim_threshold,
                max_bbs_per_comp=self.args.max_bb,
                verbose=self.verbose,
            )
            query_kwargs = {
                'n_compositions': self.args.n_compositions,
                'randomize_compositions': self.args.randomize_compositions,
                'random_seed': self.args.random_seed,
                'custom_split_sites': self.args.custom_split_sites,
                'retro_tree_depth': self.args.retro_tree_depth,
                'min_frag_size': self.args.min_frag_size,
            }
        elif self.args.healer_type == 'site':
            healer = SiteHEALER(
                bb_supplier=self.args.bb_source,
                reaction_tags=self.reaction_tags,
                max_evals_per_comp=self.args.max_evals_per_comp,
                shuffle_bbs=self.args.shuffle_bbs,
                rules=self.rules,
                struct_rules=self.struct_rules,
                verbose=self.verbose,
            )
            query_kwargs = {
                'reactive_sites': self.args.reactive_sites,
            }
        elif self.args.healer_type == 'fragment':
            healer = FragmentHEALER(
                bb_supplier=self.args.bb_source,
                reaction_tags=self.reaction_tags,
                max_evals_per_comp=self.args.max_evals_per_comp,
                shuffle_bbs=self.args.shuffle_bbs,
                sim_threshold=self.args.sim_threshold,
                max_bbs_per_comp=self.args.max_bb,
                verbose=self.verbose,
            )
            query_kwargs = {}
        else:
            raise ValueError(f"Unknown HEALER type: {self.args.healer_type}")
        
        output_columns = self.get_output_columns()
        first = True
        for smi in tqdm(self.smiles_list, desc="Enumerating", disable=self.verbose == 2, unit="molecule"):
            mol = Chem.MolFromSmiles(smi) if isinstance(smi, str) else smi
            if not mol:
                logger.warning("Skipping invalid SMILES: %s", smi)
                continue
            healer.set_query_mol(query_mol=smi, **query_kwargs)
            healer.enumerate()
            df = healer.get_results(**self.results_kwargs)
            df = df.reindex(columns=output_columns, fill_value='')
            df.to_csv(
                str(out), mode='w' if first else 'a', 
                header=first, index=False,
                columns=output_columns
            )
            first = False

        logger.info("Results saved to %s", out)

    def run_parallel(self):
        logger.info("Starting parallel enumeration with %d workers.", self.args.workers)
        init_args = (
            self.args.healer_type,
            self.args.bb_source,
            self.reaction_tags,
            self.args.max_evals_per_comp,
            self.args.shuffle_bbs,
            self.args.sim_threshold,
            self.args.max_bb,
            self.rules,
            self.struct_rules,
            self.verbose,
        )

        if self.args.healer_type == 'molecule':
            query_kwargs = {
                'n_compositions': self.args.n_compositions,
                'randomize_compositions': self.args.randomize_compositions,
                'random_seed': self.args.random_seed,
                'custom_split_sites': self.args.custom_split_sites,
                'retro_tree_depth': self.args.retro_tree_depth,
                'min_frag_size': self.args.min_frag_size,
            }
        elif self.args.healer_type == 'site':
            query_kwargs = {
                'reactive_sites': self.args.reactive_sites,
            }
        elif self.args.healer_type == 'fragment':
            query_kwargs = {}

        out = Path(self.args.output)
        if out.exists():
            logger.info("Overwriting %s", out)
            out.unlink()

        total = len(self.smiles_list)
        n_workers = self.args.workers
        batch_size = max(1, total // (n_workers * 4)) # 4 batches per worker
        batches = [
            self.smiles_list[i:i + batch_size]
            for i in range(0, total, batch_size)
        ]

        worker_fn = partial(
            _process_batch, 
            query_kwargs=query_kwargs, 
            results_kwargs=self.results_kwargs
        )
        output_columns = self.get_output_columns()
        first = True
        with Pool(processes=self.args.workers,
                  initializer=_init_worker,
                  initargs=init_args) as pool:
            for df in tqdm(pool.imap(worker_fn, batches),
                           total=len(batches), desc="Enumerating",
                           disable=self.verbose == 2, unit="batch"):
                df = df.reindex(columns=output_columns, fill_value='')
                df.to_csv(
                    str(out), mode='w' if first else 'a', 
                    header=first, index=False,
                    columns=output_columns
                )
                first = False
        
        logger.info("Results saved to %s", out)

    def run(self):
        args_dict = vars(self.args)
        args_dict['reaction_tags'] = self.reaction_tags
        args_json = json.dumps(args_dict, indent=4)
        with open('healer_args.json', 'w') as f:
            f.write(args_json)
        logger.info("HEALER arguments saved to healer_args.json")
    
        if self.args.workers > 1:
            self.run_parallel()
        else:
            self.run_sequential()

if __name__ == "__main__":
    HEALERCLI().run()

