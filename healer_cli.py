#!/usr/bin/env python3
"""
CLI for the HEALER molecule/site enumerator with optional parallel execution.
"""

import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import logging
import argparse
from pathlib import Path
import textwrap
from multiprocessing import Pool
from itertools import chain

import pandas as pd
from tqdm import tqdm

import utils
from enumerator import MoleculeEnumerator, SiteEnumerator
from rdkit.Chem.FastSDMolSupplier import FastSDMolSupplier
from rdkit import Chem

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

enumerator = None  # global for worker processes

def _init_worker(enumeration_type, bb_source, reaction_tags,
                 sim_threshold, max_bb, rules, struct_rules, verbose):
    """Initialize the global enumerator in each worker."""
    global enumerator
    if enumeration_type == 'molecule':
        enumerator = MoleculeEnumerator(
            building_blocks=bb_source,
            reaction_tags=reaction_tags,
            sim_threshold=sim_threshold,
            max_bbs_per_comp=max_bb,
            verbose=verbose,
        )
    else:
        enumerator = SiteEnumerator(
            building_blocks=bb_source,
            reaction_tags=reaction_tags,
            rules=rules,
            struct_rules=struct_rules,
            verbose=verbose,
        )

def _process_smiles(smi):
    """Enumerate a single SMILES string using the worker's enumerator."""
    if isinstance(smi, str):
        mol = Chem.MolFromSmiles(smi)
    if not mol:
        logger.warning("Skipping invalid SMILES.")
    enumerator.enumerate(molecule=smi)
    return enumerator.get_results()

class HEALERCLI:
    """Command-line interface for the HEALER enumerator."""
    def __init__(self):
        self.args = self.parse_args()
        self.reactions = self.load_reactions()
        self.reaction_tags = self.parse_reaction_tags()
        self.rules = self.parse_rules()
        self.struct_rules = self.parse_struct_rules()
        self.smiles_list = self.load_smiles()

    def parse_args(self):
        parser = argparse.ArgumentParser(
            description=textwrap.dedent('''\
                Enumerate molecules or sites using HEALER.
                Supports single-threaded or parallel execution.
            '''),
            formatter_class=argparse.RawTextHelpFormatter,
        )
        parser.add_argument('enumeration_type', choices=['molecule','site'],
                            help='Type of enumeration: molecule or site.')
        # input args
        parser.add_argument('input_smiles', help='SMILES string or input file (.smi, .csv, .sdf).')
        parser.add_argument('--header', action='store_true', help='Input CSV has a header row.')
        parser.add_argument('--column_name', default='smiles',
                            help='SMILES column in CSV if --header is set.')
        # general enumeration args
        parser.add_argument('--output', default='enumerations.csv',
                            help='Path to output CSV.')
        parser.add_argument('--bb_source', choices=['US_stock','EU_stock','Global_stock'], default='US_stock',
                            help='Building block source.')
        parser.add_argument('--reaction_tags', default='amide coupling,amide,alkylation,N-arylation,azole,amination',
                            help='Comma-separated reaction tags, or "all" for all valid tags.')
        # site enumeration args
        parser.add_argument('--MW_range', default='0:500',
                            help='Molecular weight range (min:max).')
        parser.add_argument('--HBD_range', default='0:5',
                            help='HBD count range (min:max).')
        parser.add_argument('--HBA_range', default='0:10',
                            help='HBA count range (min:max).')
        parser.add_argument('--TPSA_range', default='0:200',
                            help='TPSA range (min:max).')
        parser.add_argument('--RotB_range', default='0:10',
                            help='Rotatable bonds range (min:max).')
        parser.add_argument('--Rings_range', default='0:10',
                            help='Rings range (min:max).')
        parser.add_argument('--ArRings_range', default='0:5',
                            help='Aromatic rings range (min:max).')
        parser.add_argument('--Chiral_range', default='0:5',
                            help='Chiral centers range (min:max).')
        parser.add_argument('--structure_based_rules', default='',
                            help='Dot-separated SMARTS rules for building blocks.')
        # molecule enumeration args
        parser.add_argument('--sim_threshold', type=float, default=0.50,
                            help='Similarity threshold for building block matching.')
        parser.add_argument('--max_bb', type=int, default=100,
                            help='Max building blocks per split.')
        # other args
        parser.add_argument('--workers', type=int, default=1,
                            help='Number of parallel workers (1 for sequential).')
        parser.add_argument('--verbose', action='store_true',
                            help='Enable verbose output.')
        return parser.parse_args()

    def load_reactions(self):
        reactions = utils.load_reactions_from_json('reactions/reactions.json')
        return [r for r in reactions if r.is_valid()]

    def parse_reaction_tags(self):
        tags = self.args.reaction_tags.split(',')
        available = set(chain(*(r.tags for r in self.reactions)))
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
        return [self.args.input_smiles]

    def run_sequential(self):
        logger.info("Starting sequential enumeration for %d molecules.", len(self.smiles_list))
        out = Path(self.args.output)
        if out.exists():
            logger.info("Overwriting %s", out)
            out.unlink()

        if self.args.enumeration_type == 'molecule':
            enum = MoleculeEnumerator(
                building_blocks=self.args.bb_source,
                reaction_tags=self.reaction_tags,
                sim_threshold=self.args.sim_threshold,
                max_bbs_per_comp=self.args.max_bb,
                verbose=self.args.verbose,
            )
        else:
            enum = SiteEnumerator(
                building_blocks=self.args.bb_source,
                reaction_tags=self.reaction_tags,
                rules=self.rules,
                struct_rules=self.struct_rules,
                verbose=self.args.verbose,
            )

        first = True
        for smi in tqdm(self.smiles_list, desc="Enumerating"):
            if isinstance(smi, str):
                mol = Chem.MolFromSmiles(smi)
            if not mol:
                logger.warning("Skipping invalid SMILES.")
                continue
            enum.enumerate(molecule=mol)
            df = enum.get_results()
            df.to_csv(str(out), mode='w' if first else 'a', header=first, index=False)
            first = False

        logger.info("Results saved to %s", out)

    def run_parallel(self):
        logger.info("Starting parallel enumeration with %d workers.", self.args.workers)
        init_args = (
            self.args.enumeration_type,
            self.args.bb_source,
            self.reaction_tags,
            self.args.sim_threshold,
            self.args.max_bb,
            self.rules,
            self.struct_rules,
            self.args.verbose,
        )
        out = Path(self.args.output)
        if out.exists():
            logger.info("Overwriting %s", out)
            out.unlink()

        first = True
        with Pool(processes=self.args.workers,
                  initializer=_init_worker,
                  initargs=init_args) as pool:
            for df in tqdm(pool.imap(_process_smiles, self.smiles_list),
                           total=len(self.smiles_list),
                           desc="Enumerating"):
                df.to_csv(str(out), mode='w' if first else 'a', header=first, index=False)
                first = False

        logger.info("Results saved to %s", out)

    def run(self):
        if self.args.workers > 1:
            self.run_parallel()
        else:
            self.run_sequential()

if __name__ == "__main__":
    HEALERCLI().run()
