'''
    CLI for HEALER application.
'''
import healer.utils.rdkit_monkey_patch  # noqa: F401 - must be first

import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import json
import logging
import argparse
import tempfile
import webbrowser
import base64
from pathlib import Path
from multiprocessing import Pool
from functools import partial

import pandas as pd
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import SDMolSupplier

from healer.application.healer import MoleculeHEALER, SiteHEALER, FragmentHEALER
from healer.domain.bb_repository import get_repository
import healer.utils.utils as utils

logger = logging.getLogger(__name__)

# Global healer for worker processes
_worker_healer = None


### Input Loading ###

def load_input(input_path: str, column: str = 'smiles') -> list[str]:
    """
        Load SMILES from various input formats.
        
        Args:
            input_path: SMILES string, or path to .smi/.csv/.sdf file
            column: Column name for CSV files (default: 'smiles')
        
        Returns:
            List of SMILES strings
    """
    path = Path(input_path)
    
    # Direct SMILES string
    if not path.exists():
        mol = Chem.MolFromSmiles(input_path)
        if mol is not None:
            return [input_path]
        raise ValueError(f"Invalid SMILES or file not found: {input_path}")
    
    suffix = path.suffix.lower()
    
    if suffix == '.sdf':
        supplier = SDMolSupplier(str(path))
        return [Chem.MolToSmiles(mol) for mol in supplier if mol is not None]
    
    if suffix in ('.csv', '.smi', '.txt'):
        df = pd.read_csv(str(path))
        # Try to find SMILES column
        if column in df.columns:
            return df[column].dropna().tolist()
        # Fallback to first column
        return df.iloc[:, 0].dropna().tolist()
    
    raise ValueError(f"Unsupported file format: {suffix}")


### Config File Support ###

def load_config(config_path: str) -> dict:
    """Load configuration from JSON file."""
    with open(config_path, 'r') as f:
        return json.load(f)


def merge_args_with_config(args: argparse.Namespace, config: dict) -> argparse.Namespace:
    """Merge config file values with command-line args (CLI takes precedence)."""
    for key, value in config.items():
        # Only set if not explicitly provided on command line
        if not hasattr(args, key) or getattr(args, key) is None:
            setattr(args, key, value)
    return args


### View Command ###

def cmd_view(args: argparse.Namespace) -> None:
    """Show molecule with atom indices in browser."""
    mol = Chem.MolFromSmiles(args.smiles)
    if mol is None:
        logger.error("Invalid SMILES: %s", args.smiles)
        return
    
    # Generate SVG with atom indices
    svg_data_uri = utils.get_svg_mol(mol, show_idx=True)
    
    if args.output:
        # Save to file if requested
        svg_bytes = base64.b64decode(svg_data_uri.split(',')[1])
        with open(args.output, 'wb') as f:
            f.write(svg_bytes)
        logger.info("Saved molecule SVG to: %s", args.output)
    else:
        # Open in browser via temp file
        svg_bytes = base64.b64decode(svg_data_uri.split(',')[1])
        with tempfile.NamedTemporaryFile(suffix='.svg', delete=False) as f:
            f.write(svg_bytes)
            temp_path = f.name
        webbrowser.open(f'file://{temp_path}')


### Worker Initialization and Processing for Parallel Execution ###

def _init_worker(healer_type: str, init_kwargs: dict, bb_source: str) -> None:
    """Initialize healer in worker process."""
    global _worker_healer
    
    # Get shared repository (will use cached if available)
    bb_repo = get_repository(bb_source)
    if not bb_repo.is_loaded:
        bb_repo.load(show_progress=False)
    
    init_kwargs['bb_repository'] = bb_repo
    init_kwargs['bb_source'] = bb_source
    
    if healer_type == 'molecule':
        _worker_healer = MoleculeHEALER(**init_kwargs)
    elif healer_type == 'site':
        _worker_healer = SiteHEALER(**init_kwargs)
    elif healer_type == 'fragment':
        _worker_healer = FragmentHEALER(**init_kwargs)


def _process_molecule(
    smiles: str, 
    query_kwargs: dict, 
    enumerate_kwargs: dict,
    results_kwargs: dict
) -> pd.DataFrame:
    """Process a single molecule in worker."""
    global _worker_healer
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        logger.warning("Skipping invalid SMILES: %s", smiles)
        return pd.DataFrame()
    
    try:
        _worker_healer.set_query_mol(query_mol=smiles, **query_kwargs)
        _worker_healer.enumerate(**enumerate_kwargs)
        return _worker_healer.get_results(**results_kwargs)
    except Exception as e:
        logger.error("Error processing %s: %s", smiles, e)
        return pd.DataFrame()


### Enumeration Command Handlers ###

def get_init_kwargs(args: argparse.Namespace, healer_type: str) -> dict:
    """Build __init__ kwargs for healer class."""
    kwargs = {
        'bb_source': args.bb_source,
        'reaction_tags': args.reactions.split(',') if args.reactions != 'all' else 'all',
        'shuffle_bb_order': args.shuffle,
        'verbose': args.verbose,
    }
    
    if healer_type in ('molecule', 'fragment'):
        kwargs['sim_threshold'] = args.sim_threshold
        kwargs['max_bbs_per_frag'] = args.max_bbs_per_frag
    
    if healer_type == 'site':
        kwargs['rules'] = parse_rules(args.rules) if args.rules else {}
        kwargs['struct_rules'] = args.struct_rules.split(',') if args.struct_rules else []
    
    return kwargs


def get_query_kwargs(args: argparse.Namespace, healer_type: str) -> dict:
    """Build set_query_mol kwargs."""
    if healer_type == 'molecule':
        return {
            'n_compositions': args.n_compositions,
            'randomize_compositions': args.randomize,
            'random_seed': args.seed,
            'retro_tree_depth': args.retro_depth,
            'min_frag_size': args.min_frag_size,
        }
    elif healer_type == 'site':
        return {
            'reactive_sites': args.reactive_sites,
        }
    else:  # fragment
        return {}


def get_enumerate_kwargs(args: argparse.Namespace) -> dict:
    """Build enumerate() kwargs."""
    return {
        'max_evals_per_comp': args.max_evals,
        'max_products_per_comp': args.max_products,
        'max_total_products': args.max_total,
    }


def get_results_kwargs(args: argparse.Namespace) -> dict:
    """Build get_results() kwargs."""
    return {
        'calc_similarity': args.similarity,
        'calc_properties': args.properties,
    }


def parse_rules(rules_str: str) -> dict:
    """Parse rules string like 'MW:0:500,HBD:0:5' into dict."""
    rules = {}
    for rule in rules_str.split(','):
        parts = rule.strip().split(':')
        if len(parts) == 3:
            name, lo, hi = parts
            rules[name] = (int(lo), int(hi))
    return rules


def run_sequential(
    healer_type: str,
    smiles_list: list[str],
    init_kwargs: dict,
    query_kwargs: dict,
    enumerate_kwargs: dict,
    results_kwargs: dict,
    output_path: str,
    verbose: int
) -> None:
    """Run enumeration sequentially."""
    logger.info("Starting sequential enumeration for %d molecule(s)", len(smiles_list))
    
    # Initialize healer
    if healer_type == 'molecule':
        healer = MoleculeHEALER(**init_kwargs)
    elif healer_type == 'site':
        healer = SiteHEALER(**init_kwargs)
    elif healer_type == 'fragment':
        healer = FragmentHEALER(**init_kwargs)
    
    out = Path(output_path)
    first = True
    
    for smiles in tqdm(smiles_list, desc="Enumerating", disable=verbose >= 2):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning("Skipping invalid SMILES: %s", smiles)
            continue
        
        try:
            healer.set_query_mol(query_mol=smiles, **query_kwargs)
            healer.enumerate(**enumerate_kwargs)
            df = healer.get_results(**results_kwargs)
            
            df.to_csv(str(out), mode='w' if first else 'a', header=first, index=False)
            first = False
        except Exception as e:
            logger.error("Error processing %s: %s", smiles, e)
    
    logger.info("Results saved to %s", out)


def run_parallel(
    healer_type: str,
    smiles_list: list[str],
    init_kwargs: dict,
    query_kwargs: dict,
    enumerate_kwargs: dict,
    results_kwargs: dict,
    output_path: str,
    workers: int,
    verbose: int
) -> None:
    """Run enumeration in parallel."""
    logger.info("Starting parallel enumeration with %d workers", workers)
    
    # Extract bb_source for worker init
    bb_source = init_kwargs.get('bb_source')
    
    # Pre-load repository in main process
    bb_repo = get_repository(bb_source)
    if not bb_repo.is_loaded:
        bb_repo.load(show_progress=verbose >= 1)
    
    out = Path(output_path)
    first = True
    
    worker_fn = partial(
        _process_molecule,
        query_kwargs=query_kwargs,
        enumerate_kwargs=enumerate_kwargs,
        results_kwargs=results_kwargs
    )
    
    with Pool(
        processes=workers,
        initializer=_init_worker,
        initargs=(healer_type, init_kwargs, bb_source)
    ) as pool:
        for df in tqdm(
            pool.imap(worker_fn, smiles_list),
            total=len(smiles_list),
            desc="Enumerating",
            disable=verbose >= 2
        ):
            if df.empty:
                continue
            df.to_csv(str(out), mode='w' if first else 'a', header=first, index=False)
            first = False
    
    logger.info("Results saved to %s", out)


def cmd_enumerate(args: argparse.Namespace, healer_type: str) -> None:
    """Run enumeration for molecule/site/fragment commands."""
    # Load config if provided
    if args.config:
        config = load_config(args.config)
        args = merge_args_with_config(args, config)
    
    # Load input
    smiles_list = load_input(args.input, getattr(args, 'column', 'smiles'))
    logger.info("Loaded %d molecule(s) from input", len(smiles_list))
    
    # Build kwargs
    init_kwargs = get_init_kwargs(args, healer_type)
    query_kwargs = get_query_kwargs(args, healer_type)
    enumerate_kwargs = get_enumerate_kwargs(args)
    results_kwargs = get_results_kwargs(args)
    
    # Save args for reproducibility
    args_dict = vars(args).copy()
    args_dict.pop('func', None)  # Remove function reference
    with open('healer_args.json', 'w') as f:
        json.dump(args_dict, f, indent=2)
    
    # Run enumeration
    if args.workers > 1:
        run_parallel(
            healer_type, smiles_list, init_kwargs, query_kwargs,
            enumerate_kwargs, results_kwargs, args.output, args.workers, args.verbose
        )
    else:
        run_sequential(
            healer_type, smiles_list, init_kwargs, query_kwargs,
            enumerate_kwargs, results_kwargs, args.output, args.verbose
        )


### Argument Parser Construction ###

def add_common_args(parser: argparse.ArgumentParser) -> None:
    """Add arguments common to all enumeration commands."""
    # Input/Output
    parser.add_argument('input', help='SMILES string or input file (.smi, .csv, .sdf)')
    parser.add_argument('-o', '--output', default='healer_results.csv',
                        help='Output CSV path (default: healer_results.csv)')
    parser.add_argument('--column', default='smiles',
                        help='SMILES column name in CSV (default: smiles)')
    parser.add_argument('--config', help='JSON config file path')
    
    # Building blocks and reactions
    parser.add_argument('--bb-source', default='US_stock',
                        help='Building block source: US_stock, EU_stock, Global_stock, or path')
    parser.add_argument('--reactions', default='all',
                        help='Comma-separated reaction tags or "all" (default: all)')
    parser.add_argument('--shuffle', action='store_true',
                        help='Shuffle building block order')
    
    # Enumeration limits
    parser.add_argument('--max-evals', type=int, default=None,
                        help='Max reaction attempts per composition')
    parser.add_argument('--max-products', type=int, default=None,
                        help='Max products per composition')
    parser.add_argument('--max-total', type=int, default=None,
                        help='Max total products (stops enumeration)')
    
    # Output options
    parser.add_argument('--similarity', action='store_true',
                        help='Calculate similarity to query molecule')
    parser.add_argument('--properties', action='store_true',
                        help='Calculate molecular properties')
    
    # Execution
    parser.add_argument('--workers', type=int, default=1,
                        help='Number of parallel workers (default: 1)')
    parser.add_argument('-v', '--verbose', action='count', default=1,
                        help='Increase verbosity (-v for info, -vv for debug)')


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser with subcommands."""
    parser = argparse.ArgumentParser(
        prog='healer',
        description='HEALER: Hit Expansion by Assembling Ligands from Enumerated Reactions'
    )
    subparsers = parser.add_subparsers(dest='command', required=True)
    
    # --- molecule subcommand ---
    mol_parser = subparsers.add_parser('molecule', help='Enumerate from whole molecules')
    add_common_args(mol_parser)
    mol_parser.add_argument('--sim-threshold', type=float, default=0.5,
                            help='Similarity threshold for BB matching (default: 0.5)')
    mol_parser.add_argument('--max-bbs-per-frag', type=int, default=-1,
                            help='Max BBs per fragment, -1 for unlimited (default: -1)')
    mol_parser.add_argument('--n-compositions', type=int, default=10,
                            help='Number of compositions to consider (default: 10)')
    mol_parser.add_argument('--retro-depth', type=int, default=1,
                            help='Retrosynthesis tree depth (default: 1)')
    mol_parser.add_argument('--min-frag-size', type=int, default=3,
                            help='Minimum fragment size in heavy atoms (default: 3)')
    mol_parser.add_argument('--randomize', action='store_true',
                            help='Randomize composition order')
    mol_parser.add_argument('--seed', type=int, default=-1,
                            help='Random seed for reproducibility (default: -1)')
    mol_parser.set_defaults(func=lambda args: cmd_enumerate(args, 'molecule'))
    
    # --- site subcommand ---
    site_parser = subparsers.add_parser('site', help='Enumerate at specific reactive sites')
    add_common_args(site_parser)
    site_parser.add_argument('--reactive-sites', type=json.loads, default=None,
                             help='Atom indices for reactive sites as JSON list, e.g., "[1,2,5]"')
    site_parser.add_argument('--rules', default=None,
                             help='BB filter rules: "MW:0:500,HBD:0:5,..." (default: none)')
    site_parser.add_argument('--struct-rules', default=None,
                             help='SMARTS patterns for BB filtering, comma-separated')
    site_parser.set_defaults(func=lambda args: cmd_enumerate(args, 'site'))
    
    # --- fragment subcommand ---
    frag_parser = subparsers.add_parser('fragment', help='Enumerate from pre-split fragments')
    add_common_args(frag_parser)
    frag_parser.add_argument('--sim-threshold', type=float, default=0.5,
                             help='Similarity threshold for BB matching (default: 0.5)')
    frag_parser.add_argument('--max-bbs-per-frag', type=int, default=-1,
                             help='Max BBs per fragment, -1 for unlimited (default: -1)')
    frag_parser.set_defaults(func=lambda args: cmd_enumerate(args, 'fragment'))
    
    # --- view subcommand ---
    view_parser = subparsers.add_parser('view', help='View molecule with atom indices')
    view_parser.add_argument('smiles', help='SMILES string to visualize')
    view_parser.add_argument('-o', '--output', default=None,
                             help='Save SVG to file instead of opening in browser')
    view_parser.set_defaults(func=cmd_view)
    
    return parser


def main():
    """Main entry point."""
    parser = build_parser()
    args = parser.parse_args()
    
    # Configure logging based on verbosity
    if args.command != 'view':
        level = logging.DEBUG if args.verbose >= 2 else logging.INFO if args.verbose == 1 else logging.WARNING
        logging.basicConfig(
            level=level,
            format="%(asctime)s [%(levelname)s] %(message)s",
            datefmt="%H:%M:%S"
        )
    
    # Execute command
    args.func(args)


if __name__ == '__main__':
    main()

