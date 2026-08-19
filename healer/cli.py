'''
    CLI for HEALER application.
'''
import healer.utils.rdkit_monkey_patch  # noqa: F401 - must be first
from healer.utils.progress import progress_bar

import json
import logging
import argparse
import tempfile
import webbrowser
import base64
from pathlib import Path
from typing import Optional

import pandas as pd
from rdkit import Chem
from rdkit.Chem import SDMolSupplier

from healer.application.healer import MoleculeHEALER, SiteHEALER, FragmentHEALER
import healer.utils.utils as utils

logger = logging.getLogger(__name__)


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


### Enumeration Command Handlers ###

def get_init_kwargs(args: argparse.Namespace, healer_type: str) -> dict:
    """Build __init__ kwargs for healer class."""
    kwargs = {
        'bb_source': args.bb_source,
        'reaction_tags': args.reactions.split(',') if args.reactions != 'all' else 'all',
        'shuffle_bb_order': args.shuffle,
        'show_progress': False if args.quiet else None,
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
        'n_jobs': args.n_jobs,
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


def run_enumeration(
    healer_type: str,
    smiles_list: list[str],
    init_kwargs: dict,
    query_kwargs: dict,
    enumerate_kwargs: dict,
    results_kwargs: dict,
    output_path: str,
    show_progress: Optional[bool] = None,
) -> None:
    """Run enumeration for one or more molecules.
    
    Synthesis-level parallelism is controlled by ``n_jobs`` inside
    *enumerate_kwargs* and handled by joblib inside the synthesis loop.
    """
    n_jobs = enumerate_kwargs.get('n_jobs', 1)
    logger.info(
        "Starting enumeration for %d molecule(s) (n_jobs=%s)",
        len(smiles_list), n_jobs,
    )
    
    # Initialize healer once — reused across all molecules
    if healer_type == 'molecule':
        healer = MoleculeHEALER(**init_kwargs)
    elif healer_type == 'site':
        healer = SiteHEALER(**init_kwargs)
    elif healer_type == 'fragment':
        healer = FragmentHEALER(**init_kwargs)
    
    out = Path(output_path)
    first = True
    n_written = 0
    n_failed = 0

    with progress_bar(smiles_list, desc="Enumerating", unit="mol",
                      show_progress=show_progress) as bar:
        for i, smiles in enumerate(bar, start=1):
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                logger.warning("Skipping invalid SMILES: %s", smiles)
                n_failed += 1
                continue

            logger.debug("[%d/%d] %s", i, len(smiles_list), smiles)
            try:
                healer.set_query_mol(query_mol=smiles, **query_kwargs)
                healer.enumerate(**enumerate_kwargs)
                df = healer.get_results(**results_kwargs)

                df.to_csv(str(out), mode='w' if first else 'a', header=first, index=False)
                first = False
                n_written += len(df)
            except Exception as e:
                logger.error("Error processing %s: %s", smiles, e)
                n_failed += 1

    if n_failed:
        logger.warning("%d of %d molecule(s) failed", n_failed, len(smiles_list))
    logger.info("Wrote %d row(s) to %s", n_written, out)


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
    run_enumeration(
        healer_type, smiles_list, init_kwargs, query_kwargs,
        enumerate_kwargs, results_kwargs, args.output,
        False if args.quiet else None,
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
                        help=f'Comma-separated reaction tags or "all" (default: all). '
                             f'Available tags: {", ".join(utils.get_reaction_tags())}')
    parser.add_argument('--shuffle', action='store_true',
                        help='Shuffle building block order')
    
    # Enumeration limits (approximate; bound run time rather than give exact counts)
    parser.add_argument('--max-evals', type=int, default=None,
                        help='Max reaction attempts per composition (approximate)')
    parser.add_argument('--max-products', type=int, default=None,
                        help='Max products per composition (approximate)')
    parser.add_argument('--max-total', type=int, default=None,
                        help='Max total products, stops enumeration (approximate)')
    
    # Output options
    parser.add_argument('--similarity', action='store_true',
                        help='Calculate similarity to query molecule')
    parser.add_argument('--properties', action='store_true',
                        help='Calculate molecular properties')
    
    # Execution
    parser.add_argument('--n-jobs', type=int, default=1,
                        help='Parallel threads for synthesis loop: 1=sequential, '
                             '-1=all CPUs (default: 1)')
    parser.add_argument('-v', '--verbose', action='count', default=1,
                        help='Increase verbosity (default: info, -v for debug)')
    parser.add_argument('-q', '--quiet', action='store_true',
                        help='Only report warnings and errors, and hide progress bars')


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
        if args.quiet:
            level = logging.WARNING
        elif args.verbose >= 2:
            level = logging.DEBUG
        else:
            level = logging.INFO
        logging.basicConfig(
            level=level,
            format="%(asctime)s [%(levelname)s] %(message)s",
            datefmt="%H:%M:%S"
        )
    
    # Execute command
    args.func(args)


if __name__ == '__main__':
    main()

