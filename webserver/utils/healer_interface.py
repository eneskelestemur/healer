'''
    Interface utilities for HEALER classes to standardize web app interactions.
'''
import healer.utils.rdkit_monkey_patch as rdkit_monkey_patch

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

from typing import List, Dict, Any, Optional, Union, Tuple
import logging
from pathlib import Path
from rdkit import Chem

from healer.application.healer import MoleculeHEALER, SiteHEALER, FragmentHEALER
import healer.utils.utils as utils

logger = logging.getLogger(__name__)

# Get absolute paths for building blocks and reactions
HEALER_ROOT = Path(__file__).parent.parent.parent / 'healer'
BB_BASE_PATH = HEALER_ROOT / 'data' / 'buildingblocks'
REACTIONS_PATH = HEALER_ROOT / 'data' / 'reactions' / 'reactions.json'

# Building block paths with absolute paths
BB_PATHS = {
    "US_stock": str(BB_BASE_PATH / "Enamine_Rush-Delivery_Building_Blocks-US" / "*_processed.sdf"),
    "EU_stock": str(BB_BASE_PATH / "Enamine_Rush-Delivery_Building_Blocks-EU" / "*_processed.sdf"),
    "Global_stock": str(BB_BASE_PATH / "Enamine_Building_Blocks_Stock" / "*_processed.sdf"),
    "test": str(BB_BASE_PATH / "test_100_bb_processed.sdf"),
}


def count_molecular_fragments(smiles: str) -> int:
    """
    Count the number of molecular fragments in a SMILES string.
    
    Args:
        smiles: SMILES string to analyze
        
    Returns:
        Number of fragments (disconnected components)
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return 0
        
        # Split into fragments
        fragments = Chem.GetMolFrags(mol, asMols=True)
        return len(fragments)
    except Exception as e:
        logger.error(f"Error counting fragments in {smiles}: {e}")
        return 0


def validate_server_parameters(params: Dict[str, Any], healer_type: str = "molecule") -> Dict[str, Any]:
    """
    Validate and limit parameters when running in server mode.
    
    Args:
        params: Dictionary of parameters to validate
        healer_type: Type of healer ("molecule", "fragment", or "site")
        
    Returns:
        Dictionary of validated/limited parameters
    """
    import os
    server_mode = os.environ.get('HEALER_SERVER_MODE', 'false').lower() == 'true'
    
    if not server_mode:
        return params  # No limitations in local mode
    
    validated_params = params.copy()
    
    # Common limitations for all healers
    if 'reaction_tags' in validated_params and validated_params['reaction_tags']:
        # Remove 'all' tag in server mode as it can cause resource issues
        original_length = len(validated_params['reaction_tags'])
        if 'all' in validated_params['reaction_tags']:
            validated_params['reaction_tags'] = [tag for tag in validated_params['reaction_tags'] if tag != 'all']
            logger.warning("'all' reaction tag removed in server mode for resource management")
        
        # Limit number of reaction tags
        if len(validated_params['reaction_tags']) > 15:
            validated_params['reaction_tags'] = validated_params['reaction_tags'][:15]
            logger.warning(f"Reaction tags limited to 15 in server mode (was {original_length})")
        
        # Ensure we still have some reaction tags after filtering
        if not validated_params['reaction_tags']:
            validated_params['reaction_tags'] = ["amide coupling", "amide", "C-N bond formation", "C-N",
                                               "alkylation", "N-arylation", "azole", "amination"]
            logger.warning("No reaction tags remaining after server filtering, using default set")
    
    # Molecule/Fragment healer specific limitations
    if healer_type in ["molecule", "fragment"]:
        # Similarity threshold
        if 'sim_threshold' in validated_params:
            if validated_params['sim_threshold'] < 0.3:
                validated_params['sim_threshold'] = 0.3
                logger.warning("Similarity threshold limited to >= 0.3 in server mode")
        
        # Max BBs per composition
        if 'max_bbs_per_comp' in validated_params:
            if validated_params['max_bbs_per_comp'] > 10 or validated_params['max_bbs_per_comp'] < 1:
                validated_params['max_bbs_per_comp'] = min(10, max(1, validated_params['max_bbs_per_comp']))
                logger.warning("Max BBs/Comp limited to 1-10 in server mode")
        
        # N compositions
        if 'n_compositions' in validated_params:
            if validated_params['n_compositions'] > 50:
                validated_params['n_compositions'] = 50
                logger.warning("N compositions limited to 50 in server mode")
        
        # Max evaluations
        if 'max_evals_per_comp' in validated_params and validated_params['max_evals_per_comp']:
            if validated_params['max_evals_per_comp'] > 500 or validated_params['max_evals_per_comp'] < 1:
                validated_params['max_evals_per_comp'] = min(500, max(1, validated_params['max_evals_per_comp']))
                logger.warning("Max evaluations limited to 1-500 in server mode")
        
        # Retro depth (if applicable)
        if 'retro_depth' in validated_params:
            if validated_params['retro_depth'] > 2 or validated_params['retro_depth'] < 1:
                validated_params['retro_depth'] = min(2, max(1, validated_params['retro_depth']))
                logger.warning("Retro depth limited to 1-2 in server mode")
    
    return validated_params


def create_molecule_healer(
    bb_supplier: str = 'test',
    reaction_tags: List[str] = None,
    sim_threshold: float = 0.15,
    max_bbs_per_comp: int = -1,
    max_evals_per_comp: Optional[int] = None,
    n_compositions: int = 10,
    verbose: int = 1,
    use_fragment_healer: bool = False
) -> Union[MoleculeHEALER, FragmentHEALER]:
    '''
        Create a MoleculeHEALER or FragmentHEALER instance with the specified parameters.
        
        Args:
            bb_supplier: Building block supplier ('test', 'US_stock', 'EU_stock', 'Global_stock')
            reaction_tags: List of reaction tags to filter reactions
            sim_threshold: Similarity threshold for filtering building blocks
            max_bbs_per_comp: Maximum number of building blocks per composition
            max_evals_per_comp: Maximum number of evaluations per composition
            n_compositions: Number of compositions to consider for enumeration
            verbose: Verbosity level (0-2)
            use_fragment_healer: Whether to use FragmentHEALER instead of MoleculeHEALER
            
        Returns:
            MoleculeHEALER or FragmentHEALER instance
    '''
    if reaction_tags is None:
        reaction_tags = ["amide coupling", "amide", "C-N bond formation", "C-N",
                        "alkylation", "N-arylation", "azole", "amination"]
    
    # Validate parameters for server mode
    params = {
        'reaction_tags': reaction_tags,
        'sim_threshold': sim_threshold,
        'max_bbs_per_comp': max_bbs_per_comp,
        'max_evals_per_comp': max_evals_per_comp,
        'n_compositions': n_compositions
    }
    healer_type = "fragment" if use_fragment_healer else "molecule"
    validated_params = validate_server_parameters(params, healer_type)
    
    # Use absolute path for building blocks
    bb_path = BB_PATHS.get(bb_supplier, bb_supplier)
    
    if use_fragment_healer:
        return FragmentHEALER(
            bb_supplier=bb_path,
            reaction_tags=validated_params['reaction_tags'],
            max_evals_per_comp=validated_params['max_evals_per_comp'],
            sim_threshold=validated_params['sim_threshold'],
            max_bbs_per_comp=validated_params['max_bbs_per_comp'],
            verbose=verbose
        )
    else:
        return MoleculeHEALER(
            bb_supplier=bb_path,
            reaction_tags=validated_params['reaction_tags'],
            max_evals_per_comp=validated_params['max_evals_per_comp'],
            sim_threshold=validated_params['sim_threshold'],
            max_bbs_per_comp=validated_params['max_bbs_per_comp'],
            verbose=verbose
        )


def create_site_healer(
    bb_supplier: str = 'test',
    reaction_tags: List[str] = None,
    rules: Dict[str, Tuple[int, int]] = None,
    struct_rules: List[str] = None,
    max_evals_per_comp: Optional[int] = None,
    verbose: int = 1
) -> SiteHEALER:
    '''
        Create a SiteHEALER instance with the specified parameters.
        
        Args:
            bb_supplier: Building block supplier ('test', 'US_stock', 'EU_stock', 'Global_stock')
            reaction_tags: List of reaction tags to filter reactions
            rules: Dictionary of molecular property rules
            struct_rules: List of structural rules (SMARTS patterns)
            max_evals_per_comp: Maximum number of evaluations per composition
            verbose: Verbosity level (0-2)
            
        Returns:
            SiteHEALER instance
    '''
    if reaction_tags is None:
        reaction_tags = ["amide coupling", "amide", "C-N bond formation", "C-N",
                        "alkylation", "N-arylation", "azole", "amination"]
    
    if rules is None:
        rules = {
            'MW': (0, 500),
            'HBD': (0, 5),
            'HBA': (0, 10),
            'TPSA': (0, 200),
            'RotB': (0, 10),
            'Rings': (0, 10),
            'ArRings': (0, 5),
            'Chiral': (0, 5),
        }
    
    if struct_rules is None:
        struct_rules = []
    
    # Validate parameters for server mode
    params = {
        'reaction_tags': reaction_tags,
        'max_evals_per_comp': max_evals_per_comp
    }
    validated_params = validate_server_parameters(params, "site")
    
    # Use absolute path for building blocks
    bb_path = BB_PATHS.get(bb_supplier, bb_supplier)
    
    return SiteHEALER(
        bb_supplier=bb_path,
        reaction_tags=validated_params['reaction_tags'],
        max_evals_per_comp=validated_params['max_evals_per_comp'],
        rules=rules,
        struct_rules=struct_rules,
        verbose=verbose
    )


def run_molecule_enumeration(
    molecule: str,
    bb_supplier: str,
    reaction_tags: List[str],
    custom_sites: Optional[List[Tuple[int, int]]] = None,
    sim_threshold: float = 0.15,
    n_compositions: int = 10,
    randomize_compositions: bool = False,
    random_seed: int = -1,
    retro_tree_depth: int = 1,
    min_frag_size: int = 3,
    max_bbs_per_comp: int = -1,
    max_evals_per_comp: Optional[int] = None,
    use_fragment_healer: bool = False
) -> List[Dict[str, Any]]:
    '''
        Run molecule enumeration and return results as a list of dictionaries.
        
        Args:
            molecule: SMILES string of the query molecule
            bb_supplier: Building block supplier
            reaction_tags: List of reaction tags
            custom_sites: Custom split sites as list of bond tuples
            sim_threshold: Similarity threshold
            n_compositions: Number of compositions to generate
            randomize_compositions: Whether to randomize compositions
            random_seed: Random seed for composition generation
            retro_tree_depth: Depth of retrosynthesis tree
            min_frag_size: Minimum fragment size
            max_bbs_per_comp: Maximum building blocks per composition
            max_evals_per_comp: Maximum evaluations per composition
            use_fragment_healer: Whether to use FragmentHEALER instead of MoleculeHEALER
            
        Returns:
            List of dictionaries containing enumeration results
    '''
    try:
        # Automatically detect if we need FragmentHEALER
        num_fragments = count_molecular_fragments(molecule)
        auto_use_fragment_healer = num_fragments > 1
        
        # Use FragmentHEALER if either requested or automatically detected
        final_use_fragment_healer = use_fragment_healer or auto_use_fragment_healer
        
        healer = create_molecule_healer(
            bb_supplier=bb_supplier,
            reaction_tags=reaction_tags,
            sim_threshold=sim_threshold,
            max_bbs_per_comp=max_bbs_per_comp,
            max_evals_per_comp=max_evals_per_comp,
            n_compositions=n_compositions,
            verbose=1,
            use_fragment_healer=final_use_fragment_healer
        )
        
        # Set the query molecule with parameters
        if final_use_fragment_healer:
            # FragmentHEALER uses different parameters
            healer.set_query_mol(query_mol=molecule)
        else:
            # MoleculeHEALER parameters
            healer.set_query_mol(
                query_mol=molecule,
                n_compositions=n_compositions,
                randomize_compositions=randomize_compositions,
                random_seed=random_seed,
                custom_split_sites=[custom_sites] if custom_sites else None,
                retro_tree_depth=retro_tree_depth,
                min_frag_size=min_frag_size
            )
        
        # Run enumeration
        healer.enumerate()
        
        # Get results with similarity calculation
        results = healer.get_results(as_dict=True, calc_similarity=True)
        
        return results
        
    except Exception as e:
        logger.error(f"Error in molecule enumeration: {str(e)}")
        raise


def run_site_enumeration(
    molecule: str,
    bb_supplier: str,
    reaction_tags: List[str],
    reactive_sites: Optional[List[int]] = None,
    rules: Dict[str, Tuple[int, int]] = None,
    struct_rules: List[str] = None,
    max_evals_per_comp: Optional[int] = None
) -> List[Dict[str, Any]]:
    '''
        Run site enumeration and return results as a list of dictionaries.
        
        Args:
            molecule: SMILES string of the query molecule
            bb_supplier: Building block supplier
            reaction_tags: List of reaction tags
            reactive_sites: List of reactive site atom indices
            rules: Dictionary of molecular property rules
            struct_rules: List of structural rules (SMARTS patterns)
            max_evals_per_comp: Maximum evaluations per composition
            
        Returns:
            List of dictionaries containing enumeration results
    '''
    try:
        healer = create_site_healer(
            bb_supplier=bb_supplier,
            reaction_tags=reaction_tags,
            rules=rules,
            struct_rules=struct_rules,
            max_evals_per_comp=max_evals_per_comp,
            verbose=1
        )
        
        # Set the query molecule with reactive sites
        healer.set_query_mol(
            query_mol=molecule,
            reactive_sites=reactive_sites
        )
        
        # Run enumeration
        healer.enumerate()
        
        # Get results with similarity calculation
        results = healer.get_results(as_dict=True, calc_similarity=True)
        
        return results
        
    except Exception as e:
        logger.error(f"Error in site enumeration: {str(e)}")
        raise


def format_enumeration_results(results: List[Dict[str, Any]], app_type: str) -> List[Dict[str, Any]]:
    '''
        Format enumeration results for web app compatibility.
        
        Args:
            results: Raw enumeration results
            app_type: 'molecule' or 'site'
            
        Returns:
            Formatted results compatible with the web app
    '''
    formatted_results = []
    
    for result in results:
        formatted_result = {
            'Product': result.get('Product', ''),
            'Similarity_to_query': result.get('Similarity_to_query', 0.0)
        }
        
        # Extract building blocks
        bb_keys = [k for k in result.keys() if k.startswith('BB')]
        bb_keys.sort(key=lambda x: int(x[2:]))  # Sort BB1, BB2, etc.
        
        if app_type == 'molecule':
            # For molecule healer, we expect multiple BBs
            for i, bb_key in enumerate(bb_keys, 1):
                if result.get(bb_key):
                    formatted_result[f'BB{i}'] = result[bb_key]
        elif app_type == 'site':
            # For site healer, we typically have one BB
            if bb_keys and result.get(bb_keys[0]):
                formatted_result['BB'] = result[bb_keys[0]]
        
        # Extract reaction names
        rxn_keys = [k for k in result.keys() if k.startswith('Reaction') and k.endswith('_name')]
        if rxn_keys:
            rxn_keys.sort(key=lambda x: int(x.split('_')[0][8:]))  # Sort Reaction1_name, etc.
            reaction_names = [result.get(k, '') for k in rxn_keys if result.get(k)]
            if reaction_names:
                formatted_result['Reaction_name'] = ' -> '.join(reaction_names)
        
        formatted_results.append(formatted_result)
    
    return formatted_results


def generate_molecule_visualization(
    mol_smiles: str, 
    bb_smiles: List[str] = None, 
    legend: str = '',
    show_idx: bool = False
) -> str:
    '''
        Generate SVG visualization of a molecule with highlighted building blocks.
        
        Args:
            mol_smiles: SMILES string of the molecule
            bb_smiles: List of building block SMILES to highlight
            legend: Legend text for the image
            show_idx: Whether to show atom indices
            
        Returns:
            Base64-encoded SVG string
    '''
    try:
        if bb_smiles and len(bb_smiles) > 0:
            # Filter out empty building blocks
            bb_smiles = [bb for bb in bb_smiles if bb and bb.strip()]
            
            if bb_smiles:
                return utils.get_svg_mol_with_bbs(
                    mol=mol_smiles,
                    bbs=bb_smiles,
                    bb_colors=['purple', 'green', 'blue', 'orange', 'red'][:len(bb_smiles)],
                    legend=legend,
                    width=350,
                    height=150
                )
        
        return utils.get_svg_mol(
            mol=mol_smiles,
            legend=legend,
            show_idx=show_idx
        )
        
    except Exception as e:
        logger.error(f"Error generating molecule visualization: {str(e)}")
        raise ValueError(f"Cannot generate visualization: {str(e)}")
