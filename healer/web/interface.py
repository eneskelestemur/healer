'''
    Interface utilities for HEALER classes to standardize web app interactions.
    Adapted for the internal web package.
'''
import os
import logging
from pathlib import Path
from typing import List, Dict, Any, Optional, Union, Tuple

from rdkit import Chem

from healer.application.healer import MoleculeHEALER, SiteHEALER, FragmentHEALER


logger = logging.getLogger(__name__)

# Try to locate data directory
_env_data_dir = os.environ.get('HEALER_DATA_DIR')
if _env_data_dir:
    BB_BASE_PATH = Path(_env_data_dir) / 'buildingblocks'
    REACTIONS_PATH = Path(_env_data_dir) / 'reactions' / 'reactions.json'
else:
    HEALER_ROOT = Path(__file__).parent.parent.parent
    BB_BASE_PATH = HEALER_ROOT / 'data' / 'buildingblocks'
    REACTIONS_PATH = HEALER_ROOT / 'data' / 'reactions' / 'reactions.json'

BB_PATHS = {
    "US_stock": str(BB_BASE_PATH / "Enamine_Rush-Delivery_Building_Blocks-US" / "*_processed.sdf"),
    "EU_stock": str(BB_BASE_PATH / "Enamine_Rush-Delivery_Building_Blocks-EU" / "*_processed.sdf"),
    "Global_stock": str(BB_BASE_PATH / "Enamine_Building_Blocks_Stock" / "*_processed.sdf"),
    "test": str(BB_BASE_PATH / "test_100_bb_processed.sdf"),
}


def count_molecular_fragments(smiles: str) -> int:
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return 0
        fragments = Chem.GetMolFrags(mol, asMols=True)
        return len(fragments)
    except Exception as e:
        logger.error(f"Error counting fragments in {smiles}: {e}")
        return 0


def validate_server_parameters(params: Dict[str, Any], healer_type: str = "molecule") -> Dict[str, Any]:
    server_mode = os.environ.get('HEALER_SERVER_MODE', 'false').lower() == 'true'
    
    if not server_mode:
        return params
    
    validated_params = params.copy()
    
    if 'reaction_tags' in validated_params and validated_params['reaction_tags']:
        if 'all' in validated_params['reaction_tags']:
            validated_params['reaction_tags'] = [tag for tag in validated_params['reaction_tags'] if tag != 'all']
        
        if len(validated_params['reaction_tags']) > 15:
            validated_params['reaction_tags'] = validated_params['reaction_tags'][:15]
        
        if not validated_params['reaction_tags']:
            validated_params['reaction_tags'] = ["amide coupling", "amide", "C-N bond formation", "C-N",
                                               "alkylation", "N-arylation", "azole", "amination"]
    
    if healer_type in ["molecule", "fragment"]:
        if 'sim_threshold' in validated_params:
            if validated_params['sim_threshold'] < 0.3:
                validated_params['sim_threshold'] = 0.3
        
        if 'max_bbs_per_frag' in validated_params:
            if validated_params['max_bbs_per_frag'] > 10 or validated_params['max_bbs_per_frag'] < 1:
                validated_params['max_bbs_per_frag'] = min(10, max(1, validated_params['max_bbs_per_frag']))
        
        if 'n_compositions' in validated_params:
            if validated_params['n_compositions'] > 50:
                validated_params['n_compositions'] = 50
        
        if 'retro_depth' in validated_params:
            if validated_params['retro_depth'] > 2 or validated_params['retro_depth'] < 1:
                validated_params['retro_depth'] = min(2, max(1, validated_params['retro_depth']))

    # Common limits
    if 'max_evals_per_comp' in validated_params and validated_params['max_evals_per_comp']:
        if validated_params['max_evals_per_comp'] > 500 or validated_params['max_evals_per_comp'] < 1:
            validated_params['max_evals_per_comp'] = min(500, max(1, validated_params['max_evals_per_comp']))
            
    if 'max_products_per_comp' in validated_params and validated_params['max_products_per_comp']:
        if validated_params['max_products_per_comp'] > 100:
            validated_params['max_products_per_comp'] = 100

    if 'max_total_products' in validated_params and validated_params['max_total_products']:
        if validated_params['max_total_products'] > 1000:
            validated_params['max_total_products'] = 1000
    
    return validated_params


def create_molecule_healer(
    bb_source: str = 'test',
    reaction_tags: List[str] = None,
    sim_threshold: float = 0.15,
    max_bbs_per_frag: int = -1,
    max_evals_per_comp: Optional[int] = None,
    n_compositions: int = 10,
    verbose: int = 1,
    shuffle_bb_order: bool = False,
    use_fragment_healer: bool = False
) -> Union[MoleculeHEALER, FragmentHEALER]:
    
    if reaction_tags is None:
        reaction_tags = ["amide coupling", "amide", "C-N bond formation", "C-N",
                        "alkylation", "N-arylation", "azole", "amination"]
    
    params = {
        'reaction_tags': reaction_tags,
        'sim_threshold': sim_threshold,
        'max_bbs_per_frag': max_bbs_per_frag,
        'max_evals_per_comp': max_evals_per_comp,
        'n_compositions': n_compositions
    }
    healer_type = "fragment" if use_fragment_healer else "molecule"
    validated_params = validate_server_parameters(params, healer_type)
    
    bb_path = BB_PATHS.get(bb_source, bb_source)
    
    common_kwargs = {
        'bb_source': bb_path,
        'reaction_tags': validated_params['reaction_tags'],
        'shuffle_bb_order': shuffle_bb_order,
        'sim_threshold': validated_params['sim_threshold'],
        'max_bbs_per_frag': validated_params['max_bbs_per_frag'],
        'verbose': verbose
    }
    
    if use_fragment_healer:
        return FragmentHEALER(**common_kwargs)
    else:
        return MoleculeHEALER(**common_kwargs)


def create_site_healer(
    bb_source: str = 'test',
    reaction_tags: List[str] = None,
    rules: Dict[str, Tuple[int, int]] = None,
    struct_rules: List[str] = None,
    max_evals_per_comp: Optional[int] = None,
    verbose: int = 1,
    shuffle_bb_order: bool = False
) -> SiteHEALER:
    
    if reaction_tags is None:
        reaction_tags = ["amide coupling", "amide", "C-N bond formation", "C-N",
                        "alkylation", "N-arylation", "azole", "amination"]
    
    if rules is None:
        rules = {
            'MW': (0, 500), 'HBD': (0, 5), 'HBA': (0, 10), 'TPSA': (0, 200),
            'RotB': (0, 10), 'Rings': (0, 10), 'ArRings': (0, 5), 'Chiral': (0, 5),
        }
    
    if struct_rules is None:
        struct_rules = []
    
    params = {
        'reaction_tags': reaction_tags,
        'max_evals_per_comp': max_evals_per_comp
    }
    validated_params = validate_server_parameters(params, "site")
    bb_path = BB_PATHS.get(bb_source, bb_source)
    
    return SiteHEALER(
        bb_source=bb_path,
        reaction_tags=validated_params['reaction_tags'],
        rules=rules,
        struct_rules=struct_rules,
        shuffle_bb_order=shuffle_bb_order,
        verbose=verbose
    )


def run_molecule_enumeration(
    molecule: str,
    bb_source: str,
    reaction_tags: List[str],
    custom_sites: Optional[List[Tuple[int, int]]] = None,
    sim_threshold: float = 0.15,
    n_compositions: int = 10,
    randomize_compositions: bool = False,
    random_seed: int = -1,
    retro_tree_depth: int = 1,
    min_frag_size: int = 3,
    max_bbs_per_frag: int = -1,
    shuffle_bb_order: bool = False,
    max_evals_per_comp: Optional[int] = None,
    max_products_per_comp: Optional[int] = None,
    max_total_products: Optional[int] = None,
    use_fragment_healer: bool = False
) -> List[Dict[str, Any]]:
    
    try:
        num_fragments = count_molecular_fragments(molecule)
        auto_use_fragment_healer = num_fragments > 1
        final_use_fragment_healer = use_fragment_healer or auto_use_fragment_healer
        
        healer = create_molecule_healer(
            bb_source=bb_source,
            reaction_tags=reaction_tags,
            sim_threshold=sim_threshold,
            max_bbs_per_frag=max_bbs_per_frag,
            max_evals_per_comp=max_evals_per_comp,
            n_compositions=n_compositions,
            verbose=1,
            shuffle_bb_order=shuffle_bb_order,
            use_fragment_healer=final_use_fragment_healer
        )
        
        if final_use_fragment_healer:
            healer.set_query_mol(query_mol=molecule)
        else:
            healer.set_query_mol(
                query_mol=molecule,
                n_compositions=n_compositions,
                randomize_compositions=randomize_compositions,
                random_seed=random_seed,
                custom_split_sites=[custom_sites] if custom_sites else None,
                retro_tree_depth=retro_tree_depth,
                min_frag_size=min_frag_size
            )
        
        # Validate limit parameters
        limits = validate_server_parameters({
            'max_evals_per_comp': max_evals_per_comp,
            'max_products_per_comp': max_products_per_comp,
            'max_total_products': max_total_products
        }, "molecule")

        healer.enumerate(
            max_evals_per_comp=limits.get('max_evals_per_comp'),
            max_products_per_comp=limits.get('max_products_per_comp'),
            max_total_products=limits.get('max_total_products')
        )
        return healer.get_results(as_dict=True, calc_similarity=True, calc_properties=True)
        
    except Exception as e:
        logger.error(f"Error in molecule enumeration: {str(e)}")
        raise


def run_site_enumeration(
    molecule: str,
    bb_source: str,
    reaction_tags: List[str],
    reactive_sites: Optional[List[int]] = None,
    rules: Dict[str, Tuple[int, int]] = None,
    struct_rules: List[str] = None,
    shuffle_bb_order: bool = False,
    max_evals_per_comp: Optional[int] = None,
    max_products_per_comp: Optional[int] = None,
    max_total_products: Optional[int] = None
) -> List[Dict[str, Any]]:
    
    try:
        healer = create_site_healer(
            bb_source=bb_source,
            reaction_tags=reaction_tags,
            rules=rules,
            struct_rules=struct_rules,
            shuffle_bb_order=shuffle_bb_order,
            max_evals_per_comp=max_evals_per_comp,
            verbose=1
        )
        
        healer.set_query_mol(
            query_mol=molecule,
            reactive_sites=reactive_sites
        )
        
        # Validate limit parameters
        limits = validate_server_parameters({
            'max_evals_per_comp': max_evals_per_comp,
            'max_products_per_comp': max_products_per_comp,
            'max_total_products': max_total_products
        }, "site")

        healer.enumerate(
            max_evals_per_comp=limits.get('max_evals_per_comp'),
            max_products_per_comp=limits.get('max_products_per_comp'),
            max_total_products=limits.get('max_total_products')
        )
        return healer.get_results(as_dict=True, calc_similarity=True, calc_properties=True)
        
    except Exception as e:
        logger.error(f"Error in site enumeration: {str(e)}")
        raise

def format_enumeration_results(results: List[Dict[str, Any]], app_type: str) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    display_results = []
    complete_results = []
    
    for result in results:
        complete_result = result.copy()
        complete_results.append(complete_result)
        
        display_result = {
            'Product': result.get('Product', ''),
            'Similarity_to_query': result.get('Similarity_to_query', 0.0),
            'QED': result.get('QED', 0.0)
        }
        
        if 'stoplight_color' in result:
            display_result['stoplight_color'] = result['stoplight_color']
        
        bb_keys = [k for k in result.keys() if k.startswith('BB') and not k.endswith('_url') and not k.endswith('_id')]
        bb_keys.sort(key=lambda x: int(x[2:]))
        
        if app_type == 'molecule':
            for i, bb_key in enumerate(bb_keys, 1):
                if result.get(bb_key):
                    display_result[f'BB{i}'] = result[bb_key]
        elif app_type == 'site':
            if bb_keys and result.get(bb_keys[0]):
                display_result['BB'] = result[bb_keys[0]]
        
        rxn_keys = [k for k in result.keys() if k.startswith('Reaction') and k.endswith('_name')]
        if rxn_keys:
            rxn_keys.sort(key=lambda x: int(x.split('_')[0][8:]))
            reaction_names = [result.get(k, '') for k in rxn_keys if result.get(k)]
            if reaction_names:
                display_result['Reaction_name'] = ' -> '.join(reaction_names)
        
        display_results.append(display_result)
    
    return display_results, complete_results