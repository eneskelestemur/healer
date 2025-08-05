'''
    Shared callback functions used by both molecule and site healers.
'''
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

import pandas as pd
from pathlib import Path
from dash import Input, Output, State, callback_context, html, dcc, clientside_callback, ClientsideFunction
import healer.utils.utils as utils
from webserver.utils.healer_interface import generate_molecule_visualization

# Get absolute path to reactions file
HEALER_ROOT = Path(__file__).parent.parent.parent / 'healer'
REACTIONS_PATH = HEALER_ROOT / 'data' / 'reactions' / 'reactions.json'

# Check server mode
SERVER_MODE = os.environ.get('HEALER_SERVER_MODE', 'false').lower() == 'true'


def register_shared_callbacks(app, app_id: str):
    '''
        Register shared callbacks for both molecule and site healers.
        
        Args:
            app: Dash application instance
            app_id: Application ID ('molecule' or 'site')
    '''
    
    @app.callback(
        [Output(f"{app_id}-reaction-tags-dropdown", "options"),
         Output(f"{app_id}-reaction-tags-dropdown", "value")],
        [Input(f"{app_id}-reaction-tags-dropdown", "id")],
        [State(f"{app_id}-reaction-tags-dropdown", "value")]
    )
    def update_reaction_tags_options_and_values(_, current_value):
        '''Load reaction tags from the reaction data and preserve default values.'''
        try:
            reactions = utils.load_reactions_from_json(str(REACTIONS_PATH))
            reaction_tags = list(set(tag for r in reactions if r.is_valid() for tag in r.tags))
            reaction_tags.sort(key=lambda x: x.lower())
            reaction_tags.insert(0, 'all')
            options = [{"label": tag, "value": tag} for tag in reaction_tags]
            
            # If no current value is set, use default values
            if not current_value:
                default_value = ["amide coupling", "amide", "C-N bond formation", "C-N",
                               "alkylation", "N-arylation", "azole", "amination"]
                return options, default_value
            else:
                # Preserve current value
                return options, current_value
                
        except Exception:
            # Fallback options if loading fails
            default_tags = ["amide coupling", "amide", "C-N bond formation", "C-N",
                           "alkylation", "N-arylation", "azole", "amination"]
            options = [{"label": tag, "value": tag} for tag in default_tags]
            return options, default_tags
    
    @app.callback(
        [Output(f"{app_id}-molecule-svg", "children"),
         Output(f"{app_id}-molecule-svg", "style")],
        [Input(f"{app_id}-update-button", "n_clicks")],
        [State(f"{app_id}-molecule-input", "value")]
    )
    def update_molecule_display(n_clicks, molecule):
        '''Update the molecule visualization when the update button is clicked.'''
        ctx = callback_context
        if not ctx.triggered:
            mol_img = 'Click Update to View the Molecule'
            mol_style = {
                'width': '75%', 
                'height': '175px', 
                'text-align': 'center', 
                'fontSize': '20px', 
                'fontFamily': 'sans-serif', 
                'alignItems': 'center', 
                'justifyContent': 'center', 
                'display': 'flex', 
                'margin': '0 auto'
            }
        else:
            try:
                img = generate_molecule_visualization(molecule, show_idx=True)
                mol_img = html.Img(
                    src=img, 
                    style={'height': '100%', 'width': '80%', 'margin': '0 auto'}
                )
                mol_style = {'height': '175px', 'text-align': 'center'}
            except Exception as e:
                mol_img = f'Invalid Molecule: {str(e)}'
                mol_style = {
                    'backgroundColor': 'firebrick', 
                    'width': '80%', 
                    'height': '150px', 
                    'text-align': 'center', 
                    'fontSize': '18px', 
                    'fontFamily': 'sans-serif', 
                    'borderRadius': '20px', 
                    'alignItems': 'center', 
                    'justifyContent': 'center', 
                    'display': 'flex', 
                    'margin': '0 auto',
                    'color': 'white'
                }
        return mol_img, mol_style
    
    @app.callback(
        Output(f"{app_id}-download-enumerations", "data"),
        [Input(f"{app_id}-save-button", "n_clicks")],
        [State(f"{app_id}-enumeration-store", "data"),
         State(f"{app_id}-save-as-input", "value")],
        prevent_initial_call=True
    )
    def save_enumeration_results(n_clicks, enum_data, filename):
        '''Save enumeration results to CSV file.'''
        if enum_data and filename:
            # Use complete results for download to preserve all information
            complete_results = enum_data.get('complete_results', [])
            if complete_results:
                df = pd.DataFrame.from_records(complete_results)
                return dcc.send_data_frame(df.to_csv, filename=filename, index=False)
        return None
    
    @app.callback(
        [Output(f"{app_id}-enumeration-svg", "children", True),
         Output(f"{app_id}-enumeration-svg", "style", True)],
        [Input(f"{app_id}-enumeration-results-slider", "value")],
        [State(f"{app_id}-enumeration-store", "data")],
        prevent_initial_call=True
    )
    def update_enumeration_display(slider_val, enum_data):
        '''Update the enumeration result visualization based on slider value.'''
        if not enum_data:
            return "No data available", {'text-align': 'center'}
        
        # Use display results for visualization
        display_results = enum_data.get('display_results', [])
        if not display_results or slider_val >= len(display_results):
            return "No data available", {'text-align': 'center'}
        
        try:
            result = display_results[slider_val]
            query_mol = display_results[0]['Product']  # First result is always the query
            
            # Extract building blocks for visualization
            bb_smiles = []
            app_type = enum_data.get('app_type', app_id)  # Use stored app_type or fallback to app_id
            if app_type == "molecule":
                # For molecule healer, look for BB1, BB2, etc.
                bb_keys = [k for k in result.keys() if k.startswith('BB') and k[2:].isdigit()]
                bb_keys.sort(key=lambda x: int(x[2:]))
                bb_smiles = [result.get(k, '') for k in bb_keys if result.get(k)]
            elif app_type == "site":
                # For site healer, look for BB and include query molecule
                if result.get('BB'):
                    bb_smiles = [result['BB'], query_mol]
            
            # Generate legend
            similarity = result.get('Similarity_to_query', 0.0)
            reaction_name = result.get('Reaction_name', 'Unknown')
            legend = f"Similarity: {similarity:.3f}\n{reaction_name}"
            
            # Generate visualization
            img = generate_molecule_visualization(
                mol_smiles=result['Product'],
                bb_smiles=bb_smiles,
                legend=legend
            )
            
            enum_img = html.Img(
                src=img,
                style={'height': '100%', 'width': '80%', 'margin': '0 auto'}
            )
            enum_style = {'height': '200px', 'text-align': 'center'}
            
            return enum_img, enum_style
            
        except Exception as e:
            error_msg = f"Error displaying result: {str(e)}"
            error_style = {
                'backgroundColor': 'lightcoral',
                'width': '80%',
                'height': '150px',
                'text-align': 'center',
                'fontSize': '16px',
                'fontFamily': 'sans-serif',
                'borderRadius': '10px',
                'alignItems': 'center',
                'justifyContent': 'center',
                'display': 'flex',
                'margin': '0 auto',
                'color': 'darkred'
            }
            return error_msg, error_style
    
    # Add server mode validation callbacks if in server mode
    if SERVER_MODE and app_id == "molecule":
        @app.callback(
            [Output(f"{app_id}-sim-threshold-slider", "value"),
             Output(f"{app_id}-n-compositions-slider", "value", allow_duplicate=True),
             Output(f"{app_id}-max-bbs-input", "value"),
             Output(f"{app_id}-retro-depth-slider", "value"),
             Output(f"{app_id}-reaction-tags-dropdown", "value", allow_duplicate=True)],
            [Input(f"{app_id}-sim-threshold-slider", "value"),
             Input(f"{app_id}-n-compositions-slider", "value"),
             Input(f"{app_id}-max-bbs-input", "value"),
             Input(f"{app_id}-retro-depth-slider", "value"),
             Input(f"{app_id}-reaction-tags-dropdown", "value")],
            prevent_initial_call=True
        )
        def validate_server_limits(sim_threshold, n_compositions, max_bbs, retro_depth, reaction_tags):
            """Enforce server mode limitations."""
            # Similarity threshold must be >= 0.3
            if sim_threshold < 0.3:
                sim_threshold = 0.3
            
            # N compositions must be <= 50
            if n_compositions > 50:
                n_compositions = 50
            
            # Max BBs per comp must be 1-10 (if not -1)
            if max_bbs != -1 and (max_bbs < 1 or max_bbs > 10):
                max_bbs = min(10, max(1, max_bbs))
            
            # Retro depth must be 1-2
            if retro_depth < 1 or retro_depth > 2:
                retro_depth = min(2, max(1, retro_depth))
            
            # Reaction tags: remove 'all' and limit to 15
            if reaction_tags:
                # Remove 'all' tag in server mode
                if 'all' in reaction_tags:
                    reaction_tags = [tag for tag in reaction_tags if tag != 'all']
                    
                # Limit to 15 tags
                if len(reaction_tags) > 15:
                    reaction_tags = reaction_tags[:15]
                    
                # Ensure we still have some tags after filtering
                if not reaction_tags:
                    reaction_tags = ["amide coupling", "amide", "C-N bond formation", "C-N",
                                   "alkylation", "N-arylation", "azole", "amination"]
            
            return sim_threshold, n_compositions, max_bbs, retro_depth, reaction_tags
    
    elif SERVER_MODE and app_id == "site":
        @app.callback(
            [Output(f"{app_id}-reaction-tags-dropdown", "value", allow_duplicate=True)],
            [Input(f"{app_id}-reaction-tags-dropdown", "value")],
            prevent_initial_call=True
        )
        def validate_site_server_limits(reaction_tags):
            """Enforce server mode limitations for site healer."""
            # Reaction tags: remove 'all' and limit to 15
            if reaction_tags:
                # Remove 'all' tag in server mode
                if 'all' in reaction_tags:
                    reaction_tags = [tag for tag in reaction_tags if tag != 'all']
                    
                # Limit to 15 tags
                if len(reaction_tags) > 15:
                    reaction_tags = reaction_tags[:15]
                    
                # Ensure we still have some tags after filtering
                if not reaction_tags:
                    reaction_tags = ["amide coupling", "amide", "C-N bond formation", "C-N",
                                   "alkylation", "N-arylation", "azole", "amination"]
            
            return [reaction_tags]
