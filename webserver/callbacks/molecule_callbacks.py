'''
    Callback functions specific to the MoleculeHEALER interface.
'''
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

from typing import List, Tuple, Optional
from dash import Input, Output, State, callback_context, no_update
from webserver.utils.healer_interface import (
    run_molecule_enumeration, 
    format_enumeration_results, 
    count_molecular_fragments,
)


def register_molecule_callbacks(app, app_id: str = "molecule"):
    '''
        Register callbacks specific to the MoleculeHEALER interface.
        
        Args:
            app: Dash application instance
            app_id: Application ID for component naming
    '''
    
    @app.callback(
        [Output(f"{app_id}-random-seed-label", "style"),
         Output(f"{app_id}-random-seed-input", "style"),
         Output(f"{app_id}-random-seed-input", "disabled")],
        [Input(f"{app_id}-randomize-checkbox", "value")]
    )
    def toggle_random_seed_visibility(randomize_value):
        '''Show/hide random seed input and label based on randomize checkbox.'''
        if "randomize" in randomize_value:
            return (
                {'margin-left': 'auto', 'display': 'inline-block', 'vertical-align': 'top', 'inlineSize': '2em', 'float': 'right'},
                {'fontSize': '14px', 'display': 'inline-block', 'width': '100px', 'marginLeft': '3em',
                 'vertical-align': 'top', 'float': 'right'},
                False
            )
        else:
            return (
                {"fontSize": "14px", 'margin-right': '5px', 'display': 'none', 'vertical-align': 'top'},
                {'fontSize': '14px', 'display': 'none', 'vertical-align': 'top'},
                True
            )
    
    @app.callback(
        [Output(f"{app_id}-fragment-alert", "is_open"),
         Output(f"{app_id}-fragment-alert", "style")],
        [Input(f"{app_id}-molecule-input", "value")]
    )
    def show_fragment_alert(molecule):
        '''Show fragment healer alert when multi-component molecule is detected.'''
        if not molecule:
            return False, {'display': 'none'}
        
        num_fragments = count_molecular_fragments(molecule)
        if num_fragments > 1:
            return True, {'display': 'block'}
        else:
            return False, {'display': 'none'}
    
    @app.callback(
        [Output(f"{app_id}-enumeration-svg", "children", True),
         Output(f"{app_id}-enumeration-svg", "style", True),
         Output(f"{app_id}-enumeration-results-slider", "max"),
         Output(f"{app_id}-enumeration-results-slider", "marks"),
         Output(f"{app_id}-enumeration-store", "data")],
        [Input(f"{app_id}-enumerate-button", "n_clicks")],
        [State(f"{app_id}-molecule-input", "value"),
         State(f"{app_id}-bb-source-select", "value"),
         State(f"{app_id}-reaction-tags-dropdown", "value"),
         State(f"{app_id}-custom-sites-input", "value"),
         State(f"{app_id}-sim-threshold-slider", "value"),
         State(f"{app_id}-n-compositions-slider", "value"),
         State(f"{app_id}-randomize-checkbox", "value"),
         State(f"{app_id}-random-seed-input", "value"),
         State(f"{app_id}-retro-depth-slider", "value"),
         State(f"{app_id}-min-frag-size-slider", "value"),
         State(f"{app_id}-max-bbs-input", "value"),
         State(f"{app_id}-max-evals-input", "value")],
        background=True,
        running=[
            (Output("tabs", "disabled"), True, False),
            (Output(f"{app_id}-update-button", "disabled"), True, False),
            (Output(f"{app_id}-enumerate-button", "disabled"), True, False),
            (Output(f"{app_id}-enumeration-results-slider", "disabled"), True, False),
            (Output(f"{app_id}-save-button", "disabled"), True, False),
            (Output(f"{app_id}-cancel-button", "disabled"), False, True),
        ],
        cancel=[Input(f"{app_id}-cancel-button", "n_clicks")],
        prevent_initial_call=True
    )
    def run_molecule_enumeration_callback(
        n_clicks, molecule, bb_source, reaction_tags, custom_sites_str,
        sim_threshold, n_compositions, randomize_checkbox, random_seed,
        retro_depth, min_frag_size, max_bbs, max_evals
    ):
        '''Run molecule enumeration with the specified parameters.'''
        
        if not n_clicks or not molecule:
            return no_update, no_update, no_update, no_update, no_update
        
        try:
            # Parse custom split sites
            custom_sites = None
            if custom_sites_str and custom_sites_str.strip():
                try:
                    # Parse format like "9-10, 6-9; 12-13" -> [[(9,10), (6,9)], [(12,13)]]
                    # Semicolons separate different split groups
                    # Commas separate atom pairs within the same group
                    # Hyphens separate the two atoms in a pair
                    split_groups = []
                    for group_str in custom_sites_str.split(';'):
                        if group_str.strip():
                            group_pairs = []
                            for pair_str in group_str.split(','):
                                if pair_str.strip():
                                    # Parse atom pair (e.g., "9-10" -> (9,10))
                                    parts = [int(x.strip()) for x in pair_str.split('-')]
                                    if len(parts) == 2:
                                        group_pairs.append(tuple(parts))
                            if group_pairs:
                                split_groups.append(group_pairs)
                    custom_sites = split_groups if split_groups else None
                except Exception:
                    custom_sites = None
            
            # Parse randomization settings
            randomize_compositions = "randomize" in randomize_checkbox
            if not randomize_compositions:
                random_seed = -1
            elif random_seed is None:
                random_seed = -1
            
            # Parse max evaluations
            max_evaluations = None if max_evals is None or max_evals == "" else int(max_evals)
            
            # Parse max BBs per composition
            max_bbs_per_comp = -1 if max_bbs is None or max_bbs == "" else int(max_bbs)
            
            # Ensure reaction_tags is a list
            if not reaction_tags:
                reaction_tags = ["amide coupling", "amide", "C-N bond formation", "C-N",
                               "alkylation", "N-arylation", "azole", "amination"]
            elif isinstance(reaction_tags, str):
                reaction_tags = [reaction_tags]
            
            # Check if molecule has multiple fragments
            num_fragments = count_molecular_fragments(molecule)
            use_fragment_healer = num_fragments > 1
            
            # Run enumeration
            results = run_molecule_enumeration(
                molecule=molecule,
                bb_supplier=bb_source,
                reaction_tags=reaction_tags,
                custom_sites=custom_sites,
                sim_threshold=sim_threshold,
                n_compositions=n_compositions,
                randomize_compositions=randomize_compositions,
                random_seed=random_seed,
                retro_tree_depth=retro_depth,
                min_frag_size=min_frag_size,
                max_bbs_per_comp=max_bbs_per_comp,
                max_evals_per_comp=max_evaluations,
                use_fragment_healer=use_fragment_healer
            )
            
            # Format results for web app
            display_results, complete_results = format_enumeration_results(results, 'molecule')
            n_results = len(display_results)
            
            # Store both versions: display for UI, complete for download
            stored_data = {
                'display_results': display_results,
                'complete_results': complete_results,
                'app_type': 'molecule'
            }
            
            # Generate response based on results
            if n_results < 2:  # Only query molecule returned
                enum_img = "No Enumerations! Check Your Inputs!"
                bg_color = 'firebrick'
                slider_max = 0
                slider_marks = {0: "0"}
            else:
                enum_img = f"Slide to view {n_results} enumerated molecules!"
                bg_color = 'forestgreen'
                slider_max = n_results - 1
                # Create marks for slider
                if n_results <= 10:
                    slider_marks = {i: str(i) for i in range(n_results)}
                else:
                    step = max(1, n_results // 10)
                    slider_marks = {i: str(i) for i in range(0, n_results, step)}
                    if (n_results - 1) not in slider_marks:
                        slider_marks[n_results - 1] = str(n_results - 1)
            
            enum_style = {
                'backgroundColor': bg_color,
                'width': '80%',
                'height': '150px',
                'text-align': 'center',
                'fontSize': '20px',
                'fontFamily': 'sans-serif',
                'borderRadius': '20px',
                'alignItems': 'center',
                'justifyContent': 'center',
                'display': 'flex',
                'margin': '0 auto',
                'color': 'white'
            }
            
            return enum_img, enum_style, slider_max, slider_marks, stored_data
            
        except Exception as e:
            # Handle errors
            error_msg = f"Enumeration failed: {str(e)}"
            error_style = {
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
            
            return error_msg, error_style, 0, {0: "0"}, []
