'''
    Callback functions specific to the SiteHEALER interface.
'''
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

from typing import List, Dict, Tuple, Optional
from dash import Input, Output, State, callback_context, no_update
from webserver.utils.healer_interface import run_site_enumeration, format_enumeration_results


def register_site_callbacks(app, app_id: str = "site"):
    '''
        Register callbacks specific to the SiteHEALER interface.
        
        Args:
            app: Dash application instance
            app_id: Application ID for component naming
    '''
    
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
         State(f"{app_id}-reactive-sites-input", "value"),
         State(f"{app_id}-struct-rules-input", "value"),
         State(f"{app_id}-MW-slider", "value"),
         State(f"{app_id}-TPSA-slider", "value"),
         State(f"{app_id}-HBD-slider", "value"),
         State(f"{app_id}-HBA-slider", "value"),
         State(f"{app_id}-Rings-slider", "value"),
         State(f"{app_id}-ArRings-slider", "value"),
         State(f"{app_id}-RotB-slider", "value"),
         State(f"{app_id}-Chiral-slider", "value"),
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
    def run_site_enumeration_callback(
        n_clicks, molecule, bb_source, reaction_tags, reactive_sites_str,
        struct_rules_str, mw_range, tpsa_range, hbd_range, hba_range,
        rings_range, ar_rings_range, rotb_range, chiral_range, max_evals
    ):
        '''Run site enumeration with the specified parameters.'''
        
        if not n_clicks or not molecule:
            return no_update, no_update, no_update, no_update, no_update
        
        try:
            # Parse reactive sites
            reactive_sites = None
            if reactive_sites_str and reactive_sites_str.strip():
                try:
                    reactive_sites = [int(x.strip()) for x in reactive_sites_str.split(',') if x.strip()]
                    if not reactive_sites:
                        reactive_sites = None
                except Exception:
                    reactive_sites = None
            
            # Parse structural rules
            struct_rules = []
            if struct_rules_str and struct_rules_str.strip():
                struct_rules = [s.strip() for s in struct_rules_str.split() if s.strip()]
            
            # Parse molecular property rules
            rules = {
                'MW': tuple(mw_range) if mw_range else (0, 500),
                'TPSA': tuple(tpsa_range) if tpsa_range else (0, 200),
                'HBD': tuple(hbd_range) if hbd_range else (0, 5),
                'HBA': tuple(hba_range) if hba_range else (0, 10),
                'Rings': tuple(rings_range) if rings_range else (0, 10),
                'ArRings': tuple(ar_rings_range) if ar_rings_range else (0, 5),
                'RotB': tuple(rotb_range) if rotb_range else (0, 10),
                'Chiral': tuple(chiral_range) if chiral_range else (0, 5),
            }
            
            # Parse max evaluations
            max_evaluations = None if max_evals is None or max_evals == "" else int(max_evals)
            
            # Ensure reaction_tags is a list
            if not reaction_tags:
                reaction_tags = ["amide coupling", "amide", "C-N bond formation", "C-N",
                               "alkylation", "N-arylation", "azole", "amination"]
            elif isinstance(reaction_tags, str):
                reaction_tags = [reaction_tags]
            
            # Run enumeration
            results = run_site_enumeration(
                molecule=molecule,
                bb_source=bb_source,
                reaction_tags=reaction_tags,
                reactive_sites=reactive_sites,
                rules=rules,
                struct_rules=struct_rules,
                max_evals_per_comp=max_evaluations
            )
            
            # Format results for web app
            display_results, complete_results = format_enumeration_results(results, 'site')
            n_results = len(display_results)
            
            # Store both versions: display for UI, complete for download
            stored_data = {
                'display_results': display_results,
                'complete_results': complete_results,
                'app_type': 'site'
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
