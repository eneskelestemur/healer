'''
    Dash application for the enumerator project.
'''
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import dash_bootstrap_components as dbc
import diskcache
import pandas as pd
import healer.utils.utils as utils

from dash import dcc, html, Dash, Input, Output, State, DiskcacheManager, callback_context, no_update
from legacy.enumerator import MoleculeEnumerator, SiteEnumerator


# Diskcache for handling background jobs
cache = diskcache.Cache("./cache")
background_callback_manager = DiskcacheManager(cache)

# Initialize the Dash app
app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP], background_callback_manager=background_callback_manager)

# Shared Reaction Tags
reactions = utils.load_reactions_from_json('reactions/reactions.json')
reaction_tags = list(set(tag for r in reactions if r.is_valid() for tag in r.tags))


# Layout for Molecule and Site Enumerators
def build_layout(app_id):
    if app_id == "molecule":
        additional_inputs = [
            dbc.Row(
                [
                    dbc.Col(html.Label("Split Sites"), width=3),
                    dbc.Col(
                        dbc.Input(
                            type="text",
                            id=f"{app_id}-custom-sites-input",
                            placeholder="Sites to split molecule (e.g., 9,10; 6, 9)",
                            style={'fontSize': '14px'}
                        ),
                        width=9
                    )
                ],
                className="mb-3"
            ),
            dbc.Row(
                [
                    dbc.Col(html.Label("Sim Cutoff"), width=3),
                    dbc.Col(
                        dcc.Slider(
                            id=f"{app_id}-sim-cutoff-slider",
                            min=0, max=1, step=0.01, value=0.15,
                            marks={i / 10: str(i / 10) for i in range(11)},
                            tooltip={'placement': 'bottom', 'always_visible': False}
                        ),
                        width=9
                    )
                ],
                className="mb-3"
            )
        ]
    elif app_id == "site":
        additional_inputs = [
            dbc.Row(
                [
                    dbc.Col(html.Label("Reaction Sites"), width=3),
                    dbc.Col(dbc.Input(type="text", id=f"{app_id}-reaction-sites-input", 
                                        placeholder="Reacting atoms, e.g., 19, 20, 21",
                                        style={'fontSize': '14px'}), 
                            width=9)
                ],
                className="mb-3"
            ),
            dbc.Row(
                [
                    dbc.Col(html.Label("Structure Rules"), width=3),
                    dbc.Col(
                        dbc.Textarea(
                            id=f"{app_id}-struct-rules-input",
                            placeholder="Substructure patterns for BBs separated by whitespace, e.g., c1ccccc1 CC",
                            style={'fontSize': '14px'}
                        ),
                        width=9
                    )
                ],
                className="mb-3"
            ),
            dbc.Row(
                [
                    dbc.Col(html.Label('MW'), width=1),
                    dbc.Col(dcc.RangeSlider(id='MW-slider', min=0, max=1000, step=0.1, value=[0, 1000],
                                            marks={0: '0', 1000: '1000'},
                                            tooltip={'placement': 'bottom', 'always_visible': True}),
                            width=5),
                    dbc.Col(dcc.RangeSlider(id='TPSA-slider', min=0, max=200, step=0.1, value=[0, 200],
                                            marks={0: '0', 200: '200'},
                                            tooltip={'placement': 'bottom', 'always_visible': True}),
                            width=5),
                    dbc.Col(html.Label('TPSA'), width=1),
                ],
                className="mb-3"
            ),
            dbc.Row(
                [
                    dbc.Col(html.Label("HBD"), width=1),
                    dbc.Col(dcc.RangeSlider(id='HBD-slider', min=0, max=10, step=1, value=[0, 10],
                                            marks={0: '0', 10: '10'},
                                            tooltip={'placement': 'bottom', 'always_visible': True}),
                            width=5),
                    dbc.Col(dcc.RangeSlider(id='HBA-slider', min=0, max=10, step=1, value=[0, 10],
                                            marks={0: '0', 10: '10'},
                                            tooltip={'placement': 'bottom', 'always_visible': True}),
                            width=5),
                    dbc.Col(html.Label('HBA'), width=1),
                ],
                className="mb-3"
            ),
            dbc.Row(
                [
                    dbc.Col(html.Label('Rings'), width=1),
                    dbc.Col(dcc.RangeSlider(id='Rings-slider', min=0, max=20, step=1, value=[0, 20],
                                            marks={0: '0', 20: '20'},
                                            tooltip={'placement': 'bottom', 'always_visible': True}),
                            width=5),
                    dbc.Col(dcc.RangeSlider(id='ArRings-slider', min=0, max=10, step=1, value=[0, 10],
                                            marks={0: '0', 10: '10'},
                                            tooltip={'placement': 'bottom', 'always_visible': True}),
                            width=5),
                    dbc.Col(html.Label('Ar. Rings'), width=1),
                ],
                className="mb-3"
            ),
            dbc.Row(
                [
                    dbc.Col(html.Label('Rot. Bonds'), width=1),
                    dbc.Col(dcc.RangeSlider(id='RotBonds-slider', min=0, max=20, step=1, value=[0, 20],
                                            marks={0: '0', 20: '20'},
                                            tooltip={'placement': 'bottom', 'always_visible': True}),
                            width=5),
                    dbc.Col(dcc.RangeSlider(id='Chiral-slider', min=0, max=10, step=1, value=[0, 10],
                                            marks={0: '0', 10: '10'},
                                            tooltip={'placement': 'bottom', 'always_visible': True}),
                            width=5),
                    dbc.Col(html.Label('Chiral'), width=1),
                ],
                className="mb-3"
            ),
        ]
    else:
        additional_inputs = []

    return dbc.Container(
        [
            html.H2(f"{app_id.capitalize()} Enumerator", style={'text-align': 'center', 'margin': '10px 0px 10px 0px'}),
            dbc.Row(
                [
                    # Inputs Section
                    dbc.Col(
                        [
                            dbc.Row(
                                [
                                    dbc.Col(html.Label("Molecule"), width=3),
                                    dbc.Col(
                                        dbc.Input(
                                            type="text",
                                            id=f"{app_id}-molecule-input",
                                            placeholder="Enter molecule",
                                            value="CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C",
                                            style={'fontSize': '14px'}
                                        ),
                                        width=9
                                    )
                                ],
                                className="mb-3"
                            ),
                            dbc.Row(
                                [
                                    dbc.Col(html.Label("BB source"), width=3),
                                    dbc.Col(
                                        dbc.Select(
                                            id=f"{app_id}-bb-source-select",
                                            options=[
                                                {"label": "test", "value": "test"},
                                                {"label": "US_stock", "value": "US_stock"},
                                                {"label": "EU_stock", "value": "EU_stock"},
                                                {"label": "Global_stock", "value": "Global_stock"}
                                            ],
                                            value="test"
                                        ),
                                        width=9
                                    )
                                ],
                                className="mb-3"
                            ),
                            dbc.Row(
                                [
                                    dbc.Col(html.Label("Reaction Tags"), width=3),
                                    dbc.Col(
                                        dcc.Dropdown(
                                            id=f'{app_id}-reaction-tags-dropdown',
                                            options=[{"label": tag, "value": tag} for tag in reaction_tags],
                                            value=["amide coupling", "amide", "C-N bond formation", "C-N",
                                                   "alkylation", "N-arylation", "azole", "amination"],
                                            multi=True,
                                        ),
                                        width=9
                                    )
                                ],
                                className="mb-3"
                            ),
                            *additional_inputs,
                            dbc.Row(
                                [
                                    dbc.Col(html.Label("Save as"), width=3),
                                    dbc.Col(dbc.Input(type="text", id=f"{app_id}-save-as-input", placeholder="Filename", 
                                                    value="enumerated_molecules.csv", style={'fontSize': '14px'}),
                                            width=9)
                                ],
                                className="mb-3"
                            ),
                            dbc.Row(
                                [
                                    dbc.Col(dbc.Button("Update", id=f"{app_id}-update-button", color="primary", className="mr-2",
                                                       style={'width': '100%', 'text-align': 'center'}), width=3),
                                    dbc.Col(dbc.Button("Enumerate", id=f"{app_id}-enumerate-button", color="warning", className="mr-2",
                                                       style={'width': '100%', 'text-align': 'center'}), width=3),
                                    dbc.Col(dbc.Button("Save", id=f"{app_id}-save-button", color="success", className="mr-2",
                                                       style={'width': '100%', 'text-align': 'center'}), width=3),
                                ],
                                className="mb-3",
                                justify="end",
                            ),
                            dbc.Row(
                                [
                                    dbc.Col(dbc.Button("Cancel Job", id=f"{app_id}-cancel-button", color="danger",
                                                       style={'width': '100%', 'text-align': 'center'}), width=9),
                                ],
                                className="mb-3",
                                justify="end",
                            ),
                        ],
                        md=6,
                        style={'height': '100%', 'justify_content': 'center', 'margin': '-10px 0px 0px 0px'}
                    ),
                    # Outputs Section
                    dbc.Col(
                        [
                            html.H3("Query Molecule", style={'text-align': 'center'}),
                            html.Div(id=f'{app_id}-molecule-svg', style={'height': '150px', 'text-align': 'center'}),
                            html.H3("Enumeration Results", style={'text-align': 'center', 'margin': '20px 0px 20px 0px'}),
                            html.Div(
                                [dcc.Slider(
                                    id=f'{app_id}-enumeration-results-slider', min=0, max=0, step=1, value=0,
                                    marks={0: "0"},
                                    tooltip={'placement': 'bottom', 'always_visible': False}
                                )],
                                style={'width': '70%', 'justify-content': 'center', 'margin': '0 auto'}
                            ),
                            dcc.Loading(
                                [
                                    html.Div('Click Enumerate to View the Enumerations', id=f'{app_id}-enumeration-svg', 
                                             style={'width': '75%', 'height': '200px', 'text-align': 'center', 
                                                    'fontSize': '20px', 'fontFamily': 'sans-serif', 'alignItems': 'center', 
                                                    'justifyContent': 'center', 'display': 'flex', 'margin': '0 auto'}),
                                ],
                                id=f"{app_id}-enumeration-loading", delay_show=100, type='circle',
                                overlay_style={"visibility":"visible", "filter": "blur(2px)"}
                            ),
                            dcc.Download(id=f"{app_id}-download-enumerations"),
                            dcc.Store(id=f'{app_id}-enumeration-store')
                        ],
                        md=6,
                        style={'height': '100%', 'justify_content': 'center', 'margin': '20px 0px 0px 0px'}
                    ),
                ],
                style={'height': '100%', 'justify_content': 'center', 'margin': '40px 0px 0px 0px'}
            ),
        ],
        fluid=True,
        style={'width': '80%'}
    )

# App Layout with Tabs
app.layout = dbc.Container(
    [
        html.H1("Enumerator Dashboard", style={'text-align': 'center'}),
        html.Div(
            dcc.Tabs(
                id="tabs",
                value="molecule-tab",
                children=[
                    dcc.Tab(label="Molecule Enumerator", value="molecule-tab",
                            style={'text-align': 'center', 'justify_content': 'center'}),
                    dcc.Tab(label="Site Enumerator", value="site-tab",
                            style={'text-align': 'center', 'justify_content': 'center'}),
                ],
                style={'width': '50%', 'margin': '0 auto', 'fontSize': '25px'}  # Limit tab width and center it
            ),
            style={'text-align': 'center', 'justify_content': 'center'}
        ),
        html.Div(
            [
                html.Div(build_layout("molecule"), id="molecule-tab-content", style={'display': 'block'}),
                html.Div(build_layout("site"), id="site-tab-content", style={'display': 'none'})
            ]
        )
    ],
    fluid=True
)

@app.callback(
    [Output("molecule-tab-content", "style"),
     Output("site-tab-content", "style"),
     Output("molecule-tab-content", "children"),
     Output("site-tab-content", "children")],
    [Input("tabs", "value")]
)
def toggle_tab_content(tab_name):
    molecule_content = build_layout("molecule")
    site_content = build_layout("site")
    if tab_name == "molecule-tab":
        return {'display': 'block'}, {'display': 'none'}, molecule_content, no_update
    elif tab_name == "site-tab":
        return {'display': 'none'}, {'display': 'block'}, no_update, site_content


# Callbacks
def shared_callbacks(app_id):
    # Update Molecule Callback
    @app.callback(
        [Output(f"{app_id}-molecule-svg", "children"),
         Output(f"{app_id}-molecule-svg", "style")],
        [Input(f"{app_id}-update-button", "n_clicks"),],
        [State(f"{app_id}-molecule-input", "value")],
    )
    def update_molecule(n_clicks, molecule):
        ctx = callback_context
        if not ctx.triggered:
            mol_img = 'Click Update to View the Molecule'
            mol_style = {'width': '75%', 'height': '175px', 'text-align': 'center', 
                        'fontSize': '20px', 'fontFamily': 'sans-serif', 'alignItems': 'center', 
                        'justifyContent': 'center', 'display': 'flex', 'margin': '0 auto'}
        else:
            try:
                img = utils.get_svg_mol(molecule, show_idx=True)
                mol_img = html.Img(src=img, style={'height': '100%', 'width': '80%', 'margin': '0 auto'})
                mol_style = {'height': '175px', 'text-align': 'center'}
            except ValueError:
                mol_img = 'Invalid Molecule!'
                mol_style = {'backgroundColor': 'firebrick', 'width': '80%', 'height': '150px', 
                             'text-align': 'center', 'fontSize': '24px', 'fontFamily': 'sans-serif', 
                             'borderRadius': '20px', 'alignItems': 'center', 'justifyContent': 'center', 
                             'display': 'flex', 'margin': '0 auto'}
        return mol_img, mol_style
    
    # Save Enumerations Callback
    @app.callback(
        Output(f"{app_id}-download-enumerations", "data"),
        [Input(f"{app_id}-save-button", "n_clicks")],
        [State(f"{app_id}-enumeration-store", "data"),
         State(f"{app_id}-save-as-input", "value")],
        prevent_initial_call=True
    )
    def save_enumerations(n_clicks, enum_mols, filename):
        save_out = dcc.send_data_frame(pd.DataFrame.from_records(enum_mols).to_csv, filename=filename, index=False)
        return save_out
    
    # Results Slider Callback
    @app.callback(
        [Output(f"{app_id}-enumeration-svg", "children", True),
         Output(f"{app_id}-enumeration-svg", "style", True)],
        [Input(f"{app_id}-enumeration-results-slider", "value"),],
        [State(f"{app_id}-enumeration-store", "data")],
        prevent_initial_call=True
    )
    def slide_enumeration(slider_val, enum_mols):
        molecule = enum_mols[0]['Product']
        enum = enum_mols[slider_val]
        if app_id == "molecule":
            img = utils.get_svg_mol_with_bbs(enum['Product'], enum['BB1'], enum['BB2'],
                                             bb_colors=['purple', 'green'],
                                             legend=f"T: {enum['Similarity_to_query']:.2f}\n{enum['Reaction_name']}",)
            enum_img = html.Img(src=img,
                                style={'height': '100%', 'width': '80%', 'margin': '0 auto'})
            enum_style = {'height': '200px', 'text-align': 'center'}
        elif app_id == "site":
            img = utils.get_svg_mol_with_bbs(enum['Product'], enum['BB'], molecule,
                                             bb_colors=['purple', 'green'],
                                             legend=f"T: {enum['Similarity_to_query']:.2f}\n{enum['Reaction_name']}")
            enum_img = html.Img(src=img,
                                style={'height': '100%', 'width': '80%', 'margin': '0 auto'})
            enum_style = {'height': '200px', 'text-align': 'center'}

        return enum_img, enum_style

    # Enumeration Callback
    if app_id == "molecule":
        @app.callback(
            [Output(f"{app_id}-enumeration-svg", "children", True),
             Output(f"{app_id}-enumeration-svg", "style", True),
             Output(f"{app_id}-enumeration-results-slider", "max"),
             Output(f"{app_id}-enumeration-results-slider", "marks"),
             Output(f"{app_id}-enumeration-store", "data")],
            [Input(f"{app_id}-enumerate-button", "n_clicks"),],
            [State(f"{app_id}-molecule-input", "value"),
             State(f"{app_id}-bb-source-select", "value"),
             State(f"{app_id}-reaction-tags-dropdown", "value"),
             State(f"{app_id}-custom-sites-input", "value"),
             State(f"{app_id}-sim-cutoff-slider", "value")],
            background=True,
            running=[
                (Output("tabs", "disabled"), True, False),
                (Output(f"{app_id}-update-button", "disabled"), True, False),
                (Output(f"{app_id}-enumerate-button", "disabled"), True, False),
                (Output(f"{app_id}-enumeration-results-slider", "disabled"), True, False),
                (Output(f"{app_id}-save-button", "disabled"), True, False),
                (Output(f"{app_id}-cancel-button", "disabled"), False, True),
            ],
            cancel=[
                Input(f"{app_id}-cancel-button", "n_clicks")
            ],
            prevent_initial_call=True
        )
        def enumerate_molecule(n_clicks, molecule, bb_source, tags, custom_sites, sim_cutoff):
            custom_sites = [tuple(map(int, site.replace(' ', '').split(',')))
                            for site in custom_sites.split(';') if site] if custom_sites else []
            enumerator = MoleculeEnumerator(molecule, bb_source, tags, custom_sites, sim_threshold=sim_cutoff)
            enumerator.enumerate()
            enum_mols = enumerator.get_results(as_dict=True)
            n = len(enum_mols)

            # No enumeration returns the query molecule, so n = 1
            if n < 2: 
                enum_img = "No Enumerations! Check Your Inputs!"
                bg_color = 'firebrick'
            else:
                enum_img = f"Slide to view {n} enumerated molecules!"
                bg_color = 'forestgreen'
            enum_style = {'backgroundColor': bg_color, 'width': '80%', 'height': '150px', 
                          'text-align': 'center', 'fontSize': '24px', 'fontFamily': 'sans-serif', 
                          'borderRadius': '20px', 'alignItems': 'center', 'justifyContent': 'center', 
                          'display': 'flex', 'margin': '0 auto'}
            slider_max = n - 1
            slider_marks = {i: str(i) for i in range(0, n, max(10, round(n//10, -1)))}
            
            return enum_img, enum_style, slider_max, slider_marks, enum_mols
        
    elif app_id == "site":
        @app.callback(
            [Output(f"{app_id}-enumeration-svg", "children", True),
             Output(f"{app_id}-enumeration-svg", "style", True),
             Output(f"{app_id}-enumeration-results-slider", "max"),
             Output(f"{app_id}-enumeration-results-slider", "marks"),
             Output(f"{app_id}-enumeration-store", "data")],
            [Input(f"{app_id}-enumerate-button", "n_clicks"),],
            [State(f"{app_id}-molecule-input", "value"),
             State(f"{app_id}-bb-source-select", "value"),
             State(f"{app_id}-reaction-tags-dropdown", "value"),
             State(f"{app_id}-reaction-sites-input", "value"),
             State(f"{app_id}-struct-rules-input", "value"),
             State("MW-slider", "value"),
             State("TPSA-slider", "value"),
             State("HBD-slider", "value"),
             State("HBA-slider", "value"),
             State("Rings-slider", "value"),
             State("ArRings-slider", "value"),
             State("RotBonds-slider", "value"),
             State("Chiral-slider", "value")],
            background=True,
            running=[
                (Output("tabs", "disabled"), True, False),
                (Output(f"{app_id}-update-button", "disabled"), True, False),
                (Output(f"{app_id}-enumerate-button", "disabled"), True, False),
                (Output(f"{app_id}-enumeration-results-slider", "disabled"), True, False),
                (Output(f"{app_id}-save-button", "disabled"), True, False),
                (Output(f"{app_id}-cancel-button", "disabled"), False, True),
            ],
            cancel=[
                Input(f"{app_id}-cancel-button", "n_clicks")
            ],
            prevent_initial_call=True
        )
        def enumerate_site(n_clicks, molecule, bb_source, tags, reaction_sites, struct_rules, 
                           mw, tpsa, hbd, hba, rings, ar_rings, rot_bonds, chiral):
            reaction_sites = [int(i) for i in reaction_sites.replace(' ', '').split(',')] if reaction_sites else []
            struct_rules = [s for s in struct_rules.split()] if struct_rules else []
            rules = {'MW': mw, 'TPSA': tpsa, 'HBD': hbd, 'HBA': hba, 'Rings': rings, 'ArRings': ar_rings,
                     'RotB': rot_bonds, 'Chiral': chiral}
            enumerator = SiteEnumerator(molecule, bb_source, reaction_sites, tags, rules, struct_rules)
            enumerator.enumerate()
            enum_mols = enumerator.get_results(as_dict=True)
            n = len(enum_mols)

            # No enumeration returns the query molecule, so n = 1
            if n < 2: 
                enum_img = "No Enumerations! Check Your Inputs!"
                bg_color = 'firebrick'
            else:
                enum_img = f"Slide to view {n} enumerated molecules!"
                bg_color = 'forestgreen'
            enum_style = {'backgroundColor': bg_color, 'width': '80%', 'height': '150px', 
                        'text-align': 'center', 'fontSize': '24px', 'fontFamily': 'sans-serif', 
                        'borderRadius': '20px', 'alignItems': 'center', 'justifyContent': 'center', 
                        'display': 'flex', 'margin': '0 auto'}
            slider_max = n - 1
            slider_marks = {i: str(i) for i in range(0, n, max(10, round(n//10, -1)))}
            
            return enum_img, enum_style, slider_max, slider_marks, enum_mols


shared_callbacks("molecule")
shared_callbacks("site")

if __name__ == "__main__":
    app.run(debug=True, port=8052)
