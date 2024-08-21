'''
    Dash application for the enumerator project.
'''

# Automated Enumerator using Dash
from dash import dcc, html, Dash, Input, Output, State, callback_context, no_update

import dash_bootstrap_components as dbc
import utils

# List of unique reaction tags
reactions = utils.load_reactions_from_json('reactions/reactions.json')
reactions = [r for r in reactions if r.is_valid()]
reaction_tags = [r.tags for r in reactions if r.is_valid()]
reaction_tags = list(set([tag for tags in reaction_tags for tag in tags]))


app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

app.layout = dbc.Container(
    [
        html.H1("Automated Enumerator", style={'text-align': 'center', 'margin': '10px 0px 10px 0px'}),
        dbc.Row(
            [
                dbc.Col(
                    [
                        dbc.Row(
                            [
                                dbc.Col(html.Label("Molecule"), width=3),
                                dbc.Col(dbc.Input(type="text", id="molecule-input", placeholder="Enter molecule", 
                                                  value="CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C",
                                                  style={'fontSize': '14px'}),
                                        width=9)
                            ],
                            className="mb-3"
                        ),
                        dbc.Row(
                            [
                                dbc.Col(html.Label("BB source"), width=3),
                                dbc.Col(
                                    dbc.Select(
                                        id="bb-source-select",
                                        options=[
                                            {"label": "test", "value": "test"},
                                            {"label": "US_stocks", "value": "US_stocks"},
                                            {"label": "EU_stocks", "value": "EU_stocks"},
                                            {"label": "global_stocks", "value": "global_stocks"}
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
                                dbc.Col(html.Label("Custom Sites"), width=3),
                                dbc.Col(dbc.Input(type="text", id="custom-sites-input", 
                                                  placeholder="Custom sites to split the molecule, e.g., 9,10; 6,9",
                                                  style={'fontSize': '14px'}), 
                                        width=9)
                            ],
                            className="mb-3"
                        ),
                        dbc.Row(
                            [
                                dbc.Col(html.Label("Sim Cutoff"), width=3),
                                dbc.Col(
                                    dcc.Slider(
                                        id='sim-cutoff-slider',
                                        min=0,
                                        max=1,
                                        step=0.01,
                                        value=0.15,
                                        marks={i / 10: str(i / 10) for i in range(11)},
                                        tooltip={'placement': 'bottom', 'always_visible': False}
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
                                        id='reaction-tags-dropdown',
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
                        dbc.Row(
                            [
                                dbc.Col(html.Label("Save as"), width=3),
                                dbc.Col(dbc.Input(type="text", id="save-as-input", placeholder="Filename", 
                                                  value="enumerated_molecules.csv", style={'fontSize': '14px'}),
                                        width=9)
                            ],
                            className="mb-3"
                        ),
                        dbc.Row(
                            [
                                dbc.Col(dbc.Button("Update", id="update-button", color="primary", className="mr-2",
                                                   style={'width': '100%', 'text-align': 'center'}), width=3),
                                dbc.Col(dbc.Button("Enumerate", id="enumerate-button", color="warning", className="mr-2",
                                                   style={'width': '100%', 'text-align': 'center'}), width=3),
                                dbc.Col(dbc.Button("Save", id="save-button", color="success",
                                                   style={'width': '100%', 'text-align': 'center'}), width=3),
                            ],
                            className="mb-3",
                            justify='end',
                        ),
                        dbc.Row(
                            [
                                dbc.Col(dbc.Button(id="Cancel Enumeration Job", color="danger", 
                                                   style={'width': '100%', 'text-align': 'center'}), width=9),
                            ],
                            className="mb-3",
                            justify='end',
                        ),
                    ],
                    md=6,
                    style={'height': '100%', 'justify_content': 'center', 'margin': '20px 0px 0px 0px'}
                ),
                dbc.Col(
                    [
                        html.H3("Molecule", style={'text-align': 'center'}),
                        html.Div(id='molecule-svg', style={'height': '150px', 'text-align': 'center'}),
                        html.H3("Enumeration Results", style={'text-align': 'center', 'margin': '20px 0px 20px 0px'}),
                        html.Div(dcc.Slider(id='enumeration-results-slider', min=0, max=0, step=1, value=0,
                                            tooltip={'placement': 'bottom', 'always_visible': False}), 
                                 style={'width': '70%', 'justify-content': 'center', 'margin': '0 auto'}),
                        html.Div(id='enumeration-svg', style={'height': '150px', 'text-align': 'center'}),
                        dcc.Download(id="download-enumerations"),
                        dcc.Store(id='enumeration-store')
                    ],
                    md=6,
                ),
            ],
            style={'height': '100%', 'justify_content': 'center', 'margin': '40px 0px 0px 0px'}
        ),
    ],
    fluid=True,
    style={'width': '80%'}
)

@app.callback(
    [Output('molecule-svg', 'children'),
     Output('molecule-svg', 'style'),
     Output('enumeration-svg', 'children'),
     Output('enumeration-svg', 'style'),
     Output('enumeration-results-slider', 'max'),
     Output('enumeration-results-slider', 'marks'),
     Output('enumeration-store', 'data'),
     Output('download-enumerations', 'data')],
    [Input('update-button', 'n_clicks'),
     Input('enumerate-button', 'n_clicks'),
     Input('enumeration-results-slider', 'value'),
     Input('save-button', 'n_clicks')],
    [State('molecule-input', 'value'),
     State('bb-source-select', 'value'),
     State('custom-sites-input', 'value'),
     State('sim-cutoff-slider', 'value'),
     State('reaction-tags-dropdown', 'value'),
     State('save-as-input', 'value'),
     State('enumeration-store', 'data')]
)
def update_output(update_clicks, enumerate_clicks, enumeration_value, save_clicks, 
                  molecule, building_blocks, custom_sites, sim_threshold, reaction_tags, save_as, enumerated_molecules):
    ctx = callback_context
    if not ctx.triggered:
        mol_img = 'Click Update to View the Molecule'
        mol_style = {'width': '75%', 'height': '175px', 'text-align': 'center', 
                     'fontSize': '20px', 'fontFamily': 'sans-serif', 'alignItems': 'center', 
                     'justifyContent': 'center', 'display': 'flex', 'margin': '0 auto'}
        enum_img = 'Click Enumerate to View the Enumerations'
        enum_style = {'width': '75%', 'height': '200px', 'text-align': 'center', 
                     'fontSize': '20px', 'fontFamily': 'sans-serif', 'alignItems': 'center', 
                     'justifyContent': 'center', 'display': 'flex', 'margin': '0 auto'}
        return mol_img, mol_style, enum_img, enum_style, 0, {}, None, None
    else:
        button_id = ctx.triggered[0]['prop_id'].split('.')[0]
        if button_id == 'update-button':
            mol_img = html.Img(src=utils.get_svg_mol(molecule, show_idx=True),
                               style={'height': '100%', 'width': '80%', 'margin': '0 auto'})
            mol_style = {'height': '175px', 'text-align': 'center'}
            return mol_img, mol_style, no_update, no_update, no_update, no_update, no_update, no_update
        elif button_id == 'enumerate-button':
            custom_sites = [tuple(map(int, site.replace(' ', '').split(',')))
                            for site in custom_sites.split(';') if site] if custom_sites else []
            enumerator = utils.automated_enumerator(molecule, building_blocks, reaction_tags, custom_sites, 10, sim_threshold)
            if not enumerator.enumerated_molecules:
                mol_img = html.Img(src=utils.get_svg_mol(molecule, show_idx=True),
                                   style={'height': '100%', 'width': '80%', 'margin': '0 auto'})
                mol_style = {'height': '175px', 'text-align': 'center'}
                slider_max = 0
                slider_marks = {0: '0'}
                enum_img = "No Enumerations! Check Inputs!"
                enum_style = {'backgroundColor': 'firebrick', 'width': '75%', 'height': '150px', 
                              'text-align': 'center', 'fontSize': '24px', 'fontFamily': 'sans-serif', 
                              'borderRadius': '20px', 'alignItems': 'center', 'justifyContent': 'center', 
                              'display': 'flex', 'margin': '0 auto'}
                return mol_img, mol_style, enum_img, enum_style, slider_max, slider_marks, [], no_update
            else:
                mol_img = html.Img(src=utils.get_svg_mol(molecule, show_idx=True),
                                   style={'height': '100%', 'width': '75%', 'margin': '0 auto'})
                mol_style = {'height': '175px', 'text-align': 'center'}
                enum_img = "Slide to view enumerated molecules!"
                enum_style = {'backgroundColor': 'forestgreen', 'width': '80%', 'height': '150px', 
                              'text-align': 'center', 'fontSize': '24px', 'fontFamily': 'sans-serif', 
                              'borderRadius': '20px', 'alignItems': 'center', 'justifyContent': 'center', 
                              'display': 'flex', 'margin': '0 auto'}
                slider_max = len(enumerator.enumerated_molecules) - 1
                n = len(enumerator.enumerated_molecules)
                slider_marks = {i: str(i) for i in range(0, n, max(10, round(n//10, -1)))}
                enum_mols = enumerator.enumerated_molecules
                return mol_img, mol_style, enum_img, enum_style, slider_max, slider_marks, enum_mols, no_update
        elif button_id == 'enumeration-results-slider':
            if enumerated_molecules:
                img = utils.get_svg_mol_with_bbs(*enumerated_molecules[enumeration_value][:3],
                                                 bb_colors=['purple', 'green'],
                                                 legend=enumerated_molecules[enumeration_value][3])
                enum_img = html.Img(src=img,
                                    style={'height': '100%', 'width': '80%', 'margin': '0 auto'})
                enum_style = {'height': '200px', 'text-align': 'center'}
            else:
                enum_img = "No Enumerations! Couldn't Save!"
                enum_style = {'backgroundColor': 'firebrick', 'width': '75%', 'height': '150px', 
                              'text-align': 'center', 'fontSize': '24px', 'fontFamily': 'sans-serif', 
                              'borderRadius': '20px', 'alignItems': 'center', 'justifyContent': 'center', 
                              'display': 'flex', 'margin': '0 auto'}
            return no_update, no_update, enum_img, enum_style, no_update, no_update, no_update, no_update
        elif button_id == 'save-button':
            enum_img = f'Saved to {save_as}!'
            enum_style = {'backgroundColor': 'forestgreen', 'width': '75%', 'height': '150px',
                          'text-align': 'center', 'fontSize': '24px', 'fontFamily': 'sans-serif',
                          'borderRadius': '20px', 'alignItems': 'center', 'justifyContent': 'center', 
                          'display': 'flex', 'margin': '0 auto'}
            if enumerated_molecules:
                save_out = dcc.send_data_frame(utils.automated_enumerator_download_df(enumerated_molecules).to_csv, filename=save_as)
            else:
                save_out = None
            return no_update, no_update, enum_img, enum_style, no_update, no_update, no_update, save_out

if __name__ == '__main__':
    app.run(debug=False, jupyter_mode='tab', port=8051)
