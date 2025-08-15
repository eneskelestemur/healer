'''
    New HEALER Dash Web Application
    
    A refactored web interface for the HEALER molecular enumeration system,
    supporting both MoleculeHEALER, FragmentHEALER and SiteHEALER workflows.
'''
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from flask import request, jsonify, send_from_directory
import dash_bootstrap_components as dbc
import diskcache
from dash import Dash, dcc, html, Input, Output, DiskcacheManager, no_update, clientside_callback
from rdkit import Chem
from rdkit.Chem import rdDepictor

from layouts.molecule_layout import get_molecule_layout
from layouts.site_layout import get_site_layout

from callbacks.shared_callbacks import register_shared_callbacks
from callbacks.molecule_callbacks import register_molecule_callbacks
from callbacks.site_callbacks import register_site_callbacks
from callbacks.ketcher_callbacks import register_ketcher_callbacks

# Configuration for server vs local mode
SERVER_MODE = os.environ.get('HEALER_SERVER_MODE', 'false').lower() == 'true'


color_mode_switch = html.Span(
    [
        dbc.Label(className="fa fa-moon", html_for="color-mode-switch"),
        dbc.Switch( 
            id="color-mode-switch", 
            value=True, 
            className="d-inline-block ms-1", 
            persistence=True
        ),
        dbc.Label(className="fa fa-sun", html_for="color-mode-switch"),
    ],
    style={
        'position': 'absolute',
        'top': '20px',
        'right': '20px',
        'z-index': '1000'
    }
)

webserver_dir = os.path.dirname(os.path.abspath(__file__))
cache_dir = os.path.join(webserver_dir, 'cache')
cache = diskcache.Cache(cache_dir)
background_callback_manager = DiskcacheManager(cache)

app = Dash(
    __name__, 
    external_stylesheets=[
        dbc.themes.MORPH, 
        dbc.icons.FONT_AWESOME,
    ],
    background_callback_manager=background_callback_manager,
    suppress_callback_exceptions=True
)

app.title = "HEALER - Molecular Enumeration Dashboard"

# Configure static file serving for Ketcher standalone
@app.server.route('/ketcher/<path:filename>')
def serve_ketcher(filename):
    # Get the absolute path to the webserver directory
    webserver_dir = os.path.dirname(os.path.abspath(__file__))
    ketcher_path = os.path.join(webserver_dir, 'assets', 'ketcher-standalone-3.4.0', 'standalone')
    return send_from_directory(ketcher_path, filename)

# API endpoint to convert SMILES to Molfile for Ketcher import
@app.server.route('/api/smiles-to-molfile', methods=['POST'])
def smiles_to_molfile():
    try:
        data = request.get_json()
        smiles = data.get('smiles', '').strip()
        
        if not smiles:
            return jsonify({'error': 'No SMILES provided'}), 400
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES format'}), 400

        rdDepictor.Compute2DCoords(mol)
        
        # Convert to Molfile
        molblock = Chem.MolToMolBlock(mol)
        
        return jsonify({'molblock': molblock})
        
    except Exception as e:
        return jsonify({'error': f'Error converting SMILES: {str(e)}'}), 500

app.layout = dbc.Container(
    [
        # Color mode switch positioned in top-right corner
        color_mode_switch,
        
        html.H1(
            "HEALER Dashboard", 
            style={
                'text-align': 'center', 
                'margin': '20px 0px',
                'font-weight': 'bold'
            }
        ),
        
        html.P(
            "Molecular enumeration using retrosynthetic analysis and building block libraries",
            style={
                'text-align': 'center',
                'font-size': '16px',
                'margin-bottom': '30px'
            }
        ),
        
        # Tab selection
        html.Div(
            dcc.Tabs(
                id="tabs",
                value="molecule-tab",
                children=[
                    dcc.Tab(
                        label="Molecule HEALER", 
                        value="molecule-tab",
                        style={
                            'text-align': 'center', 
                            'padding': '12px 24px',
                            'font-size': '16px',
                            'border-radius': '8px 8px 0 0'
                        },
                        selected_style={
                            'text-align': 'center', 
                            'padding': '12px 24px',
                            'font-size': '16px',
                            'backgroundColor': "#6bb5e6",
                            # 'color': 'white',
                            'border-radius': '8px 8px 0 0',
                            'border-color': "#0099ff",
                        }
                    ),
                    dcc.Tab(
                        label="Site HEALER", 
                        value="site-tab",
                        style={
                            'text-align': 'center', 
                            'padding': '12px 24px',
                            'font-size': '16px',
                            'border-radius': '8px 8px 0 0'
                        },
                        selected_style={
                            'text-align': 'center', 
                            'padding': '12px 24px',
                            'font-size': '16px',
                            'backgroundColor': "#ca89eb",
                            # 'color': 'white',
                            'border-radius': '8px 8px 0 0',
                            'border-color': "#a200ff",
                        }
                    ),
                ],
                style={
                    'width': '60%', 
                    'margin': '0 auto',
                    'border-radius': '8px'
                }
            ),
            style={'text-align': 'center', 'margin-bottom': '20px'}
        ),
        
        # Tab content containers
        html.Div(
            [
                html.Div(
                    get_molecule_layout("molecule"), 
                    id="molecule-tab-content", 
                    style={'display': 'block'}
                ),
                html.Div(
                    get_site_layout("site"), 
                    id="site-tab-content", 
                    style={'display': 'none'}
                )
            ]
        )
    ],
    fluid=True,
    style={'min-height': '100vh', 'padding': '20px'}
)


@app.callback(
    [Output("molecule-tab-content", "style"),
     Output("site-tab-content", "style"),
     Output("molecule-tab-content", "children"),
     Output("site-tab-content", "children")],
    [Input("tabs", "value")]
)
def toggle_tab_content(tab_name):
    '''Switch between molecule and site healer interfaces.'''
    molecule_content = get_molecule_layout("molecule")
    site_content = get_site_layout("site")
    
    if tab_name == "molecule-tab":
        return (
            {'display': 'block'}, 
            {'display': 'none'}, 
            molecule_content, 
            no_update
        )
    elif tab_name == "site-tab":
        return (
            {'display': 'none'}, 
            {'display': 'block'}, 
            no_update, 
            site_content
        )
    else:
        return (
            {'display': 'block'}, 
            {'display': 'none'}, 
            molecule_content, 
            site_content
        )


# Color mode toggle clientside callback
clientside_callback(
    """
    (switchOn) => {
       document.documentElement.setAttribute("data-bs-theme", switchOn ? "light" : "dark");
       return window.dash_clientside.no_update
    }
    """,
    Output("color-mode-switch", "id"),
    Input("color-mode-switch", "value"),
)


# Register all callbacks
register_shared_callbacks(app, "molecule")
register_shared_callbacks(app, "site")
register_molecule_callbacks(app, "molecule")
register_site_callbacks(app, "site")
register_ketcher_callbacks(app, "molecule")
register_ketcher_callbacks(app, "site")


if __name__ == "__main__":
    # Run the application
    app.run(debug=True, host='0.0.0.0', port=8053)
