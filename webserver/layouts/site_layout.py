'''
    Layout components for the SiteHEALER interface.
'''
import dash_bootstrap_components as dbc
from dash import dcc, html
import os
from .base_layout import create_base_layout


def get_server_info_banner():
    """Get server mode information banner."""
    server_mode = os.environ.get('HEALER_SERVER_MODE', 'false').lower() == 'true'
    
    if server_mode:
        return dbc.Alert(
            [
                html.I(className="fas fa-info-circle me-2"),
                "Running in server mode with computational limits. ",
                html.A("Download the repository", href="https://github.com/eneskelestemur/healer", target="_blank"),
                " to run unlimited locally."
            ],
            color="info",
            className="mb-3"
        )
    else:
        return dbc.Alert(
            [
                html.I(className="fas fa-check-circle me-2"),
                "Running in local mode - no computational limits applied."
            ],
            color="success",
            className="mb-3"
        )


def get_site_layout(app_id: str = "site") -> dbc.Container:
    '''
        Create the layout for the SiteHEALER interface.
        
        Args:
            app_id: Application ID for component naming
            
        Returns:
            Dash Bootstrap Container with the complete layout
    '''
    
    # Server mode information
    server_info = get_server_info_banner()
    
    # Tool description
    tool_description = dbc.Card(
        dbc.CardBody([
            html.H5("Site HEALER", className="card-title"),
            html.P([
                "Site-specific molecular enumeration tool that focuses on modifying specific reactive sites in molecules. ",
                "The tool identifies reactive sites based on user input and generates new molecules by replacing atoms "
                "or fragments at those positions with compatible building blocks."
            ], className="card-text"),
            html.Ul([
                html.Li("Define reactive sites by atom indices"),
                html.Li("Apply molecular property filters and structural rules"),
                html.Li("Generate site-specific modifications"),
            ])
        ]),
        className="mb-4",
    )
    
    # Additional inputs specific to SiteHEALER
    additional_inputs = [
        # Reactive Sites
        dbc.Row(
            [
                dbc.Col([
                    html.Label("Reactive Sites"),
                    html.I(
                        className="fas fa-question-circle ms-1 text-muted",
                        id=f"{app_id}-reactive-sites-tooltip",
                        style={"cursor": "pointer"}
                    ),
                    dbc.Tooltip(
                        "Specify atom indices that will be modified during enumeration. "
                        "Format: comma-separated integers (e.g., '19, 20, 21').",
                        target=f"{app_id}-reactive-sites-tooltip",
                        placement="top"
                    )
                ], width=3),
                dbc.Col(
                    dbc.Input(
                        type="text", 
                        id=f"{app_id}-reactive-sites-input", 
                        placeholder="Reacting atoms, e.g., 19, 20, 21",
                        style={'fontSize': '14px'}
                    ), 
                    width=9
                )
            ],
            className="mb-3"
        ),
        
        # Structure Rules
        dbc.Row(
            [
                dbc.Col([
                    html.Label("Structure Rules"),
                    html.I(
                        className="fas fa-question-circle ms-1 text-muted",
                        id=f"{app_id}-struct-rules-tooltip",
                        style={"cursor": "pointer"}
                    ),
                    dbc.Tooltip(
                        "SMARTS patterns that building blocks must contain. "
                        "Separate multiple patterns with whitespace (e.g., 'c1ccccc1 CC').",
                        target=f"{app_id}-struct-rules-tooltip",
                        placement="top"
                    )
                ], width=3),
                dbc.Col(
                    dbc.Textarea(
                        id=f"{app_id}-struct-rules-input",
                        placeholder="Substructure patterns for BBs separated by whitespace, e.g., c1ccccc1 CC",
                        style={'fontSize': '14px', 'height': '80px'}
                    ),
                    width=9
                )
            ],
            className="mb-3"
        ),
        
        # Molecular Property Rules - Row 1: MW and TPSA
        dbc.Row(
            [
                dbc.Col(html.Label('MW'), width=1),
                dbc.Col(
                    dcc.RangeSlider(
                        id=f'{app_id}-MW-slider', 
                        min=0, max=1000, step=1, value=[0, 500],
                        marks={0: '0', 250: '250', 500: '500', 750: '750', 1000: '1000'},
                        tooltip={'placement': 'bottom', 'always_visible': True}
                    ),
                    width=5
                ),
                dbc.Col(
                    dcc.RangeSlider(
                        id=f'{app_id}-TPSA-slider', 
                        min=0, max=200, step=1, value=[0, 200],
                        marks={0: '0', 50: '50', 100: '100', 150: '150', 200: '200'},
                        tooltip={'placement': 'bottom', 'always_visible': True}
                    ),
                    width=5
                ),
                dbc.Col(html.Label('TPSA'), width=1),
            ],
            className="mb-3"
        ),
        
        # Molecular Property Rules - Row 2: HBD and HBA
        dbc.Row(
            [
                dbc.Col(html.Label("HBD"), width=1),
                dbc.Col(
                    dcc.RangeSlider(
                        id=f'{app_id}-HBD-slider', 
                        min=0, max=10, step=1, value=[0, 5],
                        marks={0: '0', 2: '2', 5: '5', 8: '8', 10: '10'},
                        tooltip={'placement': 'bottom', 'always_visible': True}
                    ),
                    width=5
                ),
                dbc.Col(
                    dcc.RangeSlider(
                        id=f'{app_id}-HBA-slider', 
                        min=0, max=10, step=1, value=[0, 10],
                        marks={0: '0', 2: '2', 5: '5', 8: '8', 10: '10'},
                        tooltip={'placement': 'bottom', 'always_visible': True}
                    ),
                    width=5
                ),
                dbc.Col(html.Label('HBA'), width=1),
            ],
            className="mb-3"
        ),
        
        # Molecular Property Rules - Row 3: Rings and Aromatic Rings
        dbc.Row(
            [
                dbc.Col(html.Label('Rings'), width=1),
                dbc.Col(
                    dcc.RangeSlider(
                        id=f'{app_id}-Rings-slider', 
                        min=0, max=20, step=1, value=[0, 10],
                        marks={0: '0', 5: '5', 10: '10', 15: '15', 20: '20'},
                        tooltip={'placement': 'bottom', 'always_visible': True}
                    ),
                    width=5
                ),
                dbc.Col(
                    dcc.RangeSlider(
                        id=f'{app_id}-ArRings-slider', 
                        min=0, max=10, step=1, value=[0, 5],
                        marks={0: '0', 2: '2', 5: '5', 8: '8', 10: '10'},
                        tooltip={'placement': 'bottom', 'always_visible': True}
                    ),
                    width=5
                ),
                dbc.Col(html.Label('Ar. Rings'), width=1),
            ],
            className="mb-3"
        ),
        
        # Molecular Property Rules - Row 4: Rotatable Bonds and Chiral Centers
        dbc.Row(
            [
                dbc.Col(html.Label('Rot. Bonds'), width=1),
                dbc.Col(
                    dcc.RangeSlider(
                        id=f'{app_id}-RotB-slider', 
                        min=0, max=20, step=1, value=[0, 10],
                        marks={0: '0', 5: '5', 10: '10', 15: '15', 20: '20'},
                        tooltip={'placement': 'bottom', 'always_visible': True}
                    ),
                    width=5
                ),
                dbc.Col(
                    dcc.RangeSlider(
                        id=f'{app_id}-Chiral-slider', 
                        min=0, max=10, step=1, value=[0, 5],
                        marks={0: '0', 2: '2', 5: '5', 8: '8', 10: '10'},
                        tooltip={'placement': 'bottom', 'always_visible': True}
                    ),
                    width=5
                ),
                dbc.Col(html.Label('Chiral'), width=1),
            ],
            className="mb-3"
        ),
    ]
    
    return create_base_layout(app_id, "Site", additional_inputs, server_info, tool_description)
