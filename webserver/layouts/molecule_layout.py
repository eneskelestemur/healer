'''
    Layout components for the MoleculeHEALER interface.
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


def get_molecule_layout(app_id: str = "molecule") -> dbc.Container:
    '''
        Create the layout for the MoleculeHEALER interface.
        
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
            html.H5("Molecule HEALER", className="card-title"),
            html.P([
                "Molecular enumeration tool that performs retrosynthetic analysis on input molecules. ",
                "Automatically switches to ", html.Strong("Fragment HEALER"), " when multi-component molecules are detected. ",
                "The tool breaks down molecules into fragments, searches for suitable building blocks, and reconstructs new molecules."
            ], className="card-text"),
            html.Ul([
                html.Li([html.Strong("Molecule HEALER:"), " Single-component molecules"]),
                html.Li([html.Strong("Fragment HEALER:"), " Multi-component molecules (fragments separated by '.')"]),
            ])
        ]),
        className="mb-4",
    )
    
    # Fragment HEALER alert (hidden by default)
    fragment_alert = dbc.Alert(
        [
            html.I(className="fas fa-info-circle me-2"),
            html.Strong("Fragment HEALER mode detected! "), 
            "Multi-component molecule detected. The system will automatically use Fragment HEALER for enumeration."
        ],
        id=f"{app_id}-fragment-alert",
        color="info",
        is_open=False,
        dismissable=True,  # Note: 'dismissable' not 'dismissible'
        className="mb-3"
    )
    
    # Additional inputs specific to MoleculeHEALER
    additional_inputs = [
        # Custom Split Sites
        dbc.Row(
            [
                dbc.Col([
                    html.Label("Split Sites"),
                    html.I(
                        className="fas fa-question-circle ms-1 text-muted",
                        id=f"{app_id}-split-sites-tooltip",
                        style={"cursor": "pointer"}
                    ),
                    dbc.Tooltip(
                        "Specify atom indices where the molecule should be split for enumeration. "
                        "You can use multiple fragment compositions separated by semicolons. "
                        "Format: 'atom1-atom2, atom3-atom4; atom5-atom6' (e.g., '9-10, 6-9; 12-13'). "
                        "Commas separate pairs within a group, semicolons separate different groups. "
                        "Used for Molecule HEALER only.",
                        target=f"{app_id}-split-sites-tooltip",
                        placement="top"
                    )
                ], width=3),
                dbc.Col(
                    dbc.Input(
                        type="text",
                        id=f"{app_id}-custom-sites-input",
                        placeholder="Sites to split molecule (e.g., 9-10, 6-9; 12-13)",
                        style={'fontSize': '14px'}
                    ),
                    width=9
                )
            ],
            className="mb-3"
        ),
        
        # Similarity Threshold
        dbc.Row(
            [
                dbc.Col([
                    html.Label("Similarity Threshold"),
                    html.I(
                        className="fas fa-question-circle ms-1 text-muted",
                        id=f"{app_id}-sim-threshold-tooltip",
                        style={"cursor": "pointer"}
                    ),
                    dbc.Tooltip(
                        "Minimum Tversky similarity for building block selection. "
                        "Lower values include more diverse building blocks. "
                        "Server limit: ≥0.3. Used for both Molecule and Fragment HEALER.",
                        target=f"{app_id}-sim-threshold-tooltip",
                        placement="top"
                    )
                ], width=3),
                dbc.Col(
                    dcc.Slider(
                        id=f"{app_id}-sim-threshold-slider",
                        min=0, max=1, step=0.01, value=0.15,
                        marks={i / 10: str(i / 10) for i in range(11)},
                        tooltip={'placement': 'bottom', 'always_visible': False}
                    ),
                    width=9
                )
            ],
            className="mb-3"
        ),
        
        # Number of Compositions
        dbc.Row(
            [
                dbc.Col([
                    html.Label("N Compositions"),
                    html.I(
                        className="fas fa-question-circle ms-1 text-muted",
                        id=f"{app_id}-n-compositions-tooltip",
                        style={"cursor": "pointer"}
                    ),
                    dbc.Tooltip(
                        "Number of molecular compositions to generate. "
                        "Use Randomize Compositions to increase diversity. "
                        "Server limit: ≤50. Used for Molecule HEALER only.",
                        target=f"{app_id}-n-compositions-tooltip",
                        placement="top"
                    )
                ], width=3),
                dbc.Col(
                    dcc.Slider(
                        id=f"{app_id}-n-compositions-slider",
                        min=1, max=1000, step=1, value=10,
                        marks={i: str(i) for i in [10, 50, 100, 500, 1000]},
                        tooltip={'placement': 'bottom', 'always_visible': False}
                    ),
                    width=9
                )
            ],
            className="mb-3"
        ),
        
        # Retro Tree Depth and Min Fragment Size
        dbc.Row(
            [
                dbc.Col([
                    html.Label("Retro Depth"),
                    html.I(
                        className="fas fa-question-circle ms-1 text-muted",
                        id=f"{app_id}-retro-depth-tooltip",
                        style={"cursor": "pointer"}
                    ),
                    dbc.Tooltip(
                        "Maximum retrosynthetic tree depth. Higher values generate more complex structures. "
                        "Server limit: 1-2. Used for Molecule HEALER only.",
                        target=f"{app_id}-retro-depth-tooltip",
                        placement="top"
                    )
                ], width=2),
                dbc.Col(
                    dcc.Slider(
                        id=f"{app_id}-retro-depth-slider",
                        min=1, max=5, step=1, value=1,
                        marks={i: str(i) for i in range(1, 6)},
                        tooltip={'placement': 'bottom', 'always_visible': False}
                    ),
                    width=4
                ),
                dbc.Col([
                    html.Label("Min Frag Size"),
                    html.I(
                        className="fas fa-question-circle ms-1 text-muted",
                        id=f"{app_id}-min-frag-size-tooltip",
                        style={"cursor": "pointer"}
                    ),
                    dbc.Tooltip(
                        "Minimum number of atoms in molecular fragments after retro-synthesis. "
                        "Used for Molecule HEALER only.",
                        target=f"{app_id}-min-frag-size-tooltip",
                        placement="top"
                    )
                ], width=2),
                dbc.Col(
                    dcc.Slider(
                        id=f"{app_id}-min-frag-size-slider",
                        min=1, max=10, step=1, value=3,
                        marks={i: str(i) for i in [1, 3, 5, 7, 10]},
                        tooltip={'placement': 'bottom', 'always_visible': False}
                    ),
                    width=4
                )
            ],
            className="mb-3"
        ),
        
        # Advanced Options
        dbc.Row(
            [
                dbc.Col([
                    html.Label("Max BBs/Comp"),
                    html.I(
                        className="fas fa-question-circle ms-1 text-muted",
                        id=f"{app_id}-max-bbs-tooltip",
                        style={"cursor": "pointer"}
                    ),
                    dbc.Tooltip(
                        "Maximum building blocks per composition. -1 for unlimited. "
                        "Server limit: 1-10. Used for both Molecule and Fragment HEALER.",
                        target=f"{app_id}-max-bbs-tooltip",
                        placement="top"
                    )
                ], width=3),
                dbc.Col(
                    dbc.Input(
                        type="number",
                        id=f"{app_id}-max-bbs-input",
                        placeholder="-1 (unlimited)",
                        value=5,
                        style={'fontSize': '14px'}
                    ),
                    width=2
                ),
                dbc.Col([
                    html.Div(
                        [
                            dbc.Checklist(
                                options=[{"label": "Randomize Compositions", "value": "randomize"}],
                                value=[],
                                id=f"{app_id}-randomize-checkbox",
                                inline=True,
                                style={'display': 'inline-block', 'margin-right': '10px', 'inline-size': '7em', 'align-items': 'center'},
                            ),
                            html.I(
                                className="fas fa-question-circle ms-1 text-muted",
                                id=f"{app_id}-randomize-tooltip",
                                style={"cursor": "pointer", 'display': 'inline-block', 'margin-right': '15px'}
                            ),
                            html.Label(
                                "Random Seed", 
                                style={
                                    'textAlign': 'right',
                                    'display': 'none',
                                    'inlineSize': '2em'
                                },
                                id=f"{app_id}-random-seed-label"
                            ),
                            dbc.Input(
                                type="number",
                                id=f"{app_id}-random-seed-input",
                                placeholder="-1 (auto)",
                                value=-1,
                                style={
                                    'fontSize': '14px', 
                                    'width': '100px',
                                    'display': 'none'
                                },
                                disabled=True
                            )
                        ],
                        id=f"{app_id}-random-seed-row",
                        style={'display': 'flex', 'align-items': 'center'}
                    ),
                    dbc.Tooltip(
                        "Randomize the order of generated compositions to increase diversity. "
                        "Used for Molecule HEALER only.",
                        target=f"{app_id}-randomize-tooltip",
                        placement="top"
                    ),
                ], width=7)
            ],
            className="mb-3"
        ),
    ]
    
    return create_base_layout(app_id, "Molecule", additional_inputs, server_info, tool_description, fragment_alert)
