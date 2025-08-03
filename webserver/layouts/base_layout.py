'''
    Base layout components shared by all HEALER interfaces.
'''
import dash_bootstrap_components as dbc
from dash import dcc, html


def create_base_layout(app_id: str, title: str, additional_inputs: list, server_info=None, tool_description=None, fragment_alert=None) -> dbc.Container:
    '''
        Create the base layout structure shared by both molecule and site healers.
        
        Args:
            app_id: Application ID for component naming
            title: Title for the interface (e.g., "Molecule", "Site")
            additional_inputs: List of additional input components
            server_info: Server mode information banner
            tool_description: Tool description card
            fragment_alert: Fragment healer alert (for molecule healer only)
            
        Returns:
            Complete layout container
    '''
    layout_components = []
    
    # Add server info banner if provided
    if server_info:
        layout_components.append(server_info)
    
    # Add tool description if provided  
    if tool_description:
        layout_components.append(tool_description)
    
    # Add fragment alert if provided
    if fragment_alert:
        layout_components.append(fragment_alert)
    
    # Add main title
    layout_components.append(html.H2(f"{title} HEALER", style={'text-align': 'center', 'margin': '10px 0px 10px 0px'}))
    
    # Add the main input/output layout
    layout_components.append(
        dbc.Row(
            [
                # Inputs Section
                dbc.Col(
                    [
                        # Molecule Input
                        dbc.Row(
                            [
                                dbc.Col(html.Label("Molecule"), width=3),
                                dbc.Col(
                                    dbc.Input(
                                        type="text",
                                        id=f"{app_id}-molecule-input",
                                        placeholder="Enter molecule SMILES",
                                        value="CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C",
                                        style={'fontSize': '14px'}
                                    ),
                                    width=9
                                )
                            ],
                            className="mb-3"
                        ),
                        
                        # Building Block Source
                        dbc.Row(
                            [
                                dbc.Col(html.Label("BB Source"), width=3),
                                dbc.Col(
                                    dbc.Select(
                                        id=f"{app_id}-bb-source-select",
                                        options=[
                                            {"label": "Test (100 BBs)", "value": "test"},
                                            {"label": "US Stock", "value": "US_stock"},
                                            {"label": "EU Stock", "value": "EU_stock"},
                                            {"label": "Global Stock", "value": "Global_stock"}
                                        ],
                                        value="test"
                                    ),
                                    width=9
                                )
                            ],
                            className="mb-3"
                        ),
                        
                        # Reaction Tags
                        dbc.Row(
                            [
                                dbc.Col([
                                    html.Label("Reaction Tags"),
                                    html.I(
                                        className="fas fa-question-circle ms-1 text-muted",
                                        id=f"{app_id}-reaction-tags-tooltip",
                                        style={"cursor": "pointer"}
                                    ),
                                    dbc.Tooltip(
                                        "Select reaction tags to synthesize molecules with. "
                                        "Server limit: 15 tags and 'all' will be removed",
                                        target=f"{app_id}-reaction-tags-tooltip",
                                        placement="top"
                                    )
                                ], width=3),
                                dbc.Col(
                                    dcc.Dropdown(
                                        id=f'{app_id}-reaction-tags-dropdown',
                                        options=[],  # Will be populated by callback
                                        value=["amide coupling", "amide", "C-N bond formation", "C-N",
                                               "alkylation", "N-arylation", "azole", "amination"],
                                        multi=True,
                                    ),
                                    width=9
                                )
                            ],
                            className="mb-3"
                        ),
                        
                        # Additional inputs specific to each healer type
                        *additional_inputs,
                        
                        # Max Evaluations
                        dbc.Row(
                            [
                                dbc.Col([
                                    html.Label("Max Evaluations"),
                                    html.I(
                                        className="fas fa-question-circle ms-1 text-muted",
                                        id=f"{app_id}-max-evals-tooltip",
                                        style={"cursor": "pointer"}
                                    ),
                                    dbc.Tooltip(
                                        "Maximum number of molecules generated for each composition. "
                                        "Higher values explore more chemical space per composition. "
                                        "Server limit: 1-500. Used for all HEALER types.",
                                        target=f"{app_id}-max-evals-tooltip",
                                        placement="top"
                                    )
                                ], width=3),
                                dbc.Col(
                                    dbc.Input(
                                        type="number",
                                        id=f"{app_id}-max-evals-input",
                                        placeholder="None (unlimited)",
                                        style={'fontSize': '14px'}
                                    ),
                                    width=9
                                )
                            ],
                            className="mb-3"
                        ),
                        
                        # Save Filename
                        dbc.Row(
                            [
                                dbc.Col(html.Label("Save as"), width=3),
                                dbc.Col(
                                    dbc.Input(
                                        type="text", 
                                        id=f"{app_id}-save-as-input", 
                                        placeholder="Filename", 
                                        value="enumerated_molecules.csv", 
                                        style={'fontSize': '14px'}
                                    ),
                                    width=9
                                )
                            ],
                            className="mb-3"
                        ),
                        
                        # Action Buttons
                        dbc.Row(
                            [
                                dbc.Col(
                                    dbc.Button(
                                        "Update", 
                                        id=f"{app_id}-update-button", 
                                        color="primary", 
                                        className="me-2",
                                        style={'width': '100%', 'text-align': 'center'}
                                    ), 
                                    width=3
                                ),
                                dbc.Col(
                                    dbc.Button(
                                        "Enumerate", 
                                        id=f"{app_id}-enumerate-button", 
                                        color="warning", 
                                        className="me-2",
                                        style={'width': '100%', 'text-align': 'center'}
                                    ), 
                                    width=3
                                ),
                                dbc.Col(
                                    dbc.Button(
                                        "Save", 
                                        id=f"{app_id}-save-button", 
                                        color="success", 
                                        className="me-2",
                                        style={'width': '100%', 'text-align': 'center'}
                                    ), 
                                    width=3
                                ),
                                dbc.Col(
                                    dbc.Button(
                                        "Cancel", 
                                        id=f"{app_id}-cancel-button", 
                                        color="danger",
                                        style={'width': '100%', 'text-align': 'center'},
                                        disabled=True
                                    ), 
                                    width=3
                                ),
                            ],
                            className="mb-3",
                            justify="center",
                        ),
                    ],
                    md=6,
                    style={'height': '100%', 'justify_content': 'center', 'margin': '-10px 0px 0px 0px'}
                ),
                
                # Outputs Section
                dbc.Col(
                    [
                        # Query Molecule Card
                        dbc.Card(
                            [
                                dbc.CardHeader(
                                    html.H4("Query Molecule", className="mb-0 text-center"),
                                    className="py-2"
                                ),
                                dbc.CardBody(
                                    html.Div(
                                        id=f'{app_id}-molecule-svg', 
                                        style={'height': '150px', 'text-align': 'center'},
                                    ),
                                    className="p-3"
                                )
                            ],
                            className="mb-4"
                        ),
                        
                        # Enumeration Results Card
                        dbc.Card(
                            [
                                dbc.CardHeader(
                                    html.H4("Enumeration Results", className="mb-0 text-center"),
                                    className="py-2"
                                ),
                                dbc.CardBody(
                                    [
                                        html.Div(
                                            [
                                                dcc.Slider(
                                                    id=f'{app_id}-enumeration-results-slider', 
                                                    min=0, max=0, step=1, value=0,
                                                    marks={0: "0"},
                                                    tooltip={'placement': 'bottom', 'always_visible': False}
                                                )
                                            ],
                                            style={'width': '70%', 'justify-content': 'center', 'margin': '0 auto 20px auto'},
                                        ),
                                        
                                        dcc.Loading(
                                            [
                                                html.Div(
                                                    'Click Enumerate to View the Enumerations', 
                                                    id=f'{app_id}-enumeration-svg', 
                                                    style={
                                                        'width': '75%', 
                                                        'height': '200px', 
                                                        'text-align': 'center', 
                                                        'fontSize': '20px', 
                                                        'fontFamily': 'sans-serif', 
                                                        'alignItems': 'center', 
                                                        'justifyContent': 'center', 
                                                        'display': 'flex', 
                                                        'margin': '0 auto'
                                                    }
                                                ),
                                            ],
                                            id=f"{app_id}-enumeration-loading", 
                                            delay_show=100, 
                                            type='circle',
                                            overlay_style={"visibility":"visible", "filter": "blur(2px)"}
                                        ),
                                    ],
                                    className="p-3"
                                )
                            ]
                        ),
                        
                        dcc.Download(id=f"{app_id}-download-enumerations"),
                        dcc.Store(id=f'{app_id}-enumeration-store')
                    ],
                    md=6,
                    style={'height': '100%', 'justify_content': 'center', 'margin': '20px 0px 0px 0px'}
                ),
            ],
            style={'height': '100%', 'justify_content': 'center', 'margin': '40px 0px 0px 0px'}
        )
    )
    
    return dbc.Container(
        layout_components,
        fluid=True,
        style={'width': '80%'}
    )
