"""
    Ketcher Molecule Sketcher Component using iframe
"""
import dash_bootstrap_components as dbc
from dash import html

def create_ketcher_modal(app_id: str) -> dbc.Modal:
    """
        Create a Ketcher molecule sketcher modal using iframe.
        
        Args:
            app_id: Application ID for component naming
            
        Returns:
            Dash Bootstrap Modal component with Ketcher iframe
    """
    return dbc.Modal([
        dbc.ModalHeader([
            dbc.ModalTitle("Draw Molecule - Ketcher Editor"),
        ]),
        dbc.ModalBody([
            # Error alert (hidden by default)
            dbc.Alert(
                "",
                id=f"{app_id}-ketcher-error",
                color="danger",
                is_open=False,
                className="mb-3"
            ),
            # Ketcher iframe container
            html.Iframe(
                id=f"{app_id}-ketcher-iframe",
                src="/ketcher/index.html",
                style={
                    "width": "100%",
                    "height": "500px",
                    "border": "1px solid #ddd",
                    "borderRadius": "4px"
                }
            ),
        ]),
        dbc.ModalFooter([
            dbc.Button(
                "Get SMILES",
                id=f"{app_id}-ketcher-get-smiles",
                color="primary",
                className="me-2"
            ),
            dbc.Button(
                "Clear",
                id=f"{app_id}-ketcher-clear",
                color="secondary",
                className="me-2"
            ),
            dbc.Button(
                "Cancel",
                id=f"{app_id}-ketcher-cancel",
                color="light"
            )
        ])
    ], 
    id=f"{app_id}-ketcher-modal", 
    size="xl",  # Extra large for better drawing space
    centered=True,
    backdrop="static"
    )


def create_ketcher_button(app_id: str) -> dbc.Button:
    """
        Create a button to open the Ketcher modal.
        
        Args:
            app_id: Application ID for component naming
            
        Returns:
            Dash Bootstrap Button component
    """
    return dbc.Button(
        "Ketcher",
        id=f"{app_id}-ketcher-btn",
        color="outline-primary",
        size="sm",
        className="h-100"
    )
