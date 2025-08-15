"""
    Callback functions for Ketcher molecule sketcher integration using iframe.
"""
from dash import Input, Output, State, clientside_callback


def register_ketcher_callbacks(app, app_id: str = "molecule"):
    """
        Register callbacks for Ketcher molecule sketcher functionality.
        
        Args:
            app: Dash application instance
            app_id: Application ID for component naming
    """
    
    # Modal open/close callback
    @app.callback(
        Output(f"{app_id}-ketcher-modal", "is_open", allow_duplicate=True),
        [Input(f"{app_id}-ketcher-btn", "n_clicks"),
         Input(f"{app_id}-ketcher-cancel", "n_clicks")],
        State(f"{app_id}-ketcher-modal", "is_open"),
        prevent_initial_call=True
    )
    def toggle_ketcher_modal(open_clicks, cancel_clicks, is_open):
        """Toggle the Ketcher modal visibility."""
        if open_clicks or cancel_clicks:
            return not is_open
        return is_open
    
    # Get SMILES from Ketcher iframe
    clientside_callback(
        """
        function(get_clicks) {
            if (get_clicks > 0) {
                try {
                    const iframe = document.getElementById('""" + app_id + """-ketcher-iframe');
                    if (iframe && iframe.contentWindow && iframe.contentWindow.ketcher) {
                        // Get SMILES from Ketcher iframe
                        return iframe.contentWindow.ketcher.getSmiles().then(function(smiles) {
                            console.log('Got SMILES from Ketcher:', smiles);
                            if (smiles && smiles.trim() !== '') {
                                // Return SMILES and close modal
                                return [smiles, false, ''];
                            } else {
                                // Show error
                                return [window.dash_clientside.no_update, window.dash_clientside.no_update, 'No molecule drawn. Please draw a molecule first.'];
                            }
                        }).catch(function(error) {
                            console.error('Error getting SMILES:', error);
                            return [window.dash_clientside.no_update, window.dash_clientside.no_update, 'Error getting molecule structure.'];
                        });
                    } else {
                        console.error('Ketcher iframe not accessible');
                        return [window.dash_clientside.no_update, window.dash_clientside.no_update, 'Ketcher not properly loaded.'];
                    }
                } catch (error) {
                    console.error('Error in get SMILES callback:', error);
                    return [window.dash_clientside.no_update, window.dash_clientside.no_update, 'Error retrieving molecule structure.'];
                }
            }
            return [window.dash_clientside.no_update, window.dash_clientside.no_update, window.dash_clientside.no_update];
        }
        """,
        [Output(f"{app_id}-molecule-input", "value", allow_duplicate=True),
         Output(f"{app_id}-ketcher-modal", "is_open", allow_duplicate=True),
         Output(f"{app_id}-ketcher-error", "children", allow_duplicate=True)],
        Input(f"{app_id}-ketcher-get-smiles", "n_clicks"),
        prevent_initial_call=True
    )
    
    # Load existing SMILES into Ketcher when modal opens (auto-import)
    clientside_callback(
        """
        function(is_open, current_smiles) {
            if (is_open && current_smiles && current_smiles.trim() !== '') {
                console.log('Modal opened with SMILES, auto-importing:', current_smiles);
                
                // Function to convert SMILES to molfile and import
                function attemptImport(retryCount = 0) {
                    const maxRetries = 5;
                    const iframe = document.getElementById('""" + app_id + """-ketcher-iframe');
                    
                    if (iframe && iframe.contentWindow && iframe.contentWindow.ketcher) {
                        console.log('Ketcher available, converting SMILES to molfile...');
                        
                        // Convert SMILES to molfile via backend (async operation, don't return promise)
                        fetch('/api/smiles-to-molfile', {
                            method: 'POST',
                            headers: {
                                'Content-Type': 'application/json',
                            },
                            body: JSON.stringify({ smiles: current_smiles })
                        })
                        .then(response => response.json())
                        .then(data => {
                            if (data.molblock) {
                                console.log('Got molblock, importing into Ketcher...');
                                return iframe.contentWindow.ketcher.setMolecule(data.molblock);
                            } else {
                                console.error('Error converting SMILES:', data.error);
                                throw new Error(data.error);
                            }
                        })
                        .then(() => {
                            console.log('Successfully auto-imported SMILES into Ketcher');
                        })
                        .catch(error => {
                            console.error('Error auto-importing SMILES:', error);
                        });
                    } else if (retryCount < maxRetries) {
                        console.log('Ketcher not yet available, retrying in', (retryCount + 1) * 300, 'ms...');
                        setTimeout(function() {
                            attemptImport(retryCount + 1);
                        }, (retryCount + 1) * 300);  // Progressive delay: 300ms, 600ms, 900ms, etc.
                    } else {
                        console.error('Failed to auto-import SMILES after maximum retries');
                    }
                }
                
                // Start the import attempt after a small initial delay
                setTimeout(function() {
                    attemptImport();
                }, 200);
            } else if (is_open) {
                console.log('Modal opened but no SMILES to import');
            }
            return window.dash_clientside.no_update;
        }
        """,
        Output(f"{app_id}-ketcher-iframe", "src", allow_duplicate=True),  # Dummy output
        Input(f"{app_id}-ketcher-modal", "is_open"),
        State(f"{app_id}-molecule-input", "value"),
        prevent_initial_call=True
    )
    
    # Clear Ketcher editor
    clientside_callback(
        """
        function(clear_clicks) {
            if (clear_clicks > 0) {
                try {
                    const iframe = document.getElementById('""" + app_id + """-ketcher-iframe');
                    if (iframe && iframe.contentWindow && iframe.contentWindow.ketcher) {
                        iframe.contentWindow.ketcher.setMolecule('').then(function() {
                            console.log('Successfully cleared Ketcher editor');
                            return '';
                        }).catch(function(error) {
                            console.error('Error clearing Ketcher:', error);
                            return 'Error clearing sketcher.';
                        });
                    } else {
                        console.error('Ketcher iframe not accessible for clearing');
                        return 'Ketcher not properly loaded.';
                    }
                } catch (error) {
                    console.error('Error clearing Ketcher:', error);
                    return 'Error clearing sketcher.';
                }
            }
            return window.dash_clientside.no_update;
        }
        """,
        Output(f"{app_id}-ketcher-error", "children", allow_duplicate=True),
        Input(f"{app_id}-ketcher-clear", "n_clicks"),
        prevent_initial_call=True
    )
