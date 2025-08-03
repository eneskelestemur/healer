# HEALER Web Application

A modern web interface for the HEALER molecular enumeration system, built with Dash and Bootstrap.

## Overview

This web application provides an intuitive interface for two main HEALER workflows:

1. **Molecule HEALER**: Enumerates molecules by splitting them into fragments and applying reactions
2. **Site HEALER**: Enumerates molecules by applying reactions at specified reactive sites

## Architecture

The application is organized into modular components:

```
webserver/
├── app.py                      # Main Dash application
├── layouts/                    # UI layout components
│   ├── molecule_layout.py      # MoleculeHEALER interface
│   └── site_layout.py          # SiteHEALER interface
├── callbacks/                  # Callback functions
│   ├── shared_callbacks.py     # Common callbacks
│   ├── molecule_callbacks.py   # MoleculeHEALER specific
│   └── site_callbacks.py       # SiteHEALER specific
├── utils/                      # Utility functions
│   └── healer_interface.py     # HEALER API wrappers
└── run_app.sh                 # Startup script
```

## Features

### Core Features
- **Automatic Fragment Detection**: Auto-switches to FragmentHEALER for multi-component molecules
- **Server/Local Mode**: Server mode enforces computational limits for shared deployment
- **Enhanced UI**: Bootstrap Minty theme with comprehensive tooltips and parameter guidance
- **Parameter Validation**: Real-time validation with server-side enforcement
- Interactive molecule visualization with highlighted building blocks
- Building block source selection (test, US stock, EU stock, Global stock)
- Reaction tag filtering
- Real-time enumeration with progress indicators
- Results export to CSV
- Background processing with cancellation support

### Molecule HEALER Specific
- Custom split site specification
- Number of compositions control
- Retrosynthesis tree depth setting
- Minimum fragment size control
- Composition randomization options
- Maximum building blocks per composition

### Site HEALER Specific
- Reactive site specification
- Molecular property filters (MW, HBD, HBA, TPSA, etc.)
- Structural rule filters (SMARTS patterns)
- Range sliders for property constraints

## Usage

### Starting the Application

```bash
# Navigate to webserver directory
cd /path/to/healer/webserver

# Local Mode (default - unlimited parameters)
python app.py
# or set explicitly
export HEALER_SERVER_MODE=false && python app.py

# Server Mode (enforced limits for shared deployment)
export HEALER_SERVER_MODE=true && python app.py
```

The application will be available at: http://localhost:8053

### Using Molecule HEALER

1. Enter a molecule SMILES string
2. Select building block source and reaction tags
3. Optionally specify custom split sites (e.g., "9,10; 6,9")
4. Adjust enumeration parameters:
   - Similarity threshold
   - Number of compositions
   - Retro tree depth
   - Minimum fragment size
5. Click "Enumerate" to start the process
6. Use the slider to browse results
7. Save results to CSV if needed

### Using Site HEALER

1. Enter a molecule SMILES string
2. Select building block source and reaction tags
3. Specify reactive sites (atom indices, e.g., "19, 20, 21")
4. Set molecular property filters using range sliders
5. Optionally add structural rules (SMARTS patterns)
6. Click "Enumerate" to start the process
7. Browse results and save as needed

## API Integration

The web app interfaces with the HEALER classes through wrapper functions in `utils/healer_interface.py`:

- `run_molecule_enumeration()`: Handles MoleculeHEALER workflow
- `run_site_enumeration()`: Handles SiteHEALER workflow
- `format_enumeration_results()`: Standardizes result format
- `generate_molecule_visualization()`: Creates SVG visualizations

## Dependencies

- dash
- dash-bootstrap-components
- diskcache
- pandas
- rdkit
- healer (main package)

## Development

### Adding New Features

1. **Layout changes**: Modify files in `layouts/`
2. **New callbacks**: Add to appropriate files in `callbacks/`
3. **API changes**: Update `utils/healer_interface.py`

### Debugging

- Set `debug=True` in `app.py` for development mode
- Check browser console for JavaScript errors
- Monitor terminal for Python exceptions

## Migration Notes

This application replaces the legacy `EnumeratorApp.py` with:

- Modular architecture for better maintainability
- Updated API calls to new HEALER classes
- Enhanced UI with additional parameter controls
- Improved error handling and user feedback
- Better code organization and documentation

## Troubleshooting

### Common Issues

1. **Reaction data not loading**: Check path to reactions.json file
2. **Building block errors**: Verify BB source file paths
3. **Enumeration failures**: Check molecule SMILES validity
4. **Visualization errors**: Ensure RDKit is properly installed

### Performance Tips

- Use "test" BB source for quick testing (100 BBs only)
- Limit number of compositions for faster enumeration
- Set max evaluations to prevent long-running jobs
- Use background processing for large enumerations
