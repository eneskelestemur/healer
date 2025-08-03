"""
Integration tests for the HEALER web application components.
Tests integration between layouts without testing the deprecated EnumeratorApp.
"""

import pytest
import os
import sys
from pathlib import Path

# Add the project root to the path
sys.path.insert(0, str(Path(__file__).parent.parent))

import dash
from dash import Dash
import dash_bootstrap_components as dbc

# We skip testing the deprecated EnumeratorApp - focus on the new modular components
from webserver.layouts.molecule_layout import get_molecule_layout
from webserver.layouts.site_layout import get_site_layout


class TestLayoutIntegration:
    """Test the integration between different layout components."""
    
    def test_layout_components_exist(self):
        """Test that all layout components can be imported and created."""
        # Test molecule layout can be created
        molecule_layout = get_molecule_layout()
        assert molecule_layout is not None
        
        # Test site layout can be created
        site_layout = get_site_layout()
        assert site_layout is not None
    
    def test_layout_components_are_dash_components(self):
        """Test that layouts return valid Dash components."""
        molecule_layout = get_molecule_layout()
        site_layout = get_site_layout()
        
        # Both should have Dash component attributes
        assert hasattr(molecule_layout, 'children') or hasattr(molecule_layout, 'id')
        assert hasattr(site_layout, 'children') or hasattr(site_layout, 'id')
    
    def test_layouts_can_be_combined(self):
        """Test that layouts can be combined in a Dash app structure."""
        molecule_layout = get_molecule_layout()
        site_layout = get_site_layout()
        
        # Should be able to create a combined layout
        combined_layout = dbc.Container([
            molecule_layout,
            site_layout
        ])
        
        assert combined_layout is not None
        assert hasattr(combined_layout, 'children')
        assert len(combined_layout.children) == 2
    
    def test_bootstrap_theme_compatibility(self):
        """Test that layouts are compatible with Bootstrap themes."""
        # Create a test app with Bootstrap theme
        app = Dash(__name__, external_stylesheets=[dbc.themes.MINTY])
        
        molecule_layout = get_molecule_layout()
        site_layout = get_site_layout()
        
        # Should be able to set layouts as app layout
        app.layout = dbc.Container([molecule_layout, site_layout])
        
        assert app.layout is not None
    
    def test_layout_ids_are_unique(self):
        """Test that layout components have unique IDs to avoid conflicts."""
        molecule_layout = get_molecule_layout()
        site_layout = get_site_layout()
        
        # Collect all IDs from both layouts
        molecule_ids = set()
        site_ids = set()
        
        def collect_ids(component, id_set):
            """Recursively collect all component IDs."""
            if hasattr(component, 'id') and component.id:
                id_set.add(component.id)
            if hasattr(component, 'children'):
                if isinstance(component.children, list):
                    for child in component.children:
                        if hasattr(child, 'id') or hasattr(child, 'children'):
                            collect_ids(child, id_set)
                elif hasattr(component.children, 'id') or hasattr(component.children, 'children'):
                    collect_ids(component.children, id_set)
        
        collect_ids(molecule_layout, molecule_ids)
        collect_ids(site_layout, site_ids)
        
        # Check for ID conflicts
        conflicts = molecule_ids.intersection(site_ids)
        assert len(conflicts) == 0, f"Found conflicting IDs between layouts: {conflicts}"


class TestEnvironmentIntegration:
    """Test integration with environment variables and configuration."""
    
    def test_server_mode_detection(self):
        """Test that server mode can be detected from environment."""
        # Test default (local) mode
        os.environ.pop('HEALER_SERVER_MODE', None)
        server_mode = os.getenv('HEALER_SERVER_MODE', 'false').lower() == 'true'
        assert server_mode == False
        
        # Test server mode
        os.environ['HEALER_SERVER_MODE'] = 'true'
        server_mode = os.getenv('HEALER_SERVER_MODE', 'false').lower() == 'true'
        assert server_mode == True
        
        # Cleanup
        os.environ.pop('HEALER_SERVER_MODE', None)
    
    def test_cache_directory_availability(self):
        """Test that cache directory is available for the app."""
        cache_dir = Path(__file__).parent.parent / "cache"
        
        # Cache directory should exist or be creatable
        if not cache_dir.exists():
            # Try to create it
            cache_dir.mkdir(exist_ok=True)
        
        assert cache_dir.exists()
        assert cache_dir.is_dir()


if __name__ == "__main__":
    pytest.main([__file__])
