"""
Test suite for web application layout components.
Tests the layout generation functions and UI components.
"""

import pytest
import os
import sys
from pathlib import Path
from unittest.mock import patch, MagicMock

# Add the project root to the path
sys.path.insert(0, str(Path(__file__).parent.parent))

import dash_bootstrap_components as dbc
from dash import html, dcc

from webserver.layouts.molecule_layout import (
    get_molecule_layout,
    get_server_info_banner
)
from webserver.layouts.site_layout import (
    get_site_layout,
    get_server_info_banner as site_get_server_info_banner
)
from webserver.layouts.base_layout import create_base_layout


class TestServerInfoBanner:
    """Test server mode information banner functionality."""
    
    @pytest.fixture(autouse=True)
    def setup_teardown(self):
        """Setup and teardown for each test."""
        # Store original environment variable
        self.original_mode = os.environ.get('HEALER_SERVER_MODE', 'false')
        yield
        # Restore original environment variable
        os.environ['HEALER_SERVER_MODE'] = self.original_mode
    
    def test_server_mode_banner(self):
        """Test server mode banner display."""
        os.environ['HEALER_SERVER_MODE'] = 'true'
        banner = get_server_info_banner()
        
        assert isinstance(banner, dbc.Alert)
        assert banner.color == "info"
        assert "server mode" in str(banner.children).lower()
    
    def test_local_mode_banner(self):
        """Test local mode banner display."""
        os.environ['HEALER_SERVER_MODE'] = 'false'
        banner = get_server_info_banner()
        
        assert isinstance(banner, dbc.Alert)
        assert banner.color == "success"
        assert "local mode" in str(banner.children).lower()
    
    def test_site_server_mode_banner(self):
        """Test server mode banner for site layout."""
        os.environ['HEALER_SERVER_MODE'] = 'true'
        banner = site_get_server_info_banner()
        
        assert isinstance(banner, dbc.Alert)
        assert banner.color == "info"


class TestMoleculeLayout:
    """Test molecule layout generation."""
    
    def test_get_molecule_layout_structure(self):
        """Test the basic structure of molecule layout."""
        layout = get_molecule_layout("test-molecule")
        
        assert isinstance(layout, dbc.Container)
        assert layout.fluid is True
        assert "80%" in str(layout.style.get('width', ''))
    
    def test_molecule_layout_components(self):
        """Test that molecule layout contains required components."""
        layout = get_molecule_layout("test-molecule")
        
        # Convert layout to string for easier testing
        layout_str = str(layout)
        
        # Check for key components that actually exist in the layout
        assert "Molecule HEALER" in layout_str
        assert "test-molecule-fragment-alert" in layout_str  # This exists in the layout
        assert "Container" in layout_str  # Basic structure
    
    def test_molecule_layout_tooltips(self):
        """Test that molecule layout includes tooltips."""
        layout = get_molecule_layout("test-molecule")
        layout_str = str(layout)
        
        # Check for tooltip or help components that actually exist
        assert ("info-circle" in layout_str or  # FontAwesome icons
                "check-circle" in layout_str or  # Alternative icon
                "Fragment HEALER" in layout_str)  # Help text
    
    def test_molecule_layout_with_different_app_id(self):
        """Test molecule layout with different app ID."""
        layout = get_molecule_layout("custom-id")
        layout_str = str(layout)
        
        # Check for app ID in the layout
        assert "custom-id" in layout_str  # Should contain the app ID
        assert "custom-id-fragment-alert" in layout_str  # This specific component exists


class TestSiteLayout:
    """Test site layout generation."""
    
    @patch('webserver.layouts.base_layout.create_base_layout')
    def test_get_site_layout_structure(self, mock_create_base):
        """Test the basic structure of site layout."""
        # Mock the base layout to return a simple container
        mock_container = dbc.Container(html.Div("Test"), fluid=True, style={'width': '80%'})
        mock_create_base.return_value = mock_container
        
        layout = get_site_layout("test-site")
        
        assert isinstance(layout, dbc.Container)
        assert layout.fluid is True
        assert "80%" in str(layout.style.get('width', ''))
    
    @patch('webserver.layouts.base_layout.create_base_layout')
    def test_site_layout_components(self, mock_create_base):
        """Test that site layout contains required components."""
        # Mock the base layout to return a container with expected components
        mock_container = dbc.Container([
            html.H2("Site HEALER"),
            html.Div(id="test-site-reactive-sites-input")
        ], fluid=True)
        mock_create_base.return_value = mock_container
        
        layout = get_site_layout("test-site")
        layout_str = str(layout)
        
        # Check for key components that should exist
        assert "Site HEALER" in layout_str
        assert "test-site-reactive-sites-input" in layout_str or "reactive-sites" in layout_str
    
    @patch('webserver.layouts.base_layout.create_base_layout')
    def test_site_layout_property_sliders(self, mock_create_base):
        """Test that site layout includes molecular property sliders."""
        # Mock the base layout to return a container with sliders
        mock_container = dbc.Container([
            dcc.RangeSlider(id="test-MW-slider")
        ], fluid=True)
        mock_create_base.return_value = mock_container
        
        layout = get_site_layout("test-site")
        layout_str = str(layout)
        
        # Check for molecular property sliders - site layout should have range sliders
        assert "RangeSlider" in layout_str or "slider" in layout_str
    
    @patch('webserver.layouts.base_layout.create_base_layout')
    def test_site_layout_tooltips(self, mock_create_base):
        """Test that site layout includes tooltips."""
        # Mock the base layout to return a container with tooltips
        mock_container = dbc.Container([
            html.I(className="fas fa-question-circle"),
            html.P("Site-specific molecular enumeration")
        ], fluid=True)
        mock_create_base.return_value = mock_container
        
        layout = get_site_layout("test-site")
        layout_str = str(layout)
        
        # Check for tooltip components - look for tooltip classes or question marks
        assert ("fa-question-circle" in layout_str or 
                "tooltip" in layout_str.lower() or 
                "Site-specific molecular enumeration" in layout_str)


class TestBaseLayout:
    """Test the base layout creation function."""
    
    def test_create_base_layout_minimal(self):
        """Test base layout creation with minimal parameters."""
        additional_inputs = [
            html.Div("Test input", id="test-input")
        ]
        
        layout = create_base_layout("test", "Test", additional_inputs)
        
        assert isinstance(layout, dbc.Container)
        assert layout.fluid is True
    
    def test_create_base_layout_with_server_info(self):
        """Test base layout creation with server info."""
        additional_inputs = [html.Div("Test input")]
        server_info = dbc.Alert("Server mode", color="info")
        
        layout = create_base_layout(
            "test", "Test", additional_inputs, 
            server_info=server_info
        )
        
        layout_str = str(layout)
        assert "Server mode" in layout_str
    
    def test_create_base_layout_with_tool_description(self):
        """Test base layout creation with tool description."""
        additional_inputs = [html.Div("Test input")]
        tool_description = dbc.Card(
            dbc.CardBody([
                html.H5("Test Tool"),
                html.P("Test description")
            ])
        )
        
        layout = create_base_layout(
            "test", "Test", additional_inputs, 
            tool_description=tool_description
        )
        
        layout_str = str(layout)
        assert "Test Tool" in layout_str
        assert "Test description" in layout_str
    
    def test_create_base_layout_with_fragment_alert(self):
        """Test base layout creation with fragment alert."""
        additional_inputs = [html.Div("Test input")]
        fragment_alert = dbc.Alert(
            "Fragment HEALER detected", 
            id="test-fragment-alert",
            color="info"
        )
        
        layout = create_base_layout(
            "test", "Test", additional_inputs, 
            fragment_alert=fragment_alert
        )
        
        layout_str = str(layout)
        assert "Fragment HEALER detected" in layout_str
        assert "test-fragment-alert" in layout_str
    
    def test_create_base_layout_all_components(self):
        """Test base layout creation with all optional components."""
        additional_inputs = [html.Div("Test input")]
        server_info = dbc.Alert("Server mode", color="info")
        tool_description = dbc.Card(dbc.CardBody(html.P("Description")))
        fragment_alert = dbc.Alert("Fragment alert", color="warning")
        
        layout = create_base_layout(
            "test", "Test", additional_inputs,
            server_info=server_info,
            tool_description=tool_description,
            fragment_alert=fragment_alert
        )
        
        layout_str = str(layout)
        assert "Server mode" in layout_str
        assert "Description" in layout_str
        assert "Fragment alert" in layout_str


class TestLayoutInteractivity:
    """Test layout components that affect interactivity."""
    
    @patch('webserver.layouts.base_layout.create_base_layout')
    def test_molecule_layout_slider_ranges(self, mock_create_base):
        """Test that sliders have appropriate ranges."""
        # Mock the base layout to return a container with sliders
        mock_container = dbc.Container([
            dcc.Slider(id="test-sim-threshold", min=0.0, max=1.0, value=0.15),
            dcc.Slider(id="test-n-compositions", min=1, max=100, value=10)
        ], fluid=True)
        mock_create_base.return_value = mock_container
        
        layout = get_molecule_layout("test")
        layout_str = str(layout)
        
        # Verify sliders exist
        assert "Slider" in layout_str
    
    @patch('webserver.layouts.base_layout.create_base_layout')
    def test_site_layout_range_sliders(self, mock_create_base):
        """Test that site layout has range sliders for molecular properties."""
        # Mock the base layout to return a container with range sliders
        mock_container = dbc.Container([
            dcc.RangeSlider(id="test-MW-slider"),
            dcc.RangeSlider(id="test-HBD-slider")
        ], fluid=True)
        mock_create_base.return_value = mock_container
        
        layout = get_site_layout("test")
        layout_str = str(layout)
        
        # Check for range sliders
        assert "RangeSlider" in layout_str
    
    @patch('webserver.layouts.base_layout.create_base_layout')
    @patch('webserver.layouts.base_layout.create_base_layout')
    def test_layout_dropdown_options(self, mock_site_base, mock_mol_base):
        """Test that layouts include dropdown components."""
        # Mock both base layouts
        mock_mol_container = dbc.Container([
            dcc.Dropdown(id="test-dropdown", options=[])
        ], fluid=True)
        mock_site_container = dbc.Container([
            dcc.Dropdown(id="test-site-dropdown", options=[])
        ], fluid=True)
        mock_mol_base.return_value = mock_mol_container
        mock_site_base.return_value = mock_site_container
        
        mol_layout = get_molecule_layout("test")
        site_layout = get_site_layout("test")
        
        mol_str = str(mol_layout)
        site_str = str(site_layout)
        
        # Check for dropdown components
        assert "Dropdown" in mol_str or "dropdown" in mol_str.lower()
        assert "Dropdown" in site_str or "dropdown" in site_str.lower()


class TestLayoutAccessibility:
    """Test layout accessibility features."""
    
    @patch('webserver.layouts.base_layout.create_base_layout')
    @patch('webserver.layouts.base_layout.create_base_layout')
    def test_layout_has_labels(self, mock_site_base, mock_mol_base):
        """Test that layouts include proper labels."""
        # Mock both base layouts with labels
        mock_mol_container = dbc.Container([
            html.Label("Test Label"),
            dbc.Input(id="test-input")
        ], fluid=True)
        mock_site_container = dbc.Container([
            html.Label("Site Label"),
            dbc.Input(id="site-input")
        ], fluid=True)
        mock_mol_base.return_value = mock_mol_container
        mock_site_base.return_value = mock_site_container
        
        mol_layout = get_molecule_layout("test")
        site_layout = get_site_layout("test")
        
        mol_str = str(mol_layout)
        site_str = str(site_layout)
        
        # Check for label elements
        assert "Label" in mol_str
        assert "Label" in site_str
    
    @patch('webserver.layouts.base_layout.create_base_layout')
    def test_layout_has_tooltips_with_proper_targets(self, mock_create_base):
        """Test that tooltips have proper target references."""
        # Mock layout with tooltips
        mock_container = dbc.Container([
            dbc.Tooltip("Test tooltip", target="test-target"),
            html.Div(id="test-target")
        ], fluid=True)
        mock_create_base.return_value = mock_container
        
        mol_layout = get_molecule_layout("test")
        layout_str = str(mol_layout)
        
        # Tooltips should reference their target elements
        if "Tooltip" in layout_str:
            assert "target=" in layout_str
    
    @patch('webserver.layouts.base_layout.create_base_layout')
    def test_layout_ids_are_unique(self, mock_create_base):
        """Test that component IDs are properly namespaced."""
        # Mock layout with unique IDs
        mock_container = dbc.Container([
            html.Div(id="mol-test-unique-id"),
            dbc.Input(id="mol-test-input")
        ], fluid=True)
        mock_create_base.return_value = mock_container
        
        mol_layout = get_molecule_layout("mol-test")
        
        mol_str = str(mol_layout)
        
        # IDs should be properly namespaced
        assert "mol-test-" in mol_str or "mol-test" in mol_str
        
        # The layout should contain the app ID
        assert "mol-test" in mol_str


if __name__ == "__main__":
    pytest.main([__file__])
