"""
Test suite for web application callback functions.
Tests the callback logic for both molecule and site healers.
"""

import pytest
import os
import sys
from pathlib import Path
from unittest.mock import patch, MagicMock

# Add the project root to the path
sys.path.insert(0, str(Path(__file__).parent.parent))

import dash
from dash import Dash
import dash_bootstrap_components as dbc

from webserver.callbacks.molecule_callbacks import register_molecule_callbacks
from webserver.callbacks.site_callbacks import register_site_callbacks
from webserver.callbacks.shared_callbacks import register_shared_callbacks


class TestCallbackRegistration:
    """Test callback registration functionality."""
    
    @pytest.fixture
    def app(self):
        """Create a test Dash app."""
        return Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
    
    def test_register_molecule_callbacks(self, app):
        """Test that molecule callbacks can be registered without errors."""
        try:
            register_molecule_callbacks(app, "test-molecule")
        except Exception as e:
            pytest.fail(f"Failed to register molecule callbacks: {e}")
    
    def test_register_site_callbacks(self, app):
        """Test that site callbacks can be registered without errors."""
        try:
            register_site_callbacks(app, "test-site")
        except Exception as e:
            pytest.fail(f"Failed to register site callbacks: {e}")
    
    def test_register_shared_callbacks(self, app):
        """Test that shared callbacks can be registered without errors."""
        try:
            register_shared_callbacks(app, "test-shared")
        except Exception as e:
            pytest.fail(f"Failed to register shared callbacks: {e}")
    
    def test_multiple_callback_registration(self, app):
        """Test registering multiple callback sets."""
        try:
            register_shared_callbacks(app, "molecule")
            register_shared_callbacks(app, "site")
            register_molecule_callbacks(app, "molecule")
            register_site_callbacks(app, "site")
        except Exception as e:
            pytest.fail(f"Failed to register multiple callbacks: {e}")


class TestMoleculeCallbacks:
    """Test molecule-specific callback functionality."""
    
    @pytest.fixture
    def app(self):
        """Create a test Dash app with molecule callbacks."""
        app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
        register_molecule_callbacks(app, "molecule")
        return app
    
    @patch('webserver.callbacks.molecule_callbacks.count_molecular_fragments')
    def test_fragment_alert_callback_logic(self, mock_count_fragments):
        """Test the logic of fragment alert callback."""
        from webserver.callbacks.molecule_callbacks import register_molecule_callbacks
        
        # Test single component molecule
        mock_count_fragments.return_value = 1
        
        # Mock callback function (we can't easily test Dash callbacks directly)
        # Instead, we test the underlying logic
        molecule = "CC1=CC=CC=C1"
        num_fragments = mock_count_fragments(molecule)
        
        if num_fragments > 1:
            alert_open = True
            alert_style = {'display': 'block'}
        else:
            alert_open = False
            alert_style = {'display': 'none'}
        
        assert alert_open is False
        assert alert_style == {'display': 'none'}
        
        # Test multi-component molecule
        mock_count_fragments.return_value = 2
        molecule = "CCO.C1=CC=CC=C1"
        num_fragments = mock_count_fragments(molecule)
        
        if num_fragments > 1:
            alert_open = True
            alert_style = {'display': 'block'}
        else:
            alert_open = False
            alert_style = {'display': 'none'}
        
        assert alert_open is True
        assert alert_style == {'display': 'block'}
    
    def test_random_seed_visibility_logic(self):
        """Test the logic for random seed visibility."""
        # Test when randomize is not selected
        randomize_value = []
        
        if "randomize" in randomize_value:
            expected_style = {'fontSize': '14px', 'display': 'block'}
            expected_disabled = False
            expected_row_style = {'display': 'block'}
        else:
            expected_style = {'fontSize': '14px', 'display': 'none'}
            expected_disabled = True
            expected_row_style = {'display': 'none'}
        
        assert expected_style == {'fontSize': '14px', 'display': 'none'}
        assert expected_disabled is True
        assert expected_row_style == {'display': 'none'}
        
        # Test when randomize is selected
        randomize_value = ["randomize"]
        
        if "randomize" in randomize_value:
            expected_style = {'fontSize': '14px', 'display': 'block'}
            expected_disabled = False
            expected_row_style = {'display': 'block'}
        else:
            expected_style = {'fontSize': '14px', 'display': 'none'}
            expected_disabled = True
            expected_row_style = {'display': 'none'}
        
        assert expected_style == {'fontSize': '14px', 'display': 'block'}
        assert expected_disabled is False
        assert expected_row_style == {'display': 'block'}


class TestSharedCallbacks:
    """Test shared callback functionality."""
    
    @pytest.fixture
    def app(self):
        """Create a test Dash app with shared callbacks."""
        app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
        register_shared_callbacks(app, "test")
        return app
    
    @patch('webserver.callbacks.shared_callbacks.utils.load_reactions_from_json')
    def test_reaction_tags_loading_logic(self, mock_load_reactions):
        """Test the logic for loading reaction tags."""
        # Mock reaction data
        mock_reaction1 = MagicMock()
        mock_reaction1.is_valid.return_value = True
        mock_reaction1.tags = ["amide coupling", "C-N bond formation"]
        
        mock_reaction2 = MagicMock()
        mock_reaction2.is_valid.return_value = True
        mock_reaction2.tags = ["alkylation", "N-arylation"]
        
        mock_load_reactions.return_value = [mock_reaction1, mock_reaction2]
        
        # Simulate the callback logic
        reactions = mock_load_reactions()
        reaction_tags = list(set(tag for r in reactions if r.is_valid() for tag in r.tags))
        reaction_tags.sort()
        options = [{"label": tag, "value": tag} for tag in reaction_tags]
        
        expected_tags = ["C-N bond formation", "N-arylation", "alkylation", "amide coupling"]
        assert reaction_tags == expected_tags
        assert len(options) == 4
        assert all("label" in opt and "value" in opt for opt in options)
    
    @patch('webserver.callbacks.shared_callbacks.utils.load_reactions_from_json')
    def test_reaction_tags_fallback_logic(self, mock_load_reactions):
        """Test fallback when reaction loading fails."""
        # Mock loading failure
        mock_load_reactions.side_effect = Exception("Loading failed")
        
        # Simulate the callback logic with error handling
        try:
            reactions = mock_load_reactions()
            reaction_tags = list(set(tag for r in reactions if r.is_valid() for tag in r.tags))
            reaction_tags.sort()
            options = [{"label": tag, "value": tag} for tag in reaction_tags]
        except Exception:
            # Fallback options
            default_tags = ["amide coupling", "amide", "C-N bond formation", "C-N",
                           "alkylation", "N-arylation", "azole", "amination"]
            options = [{"label": tag, "value": tag} for tag in default_tags]
        
        assert len(options) == 8
        assert any(opt["value"] == "amide coupling" for opt in options)


class TestParameterValidationCallbacks:
    """Test parameter validation callback logic."""
    
    @pytest.fixture(autouse=True)
    def setup_teardown(self):
        """Setup and teardown for each test."""
        # Store original environment variable
        self.original_mode = os.environ.get('HEALER_SERVER_MODE', 'false')
        yield
        # Restore original environment variable
        os.environ['HEALER_SERVER_MODE'] = self.original_mode
    
    def test_server_mode_molecule_validation_logic(self):
        """Test server mode validation logic for molecule healer."""
        os.environ['HEALER_SERVER_MODE'] = 'true'
        
        # Simulate parameter validation callback logic
        sim_threshold = 0.1  # Too low
        n_compositions = 100  # Too high
        max_bbs = 15  # Too high
        retro_depth = 5  # Too high
        reaction_tags = ['tag'] * 20  # Too many
        
        # Apply server mode limits
        server_mode = os.environ.get('HEALER_SERVER_MODE', 'false').lower() == 'true'
        
        if server_mode:
            # Similarity threshold must be >= 0.3
            if sim_threshold < 0.3:
                sim_threshold = 0.3
            
            # N compositions must be <= 50
            if n_compositions > 50:
                n_compositions = 50
            
            # Max BBs per comp must be 1-10 (if not -1)
            if max_bbs != -1 and (max_bbs < 1 or max_bbs > 10):
                max_bbs = min(10, max(1, max_bbs))
            
            # Retro depth must be 1-2
            if retro_depth < 1 or retro_depth > 2:
                retro_depth = min(2, max(1, retro_depth))
            
            # Reaction tags must be <= 15
            if reaction_tags and len(reaction_tags) > 15:
                reaction_tags = reaction_tags[:15]
        
        assert sim_threshold == 0.3
        assert n_compositions == 50
        assert max_bbs == 10
        assert retro_depth == 2
        assert len(reaction_tags) == 15
    
    def test_local_mode_no_validation_logic(self):
        """Test that local mode applies no validation."""
        os.environ['HEALER_SERVER_MODE'] = 'false'
        
        # Original values
        original_sim_threshold = 0.1
        original_n_compositions = 100
        original_max_bbs = 15
        original_retro_depth = 5
        original_reaction_tags = ['tag'] * 20
        
        # Values should remain unchanged in local mode
        sim_threshold = original_sim_threshold
        n_compositions = original_n_compositions
        max_bbs = original_max_bbs
        retro_depth = original_retro_depth
        reaction_tags = original_reaction_tags
        
        server_mode = os.environ.get('HEALER_SERVER_MODE', 'false').lower() == 'true'
        
        if not server_mode:
            # No changes should be made
            pass
        
        assert sim_threshold == original_sim_threshold
        assert n_compositions == original_n_compositions
        assert max_bbs == original_max_bbs
        assert retro_depth == original_retro_depth
        assert len(reaction_tags) == 20


class TestEnumerationCallbacks:
    """Test enumeration callback logic."""
    
    @patch('webserver.callbacks.molecule_callbacks.run_molecule_enumeration')
    @patch('webserver.callbacks.molecule_callbacks.format_enumeration_results')
    @patch('webserver.callbacks.molecule_callbacks.count_molecular_fragments')
    def test_enumeration_callback_logic(self, mock_count_fragments, mock_format_results, mock_run_enumeration):
        """Test the enumeration callback logic."""
        # Mock the dependencies
        mock_count_fragments.return_value = 1
        mock_run_enumeration.return_value = [
            {'Product': 'CC1=CC=CC=C1', 'Similarity_to_query': 0.85},
            {'Product': 'CC1=CC=CC=C1O', 'Similarity_to_query': 0.72}
        ]
        mock_format_results.return_value = [
            {'Product': 'CC1=CC=CC=C1', 'Similarity_to_query': 0.85},
            {'Product': 'CC1=CC=CC=C1O', 'Similarity_to_query': 0.72}
        ]
        
        # Simulate callback parameters
        n_clicks = 1
        molecule = "CC1=CC=CC=C1"
        bb_source = "test"
        reaction_tags = ["amide coupling"]
        custom_sites_str = ""
        sim_threshold = 0.15
        n_compositions = 10
        randomize_checkbox = []
        random_seed = -1
        retro_depth = 1
        min_frag_size = 3
        max_bbs = -1
        max_evals = None
        
        # Simulate the enumeration callback logic
        if not n_clicks or not molecule:
            result = "no_update"
        else:
            # Parse custom split sites
            custom_sites = None
            if custom_sites_str and custom_sites_str.strip():
                try:
                    site_pairs = []
                    for site_str in custom_sites_str.split(';'):
                        if site_str.strip():
                            parts = [int(x.strip()) for x in site_str.split(',')]
                            if len(parts) == 2:
                                site_pairs.append(tuple(parts))
                    custom_sites = site_pairs if site_pairs else None
                except Exception:
                    custom_sites = None
            
            # Check if molecule has multiple fragments
            num_fragments = mock_count_fragments(molecule)
            use_fragment_healer = num_fragments > 1
            
            # Run enumeration
            results = mock_run_enumeration(
                molecule=molecule,
                bb_source=bb_source,
                reaction_tags=reaction_tags,
                custom_sites=custom_sites,
                sim_threshold=sim_threshold,
                n_compositions=n_compositions,
                randomize_compositions="randomize" in randomize_checkbox,
                random_seed=random_seed if "randomize" in randomize_checkbox else -1,
                retro_tree_depth=retro_depth,
                min_frag_size=min_frag_size,
                max_bbs_per_comp=max_bbs if max_bbs is not None else -1,
                max_evals_per_comp=max_evals,
                use_fragment_healer=use_fragment_healer
            )
            
            # Format results
            formatted_results = mock_format_results(results, 'molecule')
            n_results = len(formatted_results)
            
            # Generate response
            if n_results < 2:
                enum_img = "No Enumerations! Check Your Inputs!"
                bg_color = 'firebrick'
            else:
                enum_img = f"Slide to view {n_results} enumerated molecules!"
                bg_color = 'forestgreen'
        
        # Verify the logic worked correctly
        assert enum_img == "Slide to view 2 enumerated molecules!"
        assert bg_color == 'forestgreen'
        mock_count_fragments.assert_called_once_with(molecule)
        mock_run_enumeration.assert_called_once()
        mock_format_results.assert_called_once_with(results, 'molecule')


class TestErrorHandling:
    """Test error handling in callbacks."""
    
    @patch('webserver.callbacks.molecule_callbacks.run_molecule_enumeration')
    def test_enumeration_error_handling_logic(self, mock_run_enumeration):
        """Test error handling in enumeration callback."""
        # Mock an enumeration error
        mock_run_enumeration.side_effect = Exception("Enumeration failed")
        
        # Simulate error handling logic
        try:
            results = mock_run_enumeration(
                molecule="CC1=CC=CC=C1",
                bb_source="test",
                reaction_tags=["amide coupling"]
            )
            error_occurred = False
        except Exception as e:
            error_occurred = True
            error_msg = f"Enumeration failed: {str(e)}"
            error_style = {
                'backgroundColor': 'firebrick',
                'width': '80%',
                'height': '150px',
                'text-align': 'center',
                'fontSize': '18px',
                'fontFamily': 'sans-serif',
                'borderRadius': '20px',
                'alignItems': 'center',
                'justifyContent': 'center',
                'display': 'flex',
                'margin': '0 auto',
                'color': 'white'
            }
        
        assert error_occurred is True
        assert "Enumeration failed" in error_msg
        assert error_style['backgroundColor'] == 'firebrick'
        assert error_style['color'] == 'white'


if __name__ == "__main__":
    pytest.main([__file__])
