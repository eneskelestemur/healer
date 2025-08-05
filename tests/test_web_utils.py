"""
Test suite for web application utility functions.
Tests the healer_interface.py functions including fragment detection,
parameter validation, and HEALER creation.
"""

import pytest
import os
import sys
from pathlib import Path
from unittest.mock import patch, MagicMock

# Add the project root to the path
sys.path.insert(0, str(Path(__file__).parent.parent))

from webserver.utils.healer_interface import (
    count_molecular_fragments,
    validate_server_parameters,
    create_molecule_healer,
    create_site_healer,
    run_molecule_enumeration,
    run_site_enumeration,
    format_enumeration_results,
    generate_molecule_visualization
)
from healer.application.healer import MoleculeHEALER, SiteHEALER, FragmentHEALER


class TestFragmentDetection:
    """Test fragment detection functionality."""
    
    def test_single_component_molecule(self):
        """Test fragment detection for single component molecule."""
        single_mol = "CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C"
        assert count_molecular_fragments(single_mol) == 1
    
    def test_multi_component_molecule(self):
        """Test fragment detection for multi-component molecule."""
        multi_mol = "CCO.C1=CC=CC=C1"
        assert count_molecular_fragments(multi_mol) == 2
    
    def test_three_component_molecule(self):
        """Test fragment detection for three-component molecule."""
        three_mol = "CCO.C1=CC=CC=C1.CC(C)C"
        assert count_molecular_fragments(three_mol) == 3
    
    def test_invalid_smiles(self):
        """Test fragment detection with invalid SMILES."""
        invalid_mol = "InvalidSMILES"
        assert count_molecular_fragments(invalid_mol) == 0
    
    def test_empty_smiles(self):
        """Test fragment detection with empty SMILES."""
        assert count_molecular_fragments("") == 0
    
    def test_none_input(self):
        """Test fragment detection with None input."""
        # The function should handle None gracefully and return 0
        result = count_molecular_fragments(None)
        assert result == 0


class TestParameterValidation:
    """Test server mode parameter validation."""
    
    @pytest.fixture(autouse=True)
    def setup_teardown(self):
        """Setup and teardown for each test."""
        # Store original environment variable
        self.original_mode = os.environ.get('HEALER_SERVER_MODE', 'false')
        yield
        # Restore original environment variable
        os.environ['HEALER_SERVER_MODE'] = self.original_mode
    
    def test_local_mode_no_limits(self):
        """Test that local mode applies no parameter limits."""
        os.environ['HEALER_SERVER_MODE'] = 'false'
        
        params = {
            'reaction_tags': ['tag'] * 20,
            'sim_threshold': 0.1,
            'n_compositions': 100,
            'max_bbs_per_comp': 20,
            'max_evals_per_comp': 1000,
            'retro_depth': 5
        }
        
        validated = validate_server_parameters(params, "molecule")
        
        # No limits should be applied in local mode
        assert len(validated['reaction_tags']) == 20
        assert validated['sim_threshold'] == 0.1
        assert validated['n_compositions'] == 100
        assert validated['max_bbs_per_comp'] == 20
        assert validated['max_evals_per_comp'] == 1000
        assert validated['retro_depth'] == 5
    
    def test_server_mode_molecule_limits(self):
        """Test that server mode applies limits for molecule healer."""
        os.environ['HEALER_SERVER_MODE'] = 'true'
        
        params = {
            'reaction_tags': ['tag'] * 20,
            'sim_threshold': 0.1,
            'n_compositions': 100,
            'max_bbs_per_comp': 20,
            'max_evals_per_comp': 1000,
            'retro_depth': 5
        }
        
        validated = validate_server_parameters(params, "molecule")
        
        # Server limits should be applied
        assert len(validated['reaction_tags']) == 15
        assert validated['sim_threshold'] == 0.3
        assert validated['n_compositions'] == 50
        assert validated['max_bbs_per_comp'] == 10
        assert validated['max_evals_per_comp'] == 500
        assert validated['retro_depth'] == 2
    
    def test_server_mode_site_limits(self):
        """Test that server mode applies appropriate limits for site healer."""
        os.environ['HEALER_SERVER_MODE'] = 'true'
        
        params = {
            'reaction_tags': ['tag'] * 20,
            'max_evals_per_comp': 1000
        }
        
        validated = validate_server_parameters(params, "site")
        
        # Only reaction tags should be limited for site healer
        assert len(validated['reaction_tags']) == 15
        assert validated['max_evals_per_comp'] == 1000  # No limit for site healer
    
    def test_edge_case_parameters(self):
        """Test edge cases in parameter validation."""
        os.environ['HEALER_SERVER_MODE'] = 'true'
        
        params = {
            'reaction_tags': [],
            'sim_threshold': 0.5,  # Already valid
            'n_compositions': 25,  # Already valid
            'max_bbs_per_comp': -1,  # Special case
            'max_evals_per_comp': None,  # None value
        }
        
        validated = validate_server_parameters(params, "molecule")
        
        assert validated['reaction_tags'] == []
        assert validated['sim_threshold'] == 0.5
        assert validated['n_compositions'] == 25
        assert validated['max_bbs_per_comp'] == 1  # -1 should be converted to 1
        assert validated['max_evals_per_comp'] is None


class TestHealerCreation:
    """Test HEALER instance creation."""
    
    def test_create_molecule_healer(self):
        """Test MoleculeHEALER creation."""
        healer = create_molecule_healer(
            bb_supplier='test',
            reaction_tags=['amide coupling'],
            use_fragment_healer=False
        )
        
        assert isinstance(healer, MoleculeHEALER)
        assert not isinstance(healer, FragmentHEALER)
    
    def test_create_fragment_healer(self):
        """Test FragmentHEALER creation."""
        healer = create_molecule_healer(
            bb_supplier='test',
            reaction_tags=['amide coupling'],
            use_fragment_healer=True
        )
        
        assert isinstance(healer, FragmentHEALER)
    
    def test_create_site_healer(self):
        """Test SiteHEALER creation."""
        healer = create_site_healer(
            bb_supplier='test',
            reaction_tags=['amide coupling']
        )
        
        assert isinstance(healer, SiteHEALER)
    
    def test_healer_with_custom_parameters(self):
        """Test HEALER creation with custom parameters."""
        healer = create_molecule_healer(
            bb_supplier='test',
            reaction_tags=['amide coupling', 'C-N bond formation'],
            sim_threshold=0.25,
            max_bbs_per_comp=5,
            max_evals_per_comp=10,
            n_compositions=5
        )
        
        assert isinstance(healer, MoleculeHEALER)
        assert healer.sim_threshold == 0.25
        assert healer.max_bbs_per_comp == 5
    
    def test_site_healer_with_rules(self):
        """Test SiteHEALER creation with custom rules."""
        custom_rules = {
            'MW': (100, 400),
            'HBD': (0, 3),
            'HBA': (0, 8)
        }
        
        healer = create_site_healer(
            bb_supplier='test',
            reaction_tags=['amide coupling'],
            rules=custom_rules
        )
        
        assert isinstance(healer, SiteHEALER)
        assert healer.rules['MW'] == (100, 400)
        assert healer.rules['HBD'] == (0, 3)


class TestEnumerationFunctions:
    """Test enumeration execution functions."""
    
    @pytest.fixture
    def mock_healer(self):
        """Create a mock healer for testing."""
        mock = MagicMock()
        mock.get_results.return_value = [
            {
                'Product': 'CC1=CC=CC=C1',
                'BB1': 'CCO',
                'Reaction1_name': 'amide coupling',
                'Similarity_to_query': 0.85
            },
            {
                'Product': 'CC1=CC=CC=C1O',
                'BB1': 'CCO',
                'BB2': 'C1=CC=CC=C1',
                'Reaction1_name': 'amide coupling',
                'Reaction2_name': 'C-N bond formation',
                'Similarity_to_query': 0.72
            }
        ]
        return mock
    
    @patch('webserver.utils.healer_interface.create_molecule_healer')
    def test_molecule_enumeration_auto_fragment_detection(self, mock_create, mock_healer):
        """Test automatic fragment detection in molecule enumeration."""
        mock_create.return_value = mock_healer
        
        # Test with multi-component molecule (should auto-switch to FragmentHEALER)
        multi_mol = "CCO.C1=CC=CC=C1"
        
        results = run_molecule_enumeration(
            molecule=multi_mol,
            bb_supplier='test',
            reaction_tags=['amide coupling'],
            max_evals_per_comp=2
        )
        
        # Verify FragmentHEALER was requested
        mock_create.assert_called_once()
        call_args = mock_create.call_args[1]
        assert call_args['use_fragment_healer'] is True
    
    @patch('webserver.utils.healer_interface.create_molecule_healer')
    def test_molecule_enumeration_single_component(self, mock_create, mock_healer):
        """Test molecule enumeration with single component molecule."""
        mock_create.return_value = mock_healer
        
        # Test with single component molecule (should use MoleculeHEALER)
        single_mol = "CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C"
        
        results = run_molecule_enumeration(
            molecule=single_mol,
            bb_supplier='test',
            reaction_tags=['amide coupling'],
            n_compositions=5,
            max_evals_per_comp=2
        )
        
        # Verify MoleculeHEALER was requested
        mock_create.assert_called_once()
        call_args = mock_create.call_args[1]
        assert call_args['use_fragment_healer'] is False
    
    @patch('webserver.utils.healer_interface.create_site_healer')
    def test_site_enumeration(self, mock_create, mock_healer):
        """Test site enumeration execution."""
        mock_create.return_value = mock_healer
        
        results = run_site_enumeration(
            molecule="CC1=CC=CC=C1",
            bb_supplier='test',
            reaction_tags=['amide coupling'],
            reactive_sites=[1, 2, 3],
            max_evals_per_comp=2
        )
        
        mock_create.assert_called_once()
        assert len(results) == 2
    
    def test_enumeration_error_handling(self):
        """Test error handling in enumeration functions."""
        with pytest.raises(Exception):
            run_molecule_enumeration(
                molecule="InvalidSMILES",
                bb_supplier='nonexistent',
                reaction_tags=['invalid_tag']
            )


class TestResultFormatting:
    """Test result formatting functions."""
    
    def test_format_molecule_results(self):
        """Test formatting of molecule enumeration results."""
        raw_results = [
            {
                'Product': 'CC1=CC=CC=C1',
                'BB1': 'CCO',
                'BB2': 'C1=CC=CC=C1',
                'Reaction1_name': 'amide coupling',
                'Reaction2_name': 'C-N bond formation',
                'Similarity_to_query': 0.85
            }
        ]
        
        display_results, complete_results = format_enumeration_results(raw_results, 'molecule')
        
        # Test display results
        assert len(display_results) == 1
        assert display_results[0]['Product'] == 'CC1=CC=CC=C1'
        assert display_results[0]['BB1'] == 'CCO'
        assert display_results[0]['BB2'] == 'C1=CC=CC=C1'
        assert display_results[0]['Reaction_name'] == 'amide coupling -> C-N bond formation'
        assert display_results[0]['Similarity_to_query'] == 0.85
        
        # Test complete results preserve all original data
        assert len(complete_results) == 1
        assert complete_results[0]['Product'] == 'CC1=CC=CC=C1'
        assert complete_results[0]['BB1'] == 'CCO'
        assert complete_results[0]['BB2'] == 'C1=CC=CC=C1'
        assert complete_results[0]['Reaction1_name'] == 'amide coupling'
        assert complete_results[0]['Reaction2_name'] == 'C-N bond formation'
        assert complete_results[0]['Similarity_to_query'] == 0.85
    
    def test_format_site_results(self):
        """Test formatting of site enumeration results."""
        raw_results = [
            {
                'Product': 'CC1=CC=CC=C1O',
                'BB1': 'CCO',
                'Reaction1_name': 'amide coupling',
                'Similarity_to_query': 0.72
            }
        ]
        
        formatted = format_enumeration_results(raw_results, 'site')
        
        display_results, complete_results = format_enumeration_results(raw_results, 'site')
        
        # Test display results
        assert len(display_results) == 1
        assert display_results[0]['Product'] == 'CC1=CC=CC=C1O'
        assert display_results[0]['BB'] == 'CCO'  # Site uses 'BB' not 'BB1' in display
        assert display_results[0]['Reaction_name'] == 'amide coupling'
        assert display_results[0]['Similarity_to_query'] == 0.72
        
        # Test complete results preserve all original data
        assert len(complete_results) == 1
        assert complete_results[0]['Product'] == 'CC1=CC=CC=C1O'
        assert complete_results[0]['BB1'] == 'CCO'  # Original key preserved
        assert complete_results[0]['Reaction1_name'] == 'amide coupling'
        assert complete_results[0]['Similarity_to_query'] == 0.72
    
    def test_format_empty_results(self):
        """Test formatting of empty results."""
        display_results, complete_results = format_enumeration_results([], 'molecule')
        assert display_results == []
        assert complete_results == []
    
    def test_format_results_missing_fields(self):
        """Test formatting with missing fields."""
        raw_results = [
            {
                'Product': 'CC1=CC=CC=C1',
                'Similarity_to_query': 0.85
                # Missing BB and Reaction fields
            }
        ]
        
        display_results, complete_results = format_enumeration_results(raw_results, 'molecule')
        
        # Test display results
        assert len(display_results) == 1
        assert display_results[0]['Product'] == 'CC1=CC=CC=C1'
        assert display_results[0]['Similarity_to_query'] == 0.85
        
        # Test complete results preserve all original data
        assert len(complete_results) == 1
        assert complete_results[0]['Product'] == 'CC1=CC=CC=C1'
        assert complete_results[0]['Similarity_to_query'] == 0.85


class TestVisualization:
    """Test molecule visualization functions."""
    
    @patch('webserver.utils.healer_interface.utils.get_svg_mol_with_bbs')
    def test_generate_visualization_with_bbs(self, mock_svg_with_bbs):
        """Test molecule visualization with building blocks."""
        mock_svg_with_bbs.return_value = "mock_svg_data"
        
        result = generate_molecule_visualization(
            mol_smiles="CC1=CC=CC=C1",
            bb_smiles=["CCO", "C1=CC=CC=C1"],
            legend="Test molecule"
        )
        
        mock_svg_with_bbs.assert_called_once()
        assert result == "mock_svg_data"
    
    @patch('webserver.utils.healer_interface.utils.get_svg_mol')
    def test_generate_visualization_without_bbs(self, mock_svg):
        """Test molecule visualization without building blocks."""
        mock_svg.return_value = "mock_svg_data"
        
        result = generate_molecule_visualization(
            mol_smiles="CC1=CC=CC=C1",
            bb_smiles=None,
            legend="Test molecule"
        )
        
        mock_svg.assert_called_once()
        assert result == "mock_svg_data"
    
    @patch('webserver.utils.healer_interface.utils.get_svg_mol')
    def test_generate_visualization_empty_bbs(self, mock_svg):
        """Test molecule visualization with empty building blocks list."""
        mock_svg.return_value = "mock_svg_data"
        
        result = generate_molecule_visualization(
            mol_smiles="CC1=CC=CC=C1",
            bb_smiles=["", "   "],  # Empty/whitespace only
            legend="Test molecule"
        )
        
        mock_svg.assert_called_once()
        assert result == "mock_svg_data"
    
    def test_generate_visualization_error_handling(self):
        """Test error handling in visualization generation."""
        with pytest.raises(ValueError):
            generate_molecule_visualization(
                mol_smiles="InvalidSMILES",
                bb_smiles=["CCO"]
            )


if __name__ == "__main__":
    pytest.main([__file__])
