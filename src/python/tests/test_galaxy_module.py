"""Comprehensive tests for galaxy_module"""

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import unittest
import tempfile
import shutil
from pathlib import Path
import logging
import numpy as np

from core.galaxy_module import GalaxyModel, ModelParameters
from core.constants import NP, NPGRID


class TestModelParameters(unittest.TestCase):
    """Test ModelParameters dataclass"""
    
    def test_default_initialization(self):
        """Test default initialization"""
        params = ModelParameters()
        self.assertEqual(params.name, "default")
        self.assertEqual(params.sf_mode, 0)
        self.assertEqual(params.total_mass, 1.0)
        self.assertEqual(params.sf_rate, 1.0)
        self.assertEqual(params.num_intervals, 1)
        self.assertEqual(params.exponents, [2.35])
        self.assertEqual(params.mass_limits, [1.0, 100.0])
        self.assertEqual(params.sn_cutoff, 8.0)
        self.assertEqual(params.bh_cutoff, 120.0)
        self.assertEqual(params.metallicity_id, 0)
        self.assertEqual(params.wind_id, 0)
        
    def test_custom_initialization(self):
        """Test custom initialization"""
        params = ModelParameters(
            name="Test Model",
            sf_mode=1,
            total_mass=1e6,
            sf_rate=10.0,
            num_intervals=2,
            exponents=[2.0, 2.7],
            mass_limits=[0.5, 10.0, 150.0],
            sn_cutoff=10.0,
            bh_cutoff=150.0,
            metallicity_id=14,
            wind_id=1
        )
        
        self.assertEqual(params.name, "Test Model")
        self.assertEqual(params.sf_mode, 1)
        self.assertEqual(params.total_mass, 1e6)
        self.assertEqual(params.sf_rate, 10.0)
        self.assertEqual(params.num_intervals, 2)
        self.assertEqual(params.exponents, [2.0, 2.7])
        self.assertEqual(params.mass_limits, [0.5, 10.0, 150.0])
        self.assertEqual(params.sn_cutoff, 10.0)
        self.assertEqual(params.bh_cutoff, 150.0)
        self.assertEqual(params.metallicity_id, 14)
        self.assertEqual(params.wind_id, 1)
        
    def test_time_parameters(self):
        """Test time-related parameters"""
        params = ModelParameters()
        self.assertEqual(params.time_grid, [])
        self.assertEqual(params.time_steps, 0)
        self.assertEqual(params.max_time, 0.0)
        
        # Set time parameters
        params.time_grid = [0.0, 1e6, 1e7]
        params.time_steps = 3
        params.max_time = 1e7
        
        self.assertEqual(params.time_grid, [0.0, 1e6, 1e7])
        self.assertEqual(params.time_steps, 3)
        self.assertEqual(params.max_time, 1e7)
        
    def test_spectral_parameters(self):
        """Test spectral parameters"""
        params = ModelParameters()
        self.assertEqual(params.atmosphere_models, 0)
        self.assertEqual(params.spectral_library, 0)
        self.assertEqual(params.include_lines, 0)
        self.assertEqual(params.velocity_threshold, 0)
        self.assertEqual(params.include_rsg, 0)


class TestGalaxyModel(unittest.TestCase):
    """Test GalaxyModel class comprehensively"""
    
    def setUp(self):
        """Set up test GalaxyModel instance"""
        self.galaxy = GalaxyModel()
        self.temp_dir = tempfile.mkdtemp()
        
    def tearDown(self):
        """Clean up temporary directory"""
        shutil.rmtree(self.temp_dir)
        
    def test_initialization(self):
        """Test GalaxyModel initialization"""
        self.assertIsInstance(self.galaxy.model_params, ModelParameters)
        self.assertEqual(self.galaxy.current_time, 0.0)
        self.assertEqual(self.galaxy.time_step, 0)
        self.assertEqual(self.galaxy.iteration, 0)
        self.assertEqual(len(self.galaxy.wavelength), NP)
        self.assertEqual(len(self.galaxy.spectra), NP)
        self.assertEqual(len(self.galaxy.mass_grid), NPGRID)
        self.assertEqual(len(self.galaxy.density_grid), NPGRID)
        
    def test_paths(self):
        """Test path configurations"""
        self.assertEqual(self.galaxy.data_dir, Path("data"))
        self.assertEqual(self.galaxy.output_dir, Path("output"))
        
    def test_init_module(self):
        """Test module initialization"""
        # Capture log output
        with self.assertLogs('core.galaxy_module', level='INFO') as cm:
            self.galaxy.init_module()
        
        self.assertIn('Galaxy module initialized', cm.output[0])
        
    def test_cleanup_module(self):
        """Test module cleanup"""
        with self.assertLogs('core.galaxy_module', level='INFO') as cm:
            self.galaxy.cleanup_module()
        
        self.assertIn('Galaxy module cleanup', cm.output[0])
        
    def test_open_file_read_success(self):
        """Test successful file opening for reading"""
        # Create a test file
        test_file = Path(self.temp_dir) / "test.txt"
        test_file.write_text("test content")
        
        # Open the file
        file_handle = self.galaxy.open_file(str(test_file), 'r')
        self.assertIsNotNone(file_handle)
        
        # Read content
        content = file_handle.read()
        self.assertEqual(content, "test content")
        file_handle.close()
        
    def test_open_file_read_failure(self):
        """Test file opening failure for non-existent file"""
        with self.assertLogs('core.galaxy_module', level='ERROR') as cm:
            result = self.galaxy.open_file("nonexistent.txt", 'r')
        
        self.assertIsNone(result)
        self.assertIn('File not found', cm.output[0])
        
    def test_open_file_write_success(self):
        """Test successful file opening for writing"""
        test_file = Path(self.temp_dir) / "write_test.txt"
        
        file_handle = self.galaxy.open_file(str(test_file), 'w')
        self.assertIsNotNone(file_handle)
        
        # Write content
        file_handle.write("test write")
        file_handle.close()
        
        # Verify content
        self.assertEqual(test_file.read_text(), "test write")
        
    def test_open_file_exception(self):
        """Test file opening with exception"""
        # Try to open a directory as a file
        with self.assertLogs('core.galaxy_module', level='ERROR') as cm:
            result = self.galaxy.open_file(self.temp_dir, 'r')
        
        self.assertIsNone(result)
        self.assertIn('Error opening file', cm.output[0])
        
    def test_placeholder_methods(self):
        """Test placeholder methods don't crash"""
        # These methods are not implemented yet but should not raise errors
        self.galaxy.read_parameters("dummy.txt")
        self.galaxy.read_tracks()
        self.galaxy.compute_sed()
        self.galaxy.write_output()
        
    def test_arrays_are_zero_initialized(self):
        """Test that arrays are properly zero-initialized"""
        np.testing.assert_array_equal(self.galaxy.wavelength, np.zeros(NP))
        np.testing.assert_array_equal(self.galaxy.spectra, np.zeros(NP))
        np.testing.assert_array_equal(self.galaxy.mass_grid, np.zeros(NPGRID))
        np.testing.assert_array_equal(self.galaxy.density_grid, np.zeros(NPGRID))
        
    def test_model_params_modification(self):
        """Test modifying model parameters"""
        self.galaxy.model_params.name = "Modified Model"
        self.galaxy.model_params.sf_mode = 2
        self.galaxy.model_params.total_mass = 2e6
        
        self.assertEqual(self.galaxy.model_params.name, "Modified Model")
        self.assertEqual(self.galaxy.model_params.sf_mode, 2)
        self.assertEqual(self.galaxy.model_params.total_mass, 2e6)
        
    def test_time_evolution(self):
        """Test time evolution parameters"""
        self.galaxy.current_time = 1e6
        self.galaxy.time_step = 100
        self.galaxy.iteration = 50
        
        self.assertEqual(self.galaxy.current_time, 1e6)
        self.assertEqual(self.galaxy.time_step, 100)
        self.assertEqual(self.galaxy.iteration, 50)


if __name__ == '__main__':
    unittest.main()