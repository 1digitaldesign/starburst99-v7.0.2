"""Comprehensive tests for input_parser module"""

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import unittest
import tempfile
import json
import configparser
from pathlib import Path

from file_io.input_parser import InputParser
from core.galaxy_module import ModelParameters


class TestInputParser(unittest.TestCase):
    """Test InputParser class comprehensively"""
    
    def setUp(self):
        """Set up test environment"""
        self.parser = InputParser()
        self.temp_dir = tempfile.mkdtemp()
        
    def tearDown(self):
        """Clean up test environment"""
        import shutil
        shutil.rmtree(self.temp_dir)
        
    def test_file_not_found(self):
        """Test handling of non-existent file"""
        with self.assertRaises(FileNotFoundError):
            self.parser.read_input("nonexistent_file.txt")
            
    def test_json_input(self):
        """Test JSON format input parsing"""
        json_data = {
            "name": "Test JSON Model",
            "star_formation": {
                "mode": 1,
                "total_mass": 1e6,
                "rate": 10.0
            },
            "imf": {
                "num_intervals": 2,
                "exponents": [2.0, 2.7],
                "mass_limits": [0.5, 10.0, 150.0],
                "sn_cutoff": 10.0,
                "bh_cutoff": 150.0
            },
            "model": {
                "metallicity_id": 24,
                "wind_id": 1
            }
        }
        
        json_file = Path(self.temp_dir) / "test.json"
        with open(json_file, 'w') as f:
            json.dump(json_data, f)
            
        params = self.parser.read_input(str(json_file))
        
        self.assertEqual(params.name, "Test JSON Model")
        self.assertEqual(params.sf_mode, 1)
        self.assertEqual(params.total_mass, 1e6)
        self.assertEqual(params.sf_rate, 10.0)
        self.assertEqual(params.num_intervals, 2)
        self.assertEqual(params.exponents, [2.0, 2.7])
        self.assertEqual(params.mass_limits, [0.5, 10.0, 150.0])
        self.assertEqual(params.sn_cutoff, 10.0)
        self.assertEqual(params.bh_cutoff, 150.0)
        self.assertEqual(params.metallicity_id, 24)
        self.assertEqual(params.wind_id, 1)
        
    def test_json_input_defaults(self):
        """Test JSON input with missing fields uses defaults"""
        json_data = {"name": "Minimal JSON"}
        
        json_file = Path(self.temp_dir) / "minimal.json"
        with open(json_file, 'w') as f:
            json.dump(json_data, f)
            
        params = self.parser.read_input(str(json_file))
        
        self.assertEqual(params.name, "Minimal JSON")
        self.assertEqual(params.sf_mode, 0)  # Default
        self.assertEqual(params.total_mass, 1.0)  # Default
        self.assertEqual(params.exponents, [2.35])  # Default
        
    def test_ini_input(self):
        """Test INI format input parsing"""
        ini_content = """[general]
name = Test INI Model

[star_formation]
mode = 2
total_mass = 5e5
rate = 5.0

[imf]
num_intervals = 1
exponents = 2.5
mass_limits = 2.0,80.0
sn_cutoff = 9.0
bh_cutoff = 100.0

[model]
metallicity_id = 15
wind_id = 2
"""
        
        ini_file = Path(self.temp_dir) / "test.ini"
        ini_file.write_text(ini_content)
        
        params = self.parser.read_input(str(ini_file))
        
        self.assertEqual(params.name, "Test INI Model")
        self.assertEqual(params.sf_mode, 2)
        self.assertEqual(params.total_mass, 5e5)
        self.assertEqual(params.sf_rate, 5.0)
        
    def test_ini_input_defaults(self):
        """Test INI input with missing sections uses defaults"""
        ini_content = """[general]
name = Minimal INI
"""
        
        ini_file = Path(self.temp_dir) / "minimal.ini"
        ini_file.write_text(ini_content)
        
        params = self.parser.read_input(str(ini_file))
        
        self.assertEqual(params.name, "Minimal INI")
        self.assertEqual(params.sf_mode, 0)  # Default
        
    def test_standard_input_continuous_sf(self):
        """Test standard format input with continuous star formation"""
        standard_content = """Test Standard Model
1 1000000.0 10.0
2
2.0 0.5
2.7 10.0
150.0
8.0 120.0
24 1
100 1001
1000000.0 100000000.0 1.0
"""
        
        standard_file = Path(self.temp_dir) / "test.input"
        standard_file.write_text(standard_content)
        
        params = self.parser.read_input(str(standard_file))
        
        self.assertEqual(params.name, "Test Standard Model")
        self.assertEqual(params.sf_mode, 1)
        self.assertEqual(params.total_mass, 1000000.0)
        self.assertEqual(params.sf_rate, 10.0)
        self.assertEqual(params.num_intervals, 2)
        self.assertEqual(params.exponents, [2.0, 2.7])
        self.assertEqual(params.mass_limits, [0.5, 10.0, 150.0])
        self.assertEqual(params.sn_cutoff, 8.0)
        self.assertEqual(params.bh_cutoff, 120.0)
        self.assertEqual(params.metallicity_id, 24)
        self.assertEqual(params.wind_id, 1)
        
    def test_standard_input_instantaneous_sf(self):
        """Test standard format input with instantaneous star formation"""
        standard_content = """Test Instantaneous
0 1000000.0
1
2.35 1.0
100.0
8.0 120.0
14 0
50 501
1000000.0 50000000.0 1.0
"""
        
        standard_file = Path(self.temp_dir) / "test_inst.input"
        standard_file.write_text(standard_content)
        
        params = self.parser.read_input(str(standard_file))
        
        self.assertEqual(params.name, "Test Instantaneous")
        self.assertEqual(params.sf_mode, 0)
        self.assertEqual(params.total_mass, 1000000.0)
        # No SF rate for instantaneous mode
        self.assertEqual(params.num_intervals, 1)
        
    def test_get_default_parameters(self):
        """Test default parameters generation"""
        params = self.parser.get_default_parameters()
        
        self.assertEqual(params.name, "Default Model")
        self.assertEqual(params.sf_mode, 1)
        self.assertEqual(params.total_mass, 1.0)
        self.assertEqual(params.sf_rate, 1.0)
        self.assertEqual(params.num_intervals, 1)
        self.assertEqual(params.exponents, [2.35])
        self.assertEqual(params.mass_limits, [1.0, 100.0])
        self.assertEqual(params.sn_cutoff, 8.0)
        self.assertEqual(params.bh_cutoff, 120.0)
        self.assertEqual(params.metallicity_id, 24)
        self.assertEqual(params.wind_id, 0)
        
    def test_validate_parameters_valid(self):
        """Test parameter validation with valid parameters"""
        params = self.parser.get_default_parameters()
        self.assertTrue(self.parser.validate_parameters(params))
        
    def test_validate_parameters_invalid_mass(self):
        """Test parameter validation with invalid mass"""
        params = self.parser.get_default_parameters()
        params.total_mass = -1.0
        
        with self.assertLogs('file_io.input_parser', level='ERROR') as cm:
            result = self.parser.validate_parameters(params)
        
        self.assertFalse(result)
        self.assertIn("Total mass must be positive", cm.output[0])
        
    def test_validate_parameters_invalid_sfr(self):
        """Test parameter validation with invalid star formation rate"""
        params = self.parser.get_default_parameters()
        params.sf_mode = 1
        params.sf_rate = -5.0
        
        with self.assertLogs('file_io.input_parser', level='ERROR') as cm:
            result = self.parser.validate_parameters(params)
        
        self.assertFalse(result)
        self.assertIn("Star formation rate must be positive", cm.output[0])
        
    def test_validate_parameters_imf_mismatch(self):
        """Test parameter validation with IMF parameter mismatch"""
        params = self.parser.get_default_parameters()
        params.num_intervals = 2
        params.exponents = [2.35]  # Only one exponent for 2 intervals
        
        with self.assertLogs('file_io.input_parser', level='ERROR') as cm:
            result = self.parser.validate_parameters(params)
        
        self.assertFalse(result)
        self.assertIn("doesn't match number of exponents", cm.output[0])
        
    def test_validate_parameters_mass_limits_mismatch(self):
        """Test parameter validation with mass limits mismatch"""
        params = self.parser.get_default_parameters()
        params.num_intervals = 1
        params.mass_limits = [1.0]  # Should be 2 limits for 1 interval
        
        with self.assertLogs('file_io.input_parser', level='ERROR') as cm:
            result = self.parser.validate_parameters(params)
        
        self.assertFalse(result)
        self.assertIn("Incorrect number of mass limits", cm.output[0])
        
    def test_file_format_detection(self):
        """Test correct file format detection"""
        # JSON file
        json_file = Path(self.temp_dir) / "test.json"
        json_file.write_text("{}")
        params = self.parser.read_input(str(json_file))
        self.assertIsInstance(params, ModelParameters)
        
        # INI file
        ini_file = Path(self.temp_dir) / "test.ini"
        ini_file.write_text("[general]\nname=Test")
        params = self.parser.read_input(str(ini_file))
        self.assertIsInstance(params, ModelParameters)
        
        # Standard format (other extension)
        std_file = Path(self.temp_dir) / "test.input"
        std_file.write_text("Test\n1 1.0 1.0\n1\n2.35 1.0\n100.0\n8.0 120.0\n24 0\n")
        params = self.parser.read_input(str(std_file))
        self.assertIsInstance(params, ModelParameters)


if __name__ == '__main__':
    unittest.main()