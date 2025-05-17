"""Comprehensive tests for input_parser module"""

import unittest
import tempfile
import json
import configparser
from pathlib import Path

from ..file_io.input_parser import InputParser
from ..core.galaxy_module import ModelParameters


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
        
        with self.assertLogs(level='ERROR') as cm:
            result = self.parser.validate_parameters(params)
        
        self.assertFalse(result)
        found = False
        for msg in cm.output:
            if "Total mass must be positive" in msg:
                found = True
                break
        self.assertTrue(found, f"Expected message not found in: {cm.output}")
        
    def test_validate_parameters_invalid_sfr(self):
        """Test parameter validation with invalid star formation rate"""
        params = self.parser.get_default_parameters()
        params.sf_mode = 1
        params.sf_rate = -5.0
        
        with self.assertLogs(level='ERROR') as cm:
            result = self.parser.validate_parameters(params)
        
        self.assertFalse(result)
        found = False
        for msg in cm.output:
            if "Star formation rate must be positive" in msg:
                found = True
                break
        self.assertTrue(found, f"Expected message not found in: {cm.output}")
        
    def test_validate_parameters_imf_mismatch(self):
        """Test parameter validation with IMF parameter mismatch"""
        params = self.parser.get_default_parameters()
        params.num_intervals = 2
        params.exponents = [2.35]  # Only one exponent for 2 intervals
        
        with self.assertLogs(level='ERROR') as cm:
            result = self.parser.validate_parameters(params)
        
        self.assertFalse(result)
        found = False
        for msg in cm.output:
            if "doesn't match number of exponents" in msg:
                found = True
                break
        self.assertTrue(found, f"Expected message not found in: {cm.output}")
        
    def test_validate_parameters_mass_limits_mismatch(self):
        """Test parameter validation with mass limits mismatch"""
        params = self.parser.get_default_parameters()
        params.num_intervals = 1
        params.mass_limits = [1.0]  # Should be 2 limits for 1 interval
        
        with self.assertLogs(level='ERROR') as cm:
            result = self.parser.validate_parameters(params)
        
        self.assertFalse(result)
        found = False
        for msg in cm.output:
            if "Incorrect number of mass limits" in msg:
                found = True
                break
        self.assertTrue(found, f"Expected message not found in: {cm.output}")
        
    def test_file_format_detection(self):
        """Test correct file format detection"""
        # JSON file
        json_file = Path(self.temp_dir) / "test.json"
        json_file.write_text("{}")
        params = self.parser.read_input(str(json_file))
        self.assertEqual(type(params).__name__, 'ModelParameters')
        
        # INI file
        ini_file = Path(self.temp_dir) / "test.ini"
        ini_file.write_text("[general]\nname=Test")
        params = self.parser.read_input(str(ini_file))
        self.assertEqual(type(params).__name__, 'ModelParameters')
        
        # Standard format (other extension)
        std_file = Path(self.temp_dir) / "test.input"
        # Create minimal valid standard format
        content = """MODEL DESIGNATION:                                           [NAME]
Test
CONTINUOUS STAR FORMATION (>0) OR FIXED MASS (<=0):          [ISF]
-1
TOTAL STELLAR MASS [1E6 M_SOL] IF 'FIXED MASS' IS CHOSEN:    [TOMA]
1.0
SFR [SOLAR MASSES PER YEAR] IF 'CONT. SF' IS CHOSEN:         [SFR]
1.0
NUMBER OF INTERVALS FOR THE IMF (KROUPA=2):                  [NINTERV]
1
IMF EXPONENTS (KROUPA=1.3,2.3):                              [XPONENT]
2.35
MASS BOUNDARIES FOR IMF (KROUPA=0.1,0.5,100) [SOLAR MASSES]: [XMASLIM]
1.0,100.0
SUPERNOVA CUT-OFF MASS [SOLAR MASSES]:                       [SNCUT]
8.0
BLACK HOLE CUT-OFF MASS [SOLAR MASSES]:                      [BHCUT]
120.0
METALLICITY + TRACKS:                                        [IZ]
24
WIND FLAG (0=STANDARD MASS LOSS): no standard tracks at=0.05 [IWIND]
0
TIME INTERVAL FOR WHICH SPECTRA ARE SYNTHESIZED [MYR]:       [TIME1, TIME2, LOGT]
1.,1000.,0.1
COMPUTE SYNTHETIC MAGNITUDES + OUTPUT FILES (>0=yes, <=0=no):[IO1,IO2...IO4]
+1, +1, +1, -1
"""
        std_file.write_text(content)
        try:
            params = self.parser.read_input(str(std_file))
            self.assertEqual(type(params).__name__, 'ModelParameters')
        except Exception as e:
            self.fail(f"Failed to parse standard format: {e}")


if __name__ == '__main__':
    unittest.main()