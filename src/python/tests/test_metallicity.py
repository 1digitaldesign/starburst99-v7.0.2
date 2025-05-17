"""Tests for different metallicity cases (Z020, Z034, Z040, Z100, Z200)"""

import unittest
import tempfile
import shutil
import json
from pathlib import Path

from ..starburst_main import Starburst99
from ..core.galaxy_module import GalaxyModel, ModelParameters
from ..file_io.input_parser import InputParser
from ..file_io.output_writer import OutputWriter


class TestMetallicityInputOutput(unittest.TestCase):
    """Test input/output for different metallicities"""
    
    def setUp(self):
        """Set up test environment"""
        self.temp_dir = tempfile.mkdtemp()
        self.test_data_dir = Path(self.temp_dir) / "data"
        self.test_data_dir.mkdir()
        
        # Create test track files for different metallicities
        self.create_test_track_files()
        
        # Create test atmosphere files
        self.create_test_atmosphere_files()
        
    def tearDown(self):
        """Clean up test environment"""
        shutil.rmtree(self.temp_dir)
        
    def create_test_track_files(self):
        """Create dummy track files for testing"""
        tracks_dir = self.test_data_dir / "tracks"
        tracks_dir.mkdir()
        
        # Define test metallicities and corresponding files
        track_files = {
            "Z0020v00.txt": "Z=0.020 v=0.0",
            "Z0020v40.txt": "Z=0.020 v=40.0",
            "Z0034v00.txt": "Z=0.034 v=0.0", 
            "Z0040v00.txt": "Z=0.040 v=0.0",
            "Z0100v00.txt": "Z=0.100 v=0.0",
            "Z0200v00.txt": "Z=0.200 v=0.0",
        }
        
        # Create dummy track files
        for filename, content in track_files.items():
            track_file = tracks_dir / filename
            track_file.write_text(f"# Track file: {content}\n")
            
    def create_test_atmosphere_files(self):
        """Create dummy atmosphere files for testing"""
        lejeune_dir = self.test_data_dir / "lejeune"
        lejeune_dir.mkdir()
        
        # Create atmosphere files for different metallicities
        atm_files = ["lcb97_m13.flu", "lcb97_m07.flu", "lcb97_m04.flu", 
                    "lcb97_p00.flu", "lcb97_p03.flu"]
        
        for filename in atm_files:
            atm_file = lejeune_dir / filename
            atm_file.write_text(f"# Atmosphere file: {filename}\n")
            
    def create_test_input_z020(self):
        """Create test input for Z=0.020 metallicity"""
        return {
            "name": "Test Z020 Model",
            "star_formation": {
                "mode": 1,
                "total_mass": 1.0e6,
                "rate": 10.0
            },
            "imf": {
                "num_intervals": 1,
                "exponents": [2.35],
                "mass_limits": [1.0, 100.0],
                "sn_cutoff": 8.0,
                "bh_cutoff": 120.0
            },
            "model": {
                "metallicity_id": 14,  # Z=0.020
                "wind_id": 0
            },
            "time": {
                "steps": 100,
                "min": 1.0e6,
                "max": 1.0e8
            }
        }
        
    def create_test_input_z034(self):
        """Create test input for Z=0.034 metallicity"""
        input_data = self.create_test_input_z020()
        input_data["name"] = "Test Z034 Model"
        input_data["model"]["metallicity_id"] = 24  # Assuming this maps to Z=0.034
        return input_data
        
    def create_test_input_z040(self):
        """Create test input for Z=0.040 metallicity"""
        input_data = self.create_test_input_z020()
        input_data["name"] = "Test Z040 Model"
        input_data["model"]["metallicity_id"] = 34  # Assuming this maps to Z=0.040
        return input_data
        
    def create_test_input_z100(self):
        """Create test input for Z=0.100 metallicity"""
        input_data = self.create_test_input_z020()
        input_data["name"] = "Test Z100 Model"
        input_data["model"]["metallicity_id"] = 44  # Assuming this maps to Z=0.100
        return input_data
        
    def create_test_input_z200(self):
        """Create test input for Z=0.200 metallicity"""
        input_data = self.create_test_input_z020()
        input_data["name"] = "Test Z200 Model"
        input_data["model"]["metallicity_id"] = 54  # Assuming this maps to Z=0.200
        return input_data
        
    def test_z020_input_output(self):
        """Test Z=0.020 metallicity input/output"""
        # Create input file
        input_data = self.create_test_input_z020()
        input_file = Path(self.temp_dir) / "test_z020.json"
        with open(input_file, 'w') as f:
            json.dump(input_data, f)
            
        # Parse input
        parser = InputParser()
        params = parser.read_input(str(input_file))
        
        # Check parameters
        self.assertEqual(params.name, "Test Z020 Model")
        self.assertEqual(params.metallicity_id, 14)
        
        # Test output writer
        galaxy = GalaxyModel()
        galaxy.model_params = params
        galaxy.data_dir = self.test_data_dir
        galaxy.output_dir = Path(self.temp_dir) / "output"
        galaxy.output_dir.mkdir()
        
        writer = OutputWriter()
        writer.output_dir = galaxy.output_dir
        writer.write_all_outputs(galaxy)
        
        # Check output files exist
        output_files = list(galaxy.output_dir.glob("test_z020_model_*"))
        self.assertGreater(len(output_files), 0)
        
    def test_z034_input_output(self):
        """Test Z=0.034 metallicity input/output"""
        input_data = self.create_test_input_z034()
        input_file = Path(self.temp_dir) / "test_z034.json"
        with open(input_file, 'w') as f:
            json.dump(input_data, f)
            
        parser = InputParser()
        params = parser.read_input(str(input_file))
        
        self.assertEqual(params.name, "Test Z034 Model")
        self.assertEqual(params.metallicity_id, 24)
        
    def test_z040_input_output(self):
        """Test Z=0.040 metallicity input/output"""
        input_data = self.create_test_input_z040()
        input_file = Path(self.temp_dir) / "test_z040.json"
        with open(input_file, 'w') as f:
            json.dump(input_data, f)
            
        parser = InputParser()
        params = parser.read_input(str(input_file))
        
        self.assertEqual(params.name, "Test Z040 Model")
        self.assertEqual(params.metallicity_id, 34)
        
    def test_z100_input_output(self):
        """Test Z=0.100 metallicity input/output"""
        input_data = self.create_test_input_z100()
        input_file = Path(self.temp_dir) / "test_z100.json"
        with open(input_file, 'w') as f:
            json.dump(input_data, f)
            
        parser = InputParser()
        params = parser.read_input(str(input_file))
        
        self.assertEqual(params.name, "Test Z100 Model")
        self.assertEqual(params.metallicity_id, 44)
        
    def test_z200_input_output(self):
        """Test Z=0.200 metallicity input/output"""
        input_data = self.create_test_input_z200()
        input_file = Path(self.temp_dir) / "test_z200.json"
        with open(input_file, 'w') as f:
            json.dump(input_data, f)
            
        parser = InputParser()
        params = parser.read_input(str(input_file))
        
        self.assertEqual(params.name, "Test Z200 Model")
        self.assertEqual(params.metallicity_id, 54)
        
    def test_metallicity_file_mapping(self):
        """Test correct file mapping for different metallicities"""
        starburst = Starburst99()
        starburst.galaxy.data_dir = self.test_data_dir
        
        # Test Z=0.020
        starburst.galaxy.model_params.metallicity_id = 14
        track_file = starburst._get_track_filename()
        self.assertTrue(track_file.name.startswith("Z"))
        
        # Test metallicity string mapping
        starburst._set_metallicity_string()
        self.assertEqual(starburst.namfi3, 'p00')
        self.assertEqual(starburst.nam, '020')
        
    def test_standard_format_input(self):
        """Test standard format input for different metallicities"""
        standard_input = """Test Z020 Standard
1 1000000.0 10.0
1
2.35 1.0
100.0
8.0 120.0
14 0
100 1001
1000000.0 100000000.0 1.0
"""
        
        input_file = Path(self.temp_dir) / "test_standard.input"
        input_file.write_text(standard_input)
        
        parser = InputParser()
        params = parser.read_input(str(input_file))
        
        self.assertEqual(params.name, "Test Z020 Standard")
        self.assertEqual(params.metallicity_id, 14)
        self.assertEqual(params.sf_mode, 1)
        self.assertEqual(params.total_mass, 1000000.0)
        
    def test_ini_format_input(self):
        """Test INI format input"""
        ini_content = """[general]
name = Test Z020 INI

[star_formation]
mode = 1
total_mass = 1000000.0
rate = 10.0

[imf]
num_intervals = 1
exponents = 2.35
mass_limits = 1.0,100.0
sn_cutoff = 8.0
bh_cutoff = 120.0

[model]
metallicity_id = 14
wind_id = 0
"""
        
        input_file = Path(self.temp_dir) / "test.ini"
        input_file.write_text(ini_content)
        
        parser = InputParser()
        params = parser.read_input(str(input_file))
        
        self.assertEqual(params.name, "Test Z020 INI")
        self.assertEqual(params.metallicity_id, 14)


class TestMetallicityCalculations(unittest.TestCase):
    """Test calculations for different metallicities"""
    
    def test_track_interpolation_z020(self):
        """Test track interpolation for Z=0.020"""
        # This would test actual track interpolation if implemented
        pass
        
    def test_sed_calculation_different_z(self):
        """Test SED calculation varies with metallicity"""
        # This would test that different metallicities produce different SEDs
        pass
        
    def test_output_format_consistency(self):
        """Test output format is consistent across metallicities"""
        # This would verify output format remains consistent
        pass


if __name__ == '__main__':
    unittest.main()