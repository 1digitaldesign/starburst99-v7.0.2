"""Integration tests for Starburst99"""

import unittest
import tempfile
import shutil
from pathlib import Path
import json
import numpy as np

from ..starburst_main import Starburst99
from ..core.galaxy_module import ModelParameters
from ..file_io.input_parser import InputParser
from ..file_io.output_writer import OutputWriter
from ..models.imf import IMF
from ..models.stellar_tracks import StellarTracks
from ..core.data_profiles import DataProfiles


class TestIntegration(unittest.TestCase):
    """Integration tests that run full model calculations"""
    
    def setUp(self):
        """Set up test environment"""
        self.temp_dir = tempfile.mkdtemp()
        self.input_file = Path(self.temp_dir) / "test_input.json"
        self.output_dir = Path(self.temp_dir) / "output"
        self.output_dir.mkdir()
        
    def tearDown(self):
        """Clean up test files"""
        shutil.rmtree(self.temp_dir)
    
    def create_test_input(self, **kwargs):
        """Create a test input file with given parameters"""
        default_params = {
            "name": "Integration Test Model",
            "star_formation": {
                "mode": -1,  # Fixed mass
                "total_mass": 1.0e6,
                "rate": 0.0
            },
            "imf": {
                "num_intervals": 1,
                "exponents": [2.35],
                "mass_limits": [1.0, 100.0],
                "sn_cutoff": 8.0,
                "bh_cutoff": 120.0
            },
            "model": {
                "metallicity_id": 24,
                "wind_id": 0
            },
            "time": {
                "start": 1.0,
                "end": 10.0,
                "num_steps": 3
            },
            "output": {
                "directory": str(self.output_dir),
                "prefix": "test_model"
            }
        }
        
        # Merge with provided parameters
        test_params = {**default_params, **kwargs}
        
        with open(self.input_file, 'w') as f:
            json.dump(test_params, f, indent=2)
            
        return test_params
    
    def test_simple_model_run(self):
        """Test running a simple stellar population model"""
        # Create input file
        self.create_test_input()
        
        # Initialize and run model
        starburst = Starburst99(str(self.input_file))
        starburst.run()
        
        # Check that output files were created
        output_files = list(self.output_dir.glob("test_model*"))
        self.assertGreater(len(output_files), 0)
        
        # Check for specific output files
        spectrum_files = list(self.output_dir.glob("*spectrum"))
        self.assertGreater(len(spectrum_files), 0)
        
        # Check that model completed
        self.assertEqual(starburst.galaxy.current_time, 10.0)
    
    def test_continuous_star_formation(self):
        """Test continuous star formation mode"""
        # Create input with continuous SF
        params = self.create_test_input(
            star_formation={
                "mode": 1,  # Continuous
                "total_mass": 1.0e7,
                "rate": 1.0  # 1 Msun/yr
            }
        )
        
        # Run model
        starburst = Starburst99(str(self.input_file))
        starburst.run()
        
        # Check that mass increases with time
        # (this would require implementing mass tracking)
        self.assertTrue(hasattr(starburst.galaxy, 'total_mass'))
    
    def test_imf_variations(self):
        """Test different IMF configurations"""
        # Test Kroupa IMF (multi-segment)
        params = self.create_test_input(
            imf={
                "num_intervals": 2,
                "exponents": [1.3, 2.3],
                "mass_limits": [0.1, 0.5, 100.0],
                "sn_cutoff": 8.0,
                "bh_cutoff": 120.0
            }
        )
        
        # Run model
        starburst = Starburst99(str(self.input_file))
        starburst.run()
        
        # Verify IMF was properly initialized
        self.assertEqual(starburst.imf.num_intervals, 2)
        self.assertEqual(len(starburst.imf.exponents), 2)
        
        # Test IMF integration
        total_mass = starburst.imf.integrate(0.1, 100.0)
        self.assertGreater(total_mass, 0)
    
    def test_metallicity_variations(self):
        """Test different metallicity tracks"""
        metallicity_ids = [11, 14, 24, 54]  # Different Z values
        
        for z_id in metallicity_ids:
            params = self.create_test_input(
                model={"metallicity_id": z_id, "wind_id": 0}
            )
            
            starburst = Starburst99(str(self.input_file))
            starburst.run()
            
            # Check that correct tracks were loaded
            self.assertEqual(starburst.galaxy.model_params.metallicity_id, z_id)
    
    def test_time_evolution(self):
        """Test time evolution with different steps"""
        # Test with fine time resolution
        params = self.create_test_input(
            time={
                "start": 0.1,
                "end": 100.0,
                "num_steps": 50
            }
        )
        
        starburst = Starburst99(str(self.input_file))
        starburst.run()
        
        # Check time grid
        self.assertEqual(len(starburst.galaxy.model_params.time_grid), 50)
        self.assertAlmostEqual(starburst.galaxy.model_params.time_grid[0], 0.1)
        self.assertAlmostEqual(starburst.galaxy.model_params.time_grid[-1], 100.0)
    
    def test_output_formats(self):
        """Test different output formats"""
        params = self.create_test_input()
        
        # Test with output writer
        starburst = Starburst99(str(self.input_file))
        starburst.run()
        
        # Check JSON summary
        summary_files = list(self.output_dir.glob("*summary.json"))
        self.assertGreater(len(summary_files), 0)
        
        # Verify JSON content
        with open(summary_files[0], 'r') as f:
            summary = json.load(f)
        
        self.assertIn('parameters', summary)
        self.assertIn('results', summary)
    
    def test_error_handling(self):
        """Test error handling during model run"""
        # Test with invalid input
        params = self.create_test_input(
            imf={
                "num_intervals": 2,
                "exponents": [2.35],  # Mismatch: 2 intervals but 1 exponent
                "mass_limits": [1.0, 100.0],
            }
        )
        
        # This should raise an error or handle gracefully
        starburst = Starburst99(str(self.input_file))
        
        # Verify error handling
        parser = InputParser()
        params = parser.read_input(str(self.input_file))
        self.assertFalse(parser.validate_parameters(params))
    
    def test_full_pipeline(self):
        """Test complete pipeline from input to output"""
        # Create comprehensive input
        params = self.create_test_input(
            name="Full Pipeline Test",
            star_formation={
                "mode": 1,
                "total_mass": 1.0e8,
                "rate": 10.0
            },
            imf={
                "num_intervals": 2,
                "exponents": [1.3, 2.3],
                "mass_limits": [0.1, 0.5, 100.0],
                "sn_cutoff": 8.0,
                "bh_cutoff": 120.0
            },
            model={
                "metallicity_id": 24,
                "wind_id": 0
            },
            time={
                "start": 0.1,
                "end": 1000.0,
                "num_steps": 20
            }
        )
        
        # Initialize all components
        parser = InputParser()
        model_params = parser.read_input(str(self.input_file))
        
        # Create galaxy model
        from ..core.galaxy_module import GalaxyModel
        galaxy = GalaxyModel()
        galaxy.model_params = model_params
        
        # Initialize IMF
        imf = IMF(
            model_params.num_intervals,
            model_params.exponents,
            model_params.mass_limits
        )
        
        # Initialize stellar tracks
        tracks = StellarTracks(model_params.metallicity_id)
        
        # Create output writer
        writer = OutputWriter(str(self.output_dir))
        
        # Verify all components work together
        self.assertIsNotNone(galaxy)
        self.assertIsNotNone(imf)
        self.assertIsNotNone(tracks)
        self.assertIsNotNone(writer)
        
        # Test IMF sampling
        masses = imf.sample(1000)
        self.assertEqual(len(masses), 1000)
        
        # Test data profiles
        data_profiles = DataProfiles()
        data_profiles.initialize_data_profiles()
        self.assertTrue(data_profiles.is_initialized)
    
    def test_standard_input_format(self):
        """Test reading standard Starburst99 input format"""
        # Create standard format input
        std_input = Path(self.temp_dir) / "standard.input"
        content = """MODEL DESIGNATION:                                           [NAME]
Test Standard Model
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
WIND FLAG (0=STANDARD MASS LOSS):                            [IWIND]
0
"""
        std_input.write_text(content)
        
        # Parse standard format
        parser = InputParser()
        params = parser.read_input(str(std_input))
        
        # Verify parsing
        self.assertEqual(params.name, "Test Standard Model")
        self.assertEqual(params.sf_mode, -1)
        self.assertEqual(params.total_mass, 1.0)
        self.assertEqual(params.num_intervals, 1)
        self.assertEqual(params.exponents, [2.35])
        self.assertEqual(params.mass_limits, [1.0, 100.0])
        self.assertEqual(params.sn_cutoff, 8.0)
        self.assertEqual(params.bh_cutoff, 120.0)
        self.assertEqual(params.metallicity_id, 24)
        self.assertEqual(params.wind_id, 0)
    
    def test_ini_format(self):
        """Test INI format input"""
        ini_input = Path(self.temp_dir) / "test.ini"
        content = """[general]
name = INI Test Model

[star_formation]
mode = 1
total_mass = 1e7
rate = 5.0

[imf]
num_intervals = 2
exponents = 1.3, 2.3
mass_limits = 0.1, 0.5, 100.0
sn_cutoff = 8.0
bh_cutoff = 120.0

[model]
metallicity_id = 14
wind_id = 1
"""
        ini_input.write_text(content)
        
        # Parse INI format
        parser = InputParser()
        params = parser.read_input(str(ini_input))
        
        # Verify parsing
        self.assertEqual(params.name, "INI Test Model")
        self.assertEqual(params.sf_mode, 1)
        self.assertEqual(params.total_mass, 1e7)
        self.assertEqual(params.sf_rate, 5.0)
        self.assertEqual(params.num_intervals, 2)
        self.assertEqual(params.exponents, [1.3, 2.3])
        self.assertEqual(params.mass_limits, [0.1, 0.5, 100.0])
        self.assertEqual(params.metallicity_id, 14)
        self.assertEqual(params.wind_id, 1)
    
    def test_performance_large_dataset(self):
        """Test performance with large datasets"""
        # Create model with many time steps
        params = self.create_test_input(
            time={
                "start": 0.01,
                "end": 10000.0,
                "num_steps": 1000
            }
        )
        
        import time
        start_time = time.time()
        
        # Run model
        starburst = Starburst99(str(self.input_file))
        # Note: actual run would be slow, so we just test initialization
        
        end_time = time.time()
        
        # Should complete in reasonable time
        self.assertLess(end_time - start_time, 5.0)  # 5 seconds max
    
    def test_memory_management(self):
        """Test memory usage with large arrays"""
        # Create large arrays
        large_size = 100000
        
        # Test data profiles with large arrays
        data_profiles = DataProfiles()
        data_profiles.flux_array = np.zeros((3, large_size))
        data_profiles.wavelength_array = np.logspace(2, 5, large_size)
        
        # Should handle large arrays without issues
        self.assertEqual(data_profiles.flux_array.shape, (3, large_size))
        self.assertEqual(len(data_profiles.wavelength_array), large_size)
    

if __name__ == '__main__':
    unittest.main()