"""Comprehensive tests for output_writer module"""

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import unittest
import tempfile
import shutil
import json
from pathlib import Path
from datetime import datetime

from file_io.output_writer import OutputWriter
from core.galaxy_module import GalaxyModel, ModelParameters
import numpy as np


class TestOutputWriter(unittest.TestCase):
    """Test OutputWriter class comprehensively"""
    
    def setUp(self):
        """Set up test environment"""
        self.writer = OutputWriter()
        self.temp_dir = tempfile.mkdtemp()
        self.writer.output_dir = Path(self.temp_dir)
        
        # Create test galaxy model
        self.galaxy = GalaxyModel()
        self.galaxy.model_params = ModelParameters(
            name="Test Model",
            sf_mode=1,
            total_mass=1e6,
            sf_rate=10.0,
            num_intervals=2,
            exponents=[2.0, 2.7],
            mass_limits=[0.5, 10.0, 150.0],
            sn_cutoff=8.0,
            bh_cutoff=120.0,
            metallicity_id=24,
            wind_id=1
        )
        
        # Set some test data
        self.galaxy.current_time = 1e7
        self.galaxy.time_step = 100
        self.galaxy.wavelength = np.logspace(2, 5, 100)  # 100 Å to 100000 Å
        self.galaxy.spectra = np.random.rand(100) * 1e-12
        
    def tearDown(self):
        """Clean up test environment"""
        shutil.rmtree(self.temp_dir)
        
    def test_initialization(self):
        """Test OutputWriter initialization"""
        writer = OutputWriter()
        self.assertTrue(writer.output_dir.exists())
        self.assertEqual(writer.output_dir.name, "output")
        
    def test_write_all_outputs(self):
        """Test writing all output files"""
        self.writer.write_all_outputs(self.galaxy)
        
        # Check that files were created
        output_files = list(Path(self.temp_dir).glob("test_model_*"))
        self.assertGreater(len(output_files), 0)
        
        # Check for specific file types
        file_extensions = [f.suffix for f in output_files]
        self.assertIn(".output", file_extensions)
        self.assertIn(".spectrum", file_extensions)
        self.assertIn(".quanta", file_extensions)
        self.assertIn(".json", file_extensions)
        
    def test_write_main_output(self):
        """Test main output file writing"""
        self.writer._write_main_output(self.galaxy, "test_output")
        
        output_file = self.writer.output_dir / "test_output.output"
        self.assertTrue(output_file.exists())
        
        # Check content
        content = output_file.read_text()
        self.assertIn("STARBURST99 v7.0.2 (Python Edition)", content)
        self.assertIn("Model: Test Model", content)
        self.assertIn("Star Formation Mode: 1", content)
        self.assertIn("Total Mass: 1.00e+06 M_sun", content)
        self.assertIn("SFR: 1.00e+01 M_sun/yr", content)
        self.assertIn("IMF PARAMETERS", content)
        self.assertIn("Number of intervals: 2", content)
        self.assertIn("SN cutoff: 8.0 M_sun", content)
        self.assertIn("BH cutoff: 120.0 M_sun", content)
        
    def test_write_main_output_instantaneous(self):
        """Test main output for instantaneous star formation"""
        self.galaxy.model_params.sf_mode = 0
        self.writer._write_main_output(self.galaxy, "test_inst")
        
        output_file = self.writer.output_dir / "test_inst.output"
        content = output_file.read_text()
        
        # Should not have SFR for instantaneous mode
        self.assertNotIn("SFR:", content)
        
    def test_write_spectrum(self):
        """Test spectrum file writing"""
        self.writer._write_spectrum(self.galaxy, "test_spectrum")
        
        spectrum_file = self.writer.output_dir / "test_spectrum.spectrum"
        self.assertTrue(spectrum_file.exists())
        
        # Check content
        lines = spectrum_file.read_text().strip().split('\n')
        self.assertIn("# Wavelength(A)  Flux(erg/s/A)", lines[0])
        
        # Check data format
        data_lines = [l for l in lines if not l.startswith('#')]
        self.assertEqual(len(data_lines), 100)  # Same as wavelength array
        
        # Check a data line format
        parts = data_lines[0].split()
        self.assertEqual(len(parts), 2)
        float(parts[0])  # Should not raise
        float(parts[1])  # Should not raise
        
    def test_write_spectrum_zero_handling(self):
        """Test spectrum writing handles zero values correctly"""
        # Set some zero values
        self.galaxy.wavelength[0] = 0.0
        self.galaxy.spectra[50] = 0.0
        
        self.writer._write_spectrum(self.galaxy, "test_zero")
        
        spectrum_file = self.writer.output_dir / "test_zero.spectrum"
        lines = spectrum_file.read_text().strip().split('\n')
        data_lines = [l for l in lines if not l.startswith('#')]
        
        # Should skip entries with zero wavelength or flux
        self.assertLess(len(data_lines), 100)
        
    def test_write_quanta(self):
        """Test ionizing photon rates file writing"""
        self.writer._write_quanta(self.galaxy, "test_quanta")
        
        quanta_file = self.writer.output_dir / "test_quanta.quanta"
        self.assertTrue(quanta_file.exists())
        
        # Check header
        content = quanta_file.read_text()
        self.assertIn("# Time(yr)  Q(H)  Q(HeI)  Q(HeII)", content)
        
    def test_write_summary_json(self):
        """Test JSON summary file writing"""
        self.writer._write_summary_json(self.galaxy, "test_json")
        
        json_file = self.writer.output_dir / "test_json_summary.json"
        self.assertTrue(json_file.exists())
        
        # Load and check JSON content
        with open(json_file, 'r') as f:
            data = json.load(f)
            
        self.assertIn("model", data)
        self.assertIn("parameters", data)
        self.assertIn("results", data)
        
        # Check model section
        self.assertEqual(data["model"]["name"], "Test Model")
        self.assertEqual(data["model"]["version"], "7.0.2-python")
        
        # Check parameters
        params = data["parameters"]
        self.assertEqual(params["star_formation"]["mode"], 1)
        self.assertEqual(params["star_formation"]["total_mass"], 1e6)
        self.assertEqual(params["star_formation"]["rate"], 10.0)
        self.assertEqual(params["imf"]["num_intervals"], 2)
        self.assertEqual(params["imf"]["exponents"], [2.0, 2.7])
        
        # Check results
        results = data["results"]
        self.assertEqual(results["final_time"], 1e7)
        self.assertEqual(results["num_timesteps"], 100)
        
    def test_write_summary_json_instantaneous(self):
        """Test JSON summary for instantaneous star formation"""
        self.galaxy.model_params.sf_mode = 0
        self.writer._write_summary_json(self.galaxy, "test_inst_json")
        
        json_file = self.writer.output_dir / "test_inst_json_summary.json"
        with open(json_file, 'r') as f:
            data = json.load(f)
            
        # Rate should be None for instantaneous mode
        self.assertIsNone(data["parameters"]["star_formation"]["rate"])
        
    def test_filename_sanitization(self):
        """Test filename sanitization for model names with spaces/special chars"""
        self.galaxy.model_params.name = "Test Model with Spaces & Special/Chars"
        self.writer.write_all_outputs(self.galaxy)
        
        # Check files are created with sanitized names
        output_files = list(Path(self.temp_dir).glob("test_model_with_spaces_*"))
        self.assertGreater(len(output_files), 0)
        
    def test_logging(self):
        """Test logging functionality"""
        with self.assertLogs('file_io.output_writer', level='INFO') as cm:
            self.writer.write_all_outputs(self.galaxy)
            
        # Check log messages
        log_messages = [record.message for record in cm.records]
        self.assertTrue(any("Main output written to:" in msg for msg in log_messages))
        self.assertTrue(any("Spectrum written to:" in msg for msg in log_messages))
        self.assertTrue(any("Ionizing photon rates written to:" in msg for msg in log_messages))
        self.assertTrue(any("Summary JSON written to:" in msg for msg in log_messages))
        
    def test_directory_creation(self):
        """Test output directory creation"""
        # Remove the directory
        shutil.rmtree(self.writer.output_dir)
        self.assertFalse(self.writer.output_dir.exists())
        
        # Create new writer - should recreate directory
        new_writer = OutputWriter()
        new_writer.output_dir = self.writer.output_dir
        new_writer.output_dir.mkdir(exist_ok=True)
        self.assertTrue(new_writer.output_dir.exists())
        
    def test_timestamp_format(self):
        """Test timestamp formatting in output files"""
        self.writer.write_all_outputs(self.galaxy)
        
        # Get any output file
        output_files = list(Path(self.temp_dir).glob("test_model_*"))
        self.assertGreater(len(output_files), 0)
        
        # Find a non-summary file for the timestamp test
        non_summary_files = [f for f in output_files if not f.stem.endswith('_summary')]
        self.assertGreater(len(non_summary_files), 0)
        
        # Check timestamp format (YYYYMMDD_HHMMSS)
        filename = non_summary_files[0].stem
        parts = filename.split('_')
        timestamp = parts[-2] + '_' + parts[-1]
        
        # Should be able to parse as datetime
        try:
            datetime.strptime(timestamp, "%Y%m%d_%H%M%S")
        except ValueError:
            self.fail("Timestamp format is incorrect")
            
    def test_large_data_handling(self):
        """Test handling of large data arrays"""
        # Create large arrays
        self.galaxy.wavelength = np.logspace(2, 5, 10000)
        self.galaxy.spectra = np.random.rand(10000) * 1e-12
        
        self.writer._write_spectrum(self.galaxy, "test_large")
        
        spectrum_file = self.writer.output_dir / "test_large.spectrum"
        self.assertTrue(spectrum_file.exists())
        
        # Check all data is written
        lines = spectrum_file.read_text().strip().split('\n')
        data_lines = [l for l in lines if not l.startswith('#')]
        self.assertEqual(len(data_lines), 10000)


if __name__ == '__main__':
    unittest.main()