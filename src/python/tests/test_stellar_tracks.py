"""Comprehensive tests for stellar_tracks module"""

import unittest
import tempfile
import shutil
import numpy as np
from pathlib import Path

from ..models.stellar_tracks import StellarTracks


class TestStellarTracks(unittest.TestCase):
    """Test StellarTracks class comprehensively"""
    
    def setUp(self):
        """Set up test environment"""
        self.temp_dir = tempfile.mkdtemp()
        self.test_data_dir = Path(self.temp_dir) / "test_tracks"
        self.test_data_dir.mkdir()
        
        self.stellar_tracks = StellarTracks(data_dir=self.test_data_dir)
        
        # Create test track files
        self.create_test_track_files()
        
    def tearDown(self):
        """Clean up test environment"""
        shutil.rmtree(self.temp_dir)
        
    def create_test_track_files(self):
        """Create dummy track files for testing"""
        # Simple track file for Z=0.020
        track_content = """# Test track file for Z=0.020
M=1.0
1.0e6  3.0  5700  1.0  1.0e-7  1
2.0e6  3.2  5600  1.1  1.1e-7  1
5.0e6  3.5  5500  1.2  1.2e-7  1
1.0e7  3.8  5400  1.3  1.5e-7  2

M=5.0
1.0e6  3.8  10000  5.0  1.0e-6  1
2.0e6  4.0  9500   5.5  1.5e-6  1
5.0e6  4.2  9000   6.0  2.0e-6  2
1.0e7  4.5  8500   7.0  3.0e-6  3

M=20.0
1.0e6  4.8  20000  20.0  1.0e-5  1
2.0e6  5.0  19000  22.0  2.0e-5  1
5.0e6  5.2  18000  25.0  5.0e-5  2
1.0e7  5.5  17000  30.0  1.0e-4  3
"""
        
        z020_file = self.test_data_dir / "Z0020v00.txt"
        z020_file.write_text(track_content)
        
        # Create another metallicity file
        z040_content = track_content.replace("Z=0.020", "Z=0.040")
        z040_file = self.test_data_dir / "Z0040v00.txt"
        z040_file.write_text(z040_content)
        
    def test_initialization(self):
        """Test StellarTracks initialization"""
        self.assertEqual(self.stellar_tracks.data_dir, self.test_data_dir)
        self.assertEqual(len(self.stellar_tracks.tracks), 0)
        self.assertEqual(len(self.stellar_tracks.metallicities), 0)
        
    def test_load_tracks_success(self):
        """Test successful track loading"""
        self.stellar_tracks.load_tracks("Z0020v00")
        
        self.assertIn("Z0020v00", self.stellar_tracks.tracks)
        self.assertIn("Z0020v00", self.stellar_tracks.metallicities)
        
        # Check track data structure
        track_data = self.stellar_tracks.tracks["Z0020v00"]
        self.assertIn(1.0, track_data)
        self.assertIn(5.0, track_data)
        self.assertIn(20.0, track_data)
        
        # Check data arrays
        mass1_data = track_data[1.0]
        self.assertIn('time', mass1_data)
        self.assertIn('luminosity', mass1_data)
        self.assertIn('temperature', mass1_data)
        self.assertIn('radius', mass1_data)
        self.assertIn('mass_loss_rate', mass1_data)
        self.assertIn('stellar_type', mass1_data)
        
        # Check array lengths
        self.assertEqual(len(mass1_data['time']), 4)
        self.assertEqual(len(mass1_data['luminosity']), 4)
        
        # Check values
        np.testing.assert_array_equal(mass1_data['time'], [1.0e6, 2.0e6, 5.0e6, 1.0e7])
        np.testing.assert_array_equal(mass1_data['luminosity'], [3.0, 3.2, 3.5, 3.8])
        
    def test_load_tracks_file_not_found(self):
        """Test track loading with non-existent file"""
        with self.assertRaises(FileNotFoundError):
            self.stellar_tracks.load_tracks("nonexistent")
            
    def test_load_tracks_empty_file(self):
        """Test loading empty track file"""
        empty_file = self.test_data_dir / "empty.txt"
        empty_file.write_text("")
        
        self.stellar_tracks.load_tracks("empty")
        
        # Should load but have no masses
        self.assertIn("empty", self.stellar_tracks.tracks)
        self.assertEqual(len(self.stellar_tracks.tracks["empty"]), 0)
        
    def test_load_tracks_with_comments(self):
        """Test loading tracks with comment lines"""
        comment_content = """# This is a comment
# Another comment

M=2.0
# Comment in between
1.0e6  3.5  8000  2.0  5.0e-7  1
2.0e6  3.6  7900  2.1  6.0e-7  1

# Final comment
"""
        
        comment_file = self.test_data_dir / "comment_test.txt"
        comment_file.write_text(comment_content)
        
        self.stellar_tracks.load_tracks("comment_test")
        
        track_data = self.stellar_tracks.tracks["comment_test"]
        self.assertIn(2.0, track_data)
        self.assertEqual(len(track_data[2.0]['time']), 2)
        
    def test_interpolate_track_exact_mass(self):
        """Test track interpolation for exact mass match"""
        self.stellar_tracks.load_tracks("Z0020v00")
        
        # Exact mass and time
        props = self.stellar_tracks.interpolate_track(5.0, 2.0e6, "Z0020v00")
        
        self.assertAlmostEqual(props['luminosity'], 4.0)
        self.assertAlmostEqual(props['temperature'], 9500)
        self.assertAlmostEqual(props['radius'], 5.5)
        self.assertAlmostEqual(props['mass_loss_rate'], 1.5e-6)
        self.assertEqual(props['stellar_type'], 1)
        
    def test_interpolate_track_time_interpolation(self):
        """Test time interpolation within a track"""
        self.stellar_tracks.load_tracks("Z0020v00")
        
        # Interpolate at intermediate time
        props = self.stellar_tracks.interpolate_track(5.0, 3.5e6, "Z0020v00")
        
        # Should be between values at 2e6 and 5e6
        self.assertTrue(4.0 < props['luminosity'] < 4.2)
        self.assertTrue(9000 < props['temperature'] < 9500)
        
    def test_interpolate_track_mass_interpolation(self):
        """Test mass interpolation between tracks"""
        self.stellar_tracks.load_tracks("Z0020v00")
        
        # Mass between 5.0 and 20.0
        props = self.stellar_tracks.interpolate_track(10.0, 2.0e6, "Z0020v00")
        
        # Should be between values for M=5 and M=20
        self.assertTrue(4.0 < props['luminosity'] < 5.0)
        self.assertTrue(9500 < props['temperature'] < 19000)
        
    def test_interpolate_track_below_minimum_mass(self):
        """Test interpolation below minimum mass"""
        self.stellar_tracks.load_tracks("Z0020v00")
        
        # Mass below all tracks
        props = self.stellar_tracks.interpolate_track(0.5, 2.0e6, "Z0020v00")
        
        # Should use lowest mass track (1.0)
        self.assertAlmostEqual(props['luminosity'], 3.2)
        self.assertAlmostEqual(props['temperature'], 5600)
        
    def test_interpolate_track_above_maximum_mass(self):
        """Test interpolation above maximum mass"""
        self.stellar_tracks.load_tracks("Z0020v00")
        
        # Mass above all tracks
        props = self.stellar_tracks.interpolate_track(50.0, 2.0e6, "Z0020v00")
        
        # Should use highest mass track (20.0)
        self.assertAlmostEqual(props['luminosity'], 5.0)
        self.assertAlmostEqual(props['temperature'], 19000)
        
    def test_interpolate_track_stellar_type(self):
        """Test stellar type interpolation (discrete value)"""
        self.stellar_tracks.load_tracks("Z0020v00")
        
        # Mass interpolation with different stellar types
        props1 = self.stellar_tracks.interpolate_track(3.0, 2.0e6, "Z0020v00")
        props2 = self.stellar_tracks.interpolate_track(15.0, 2.0e6, "Z0020v00")
        
        # Should use type from dominant mass
        self.assertEqual(props1['stellar_type'], 1)  # Closer to M=1
        self.assertEqual(props2['stellar_type'], 1)  # Closer to M=20
        
    def test_interpolate_track_not_loaded(self):
        """Test interpolation triggers loading if not loaded"""
        # Track not loaded yet
        self.assertNotIn("Z0020v00", self.stellar_tracks.tracks)
        
        # Should load automatically
        props = self.stellar_tracks.interpolate_track(5.0, 2.0e6, "Z0020v00")
        
        self.assertIn("Z0020v00", self.stellar_tracks.tracks)
        self.assertIsNotNone(props)
        
    def test_get_lifetime_exact_mass(self):
        """Test getting lifetime for exact mass"""
        self.stellar_tracks.load_tracks("Z0020v00")
        
        lifetime = self.stellar_tracks.get_lifetime(5.0, "Z0020v00")
        
        # Should be the last time point
        self.assertEqual(lifetime, 1.0e7)
        
    def test_get_lifetime_interpolated_mass(self):
        """Test getting lifetime for interpolated mass"""
        self.stellar_tracks.load_tracks("Z0020v00")
        
        # Mass between 1.0 and 5.0
        lifetime = self.stellar_tracks.get_lifetime(2.5, "Z0020v00")
        
        # Should be between lifetimes of 1.0 and 5.0
        self.assertTrue(1.0e7 <= lifetime <= 1.0e7)  # Both have same lifetime in test data
        
        # For log interpolation test, create tracks with different lifetimes
        varied_content = """M=1.0
1.0e6  3.0  5700  1.0  1.0e-7  1
1.0e8  3.0  5700  1.0  1.0e-7  1

M=10.0
1.0e6  4.0  10000  10.0  1.0e-6  1
1.0e7  4.0  10000  10.0  1.0e-6  1
"""
        
        varied_file = self.test_data_dir / "varied_lifetime.txt"
        varied_file.write_text(varied_content)
        
        varied_tracks = StellarTracks(data_dir=self.test_data_dir)
        varied_tracks.load_tracks("varied_lifetime")
        
        # Test log interpolation
        lifetime = varied_tracks.get_lifetime(3.16, "varied_lifetime")  # sqrt(10)
        
        # Should be approximately geometric mean of 1e8 and 1e7
        expected = np.sqrt(1.0e8 * 1.0e7)
        # Allow for some tolerance in the log interpolation
        self.assertAlmostEqual(lifetime, expected, delta=expected*0.001)
        
    def test_get_lifetime_edge_cases(self):
        """Test lifetime for edge cases"""
        self.stellar_tracks.load_tracks("Z0020v00")
        
        # Below minimum mass
        lifetime = self.stellar_tracks.get_lifetime(0.5, "Z0020v00")
        self.assertEqual(lifetime, 1.0e7)  # Uses M=1.0
        
        # Above maximum mass
        lifetime = self.stellar_tracks.get_lifetime(50.0, "Z0020v00")
        self.assertEqual(lifetime, 1.0e7)  # Uses M=20.0
        
    def test_multiple_metallicities(self):
        """Test handling multiple metallicities"""
        self.stellar_tracks.load_tracks("Z0020v00")
        self.stellar_tracks.load_tracks("Z0040v00")
        
        self.assertEqual(len(self.stellar_tracks.metallicities), 2)
        self.assertIn("Z0020v00", self.stellar_tracks.metallicities)
        self.assertIn("Z0040v00", self.stellar_tracks.metallicities)
        
        # Interpolate from different metallicities
        props1 = self.stellar_tracks.interpolate_track(5.0, 2.0e6, "Z0020v00")
        props2 = self.stellar_tracks.interpolate_track(5.0, 2.0e6, "Z0040v00")
        
        # Should give same values (test files are identical except for header)
        self.assertAlmostEqual(props1['luminosity'], props2['luminosity'])
        
    def test_array_properties(self):
        """Test that all arrays are numpy arrays"""
        self.stellar_tracks.load_tracks("Z0020v00")
        
        track_data = self.stellar_tracks.tracks["Z0020v00"][1.0]
        
        for key in ['time', 'luminosity', 'temperature', 'radius', 
                   'mass_loss_rate', 'stellar_type']:
            self.assertIsInstance(track_data[key], np.ndarray)


if __name__ == '__main__':
    unittest.main()