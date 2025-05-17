"""Tests for edge cases in stellar tracks to achieve 100% coverage"""

import unittest
import tempfile
import shutil
from pathlib import Path
import numpy as np

from ..models.stellar_tracks import StellarTracks


class TestStellarTracksEdgeCases(unittest.TestCase):
    """Test edge cases in stellar tracks for 100% coverage"""
    
    def setUp(self):
        """Set up test environment"""
        self.temp_dir = tempfile.mkdtemp()
        self.test_data_dir = Path(self.temp_dir) / "test_tracks"
        self.test_data_dir.mkdir()
        
    def tearDown(self):
        """Clean up test environment"""
        shutil.rmtree(self.temp_dir)
    
    def test_interpolate_single_track_edge_case(self):
        """Test single track interpolation edge case (line 184)"""
        stellar_tracks = StellarTracks(data_dir=self.test_data_dir)
        
        # Create a track file with different stellar types
        track_content = """M=10.0
1.0e6  4.0  15000  10.0  1.0e-6  1
2.0e6  4.1  14500  10.5  1.5e-6  1
3.0e6  4.2  14000  11.0  2.0e-6  2
4.0e6  4.3  13500  11.5  2.5e-6  2
5.0e6  4.4  13000  12.0  3.0e-6  3
"""
        
        test_file = self.test_data_dir / "test_single.txt"
        test_file.write_text(track_content)
        
        stellar_tracks.load_tracks("test_single")
        
        # Test interpolation that hits the mass match case
        props = stellar_tracks.interpolate_track(10.0, 2.5e6, "test_single")
        
        # Properties should be interpolated between times
        self.assertTrue(4.1 < props['luminosity'] < 4.2)
        # Stellar type is from nearest time point
        
        # Test with time that gives stellar type 1 or 2 depending on interpolation
        props = stellar_tracks.interpolate_track(10.0, 1.5e6, "test_single")
        self.assertEqual(props['stellar_type'], 1)  # Closer to 2e6 which has type 1


if __name__ == '__main__':
    unittest.main()