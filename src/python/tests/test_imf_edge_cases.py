"""Tests for edge cases in IMF to achieve 100% coverage"""

import unittest
import numpy as np
from ..models.imf import IMF


class TestIMFEdgeCases(unittest.TestCase):
    """Test edge cases in IMF for 100% coverage"""
    
    def test_integrate_with_alpha_equals_1(self):
        """Test integration when alpha equals 1 (line 114)"""
        # Create IMF where alpha = 1
        imf = IMF(
            intervals=1,
            exponents=[1.0],  # alpha = 1
            mass_limits=[1.0, 10.0]
        )
        
        # Test integration
        result = imf.integrate(2.0, 5.0)
        
        # When alpha = 1, we should use the log formula
        expected = np.log(5.0 / 2.0)
        self.assertAlmostEqual(result, expected, places=10)
    
    def test_sample_edge_case(self):
        """Test sampling edge case (line 160)"""
        # Create a simple IMF
        imf = IMF(
            intervals=1,
            exponents=[2.35],
            mass_limits=[1.0, 100.0]
        )
        
        # Sample with specific seed to ensure consistent results
        masses = imf.sample(10, seed=12345)
        
        # Check that all masses are in valid range
        self.assertTrue(np.all(masses >= 1.0))
        self.assertTrue(np.all(masses <= 100.0))
        self.assertEqual(len(masses), 10)


if __name__ == '__main__':
    unittest.main()