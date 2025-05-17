"""Comprehensive tests for utilities module"""

import unittest
import numpy as np
from ..utils.utilities import exp10, linear_interp, integer_to_string


class TestUtilities(unittest.TestCase):
    """Test utility functions comprehensively"""
    
    def test_exp10_scalar(self):
        """Test exp10 with scalar values"""
        self.assertAlmostEqual(exp10(0), 1.0)
        self.assertAlmostEqual(exp10(1), 10.0)
        self.assertAlmostEqual(exp10(2), 100.0)
        self.assertAlmostEqual(exp10(-1), 0.1)
        self.assertAlmostEqual(exp10(-2), 0.01)
        self.assertAlmostEqual(exp10(0.5), np.sqrt(10))
        
    def test_exp10_array(self):
        """Test exp10 with numpy arrays"""
        x = np.array([0, 1, 2, -1])
        expected = np.array([1.0, 10.0, 100.0, 0.1])
        result = exp10(x)
        np.testing.assert_array_almost_equal(result, expected)
        
    def test_exp10_edge_cases(self):
        """Test exp10 edge cases"""
        # Very large values
        self.assertTrue(np.isfinite(exp10(300)))
        # Very small values
        self.assertAlmostEqual(exp10(-300), 0.0)
        
    def test_linear_interp_normal(self):
        """Test linear interpolation normal cases"""
        x_arr = np.array([0, 1, 2, 3, 4])
        y_arr = np.array([0, 2, 4, 6, 8])
        
        # Exact points
        self.assertAlmostEqual(linear_interp(0, x_arr, y_arr), 0.0)
        self.assertAlmostEqual(linear_interp(1, x_arr, y_arr), 2.0)
        self.assertAlmostEqual(linear_interp(2, x_arr, y_arr), 4.0)
        
        # Interpolated points
        self.assertAlmostEqual(linear_interp(0.5, x_arr, y_arr), 1.0)
        self.assertAlmostEqual(linear_interp(1.5, x_arr, y_arr), 3.0)
        self.assertAlmostEqual(linear_interp(2.5, x_arr, y_arr), 5.0)
        self.assertAlmostEqual(linear_interp(3.75, x_arr, y_arr), 7.5)
        
    def test_linear_interp_extrapolation(self):
        """Test linear interpolation extrapolation"""
        x_arr = np.array([0, 1, 2, 3])
        y_arr = np.array([0, 2, 4, 6])
        
        # Below range - should return first value
        self.assertAlmostEqual(linear_interp(-1, x_arr, y_arr), 0.0)
        self.assertAlmostEqual(linear_interp(-10, x_arr, y_arr), 0.0)
        
        # Above range - should return last value
        self.assertAlmostEqual(linear_interp(4, x_arr, y_arr), 6.0)
        self.assertAlmostEqual(linear_interp(10, x_arr, y_arr), 6.0)
        
    def test_linear_interp_edge_cases(self):
        """Test linear interpolation edge cases"""
        # Single point
        x_arr = np.array([1.0])
        y_arr = np.array([5.0])
        self.assertAlmostEqual(linear_interp(0.5, x_arr, y_arr), 5.0)
        self.assertAlmostEqual(linear_interp(1.0, x_arr, y_arr), 5.0)
        self.assertAlmostEqual(linear_interp(2.0, x_arr, y_arr), 5.0)
        
        # Two points
        x_arr = np.array([0.0, 1.0])
        y_arr = np.array([0.0, 10.0])
        self.assertAlmostEqual(linear_interp(0.5, x_arr, y_arr), 5.0)
        
        # Non-uniform spacing
        x_arr = np.array([0, 1, 10, 100])
        y_arr = np.array([0, 1, 10, 100])
        self.assertAlmostEqual(linear_interp(5.5, x_arr, y_arr), 5.5)
        self.assertAlmostEqual(linear_interp(55, x_arr, y_arr), 55)
        
    def test_integer_to_string_no_padding(self):
        """Test integer to string conversion without padding"""
        self.assertEqual(integer_to_string(0), "0")
        self.assertEqual(integer_to_string(1), "1")
        self.assertEqual(integer_to_string(42), "42")
        self.assertEqual(integer_to_string(999), "999")
        self.assertEqual(integer_to_string(-5), "-5")
        
    def test_integer_to_string_with_padding(self):
        """Test integer to string conversion with padding"""
        self.assertEqual(integer_to_string(0, 3), "000")
        self.assertEqual(integer_to_string(1, 3), "001")
        self.assertEqual(integer_to_string(42, 3), "042")
        self.assertEqual(integer_to_string(999, 3), "999")
        self.assertEqual(integer_to_string(42, 5), "00042")
        
    def test_integer_to_string_edge_cases(self):
        """Test integer to string edge cases"""
        # Width smaller than number
        self.assertEqual(integer_to_string(1234, 2), "1234")
        
        # Width of 0 (same as no padding)
        self.assertEqual(integer_to_string(42, 0), "42")
        
        # Negative numbers with padding
        self.assertEqual(integer_to_string(-5, 3), "-05")
        self.assertEqual(integer_to_string(-42, 5), "-0042")


if __name__ == '__main__':
    unittest.main()