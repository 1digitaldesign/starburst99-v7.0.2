"""Comprehensive tests for IMF module"""

import unittest
import numpy as np
from ..models.imf import IMF


class TestIMF(unittest.TestCase):
    """Test Initial Mass Function class comprehensively"""
    
    def setUp(self):
        """Set up test IMF instances"""
        # Standard Salpeter IMF
        self.salpeter_imf = IMF(
            intervals=1,
            exponents=[2.35],
            mass_limits=[1.0, 100.0]
        )
        
        # Multi-segment IMF
        self.multi_imf = IMF(
            intervals=3,
            exponents=[1.3, 2.3, 2.7],
            mass_limits=[0.1, 0.5, 8.0, 120.0]
        )
        
        # Flat IMF (exponent = -1)
        self.flat_imf = IMF(
            intervals=1,
            exponents=[-1.0],
            mass_limits=[1.0, 10.0]
        )
        
    def test_initialization(self):
        """Test IMF initialization"""
        self.assertEqual(self.salpeter_imf.intervals, 1)
        np.testing.assert_array_equal(self.salpeter_imf.exponents, [2.35])
        np.testing.assert_array_equal(self.salpeter_imf.mass_limits, [1.0, 100.0])
        
    def test_validation_errors(self):
        """Test parameter validation"""
        # Wrong number of exponents
        with self.assertRaises(ValueError) as context:
            IMF(intervals=2, exponents=[2.35], mass_limits=[1.0, 10.0, 100.0])
        self.assertIn("Number of exponents must match", str(context.exception))
        
        # Wrong number of mass limits
        with self.assertRaises(ValueError) as context:
            IMF(intervals=1, exponents=[2.35], mass_limits=[1.0])
        self.assertIn("Number of mass limits must be intervals + 1", str(context.exception))
        
        # Non-increasing mass limits
        with self.assertRaises(ValueError) as context:
            IMF(intervals=1, exponents=[2.35], mass_limits=[10.0, 1.0])
        self.assertIn("Mass limits must be in increasing order", str(context.exception))
        
    def test_normalizations(self):
        """Test normalization calculations"""
        # Check that normalizations are calculated
        self.assertEqual(len(self.salpeter_imf.normalizations), 1)
        self.assertEqual(len(self.multi_imf.normalizations), 3)
        
        # For flat IMF, normalization should use log formula
        self.assertAlmostEqual(
            self.flat_imf.normalizations[0],
            np.log(10.0 / 1.0)
        )
        
    def test_xi_salpeter(self):
        """Test IMF value calculation for Salpeter IMF"""
        # Inside range
        self.assertGreater(self.salpeter_imf.xi(10.0), 0.0)
        self.assertGreater(self.salpeter_imf.xi(1.0), 0.0)
        self.assertGreater(self.salpeter_imf.xi(100.0), 0.0)
        
        # Outside range
        self.assertEqual(self.salpeter_imf.xi(0.5), 0.0)
        self.assertEqual(self.salpeter_imf.xi(200.0), 0.0)
        
        # Check relative values (higher mass = lower IMF value for positive exponent)
        self.assertGreater(self.salpeter_imf.xi(1.0), self.salpeter_imf.xi(10.0))
        self.assertGreater(self.salpeter_imf.xi(10.0), self.salpeter_imf.xi(100.0))
        
    def test_xi_multi_segment(self):
        """Test IMF value calculation for multi-segment IMF"""
        # Test values in different segments
        self.assertGreater(self.multi_imf.xi(0.2), 0.0)  # First segment
        self.assertGreater(self.multi_imf.xi(1.0), 0.0)  # Second segment
        self.assertGreater(self.multi_imf.xi(50.0), 0.0)  # Third segment
        
        # Test boundaries
        self.assertGreater(self.multi_imf.xi(0.5), 0.0)  # Boundary value
        self.assertGreater(self.multi_imf.xi(8.0), 0.0)  # Boundary value
        
    def test_xi_flat(self):
        """Test IMF value calculation for flat IMF"""
        # For flat IMF (exponent = -1), xi(m) = 1/m
        val1 = self.flat_imf.xi(2.0)
        val2 = self.flat_imf.xi(5.0)
        val3 = self.flat_imf.xi(8.0)
        
        # Check that xi(m) = 1/m
        self.assertAlmostEqual(val1, 0.5)
        self.assertAlmostEqual(val2, 0.2)
        self.assertAlmostEqual(val3, 0.125)
        
    def test_integrate_salpeter(self):
        """Test integration for Salpeter IMF"""
        # Full range
        total = self.salpeter_imf.integrate(1.0, 100.0)
        self.assertGreater(total, 0.0)
        
        # Partial ranges
        low_mass = self.salpeter_imf.integrate(1.0, 10.0)
        high_mass = self.salpeter_imf.integrate(10.0, 100.0)
        
        self.assertGreater(low_mass, 0.0)
        self.assertGreater(high_mass, 0.0)
        self.assertAlmostEqual(total, low_mass + high_mass, places=10)
        
        # Outside range
        self.assertEqual(self.salpeter_imf.integrate(0.1, 0.5), 0.0)
        self.assertEqual(self.salpeter_imf.integrate(200.0, 300.0), 0.0)
        
        # Mixed range
        mixed = self.salpeter_imf.integrate(0.5, 50.0)
        self.assertAlmostEqual(mixed, self.salpeter_imf.integrate(1.0, 50.0))
        
    def test_integrate_multi_segment(self):
        """Test integration for multi-segment IMF"""
        # Each segment
        seg1 = self.multi_imf.integrate(0.1, 0.5)
        seg2 = self.multi_imf.integrate(0.5, 8.0)
        seg3 = self.multi_imf.integrate(8.0, 120.0)
        
        self.assertGreater(seg1, 0.0)
        self.assertGreater(seg2, 0.0)
        self.assertGreater(seg3, 0.0)
        
        # Total should equal sum of segments
        total = self.multi_imf.integrate(0.1, 120.0)
        self.assertAlmostEqual(total, seg1 + seg2 + seg3, places=10)
        
        # Across segments
        cross = self.multi_imf.integrate(0.3, 20.0)
        self.assertGreater(cross, 0.0)
        
    def test_integrate_flat(self):
        """Test integration for flat IMF"""
        # For flat IMF with xi(m) = 1/m, integral is ln(m_high/m_low)
        int1 = self.flat_imf.integrate(1.0, 2.0)
        int2 = self.flat_imf.integrate(2.0, 4.0)
        
        # Check values
        self.assertAlmostEqual(int1, np.log(2.0), places=10)
        self.assertAlmostEqual(int2, np.log(2.0), places=10)
        
        # They should be equal (both are ln(2))
        self.assertAlmostEqual(int1, int2, places=10)
        
    def test_integrate_edge_cases(self):
        """Test integration edge cases"""
        # Same limits
        self.assertEqual(self.salpeter_imf.integrate(10.0, 10.0), 0.0)
        
        # Reversed limits
        self.assertEqual(self.salpeter_imf.integrate(50.0, 10.0), 0.0)
        
        # Very small range
        small = self.salpeter_imf.integrate(10.0, 10.001)
        self.assertGreater(small, 0.0)
        self.assertLess(small, 0.1)
        
    def test_sample(self):
        """Test mass sampling from IMF"""
        # Sample masses
        n_stars = 1000
        masses = self.salpeter_imf.sample(n_stars, seed=42)
        
        # Check properties
        self.assertEqual(len(masses), n_stars)
        self.assertTrue(np.all(masses >= 1.0))
        self.assertTrue(np.all(masses <= 100.0))
        
        # Check reproducibility
        masses2 = self.salpeter_imf.sample(n_stars, seed=42)
        np.testing.assert_array_equal(masses, masses2)
        
        # Check distribution roughly follows IMF
        # More low-mass stars than high-mass
        low_mass_count = np.sum(masses < 10.0)
        high_mass_count = np.sum(masses > 10.0)
        self.assertGreater(low_mass_count, high_mass_count)
        
    def test_sample_multi_segment(self):
        """Test sampling from multi-segment IMF"""
        n_stars = 1000
        masses = self.multi_imf.sample(n_stars, seed=42)
        
        # Check all masses are in valid range
        self.assertTrue(np.all(masses >= 0.1))
        self.assertTrue(np.all(masses <= 120.0))
        
        # Check masses in each segment
        seg1_count = np.sum((masses >= 0.1) & (masses < 0.5))
        seg2_count = np.sum((masses >= 0.5) & (masses < 8.0))
        seg3_count = np.sum((masses >= 8.0) & (masses <= 120.0))
        
        # All segments should have some stars
        self.assertGreater(seg1_count, 0)
        self.assertGreater(seg2_count, 0)
        self.assertGreater(seg3_count, 0)
        
        # Total should equal n_stars
        self.assertEqual(seg1_count + seg2_count + seg3_count, n_stars)
        
    def test_sample_edge_cases(self):
        """Test sampling edge cases"""
        # Sample zero stars
        masses = self.salpeter_imf.sample(0, seed=42)
        self.assertEqual(len(masses), 0)
        
        # Sample one star
        masses = self.salpeter_imf.sample(1, seed=42)
        self.assertEqual(len(masses), 1)
        self.assertGreaterEqual(masses[0], 1.0)
        self.assertLessEqual(masses[0], 100.0)


if __name__ == '__main__':
    unittest.main()