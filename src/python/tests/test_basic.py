"""Basic tests for Starburst99 Python implementation"""

import unittest
import numpy as np
from ..core.constants import *
from ..core.galaxy_module import GalaxyModel, ModelParameters
from ..models.imf import IMF
from ..utils.utilities import exp10, linear_interp


class TestConstants(unittest.TestCase):
    """Test physical constants"""
    
    def test_constants_values(self):
        """Test that constants have expected values"""
        self.assertAlmostEqual(PI, 3.14159265, places=6)
        self.assertAlmostEqual(SOLAR_MASS, 1.989e33)
        self.assertAlmostEqual(C_LIGHT, 2.99792458e10)


class TestUtilities(unittest.TestCase):
    """Test utility functions"""
    
    def test_exp10(self):
        """Test exp10 function"""
        self.assertAlmostEqual(exp10(0), 1.0)
        self.assertAlmostEqual(exp10(1), 10.0)
        self.assertAlmostEqual(exp10(2), 100.0)
        self.assertAlmostEqual(exp10(-1), 0.1)
    
    def test_linear_interp(self):
        """Test linear interpolation"""
        x_arr = np.array([0, 1, 2, 3])
        y_arr = np.array([0, 2, 4, 6])
        
        self.assertAlmostEqual(linear_interp(0.5, x_arr, y_arr), 1.0)
        self.assertAlmostEqual(linear_interp(1.5, x_arr, y_arr), 3.0)
        self.assertAlmostEqual(linear_interp(2.5, x_arr, y_arr), 5.0)
        
        # Test extrapolation
        self.assertAlmostEqual(linear_interp(-1, x_arr, y_arr), 0.0)
        self.assertAlmostEqual(linear_interp(4, x_arr, y_arr), 6.0)


class TestIMF(unittest.TestCase):
    """Test Initial Mass Function"""
    
    def setUp(self):
        """Set up test IMF"""
        self.imf = IMF(
            intervals=1,
            exponents=[2.35],
            mass_limits=[1.0, 100.0]
        )
    
    def test_imf_creation(self):
        """Test IMF object creation"""
        self.assertEqual(self.imf.intervals, 1)
        self.assertEqual(len(self.imf.exponents), 1)
        self.assertEqual(len(self.imf.mass_limits), 2)
    
    def test_imf_value(self):
        """Test IMF value calculation"""
        # Should return non-zero for valid masses
        self.assertGreater(self.imf.xi(10.0), 0.0)
        
        # Should return zero outside mass range
        self.assertEqual(self.imf.xi(0.5), 0.0)
        self.assertEqual(self.imf.xi(200.0), 0.0)
    
    def test_imf_integration(self):
        """Test IMF integration"""
        # Integration over full range should give total number
        total = self.imf.integrate(1.0, 100.0)
        self.assertGreater(total, 0.0)
        
        # Integration over partial range
        partial = self.imf.integrate(10.0, 50.0)
        self.assertGreater(partial, 0.0)
        self.assertLess(partial, total)


class TestGalaxyModel(unittest.TestCase):
    """Test GalaxyModel class"""
    
    def setUp(self):
        """Set up test galaxy model"""
        self.galaxy = GalaxyModel()
    
    def test_galaxy_creation(self):
        """Test galaxy model creation"""
        self.assertIsInstance(self.galaxy.model_params, ModelParameters)
        self.assertEqual(self.galaxy.current_time, 0.0)
        self.assertEqual(len(self.galaxy.wavelength), NP)
        self.assertEqual(len(self.galaxy.spectra), NP)
    
    def test_model_parameters(self):
        """Test model parameters"""
        params = ModelParameters(
            name="Test Model",
            sf_mode=1,
            total_mass=1e6,
            sf_rate=10.0
        )
        
        self.assertEqual(params.name, "Test Model")
        self.assertEqual(params.sf_mode, 1)
        self.assertEqual(params.total_mass, 1e6)
        self.assertEqual(params.sf_rate, 10.0)


if __name__ == '__main__':
    unittest.main()