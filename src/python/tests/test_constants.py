"""Comprehensive tests for constants module"""

import unittest
import numpy as np
from ..core.constants import *


class TestConstants(unittest.TestCase):
    """Test all physical constants"""
    
    def test_mathematical_constants(self):
        """Test mathematical constants"""
        self.assertAlmostEqual(PI, 3.14159265358979323846, places=10)
        
    def test_astronomical_constants(self):
        """Test astronomical constants"""
        self.assertAlmostEqual(SOLAR_MASS, 1.989e33, places=30)
        self.assertAlmostEqual(SOLAR_LUM, 3.826e33, places=30)
        self.assertAlmostEqual(YEAR_IN_SEC, 3.1557e7, places=3)
        self.assertAlmostEqual(PARSEC, 3.0857e18, places=16)
        self.assertAlmostEqual(LSUN_MW, 4.0e10, places=8)
        
    def test_atomic_constants(self):
        """Test atomic constants"""
        self.assertAlmostEqual(K_BOLTZ, 1.380649e-16, places=20)
        self.assertAlmostEqual(H_PLANCK, 6.62607015e-27, places=33)
        self.assertAlmostEqual(C_LIGHT, 2.99792458e10, places=6)
        self.assertAlmostEqual(SIGMA_SB, 5.670374419e-5, places=12)
        
    def test_module_parameters(self):
        """Test module-level parameters"""
        self.assertEqual(NMAXINT, 10)
        self.assertEqual(NMAXINT1, 11)
        self.assertEqual(NP, 860)
        self.assertEqual(NP1, 1415)
        self.assertEqual(NPGRID, 3000)
        
    def test_file_units(self):
        """Test file unit numbers"""
        self.assertEqual(UN_INPUT, 10)
        self.assertEqual(UN_OUTPUT, 11)
        self.assertEqual(UN_SPECTRUM, 12)
        self.assertEqual(UN_QUANTA, 13)
        self.assertEqual(UN_SNR, 14)
        self.assertEqual(UN_POWER, 15)
        self.assertEqual(UN_SPTYP, 16)
        self.assertEqual(UN_YIELD, 17)
        self.assertEqual(UN_UVLINE, 18)
        self.assertEqual(UN_COLOR, 19)
        self.assertEqual(UN_ATM, 20)
        self.assertEqual(UN_DEBUG, 21)
        self.assertEqual(UN_WRLINE, 22)
        self.assertEqual(UN_HIRES, 23)
        self.assertEqual(UN_FUSE, 24)


if __name__ == '__main__':
    unittest.main()