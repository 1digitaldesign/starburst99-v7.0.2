"""Comprehensive tests for data_profiles module"""

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import unittest
from core.data_profiles import DataProfiles


class TestDataProfiles(unittest.TestCase):
    """Test DataProfiles class comprehensively"""
    
    def setUp(self):
        """Set up test DataProfiles instance"""
        self.data_profiles = DataProfiles()
    
    def test_initialization(self):
        """Test DataProfiles initialization"""
        self.assertFalse(self.data_profiles.is_initialized)
        
    def test_initialize_data_profiles(self):
        """Test data profiles initialization"""
        self.data_profiles.initialize_data_profiles()
        self.assertTrue(self.data_profiles.is_initialized)
        
    def test_multiple_initializations(self):
        """Test multiple initializations"""
        self.data_profiles.initialize_data_profiles()
        self.assertTrue(self.data_profiles.is_initialized)
        
        # Initialize again - should still work
        self.data_profiles.initialize_data_profiles()
        self.assertTrue(self.data_profiles.is_initialized)
        
    def test_is_initialized_method(self):
        """Test is_initialized method"""
        # Before initialization
        self.assertFalse(self.data_profiles.is_initialized_check())
        
        # After initialization
        self.data_profiles.initialize_data_profiles()
        self.assertTrue(self.data_profiles.is_initialized_check())
        
    def test_state_persistence(self):
        """Test that initialization state persists"""
        # Create new instance
        dp1 = DataProfiles()
        self.assertFalse(dp1.is_initialized)
        
        # Initialize it
        dp1.initialize_data_profiles()
        self.assertTrue(dp1.is_initialized)
        
        # Create another instance - should be independent
        dp2 = DataProfiles()
        self.assertFalse(dp2.is_initialized)


if __name__ == '__main__':
    unittest.main()