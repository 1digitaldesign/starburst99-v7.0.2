"""
Starburst99 Python Implementation
================================

A modern Python implementation of the Starburst99 stellar population 
synthesis code, originally written in Fortran.

This package provides modules for:
- Stellar population synthesis calculations
- Spectral energy distribution modeling
- Stellar feedback and chemical yield computations
"""

__version__ = '7.0.2-python'
__author__ = 'Starburst99 Team'

from . import core
from . import file_io
from . import models
from . import utils