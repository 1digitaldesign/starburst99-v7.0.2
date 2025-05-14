"""Physical constants used in Starburst99 calculations"""

import numpy as np

# Mathematical constants
PI = np.pi

# Astronomical constants
SOLAR_MASS = 1.989e33      # g
SOLAR_LUM = 3.826e33       # erg/s
YEAR_IN_SEC = 3.1557e7     # seconds
PARSEC = 3.0857e18         # cm
LSUN_MW = 4.0e10          # L_sun

# Atomic constants
K_BOLTZ = 1.380649e-16     # erg/K
H_PLANCK = 6.62607015e-27  # erg-s
C_LIGHT = 2.99792458e10    # cm/s
SIGMA_SB = 5.670374419e-5  # erg cm^-2 s^-1 K^-4

# Module-level parameters
NMAXINT = 10               # Maximum number of IMF intervals
NMAXINT1 = 11              # Maximum number of IMF mass boundaries
NP = 860                   # Size of standard wavelength grid
NP1 = 1415                 # Size of extended wavelength grid
NPGRID = 3000              # Size of mass grid

# File unit mapping (for compatibility with Fortran code)
UN_INPUT = 10              # Input parameter file
UN_OUTPUT = 11             # Main output file  
UN_SPECTRUM = 12           # Spectrum output
UN_QUANTA = 13             # Ionizing photon output
UN_SNR = 14                # Supernova rate output
UN_POWER = 15              # Wind power output
UN_SPTYP = 16              # Spectral types output
UN_YIELD = 17              # Chemical yields output
UN_UVLINE = 18             # UV line spectrum output
UN_COLOR = 19              # Color index output
UN_ATM = 20                # Atmosphere model data
UN_DEBUG = 21              # Debug output
UN_WRLINE = 22             # WR line output
UN_HIRES = 23              # High-resolution output
UN_FUSE = 24               # FUSE output