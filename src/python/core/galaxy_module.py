"""Galaxy module for shared data and constants"""

import numpy as np
from typing import Optional, List, Tuple
from dataclasses import dataclass, field
import logging
from pathlib import Path

from core.constants import *
from utils.utilities import exp10, linear_interp


@dataclass
class ModelParameters:
    """Type containing all model configuration parameters"""
    name: str = "default"
    sf_mode: int = 0                    # Star formation mode: >0 continuous, <=0 fixed mass
    total_mass: float = 1.0             # Total stellar mass (10^6 solar masses)
    sf_rate: float = 1.0                # Star formation rate (solar masses per year)
    
    # IMF parameters
    num_intervals: int = 1              # Number of intervals for the IMF
    exponents: List[float] = field(default_factory=lambda: [2.35])
    mass_limits: List[float] = field(default_factory=lambda: [1.0, 100.0])
    sn_cutoff: float = 8.0              # Supernova cut-off mass (solar masses)
    bh_cutoff: float = 120.0            # Black hole cut-off mass (solar masses)
    
    # Model selection parameters
    metallicity_id: int = 0             # Metallicity + tracks identifier
    wind_id: int = 0                    # Mass loss rate selection
    time_grid: List[float] = field(default_factory=list)
    time_steps: int = 0
    max_time: float = 0.0
    
    # Spectral parameters
    atmosphere_models: int = 0          # Stellar atmosphere models
    spectral_library: int = 0           # Spectral library selection
    include_lines: int = 0              # Include spectral lines
    velocity_threshold: int = 0         # Velocity threshold
    include_rsg: int = 0                # Include red supergiants


@dataclass
class TrackData:
    """Type for storing evolutionary track information"""
    num_tracks: int = 0
    num_time_steps: List[int] = field(default_factory=list)
    num_wr_types: List[int] = field(default_factory=list)
    masses: List[float] = field(default_factory=list)
    time_grid: List[List[float]] = field(default_factory=list)
    luminosities: List[List[float]] = field(default_factory=list)
    temperatures: List[List[float]] = field(default_factory=list)
    stellar_types: List[List[int]] = field(default_factory=list)
    mass_loss_rates: List[List[float]] = field(default_factory=list)


class GalaxyModel:
    """
    Galaxy code module for shared data and constants.
    
    This class provides shared variables, constants, and utility functions
    for the galaxy stellar population synthesis code, using modern Python
    features to replace outdated Fortran COMMON blocks and implicit global variables.
    """
    
    def __init__(self):
        """Initialize the GalaxyModel with default values"""
        # Model parameters
        self.model_params = ModelParameters()
        
        # Track data
        self.tracks = TrackData()
        
        # Output arrays
        self.wavelength = np.zeros(NP)
        self.spectra = np.zeros(NP)
        self.mass_grid = np.zeros(NPGRID)
        self.density_grid = np.zeros(NPGRID)
        
        # Time and iteration parameters
        self.current_time = 0.0
        self.time_step = 0
        self.iteration = 0
        
        # Path configurations
        self.data_dir = Path("data")
        self.output_dir = Path("output")
        
        # Logger
        self.logger = logging.getLogger(__name__)
        
    def init_module(self):
        """Initialize module variables and set up logging"""
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
        self.logger.info("Galaxy module initialized")
        
    def cleanup_module(self):
        """Clean up resources and close any open files"""
        self.logger.info("Galaxy module cleanup")
        
    def open_file(self, filename: str, mode: str = 'r') -> Optional[object]:
        """
        Open a file with error handling.
        
        Args:
            filename: Path to the file
            mode: File opening mode ('r', 'w', etc.)
            
        Returns:
            File handle or None if error
        """
        try:
            file_path = Path(filename)
            if mode == 'r' and not file_path.exists():
                self.logger.error(f"File not found: {filename}")
                return None
            
            return open(file_path, mode)
        except Exception as e:
            self.logger.error(f"Error opening file {filename}: {e}")
            return None
    
    def read_parameters(self, filename: str):
        """Read model parameters from input file"""
        # Implementation will read the input parameter file
        # and populate self.model_params
        pass
    
    def read_tracks(self):
        """Read evolutionary tracks"""
        # Implementation will read track data and populate self.tracks
        pass
    
    def compute_sed(self):
        """Compute spectral energy distribution"""
        # Main calculation routine
        pass
    
    def write_output(self):
        """Write results to output files"""
        # Output writing routine
        pass