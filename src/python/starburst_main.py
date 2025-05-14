#!/usr/bin/env python3
"""
GALAXY - Stellar Population Synthesis Code
=========================================

Main program for galaxy/starburst population synthesis code (Starburst99).

This program computes observable parameters for populations of massive stars, 
including spectral energy distributions, stellar feedback, and chemical yields.

Original version: Claus Leitherer (August 1998)
Last major update: August 2014
Python version: [2024]
"""

import sys
import argparse
import logging
from datetime import datetime
from pathlib import Path

import numpy as np

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from core.galaxy_module import GalaxyModel
from core.data_profiles import DataProfiles
from file_io.input_parser import InputParser
from file_io.output_writer import OutputWriter


class Starburst99:
    """Main class for Starburst99 stellar population synthesis calculations"""
    
    def __init__(self, input_file: str = None):
        """
        Initialize Starburst99 with optional input file.
        
        Args:
            input_file: Path to input parameter file
        """
        self.galaxy = GalaxyModel()
        self.data_profiles = DataProfiles()
        self.input_parser = InputParser()
        self.output_writer = OutputWriter()
        
        self.input_file = input_file
        self.start_time = datetime.now()
        
        # Set up logging
        self._setup_logging()
        
    def _setup_logging(self):
        """Configure logging for the application"""
        log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        logging.basicConfig(
            level=logging.INFO,
            format=log_format,
            handlers=[
                logging.FileHandler('starburst99.log'),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger('Starburst99')
        
    def run(self):
        """Main execution routine"""
        self.logger.info("GALAXY - Stellar Population Synthesis Code (Python Edition)")
        self.logger.info("=" * 56)
        self.logger.info(f"Run started at: {self.start_time}")
        self.logger.info("")
        
        try:
            # Initialize modules
            self.logger.info("Initializing modules...")
            self.galaxy.init_module()
            self.data_profiles.initialize_data_profiles()
            
            # Read input parameters
            self.logger.info("Reading input parameters...")
            if self.input_file:
                self.galaxy.model_params = self.input_parser.read_input(self.input_file)
            else:
                # Use default parameters or prompt for input
                self.galaxy.model_params = self.input_parser.get_default_parameters()
            
            # Read evolutionary tracks
            self.logger.info("Reading evolutionary tracks...")
            self._read_tracks()
            
            # Set metallicity string for filenames
            self._set_metallicity_string()
            
            # Read atmospheric and opacity data
            self.logger.info("Reading atmospheric and opacity data...")
            self._read_atmosphere_data()
            
            # Main calculation loop
            self.logger.info("Starting main calculations...")
            self._main_calculation_loop()
            
            # Write output files
            self.logger.info("Writing output files...")
            self._write_output()
            
            # Cleanup
            self.logger.info("Cleaning up...")
            self.galaxy.cleanup_module()
            
            end_time = datetime.now()
            duration = end_time - self.start_time
            self.logger.info(f"Run completed at: {end_time}")
            self.logger.info(f"Total runtime: {duration}")
            
        except Exception as e:
            self.logger.error(f"Error during execution: {e}", exc_info=True)
            sys.exit(1)
    
    def _read_tracks(self):
        """Read evolutionary tracks based on metallicity"""
        # Implementation will load track data files
        track_file = self._get_track_filename()
        self.galaxy.read_tracks(track_file)
    
    def _get_track_filename(self):
        """Determine track filename based on metallicity ID"""
        # Map metallicity_id to track file
        metallicity_map = {
            11: "Z0020v00.txt", 21: "Z0020v40.txt",
            12: "Z0020v00.txt", 22: "Z0020v40.txt",
            13: "Z0140v00.txt", 23: "Z0140v40.txt",
            14: "Z0140v00.txt", 24: "Z0140v40.txt",
            15: "Z0140v00.txt", 25: "Z0140v40.txt",
        }
        
        filename = metallicity_map.get(self.galaxy.model_params.metallicity_id, "Z0140v00.txt")
        return self.galaxy.data_dir / "tracks" / filename
    
    def _set_metallicity_string(self):
        """Set metallicity string for filenames based on selected tracks"""
        z_id = self.galaxy.model_params.metallicity_id
        
        metallicity_strings = {
            (11, 21, 31, 41): ('m13', '001'),
            (12, 22, 32, 42): ('m07', '004'),
            (13, 23, 33, 43): ('m04', '008'),
            (14, 24, 34, 44): ('p00', '020'),
            (15, 25, 35, 45): ('p03', '040'),
            (51, 61): ('m13', '001'),
            (52, 62): ('m07', '004'),
            (53, 63): ('m04', '008'),
            (54, 64): ('p00', '020'),
            (55, 65): ('p03', '040'),
        }
        
        for ids, strings in metallicity_strings.items():
            if z_id in ids:
                self.namfi3, self.nam = strings
                break
        else:
            self.namfi3, self.nam = 'p00', '020'
            
        # Override if needed based on wind scale
        if self.galaxy.model_params.wind_id < 0:
            self.nam = '020'
    
    def _read_atmosphere_data(self):
        """Read atmospheric and opacity data"""
        # Read Lejeune atmospheres
        atm_file = self.galaxy.data_dir / "lejeune" / f"lcb97_{self.namfi3}.flu"
        
        if not atm_file.exists():
            self.logger.error(f"Cannot find Lejeune atmosphere file: {atm_file}")
            sys.exit(1)
            
        # Read atmosphere data
        # Implementation will parse the atmosphere file
        
    def _main_calculation_loop(self):
        """Main calculation loop for population synthesis"""
        # Implementation of the main computational loop
        time_steps = self.galaxy.model_params.time_steps
        
        for step in range(time_steps):
            self.logger.info(f"Processing time step {step + 1}/{time_steps}")
            
            # Compute stellar populations
            self._compute_stellar_population(step)
            
            # Compute spectra
            self._compute_spectra(step)
            
            # Compute feedback quantities
            self._compute_feedback(step)
            
            # Update time
            self.galaxy.current_time = self.galaxy.model_params.time_grid[step]
    
    def _compute_stellar_population(self, step: int):
        """Compute stellar population for current time step"""
        # Implementation of stellar population calculations
        pass
    
    def _compute_spectra(self, step: int):
        """Compute synthetic spectra for current time step"""
        # Implementation of spectral synthesis
        pass
    
    def _compute_feedback(self, step: int):
        """Compute stellar feedback quantities"""
        # Implementation of feedback calculations
        pass
    
    def _write_output(self):
        """Write all output files"""
        self.output_writer.write_all_outputs(self.galaxy)


def main():
    """Main entry point for the application"""
    parser = argparse.ArgumentParser(
        description='Starburst99 - Stellar Population Synthesis Code',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        'input_file',
        nargs='?',
        help='Input parameter file'
    )
    
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose output'
    )
    
    parser.add_argument(
        '--version',
        action='version',
        version='Starburst99 v7.0.2 (Python Edition)'
    )
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Create and run Starburst99
    starburst = Starburst99(input_file=args.input_file)
    starburst.run()


if __name__ == '__main__':
    main()