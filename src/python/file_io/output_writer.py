"""Output writer for Starburst99 results"""

import logging
from pathlib import Path
from typing import Dict, List, Any
import json
import numpy as np
from datetime import datetime

from core.galaxy_module import GalaxyModel


class OutputWriter:
    """Writer for Starburst99 output files"""
    
    def __init__(self):
        self.logger = logging.getLogger(__name__)
        self.output_dir = Path("output")
        self.output_dir.mkdir(exist_ok=True)
    
    def write_all_outputs(self, galaxy: GalaxyModel):
        """
        Write all output files.
        
        Args:
            galaxy: GalaxyModel instance with results
        """
        # Sanitize the model name for use in filenames
        sanitized_name = galaxy.model_params.name
        for char in [' ', '/', '&', '\\', ':', '*', '?', '"', '<', '>', '|']:
            sanitized_name = sanitized_name.replace(char, '_')
        base_name = sanitized_name.lower()
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Write different output files
        self._write_main_output(galaxy, f"{base_name}_{timestamp}")
        self._write_spectrum(galaxy, f"{base_name}_{timestamp}")
        self._write_quanta(galaxy, f"{base_name}_{timestamp}")
        self._write_summary_json(galaxy, f"{base_name}_{timestamp}")
    
    def _write_main_output(self, galaxy: GalaxyModel, base_name: str):
        """Write main output file with model parameters and results"""
        output_file = self.output_dir / f"{base_name}.output"
        
        with open(output_file, 'w') as f:
            # Header
            f.write("STARBURST99 v7.0.2 (Python Edition)\n")
            f.write("=" * 50 + "\n")
            f.write(f"Model: {galaxy.model_params.name}\n")
            f.write(f"Date: {datetime.now()}\n")
            f.write("\n")
            
            # Model parameters
            f.write("MODEL PARAMETERS\n")
            f.write("-" * 20 + "\n")
            f.write(f"Star Formation Mode: {galaxy.model_params.sf_mode}\n")
            f.write(f"Total Mass: {galaxy.model_params.total_mass:.2e} M_sun\n")
            if galaxy.model_params.sf_mode > 0:
                f.write(f"SFR: {galaxy.model_params.sf_rate:.2e} M_sun/yr\n")
            f.write("\n")
            
            # IMF parameters
            f.write("IMF PARAMETERS\n")
            f.write("-" * 20 + "\n")
            f.write(f"Number of intervals: {galaxy.model_params.num_intervals}\n")
            for i in range(galaxy.model_params.num_intervals):
                f.write(f"  Interval {i+1}: α={galaxy.model_params.exponents[i]:.2f}, "
                       f"M=[{galaxy.model_params.mass_limits[i]:.1f}, "
                       f"{galaxy.model_params.mass_limits[i+1]:.1f}] M_sun\n")
            f.write(f"SN cutoff: {galaxy.model_params.sn_cutoff:.1f} M_sun\n")
            f.write(f"BH cutoff: {galaxy.model_params.bh_cutoff:.1f} M_sun\n")
            f.write("\n")
            
            # Results summary
            f.write("RESULTS SUMMARY\n")
            f.write("-" * 20 + "\n")
            # Write key results here
            
        self.logger.info(f"Main output written to: {output_file}")
    
    def _write_spectrum(self, galaxy: GalaxyModel, base_name: str):
        """Write spectral energy distribution"""
        output_file = self.output_dir / f"{base_name}.spectrum"
        
        with open(output_file, 'w') as f:
            # Header
            f.write("# Wavelength(A)  Flux(erg/s/A)\n")
            
            # Write wavelength and flux data
            for i in range(len(galaxy.wavelength)):
                if galaxy.wavelength[i] > 0 and galaxy.spectra[i] > 0:
                    f.write(f"{galaxy.wavelength[i]:<12.3f}  {galaxy.spectra[i]:<12.5e}\n")
        
        self.logger.info(f"Spectrum written to: {output_file}")
    
    def _write_quanta(self, galaxy: GalaxyModel, base_name: str):
        """Write ionizing photon rates"""
        output_file = self.output_dir / f"{base_name}.quanta"
        
        with open(output_file, 'w') as f:
            # Header
            f.write("# Time(yr)  Q(H)  Q(HeI)  Q(HeII)\n")
            
            # Write ionizing photon data
            # Implementation would write actual calculated values
            
        self.logger.info(f"Ionizing photon rates written to: {output_file}")
    
    def _write_summary_json(self, galaxy: GalaxyModel, base_name: str):
        """Write summary in JSON format for easy parsing"""
        output_file = self.output_dir / f"{base_name}_summary.json"
        
        summary = {
            "model": {
                "name": galaxy.model_params.name,
                "timestamp": datetime.now().isoformat(),
                "version": "7.0.2-python"
            },
            "parameters": {
                "star_formation": {
                    "mode": galaxy.model_params.sf_mode,
                    "total_mass": galaxy.model_params.total_mass,
                    "rate": galaxy.model_params.sf_rate if galaxy.model_params.sf_mode > 0 else None
                },
                "imf": {
                    "num_intervals": galaxy.model_params.num_intervals,
                    "exponents": galaxy.model_params.exponents,
                    "mass_limits": galaxy.model_params.mass_limits,
                    "sn_cutoff": galaxy.model_params.sn_cutoff,
                    "bh_cutoff": galaxy.model_params.bh_cutoff
                },
                "tracks": {
                    "metallicity_id": galaxy.model_params.metallicity_id,
                    "wind_id": galaxy.model_params.wind_id
                }
            },
            "results": {
                # Add computed results here
                "final_time": galaxy.current_time,
                "num_timesteps": galaxy.time_step
            }
        }
        
        with open(output_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        self.logger.info(f"Summary JSON written to: {output_file}")