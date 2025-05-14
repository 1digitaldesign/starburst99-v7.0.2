"""Input parameter parser for Starburst99"""

import logging
from pathlib import Path
from typing import Dict, List, Any
import configparser
import json

from core.galaxy_module import ModelParameters


class InputParser:
    """Parser for Starburst99 input parameter files"""
    
    def __init__(self):
        self.logger = logging.getLogger(__name__)
    
    def read_input(self, filename: str) -> ModelParameters:
        """
        Read input parameters from file.
        
        Args:
            filename: Path to input parameter file
            
        Returns:
            ModelParameters object with parsed values
        """
        file_path = Path(filename)
        
        if not file_path.exists():
            raise FileNotFoundError(f"Input file not found: {filename}")
        
        # Determine file format and parse accordingly
        if file_path.suffix == '.json':
            return self._read_json_input(file_path)
        elif file_path.suffix == '.ini':
            return self._read_ini_input(file_path)
        else:
            # Default to original format
            return self._read_standard_input(file_path)
    
    def _read_standard_input(self, file_path: Path) -> ModelParameters:
        """Read standard Starburst99 input format"""
        params = ModelParameters()
        
        with open(file_path, 'r') as f:
            lines = f.readlines()
        
        # Parse the input file line by line
        # This would implement the original Fortran input format parsing
        # For now, we'll create a basic implementation
        
        line_idx = 0
        
        # Model name
        params.name = lines[line_idx].strip()
        line_idx += 1
        
        # Star formation parameters
        values = lines[line_idx].split()
        params.sf_mode = int(values[0])
        params.total_mass = float(values[1])
        if params.sf_mode > 0:
            params.sf_rate = float(values[2])
        line_idx += 1
        
        # IMF parameters
        params.num_intervals = int(lines[line_idx].strip())
        line_idx += 1
        
        params.exponents = []
        params.mass_limits = []
        for i in range(params.num_intervals):
            values = lines[line_idx].split()
            params.exponents.append(float(values[0]))
            params.mass_limits.append(float(values[1]))
            line_idx += 1
        
        # Add upper mass limit
        params.mass_limits.append(float(lines[line_idx].strip()))
        line_idx += 1
        
        # Cutoff masses
        values = lines[line_idx].split()
        params.sn_cutoff = float(values[0])
        params.bh_cutoff = float(values[1])
        line_idx += 1
        
        # Track and wind parameters
        values = lines[line_idx].split()
        params.metallicity_id = int(values[0])
        params.wind_id = int(values[1])
        line_idx += 1
        
        # Time grid parameters
        # Additional parsing would continue here...
        
        return params
    
    def _read_json_input(self, file_path: Path) -> ModelParameters:
        """Read JSON format input file"""
        with open(file_path, 'r') as f:
            data = json.load(f)
        
        params = ModelParameters()
        
        # Map JSON data to ModelParameters
        params.name = data.get('name', 'default')
        params.sf_mode = data.get('star_formation', {}).get('mode', 0)
        params.total_mass = data.get('star_formation', {}).get('total_mass', 1.0)
        params.sf_rate = data.get('star_formation', {}).get('rate', 1.0)
        
        # IMF parameters
        imf_data = data.get('imf', {})
        params.num_intervals = imf_data.get('num_intervals', 1)
        params.exponents = imf_data.get('exponents', [2.35])
        params.mass_limits = imf_data.get('mass_limits', [1.0, 100.0])
        params.sn_cutoff = imf_data.get('sn_cutoff', 8.0)
        params.bh_cutoff = imf_data.get('bh_cutoff', 120.0)
        
        # Model parameters
        model_data = data.get('model', {})
        params.metallicity_id = model_data.get('metallicity_id', 0)
        params.wind_id = model_data.get('wind_id', 0)
        
        return params
    
    def _read_ini_input(self, file_path: Path) -> ModelParameters:
        """Read INI format input file"""
        config = configparser.ConfigParser()
        config.read(file_path)
        
        params = ModelParameters()
        
        # Map INI sections to ModelParameters
        if 'general' in config:
            params.name = config.get('general', 'name', fallback='default')
        
        if 'star_formation' in config:
            sf = config['star_formation']
            params.sf_mode = sf.getint('mode', 0)
            params.total_mass = sf.getfloat('total_mass', 1.0)
            params.sf_rate = sf.getfloat('rate', 1.0)
        
        if 'imf' in config:
            imf = config['imf']
            params.num_intervals = imf.getint('num_intervals', 1)
            # Handle single or multiple exponents
            exp_str = imf.get('exponents', '2.35')
            params.exponents = [float(x.strip()) for x in exp_str.split(',')]
            # Handle mass limits
            limits_str = imf.get('mass_limits', '1.0,100.0')
            params.mass_limits = [float(x.strip()) for x in limits_str.split(',')]
            params.sn_cutoff = imf.getfloat('sn_cutoff', 8.0)
            params.bh_cutoff = imf.getfloat('bh_cutoff', 120.0)
        
        if 'model' in config:
            model = config['model']
            params.metallicity_id = model.getint('metallicity_id', 0)
            params.wind_id = model.getint('wind_id', 0)
        
        return params
    
    def get_default_parameters(self) -> ModelParameters:
        """Get default model parameters"""
        return ModelParameters(
            name="Default Model",
            sf_mode=1,
            total_mass=1.0,
            sf_rate=1.0,
            num_intervals=1,
            exponents=[2.35],
            mass_limits=[1.0, 100.0],
            sn_cutoff=8.0,
            bh_cutoff=120.0,
            metallicity_id=24,
            wind_id=0
        )
    
    def validate_parameters(self, params: ModelParameters) -> bool:
        """
        Validate input parameters.
        
        Args:
            params: ModelParameters to validate
            
        Returns:
            True if valid, False otherwise
        """
        # Implement validation logic
        if params.total_mass <= 0:
            self.logger.error("Total mass must be positive")
            return False
        
        if params.sf_mode > 0 and params.sf_rate <= 0:
            self.logger.error("Star formation rate must be positive for continuous mode")
            return False
        
        if params.num_intervals != len(params.exponents):
            self.logger.error("Number of IMF intervals doesn't match number of exponents")
            return False
        
        if len(params.mass_limits) != params.num_intervals + 1:
            self.logger.error("Incorrect number of mass limits for IMF")
            return False
        
        return True