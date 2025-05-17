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
        # This implements the original Fortran input format parsing
        
        line_idx = 0
        
        # When the format includes label lines, skip them
        if len(lines) > 0 and "MODEL DESIGNATION:" in lines[0]:
            # Skip label line and read model name
            line_idx += 1  # Skip "MODEL DESIGNATION:" line
            if line_idx >= len(lines):
                self.logger.error("Unexpected end of file while reading model name")
                raise ValueError("Invalid input file format")
            params.name = lines[line_idx].strip()
            line_idx += 1
            
            # Skip SF mode label and read value
            line_idx += 1  # Skip "CONTINUOUS STAR FORMATION..." line
            params.sf_mode = int(lines[line_idx].strip())
            line_idx += 1
        else:
            # Direct format without labels - for testing
            params.name = lines[line_idx].strip()
            line_idx += 1
            
            # Parse SF mode, total mass, and SFR from a single line
            parts = lines[line_idx].strip().split()
            params.sf_mode = int(parts[0])
            params.total_mass = float(parts[1])
            # For instantaneous SF (mode 0), SFR is not provided
            if len(parts) > 2:
                params.sf_rate = float(parts[2])
            else:
                params.sf_rate = 0.0
            line_idx += 1
        
        # If we used labeled format above, continue with that format
        if "MODEL DESIGNATION:" in lines[0]:
            # Skip total mass label and read value
            line_idx += 1  # Skip "TOTAL STELLAR MASS..." line
            params.total_mass = float(lines[line_idx].strip())
            line_idx += 1
            
            # Skip SFR label and read value
            line_idx += 1  # Skip "SFR..." line
            params.sf_rate = float(lines[line_idx].strip())
            line_idx += 1
            
            # Skip IMF intervals label and read value
            line_idx += 1  # Skip "NUMBER OF INTERVALS..." line
            params.num_intervals = int(lines[line_idx].strip())
            line_idx += 1
        else:
            # Direct format without labels
            params.num_intervals = int(lines[line_idx].strip())
            line_idx += 1
        
        # Read exponents
        if "MODEL DESIGNATION:" in lines[0]:
            # Skip exponents label and read values
            line_idx += 1  # Skip "IMF EXPONENTS..." line
            exp_values = lines[line_idx].strip().split(',')
            params.exponents = [float(x.strip()) for x in exp_values]
            line_idx += 1
            
            # Skip mass limits label and read values
            line_idx += 1  # Skip "MASS BOUNDARIES..." line
            mass_values = lines[line_idx].strip().split(',')
            params.mass_limits = [float(x.strip()) for x in mass_values]
            line_idx += 1
        else:
            # Direct format - parse alternating exponents and mass limits
            params.exponents = []
            params.mass_limits = []
            
            # Read pairs of exponent and mass limit
            for i in range(params.num_intervals):
                exp_mass = lines[line_idx].strip().split()
                params.exponents.append(float(exp_mass[0]))
                params.mass_limits.append(float(exp_mass[1]))
                line_idx += 1
            
            # Read the final mass limit
            final_mass = float(lines[line_idx].strip())
            params.mass_limits.append(final_mass)
            line_idx += 1
        
        # Read SN and BH cutoffs
        if "MODEL DESIGNATION:" in lines[0]:
            # Skip SN cutoff label and read value
            line_idx += 1  # Skip "SUPERNOVA CUT-OFF..." line
            params.sn_cutoff = float(lines[line_idx].strip())
            line_idx += 1
            
            # Skip BH cutoff label and read value
            line_idx += 1  # Skip "BLACK HOLE CUT-OFF..." line
            params.bh_cutoff = float(lines[line_idx].strip())
            line_idx += 1
        else:
            # Direct format - parse cutoffs
            cutoff_parts = lines[line_idx].strip().split()
            params.sn_cutoff = float(cutoff_parts[0])
            params.bh_cutoff = float(cutoff_parts[1])
            line_idx += 1
        
        # Read metallicity and wind model
        if "MODEL DESIGNATION:" in lines[0]:
            # Skip metallicity label and read value (may include text after number)
            line_idx += 1  # Skip "METALLICITY + TRACKS:" line
            while line_idx < len(lines) and not lines[line_idx].strip()[0].isdigit():
                line_idx += 1  # Skip additional metallicity description lines
            if line_idx < len(lines):
                params.metallicity_id = int(lines[line_idx].strip().split()[0])
            line_idx += 1
            
            # Skip wind model label and read value
            line_idx += 1  # Skip "WIND MODEL..." line
            if line_idx < len(lines):
                params.wind_id = int(lines[line_idx].strip())
            line_idx += 1
        else:
            # Direct format - parse metallicity and wind model
            if line_idx < len(lines):
                metal_parts = lines[line_idx].strip().split()
                params.metallicity_id = int(metal_parts[0])
                params.wind_id = int(metal_parts[1])
                line_idx += 1
            
            # Parse time grid parameters if present
            if line_idx + 2 < len(lines):
                time_parts = lines[line_idx].strip().split()
                if len(time_parts) >= 3:  # For example: 100 1001
                    params.time_steps = int(time_parts[0])
                    # params.max_time = float(time_parts[1])
                    line_idx += 1
                    
                    # Parse more time parameters if present
                    if line_idx < len(lines):
                        more_time_parts = lines[line_idx].strip().split()
                        if len(more_time_parts) >= 3:
                            # Start time, end time, num steps
                            start_time = float(more_time_parts[0])
                            end_time = float(more_time_parts[1])
                            num_steps = float(more_time_parts[2])
                            params.time_grid = self._generate_time_grid(start_time, end_time, int(num_steps))
                            params.time_steps = int(num_steps)
                            params.max_time = end_time
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
        
        # Output parameters  
        output_data = data.get('output', {})
        params.output_directory = output_data.get('directory')
        params.output_prefix = output_data.get('prefix', 'model')
        
        # Time grid parameters
        time_data = data.get('time', {})
        if 'start' in time_data and 'end' in time_data and 'num_steps' in time_data:
            start = time_data['start']
            end = time_data['end']
            num_steps = time_data['num_steps']
            params.time_grid = self._generate_time_grid(start, end, num_steps)
            params.time_steps = num_steps
            params.max_time = end
        
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
    
    def _generate_time_grid(self, start: float, end: float, num_steps: int) -> List[float]:
        """Generate time grid with logarithmic spacing"""
        if num_steps <= 1:
            return [end]
        
        # Use logarithmic spacing for time grid
        import numpy as np
        return list(np.logspace(np.log10(start), np.log10(end), num_steps))
    
    def validate_parameters(self, params: ModelParameters) -> bool:
        """
        Validate model parameters for consistency.
        
        Args:
            params: ModelParameters to validate
            
        Returns:
            True if parameters are valid, False otherwise
        """
        valid = True
        
        # Check total mass
        if params.total_mass <= 0:
            self.logger.error("Total mass must be positive")
            valid = False
        
        # Check star formation rate for continuous mode
        if params.sf_mode > 0 and params.sf_rate <= 0:
            self.logger.error("Star formation rate must be positive for continuous mode")
            valid = False
        
        # Check IMF parameters
        if len(params.exponents) != params.num_intervals:
            self.logger.error("Number of IMF intervals doesn't match number of exponents")
            valid = False
        
        if len(params.mass_limits) != params.num_intervals + 1:
            self.logger.error("Incorrect number of mass limits for IMF intervals")
            valid = False
        
        return valid
    
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