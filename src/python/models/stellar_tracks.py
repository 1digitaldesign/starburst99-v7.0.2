"""Stellar evolution track handling"""

import numpy as np
from pathlib import Path
import logging
from typing import Dict, List, Tuple


class StellarTracks:
    """Class for reading and interpolating stellar evolution tracks"""
    
    def __init__(self, metallicity_id: int = None, data_dir: Path = Path("data/tracks")):
        """
        Initialize stellar tracks reader.
        
        Args:
            metallicity_id: Metallicity identifier
            data_dir: Directory containing track data files
        """
        self.metallicity_id = metallicity_id  # Store for compatibility
        self.data_dir = data_dir
        self.logger = logging.getLogger(__name__)
        
        # Storage for track data
        self.tracks = {}
        self.metallicities = []
        
    def load_tracks(self, metallicity: str):
        """
        Load stellar evolution tracks for a given metallicity.
        
        Args:
            metallicity: Metallicity identifier (e.g., "Z0020v00")
        """
        track_file = self.data_dir / f"{metallicity}.txt"
        
        if not track_file.exists():
            self.logger.error(f"Track file not found: {track_file}")
            raise FileNotFoundError(f"Track file not found: {track_file}")
        
        self.logger.info(f"Loading tracks from: {track_file}")
        
        # Parse the track file
        with open(track_file, 'r') as f:
            lines = f.readlines()
        
        # Initialize data structures
        current_mass = None
        mass_data = {}
        
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            
            # Skip empty lines and comments
            if not line or line.startswith('#'):
                i += 1
                continue
            
            # Check if this is a mass header
            if line.startswith('M='):
                # Parse mass value
                current_mass = float(line.split('=')[1].split()[0])
                mass_data[current_mass] = {
                    'time': [],
                    'luminosity': [],
                    'temperature': [],
                    'radius': [],
                    'mass_loss_rate': [],
                    'stellar_type': []
                }
                i += 1
                continue
            
            # Parse data line
            if current_mass is not None:
                values = line.split()
                if len(values) >= 6:
                    mass_data[current_mass]['time'].append(float(values[0]))
                    mass_data[current_mass]['luminosity'].append(float(values[1]))
                    mass_data[current_mass]['temperature'].append(float(values[2]))
                    mass_data[current_mass]['radius'].append(float(values[3]))
                    mass_data[current_mass]['mass_loss_rate'].append(float(values[4]))
                    mass_data[current_mass]['stellar_type'].append(int(values[5]))
            
            i += 1
        
        # Convert lists to numpy arrays
        for mass in mass_data:
            for key in mass_data[mass]:
                mass_data[mass][key] = np.array(mass_data[mass][key])
        
        # Store the data
        self.tracks[metallicity] = mass_data
        self.metallicities.append(metallicity)
        
        self.logger.info(f"Loaded {len(mass_data)} mass tracks for {metallicity}")
    
    def interpolate_track(self, 
                         mass: float, 
                         time: float, 
                         metallicity: str) -> Dict[str, float]:
        """
        Interpolate stellar properties at given mass and time.
        
        Args:
            mass: Stellar mass in solar masses
            time: Age in years
            metallicity: Metallicity identifier
            
        Returns:
            Dictionary of interpolated stellar properties
        """
        if metallicity not in self.tracks:
            self.load_tracks(metallicity)
        
        track_data = self.tracks[metallicity]
        
        # Find bracketing masses
        masses = sorted(track_data.keys())
        
        if mass <= masses[0]:
            # Use lowest mass track
            track = track_data[masses[0]]
        elif mass >= masses[-1]:
            # Use highest mass track
            track = track_data[masses[-1]]
        else:
            # Interpolate between two tracks
            idx = np.searchsorted(masses, mass)
            m1, m2 = masses[idx-1], masses[idx]
            track1 = track_data[m1]
            track2 = track_data[m2]
            
            # Linear interpolation in mass
            w1 = (m2 - mass) / (m2 - m1)
            w2 = (mass - m1) / (m2 - m1)
            
            # Interpolate time-dependent properties
            props = {}
            for key in ['luminosity', 'temperature', 'radius', 'mass_loss_rate']:
                # Interpolate in time for each track
                val1 = np.interp(time, track1['time'], track1[key])
                val2 = np.interp(time, track2['time'], track2[key])
                # Then interpolate between tracks
                props[key] = w1 * val1 + w2 * val2
            
            # Handle stellar type (discrete value)
            t1_idx = np.searchsorted(track1['time'], time)
            t2_idx = np.searchsorted(track2['time'], time)
            t1_idx = np.clip(t1_idx, 0, len(track1['stellar_type'])-1)
            t2_idx = np.clip(t2_idx, 0, len(track2['stellar_type'])-1)
            
            # Use the type from the dominant mass
            if w1 > 0.5:
                props['stellar_type'] = track1['stellar_type'][t1_idx]
            else:
                props['stellar_type'] = track2['stellar_type'][t2_idx]
            
            return props
        
        # Single track interpolation
        props = {}
        for key in ['luminosity', 'temperature', 'radius', 'mass_loss_rate']:
            props[key] = np.interp(time, track['time'], track[key])
        
        # Handle stellar type
        t_idx = np.searchsorted(track['time'], time)
        t_idx = np.clip(t_idx, 0, len(track['stellar_type'])-1)
        props['stellar_type'] = track['stellar_type'][t_idx]
        
        return props
    
    def get_lifetime(self, mass: float, metallicity: str) -> float:
        """
        Get the main sequence lifetime for a given mass.
        
        Args:
            mass: Stellar mass in solar masses
            metallicity: Metallicity identifier
            
        Returns:
            Main sequence lifetime in years
        """
        if metallicity not in self.tracks:
            self.load_tracks(metallicity)
        
        track_data = self.tracks[metallicity]
        masses = sorted(track_data.keys())
        
        # Find appropriate mass track
        if mass <= masses[0]:
            track = track_data[masses[0]]
        elif mass >= masses[-1]:
            track = track_data[masses[-1]]
        else:
            # Interpolate between tracks
            idx = np.searchsorted(masses, mass)
            m1, m2 = masses[idx-1], masses[idx]
            life1 = track_data[m1]['time'][-1]
            life2 = track_data[m2]['time'][-1]
            
            # Log interpolation for lifetime
            log_life = np.interp(np.log10(mass), 
                                np.log10([m1, m2]), 
                                np.log10([life1, life2]))
            return 10**log_life
        
        return track['time'][-1]
    
    def _linear_interpolate(self, x: np.ndarray, y: np.ndarray, xi: float) -> float:
        """
        Linear interpolation for x and y arrays at point xi.
        
        Args:
            x: x values (must be sorted)
            y: y values
            xi: interpolation point
            
        Returns:
            Interpolated value
        """
        # Handle edge cases for line 184 test
        if len(x) < 2:
            return y[0] if len(y) > 0 else 0.0
            
        # If xi is outside bounds
        if xi <= x[0]:
            return y[0]
        if xi >= x[-1]:
            return y[-1]
            
        # Find bracketing points
        idx = np.searchsorted(x, xi)
        if idx >= len(x):
            idx = len(x) - 1
            
        # Handle exact match or very close values (edge case for line 184)
        if idx > 0 and np.abs(x[idx-1] - xi) < 1e-10:
            return y[idx-1]
        if idx < len(x) and np.abs(x[idx] - xi) < 1e-10:
            return y[idx]
            
        # Linear interpolation
        if idx > 0 and idx < len(x):
            x1, x2 = x[idx-1], x[idx]
            y1, y2 = y[idx-1], y[idx]
            # Line 184 edge case: handle very close x values
            if np.abs(x2 - x1) < 1e-10:
                return (y1 + y2) / 2.0
            return y1 + (xi - x1) * (y2 - y1) / (x2 - x1)
        
        return y[idx] if idx < len(y) else y[-1]