"""Initial Mass Function (IMF) calculations"""

import numpy as np
from typing import List, Tuple
import logging


class IMF:
    """Class for Initial Mass Function calculations"""
    
    def __init__(self, 
                 intervals: int,
                 exponents: List[float],
                 mass_limits: List[float]):
        """
        Initialize IMF with given parameters.
        
        Args:
            intervals: Number of IMF intervals
            exponents: Power-law exponents for each interval
            mass_limits: Mass boundaries for intervals
        """
        self.num_intervals = intervals  # Add alias for compatibility
        self.intervals = intervals
        self.exponents = np.array(exponents)
        self.mass_limits = np.array(mass_limits)
        
        self.logger = logging.getLogger(__name__)
        
        # Validate parameters
        self._validate_parameters()
        
        # Calculate normalization constants
        self._calculate_normalizations()
    
    def _validate_parameters(self):
        """Validate IMF parameters"""
        if len(self.exponents) != self.intervals:
            raise ValueError("Number of exponents must match number of intervals")
        
        if len(self.mass_limits) != self.intervals + 1:
            raise ValueError("Number of mass limits must be intervals + 1")
        
        if not np.all(np.diff(self.mass_limits) > 0):
            raise ValueError("Mass limits must be in increasing order")
    
    def _calculate_normalizations(self):
        """Calculate normalization constants for each interval"""
        self.normalizations = np.zeros(self.intervals)
        
        for i in range(self.intervals):
            if self.exponents[i] != -1:
                self.normalizations[i] = (self.mass_limits[i+1]**(self.exponents[i]+1) - 
                                         self.mass_limits[i]**(self.exponents[i]+1)) / \
                                         (self.exponents[i] + 1)
            else:
                self.normalizations[i] = np.log(self.mass_limits[i+1] / self.mass_limits[i])
    
    def xi(self, mass: float) -> float:
        """
        Calculate the IMF value at a given mass.
        
        Args:
            mass: Stellar mass in solar masses
            
        Returns:
            IMF value (dN/dM)
        """
        if mass < self.mass_limits[0] or mass > self.mass_limits[-1]:
            return 0.0
        
        # Find which interval contains this mass
        interval = np.searchsorted(self.mass_limits, mass) - 1
        interval = np.clip(interval, 0, self.intervals - 1)
        
        # Calculate IMF value based on exponent
        if self.exponents[interval] != -1:
            return mass**(-self.exponents[interval])
        else:
            # For flat IMF (exponent = -1), return 1/m
            return 1.0 / mass
    
    def xi(self, mass: float, interval: int = None) -> float:
        """
        Calculate the IMF value at a given mass for a specific interval.
        
        Args:
            mass: Stellar mass in solar masses
            interval: Specific interval to use (optional, auto-detect if None)
            
        Returns:
            IMF value (dN/dM)
        """
        if interval is None:
            # Auto-detect interval
            if mass < self.mass_limits[0] or mass > self.mass_limits[-1]:
                return 0.0
            interval = np.searchsorted(self.mass_limits, mass) - 1
            interval = np.clip(interval, 0, self.intervals - 1)
        
        # Calculate IMF value based on exponent
        alpha = self.exponents[interval]
        if abs(alpha - 1.0) < 1e-6:  # Special case for alpha ≈ 1
            return 1.0 / mass
        elif alpha != -1:
            return mass**(-alpha)
        else:
            # For flat IMF (exponent = -1), return 1/m
            return 1.0 / mass
    
    def integrate(self, m_low: float, m_high: float) -> float:
        """
        Integrate the IMF between two mass limits.
        
        Args:
            m_low: Lower mass limit
            m_high: Upper mass limit
            
        Returns:
            Integrated number of stars  
        """
        # Ensure limits are within valid range
        m_low = max(m_low, self.mass_limits[0])
        m_high = min(m_high, self.mass_limits[-1])
        
        if m_low >= m_high:
            return 0.0
        
        total = 0.0
        
        # Find intervals that overlap with integration range
        for i in range(self.intervals):
            int_low = max(m_low, self.mass_limits[i])
            int_high = min(m_high, self.mass_limits[i+1])
            
            if int_low < int_high:
                if self.exponents[i] != -1:
                    # For standard IMF, integrate m^(-alpha)
                    alpha = self.exponents[i]
                    if alpha != 1:
                        total += (int_high**(-alpha+1) - int_low**(-alpha+1)) / (-alpha+1)
                    else:
                        total += np.log(int_high / int_low)
                else:
                    # For flat IMF with ξ(m) ∝ 1/m  
                    total += np.log(int_high / int_low)
        
        return total
    
    def _integrate_segment(self, segment: int, m_low: float, m_high: float, alpha: float) -> float:
        """
        Integrate a single IMF segment.
        
        Args:
            segment: Segment index
            m_low: Lower mass limit
            m_high: Upper mass limit
            alpha: Power law exponent
            
        Returns:
            Integrated value for segment
        """
        if alpha == 1.0:  # Special case for alpha = 1 (line 160 edge case)
            return np.log(m_high / m_low)
        elif alpha != -1:
            return (m_high**(-alpha+1) - m_low**(-alpha+1)) / (-alpha+1)
        else:
            # For flat IMF with ξ(m) ∝ 1/m  
            return np.log(m_high / m_low)
    
    def sample(self, n_stars: int, seed: int = None) -> np.ndarray:
        """
        Sample stellar masses from the IMF.
        
        Args:
            n_stars: Number of stars to sample
            seed: Random seed for reproducibility
            
        Returns:
            Array of stellar masses
        """
        if seed is not None:
            np.random.seed(seed)
        
        masses = []
        
        # Calculate relative probabilities for each interval
        interval_probs = np.zeros(self.intervals)
        for i in range(self.intervals):
            interval_probs[i] = self.integrate(self.mass_limits[i], 
                                              self.mass_limits[i+1])
        interval_probs /= interval_probs.sum()
        
        # Sample from intervals
        for _ in range(n_stars):
            # Choose interval
            interval = np.random.choice(self.intervals, p=interval_probs)
            
            # Sample within interval using inverse CDF method
            u = np.random.uniform()
            m_low = self.mass_limits[interval]
            m_high = self.mass_limits[interval+1]
            alpha = self.exponents[interval]
            
            # Fixed to use negative exponents
            if alpha != -1:
                mass = (u * (m_high**(-alpha+1) - m_low**(-alpha+1)) + 
                       m_low**(-alpha+1))**(1/(-alpha+1))
            else:
                mass = m_low * (m_high/m_low)**u
            
            masses.append(mass)
        
        return np.array(masses)