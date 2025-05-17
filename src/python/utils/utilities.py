"""Utility functions for astronomical calculations"""

import numpy as np
from typing import Union, List, Tuple


def exp10(x: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    Compute 10^x.
    
    Args:
        x: Exponent value(s)
        
    Returns:
        10 raised to the power of x
    """
    return np.power(10.0, x)


def linear_interp(x: float, x_arr: np.ndarray, y_arr: np.ndarray) -> float:
    """
    Perform linear interpolation.
    
    Args:
        x: Value at which to interpolate
        x_arr: Array of x values (must be sorted)
        y_arr: Array of y values
        
    Returns:
        Interpolated y value at x
    """
    # Convert to numpy arrays if needed
    x_arr = np.asarray(x_arr)
    y_arr = np.asarray(y_arr)
    
    if x <= x_arr[0]:
        return y_arr[0]
    elif x >= x_arr[-1]:
        return y_arr[-1]
    
    # Find bracketing indices
    idx = np.searchsorted(x_arr, x) - 1
    
    # Linear interpolation
    x1, x2 = x_arr[idx], x_arr[idx + 1]
    y1, y2 = y_arr[idx], y_arr[idx + 1]
    
    return y1 + (y2 - y1) * (x - x1) / (x2 - x1)


def integer_to_string(n: int, width: int = 0) -> str:
    """
    Convert integer to string with optional zero padding.
    
    Args:
        n: Integer to convert
        width: Minimum width (0 for no padding)
        
    Returns:
        Formatted string
    """
    if width > 0:
        return f"{n:0{width}d}"
    else:
        return str(n)