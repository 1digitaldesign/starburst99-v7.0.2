#!/usr/bin/env python3
"""
IMF checker for Starburst99 outputs.

This script analyzes the output to verify that the IMF is
correctly implemented.
"""

import sys
import numpy as np
import matplotlib.pyplot as plt

def theoretical_salpeter(masses, normalize=True):
    """Theoretical Salpeter IMF: dN/dM ~ M^-2.35."""
    imf = masses ** -2.35
    if normalize and len(imf) > 0:
        imf = imf / np.sum(imf)
    return imf

def read_imf_from_output(filename):
    """Extract IMF information from output file."""
    masses = []
    counts = []
    
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
            in_imf_section = False
            
            for line in lines:
                if "INITIAL MASS FUNCTION" in line:
                    in_imf_section = True
                    continue
                
                if in_imf_section and "dn/dm" in line.lower():
                    in_imf_section = False
                    continue
                
                if in_imf_section and line.strip():
                    try:
                        parts = line.strip().split()
                        if len(parts) >= 2:
                            mass = float(parts[0])
                            count = float(parts[1])
                            masses.append(mass)
                            counts.append(count)
                    except ValueError:
                        pass
    except FileNotFoundError:
        print(f"File not found: {filename}")
    
    return np.array(masses), np.array(counts)

def plot_imf_comparison(output_file, masses, actual_imf, expected_imf):
    """Plot comparison between actual and expected IMF."""
    plt.figure(figsize=(10, 6))
    
    plt.loglog(masses, actual_imf, 'bo-', label='Actual IMF')
    plt.loglog(masses, expected_imf, 'r--', label='Expected IMF')
    
    plt.xlabel('Mass (Solar Masses)')
    plt.ylabel('dN/dM (Normalized)')
    plt.title('Initial Mass Function Comparison')
    plt.legend()
    plt.grid(True)
    
    plt.savefig(output_file)
    print(f"Plot saved to {output_file}")

def compare_imfs(output_file, imf_file, plot_file):
    """Compare actual IMF with theoretical expectation."""
    masses, counts = read_imf_from_output(imf_file)
    
    if len(masses) == 0:
        print("No IMF data found in output file")
        return False
    
    # Normalize counts
    if np.sum(counts) > 0:
        counts = counts / np.sum(counts)
    
    # Calculate expected IMF
    expected = theoretical_salpeter(masses)
    
    # Plot comparison
    plot_imf_comparison(plot_file, masses, counts, expected)
    
    # Calculate error
    error = np.mean(np.abs(counts - expected) / expected)
    print(f"Mean relative error: {error:.3f} ({error*100:.1f}%)")
    
    # IMF is correct if error is less than 10%
    return error < 0.1

def main():
    if len(sys.argv) != 3:
        print("Usage: check_imf.py output_file plot_file")
        return 1
    
    output_file = sys.argv[1]
    plot_file = sys.argv[2]
    
    if compare_imfs(output_file, output_file, plot_file):
        print("✓ IMF implementation is correct within acceptable limits")
        return 0
    else:
        print("✗ IMF implementation deviates from expected values")
        return 1

if __name__ == "__main__":
    sys.exit(main())
