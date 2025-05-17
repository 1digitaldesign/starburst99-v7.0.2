#!/usr/bin/env python3
"""
Energy conservation checker for Starburst99 outputs.

This script analyzes the spectrum, power, and SNR outputs to verify
energy conservation within expected tolerances.
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

def read_spectrum(filename):
    """Read spectrum file and calculate total energy."""
    data = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        # Skip header line
        for line in lines[1:]:
            try:
                values = line.strip().split()
                wavelength = float(values[0])  # Angstroms
                flux = float(values[1])        # erg/s/A
                data.append((wavelength, flux))
            except (ValueError, IndexError):
                continue
    
    if not data:
        return 0.0
        
    data = np.array(data)
    
    # Calculate integrated luminosity (simple trapezoidal rule)
    total_energy = 0.0
    for i in range(len(data)-1):
        dw = data[i+1][0] - data[i][0]
        avg_flux = (data[i][1] + data[i+1][1]) / 2.0
        total_energy += avg_flux * dw
    
    return total_energy

def read_mechanical_energy(power_file):
    """Read mechanical energy from power file."""
    try:
        with open(power_file, 'r') as f:
            lines = f.readlines()
            # Skip header line
            data = []
            for line in lines[1:]:
                try:
                    values = line.strip().split()
                    time = float(values[0])       # Myr
                    power = float(values[1])      # erg/s
                    data.append((time, power))
                except (ValueError, IndexError):
                    continue
            
            if not data:
                return 0.0
            
            # Last value is the current power
            return data[-1][1]
    except FileNotFoundError:
        return 0.0

def read_snr_energy(snr_file):
    """Read supernova energy from SNR file."""
    try:
        with open(snr_file, 'r') as f:
            lines = f.readlines()
            # Skip header line
            data = []
            for line in lines[1:]:
                try:
                    values = line.strip().split()
                    time = float(values[0])       # Myr
                    rate = float(values[1])       # SN/yr
                    data.append((time, rate))
                except (ValueError, IndexError):
                    continue
            
            if not data:
                return 0.0
            
            # Last value is the current SN rate
            # Multiply by typical SN energy (10^51 erg)
            return data[-1][1] * 1.0e51
    except FileNotFoundError:
        return 0.0

def plot_energy_budget(spectrum_files, output_file):
    """Plot energy budget evolution over time."""
    times = []
    rad_energies = []
    mech_energies = []
    sn_energies = []
    
    base_prefix = os.path.commonprefix(spectrum_files)
    
    for spec_file in sorted(spectrum_files):
        # Extract time from filename (assuming pattern like "conservation.spectrum1.50")
        try:
            time_str = spec_file.split('.')[-1]
            time = float(time_str)
        except (IndexError, ValueError):
            time = 0.0
        
        power_file = spec_file.replace('spectrum', 'power')
        snr_file = spec_file.replace('spectrum', 'snr')
        
        rad_energy = read_spectrum(spec_file)
        mech_energy = read_mechanical_energy(power_file)
        sn_energy = read_snr_energy(snr_file)
        
        times.append(time)
        rad_energies.append(rad_energy)
        mech_energies.append(mech_energy)
        sn_energies.append(sn_energy)
    
    # Convert to numpy arrays
    times = np.array(times)
    rad_energies = np.array(rad_energies)
    mech_energies = np.array(mech_energies)
    sn_energies = np.array(sn_energies)
    
    # Calculate total energy
    total_energies = rad_energies + mech_energies + sn_energies
    
    # Normalize to initial energy
    if len(total_energies) > 0 and total_energies[0] > 0:
        rad_energies = rad_energies / total_energies[0]
        mech_energies = mech_energies / total_energies[0]
        sn_energies = sn_energies / total_energies[0]
        total_energies = total_energies / total_energies[0]
    
    # Create plot
    plt.figure(figsize=(12, 8))
    plt.plot(times, rad_energies, 'b-', label='Radiative')
    plt.plot(times, mech_energies, 'g-', label='Mechanical')
    plt.plot(times, sn_energies, 'r-', label='Supernova')
    plt.plot(times, total_energies, 'k--', label='Total')
    
    plt.xlabel('Time (Myr)')
    plt.ylabel('Normalized Energy')
    plt.title('Energy Budget Evolution')
    plt.legend()
    plt.grid(True)
    
    plt.savefig(output_file)
    print(f"Plot saved to {output_file}")
    
    # Calculate energy conservation
    if len(total_energies) > 1:
        max_deviation = np.max(np.abs(total_energies - 1.0))
        print(f"Maximum energy deviation: {max_deviation:.3f} ({max_deviation*100:.1f}%)")
        
        # Return True if energy is conserved within 10%
        return max_deviation < 0.1
    return False

def main():
    if len(sys.argv) < 3:
        print("Usage: check_energy.py output_plot.png spectrum_file1 [spectrum_file2 ...]")
        return 1
    
    output_file = sys.argv[1]
    spectrum_files = sys.argv[2:]
    
    if not spectrum_files:
        print("No spectrum files provided")
        return 1
    
    conserved = plot_energy_budget(spectrum_files, output_file)
    
    if conserved:
        print("✓ Energy is conserved within acceptable limits")
        return 0
    else:
        print("✗ Energy conservation exceeds acceptable limits")
        return 1

if __name__ == "__main__":
    sys.exit(main())
