#!/bin/bash
#==============================================================================
# Starburst99 Numerical Testing
#==============================================================================
# This script runs tests specifically to verify the numerical behavior
# of the Starburst99 code, checking for conservation laws and
# comparing outputs to analytical solutions where possible.
#
# Usage: ./numerical_tests.sh
#==============================================================================

# Set up environment
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
TEST_OUTPUT_DIR="$SCRIPT_DIR/output/numerical"
TOOLS_DIR="$SCRIPT_DIR/tools"

# Set terminal colors
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Create necessary directories
mkdir -p "$TEST_OUTPUT_DIR" "$TOOLS_DIR"

# Function to check if galaxy executable exists
check_executable() {
  if [ ! -e "$ROOT_DIR/galaxy" ] && [ ! -e "$ROOT_DIR/bin/galaxy" ]; then
    echo -e "${RED}ERROR: Galaxy executable not found.${NC}"
    echo "Build the code first with: cd $ROOT_DIR && make"
    exit 1
  fi
}

# Create input file for conservation test
create_conservation_input() {
  local output_file="$TEST_OUTPUT_DIR/conservation.input"
  
  cat > "$output_file" << EOF
conservation     Model designation
-1     Star formation mode: >0 continuous, -1 instantaneous
1.0     Total stellar mass (10^6 Msun)
1.0     Star formation rate (Msun/yr)
1     Number of intervals for the IMF
2.35     Exponent for the IMF
0.1 100.     Mass boundaries for the IMF (Msun)
8.0     SNII cutoff mass (Msun)
120.0     Black hole cutoff mass (Msun)
4     Metallicity + tracks identifier (solar)
0     Wind model
0.0    Initial time (10^6 yr)
0     Time scale: 0=linear, 1=logarithmic
0.5     Time step (10^6 yr) if jtime=0
10     Number of steps if jtime=1
50.0     Last grid point (10^6 yr)
0     Atmosphere model
0     Metalicity of high-resolution models
1     UV line spectrum (1=solar, 2=LMC/SMC)
0 0     RSG parameters for individual features
1 1 1 1 1     Output files 1-5
1 1 1 1 1     Output files 6-10
1 1 1 1 1     Output files 11-15
EOF
  
  echo "$output_file"
}

# Create Python utility to check energy conservation
create_energy_checker() {
  local output_file="$TOOLS_DIR/check_energy.py"
  
  cat > "$output_file" << EOF
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
EOF

  chmod +x "$output_file"
  echo "$output_file"
}

# Function to create a simple Python mass function checker
create_imf_checker() {
  local output_file="$TOOLS_DIR/check_imf.py"
  
  cat > "$output_file" << EOF
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
EOF

  chmod +x "$output_file"
  echo "$output_file"
}

# Run the conservation test
run_conservation_test() {
  echo -e "${BLUE}Running energy conservation test...${NC}"
  
  # Prepare test environment
  cd "$TEST_OUTPUT_DIR"
  mkdir -p plots
  
  # Create symbolic links to data directories if they don't exist
  if [ ! -e tracks ]; then ln -sf "$ROOT_DIR/data/tracks" tracks; fi
  if [ ! -e lejeune ]; then ln -sf "$ROOT_DIR/data/lejeune" lejeune; fi
  if [ ! -e auxil ]; then ln -sf "$ROOT_DIR/data/auxil" auxil; fi
  
  # Create input file
  local input_file=$(create_conservation_input)
  cp "$input_file" fort.1
  
  # Run galaxy code
  echo "  Executing galaxy for conservation test..."
  "$ROOT_DIR/galaxy" > conservation.log 2>&1
  local result=$?
  
  if [ $result -ne 0 ]; then
    echo -e "${RED}  ✗ Test failed: Galaxy execution error${NC}"
    echo "  See conservation.log for details"
    return 1
  fi
  
  # Save output files
  "$ROOT_DIR/scripts/save_output" "conservation" 1 >> conservation.log
  
  # Run multiple timesteps to check conservation over time
  local time_steps=(0 5 10 15 20 25 30 35 40 45 50)
  local spectrum_files=()
  
  for time in "${time_steps[@]}"; do
    # Update input file with new initial time
    sed -i.bak "s/^0.0    Initial time (10^6 yr)/$time    Initial time (10^6 yr)/" fort.1
    
    # Run galaxy code for this timestep
    echo "  Running for time: $time Myr"
    "$ROOT_DIR/galaxy" > "conservation.log.$time" 2>&1
    
    # Save output with timestep suffix
    "$ROOT_DIR/scripts/save_output" "conservation" "1.$time" >> "conservation.log.$time"
    
    # Add spectrum file to list for analysis
    spectrum_files+=("conservation.spectrum1.$time")
  done
  
  # Check energy conservation using Python tool
  local energy_checker=$(create_energy_checker)
  
  if command -v python3 >/dev/null 2>&1; then
    echo "  Analyzing energy conservation..."
    python3 "$energy_checker" "plots/energy_conservation.png" "${spectrum_files[@]}"
    local result=$?
    
    if [ $result -eq 0 ]; then
      echo -e "${GREEN}  ✓ Energy conservation test passed${NC}"
    else
      echo -e "${RED}  ✗ Energy conservation test failed${NC}"
    fi
    
    return $result
  else
    echo -e "${YELLOW}  ⚠ Cannot check energy conservation: Python 3 not available${NC}"
    echo "  Install Python 3 with numpy and matplotlib to run this test"
    return 0  # Don't fail the test if Python is not available
  fi
}

# Run the IMF test
run_imf_test() {
  echo -e "\n${BLUE}Running IMF implementation test...${NC}"
  
  # Prepare test environment
  cd "$TEST_OUTPUT_DIR"
  
  # Create symbolic links to data directories if they don't exist
  if [ ! -e tracks ]; then ln -sf "$ROOT_DIR/data/tracks" tracks; fi
  if [ ! -e lejeune ]; then ln -sf "$ROOT_DIR/data/lejeune" lejeune; fi
  if [ ! -e auxil ]; then ln -sf "$ROOT_DIR/data/auxil" auxil; fi
  
  # Create input file - using a simple Salpeter IMF for clear comparison
  cat > imf_test.input << EOF
imf_test     Model designation
-1     Star formation mode: >0 continuous, -1 instantaneous
1.0     Total stellar mass (10^6 Msun)
1.0     Star formation rate (Msun/yr)
1     Number of intervals for the IMF
2.35     Exponent for the IMF
1.0 100.     Mass boundaries for the IMF (Msun)
8.0     SNII cutoff mass (Msun)
120.0     Black hole cutoff mass (Msun)
4     Metallicity + tracks identifier (solar)
0     Wind model
0.0    Initial time (10^6 yr)
0     Time scale: 0=linear, 1=logarithmic
1.0     Time step (10^6 yr) if jtime=0
10     Number of steps if jtime=1
10.0     Last grid point (10^6 yr)
0     Atmosphere model
0     Metalicity of high-resolution models
1     UV line spectrum (1=solar, 2=LMC/SMC)
0 0     RSG parameters for individual features
1 1 1 1 1     Output files 1-5
1 1 1 1 1     Output files 6-10
1 1 1 1 1     Output files 11-15
EOF
  
  # Link fort.1
  cp imf_test.input fort.1
  
  # Run galaxy code
  echo "  Executing galaxy for IMF test..."
  "$ROOT_DIR/galaxy" > imf_test.log 2>&1
  local result=$?
  
  if [ $result -ne 0 ]; then
    echo -e "${RED}  ✗ Test failed: Galaxy execution error${NC}"
    echo "  See imf_test.log for details"
    return 1
  fi
  
  # Save output files
  "$ROOT_DIR/scripts/save_output" "imf_test" 1 >> imf_test.log
  
  # Check IMF using Python tool
  local imf_checker=$(create_imf_checker)
  
  if command -v python3 >/dev/null 2>&1; then
    echo "  Analyzing IMF implementation..."
    python3 "$imf_checker" "imf_test.output1" "plots/imf_comparison.png"
    local result=$?
    
    if [ $result -eq 0 ]; then
      echo -e "${GREEN}  ✓ IMF implementation test passed${NC}"
    else
      echo -e "${RED}  ✗ IMF implementation test failed${NC}"
    fi
    
    return $result
  else
    echo -e "${YELLOW}  ⚠ Cannot check IMF implementation: Python 3 not available${NC}"
    echo "  Install Python 3 with numpy and matplotlib to run this test"
    return 0  # Don't fail the test if Python is not available
  fi
}

# Main function to run numerical tests
run_numerical_tests() {
  local total_failures=0
  
  # Check if executable exists
  check_executable
  
  # Run energy conservation test
  run_conservation_test
  total_failures=$((total_failures + $?))
  
  # Run IMF implementation test
  run_imf_test
  total_failures=$((total_failures + $?))
  
  # Summary
  echo -e "\n${BLUE}===== Numerical Test Summary =====${NC}"
  if [ $total_failures -eq 0 ]; then
    echo -e "${GREEN}All numerical tests passed!${NC}"
    return 0
  else
    echo -e "${RED}$total_failures numerical test(s) failed.${NC}"
    return 1
  fi
}

# Execute the numerical tests
run_numerical_tests
exit $?