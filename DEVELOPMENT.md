# Development Guide

This file provides guidance for working with code in this repository.

## Galaxy - Starburst99 Stellar Population Synthesis Code

Galaxy (also known as Starburst99) is a Fortran code for modeling the spectrophotometric and related properties of star-forming galaxies and stellar populations. The code was developed by Claus Leitherer and colleagues, with the first version released in 1998.

## Key Commands

### Building the Code

```bash
# Compile the code with optimizations
make

# If debugging is needed, edit the Makefile to use these flags:
# FFLAGS = -c -g -C -sb -Nl350
# LFLAGS = -g
```

### Running the Code

```bash
# Navigate to the output directory
cd output

# Run the code (may need to update paths in go_galaxy first)
./go_galaxy

# Save output files with custom name and extension
../save_output custom_model 1
```

## File Structure

- `galaxy.f90`: Main source code in Fortran 90
- `Makefile`: Compilation configuration
- `output/go_galaxy`: Shell script to run the code
- `save_output`: Script to save and rename output files
- `output/standard.input1`: Default input parameters
- Auxiliary directories:
  - `tracks/`: Stellar evolutionary tracks
  - `lejeune/`: Atmosphere models
  - `auxil/`: Spectral libraries and calibration data

## Configuration

### Input Parameters

The model parameters are configured in the input file (default: `output/standard.input1`). Key parameters include:

- Star formation mode (instantaneous or continuous)
- Initial mass function (IMF) specification
- Metallicity and stellar evolution tracks
- Time steps and age range
- Output file selection

To create a new model:
1. Copy and modify the `standard.input1` file
2. Update the `ninput` variable in `go_galaxy` script
3. Run the code
4. Save outputs with a unique name

### Modifying Path Configuration

Before running the code, update the directory paths in the `go_galaxy` script:
- `drun`: Directory where output files will be written
- `dcode`: Directory where the executable is located
- `dlib`: Directory containing auxiliary data (tracks, atmosphere models, etc.)

## Code Architecture

The Galaxy/Starburst99 code implements evolutionary synthesis modeling of stellar populations, calculating:

- Spectral energy distributions from UV to near-IR
- Stellar wind and supernova feedback
- Chemical yields
- Photometric properties
- High-resolution spectral features

The code uses stellar evolutionary tracks from various sources (Geneva, Padova) and combines them with atmospheric models (Lejeune, Schmutz, Hillier, Pauldrach) to synthesize spectral properties of composite stellar populations.

Output files include numerous physical properties with naming conventions like:
- `*.spectrum*`: Spectral energy distributions
- `*.color*`: Photometric properties
- `*.quanta*`: Ionizing photon production
- `*.power*`: Mechanical luminosity
- `*.yield*`: Chemical yields