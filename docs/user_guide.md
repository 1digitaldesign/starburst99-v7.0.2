# Starburst99 Population Synthesis Code

## Overview

Starburst99 (formerly "Galaxy") is a stellar population synthesis code for modeling the spectrophotometric and related properties of star-forming galaxies. This code calculates observable parameters of stellar populations including:

- Spectral energy distributions from UV to near-IR
- Stellar wind feedback and supernova rates
- Chemical yields from stellar evolution
- Line spectra at various wavelengths and resolutions
- Photometric properties (colors, magnitudes)
- Ionizing photon production

This version has been modernized to use Fortran 2018 standards while maintaining full compatibility with the original functionality.

**Citation:** Starburst99 should be cited as Leitherer et al. (1999, ApJS, 123, 3), Vazquez & Leitherer (2005; ApJ, 621, 695), Leitherer et al. (2010; ApJS, 189, 309), and Leitherer et al. (2014, ApJS, 213, 1).

## Version History

- **First version:** August 12, 1998 (Claus Leitherer)
- **Version 4.0:** July 2002 - Added blanketed WR models
- **Version 5.0:** December 2004 - Added Padova tracks and high-resolution optical library
- **Version 6.0:** August 25, 2010 - Added theoretical UV spectra
- **Version 7.0.0:** March 2014 - Added rotating tracks and Wolf-Rayet library
- **Version 7.1.0:** May 2024 - Complete Fortran 2018 modernization

## Directory Structure

The package follows a standard code repository structure:

```
starburst99/
├── src/               # Source code files
│   ├── galaxy_module.f90           # Core module with data structures
│   ├── galaxy_module_error.f90     # Error handling submodule
│   └── starburst_main.f90          # Main program code
├── bin/               # Executable files
├── data/              # Data files
│   ├── tracks/        # Stellar evolutionary tracks
│   ├── lejeune/       # Model atmospheres
│   └── auxil/         # Auxiliary data files
├── json_data/         # JSON versions of data files
├── tools/             # Utility tools
│   └── converters/    # Data conversion utilities
├── scripts/           # Runtime scripts
│   ├── go_galaxy      # Main execution script
│   └── save_output    # Output management script
├── output/            # Runtime output directory
├── docs/              # Documentation
└── tests/             # Test files
```

## Building the Code

### Prerequisites
- A Fortran 2018 compatible compiler (e.g., gfortran 8.0+)
- Make build utility

### Compilation
```bash
# Standard build (optimized)
make

# Debug build with additional validation
make debug

# Clean build artifacts
make clean

# Setup runtime environment
make setup

# Build JSON converter only
make json

# Convert all data files to JSON (Fortran version)
make convert-json

# Convert all data files to JSON (Python version - handles binary files)
make convert-json-py
```

### Optional Installation
```bash
# Install executable to ~/bin
make install
```

## Running the Code

1. Ensure directory structure is set up correctly with `make setup`
2. Modify input parameters in `output/standard.input1` as needed
3. Run the code:
```bash
# Using the Makefile
make run

# Or directly:
cd output
./go_galaxy [input_file] [output_prefix] [ext_number]
```

Example:
```bash
cd output
./go_galaxy standard.input1 mymodel 1
```

This will generate output files with names like `mymodel.spectrum1`, `mymodel.color1`, etc.

## Input Parameters

The file `standard.input1` contains all model parameters. Key parameters include:

- **Model Designation**: Identifier for the model
- **Star Formation Mode**: Instantaneous (-1) or continuous (>0)
- **Stellar Mass Parameters**: Total mass or SFR
- **IMF Settings**: Exponents and mass boundaries
- **Metallicity Selection**: Choose from Geneva or Padova tracks at various metallicities
- **Time Range**: Initial time, time step, and maximum age
- **Model Atmosphere**: Choice of atmospheric models
- **Output Selection**: Flags for different output products

## Modernization Features

The code has been modernized with:

1. **Modular Structure**
   - Properly encapsulated modules with explicit interfaces
   - Derived types with type-bound procedures
   - Submodules for implementation details

2. **Fortran 2018 Features**
   - Enhanced error handling
   - Allocatable character strings
   - Block constructs
   - Pure and elemental procedures
   - Associate constructs for improved readability
   - Abstract interfaces

3. **Improved Build System**
   - Modern Makefile with proper dependencies
   - Multiple build configurations (debug, release)
   - Better error reporting during compilation

4. **Enhanced Runtime Experience**
   - Better error messages
   - Progress reporting
   - Improved script usability
   
5. **Data Format Modernization**
   - JSON conversion utilities for all data files
   - Machine-readable format compatible with modern tools
   - Preserves all scientific data while adding metadata

## Documentation

Additional documentation files:
- `docs/README.cleanup`: Details of modernization changes
- `docs/README`: Original documentation from previous versions

## Support

While there is no official help desk, the authors may assist with questions and issues on a time-available basis. The code is distributed freely, and users accept sole responsibility for results produced by the code.

---

*Original authors: Claus Leitherer, Carmelle Robert, Daniel Schaerer, Jeff Goldader, Rosa Gonzalez-Delgado, & Duilia de Mello*

*Disclaimer: This code is distributed freely to the community. The user accepts sole responsibility for the results produced by the code. Although every effort has been made to identify and eliminate errors, we accept no responsibility for erroneous model predictions.*
