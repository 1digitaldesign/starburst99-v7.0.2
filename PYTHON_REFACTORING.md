# Python Refactoring of Starburst99

This document summarizes the Python refactoring of the Starburst99 stellar population synthesis code.

## Overview

The original Fortran code has been refactored into a modern Python implementation using:
- Python 3.8+ features
- Object-oriented design
- Modern module architecture
- Type hints and documentation

## Directory Structure

```
src/python/
├── core/               # Core modules
│   ├── constants.py    # Physical constants
│   ├── galaxy_module.py # Main galaxy model
│   └── data_profiles.py # Data profile handling
├── file_io/            # Input/Output modules (renamed from 'io' to avoid conflicts)
│   ├── input_parser.py # Input file parsing
│   └── output_writer.py # Output file writing
├── models/             # Physical models
│   ├── imf.py         # Initial Mass Function
│   └── stellar_tracks.py # Stellar evolution tracks
├── utils/              # Utility functions
│   └── utilities.py    # Helper functions
├── tests/              # Unit tests
│   └── test_basic.py  # Basic tests
└── starburst_main.py  # Main program
```

## Key Improvements

1. **Modular Architecture**: Code is organized into logical modules instead of monolithic files
2. **Better Error Handling**: Comprehensive exception handling and logging
3. **Multiple Input Formats**: Support for standard, JSON, and INI input formats
4. **Type Safety**: Type hints throughout the codebase
5. **Unit Tests**: Test suite for validation
6. **Modern Python**: Uses dataclasses, pathlib, and other modern features

## Components Converted

### Core Components
- ✅ data_profiles.f90 → data_profiles.py
- ✅ galaxy_module.f90 → galaxy_module.py
- ✅ starburst_main.f90 → starburst_main.py

### Supporting Modules
- ✅ Physical constants module
- ✅ Initial Mass Function (IMF) calculations
- ✅ Stellar evolution track handling
- ✅ Input/Output parsers and writers
- ✅ Utility functions

## Testing

All tests pass successfully:
```
tests/test_basic.py::TestConstants::test_constants_values PASSED
tests/test_basic.py::TestUtilities::test_exp10 PASSED
tests/test_basic.py::TestUtilities::test_linear_interp PASSED
tests/test_basic.py::TestIMF::test_imf_creation PASSED
tests/test_basic.py::TestIMF::test_imf_integration PASSED
tests/test_basic.py::TestIMF::test_imf_value PASSED
tests/test_basic.py::TestGalaxyModel::test_galaxy_creation PASSED
tests/test_basic.py::TestGalaxyModel::test_model_parameters PASSED
```

## Usage

### Command Line
```bash
python starburst_main.py input_file.txt
```

### Python API
```python
from core.galaxy_module import GalaxyModel
from file_io.input_parser import InputParser

galaxy = GalaxyModel()
parser = InputParser()
galaxy.model_params = parser.read_input('input_file.txt')
```

## Installation

```bash
cd src/python
pip install -r requirements.txt
python -m pytest tests/  # Run tests
```

## Next Steps

1. Complete implementation of calculation routines
2. Add more comprehensive tests
3. Validate against Fortran output
4. Add documentation and examples
5. Optimize performance-critical sections

## Notes

- The `io` module was renamed to `file_io` to avoid conflicts with Python's built-in `io` module
- Import statements use absolute imports for better compatibility
- The implementation follows modern Python best practices