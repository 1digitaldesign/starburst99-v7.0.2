# Starburst99 Python Implementation

This is a modern Python implementation of the Starburst99 stellar population synthesis code, originally written in Fortran.

## Installation

1. Install dependencies:
```bash
pip install -r requirements.txt
```

2. Install the package in development mode:
```bash
pip install -e .
```

## Usage

### Command Line Interface

Run the main program with:
```bash
python starburst_main.py input_file.txt
```

Or after installation:
```bash
starburst99 input_file.txt
```

### Python API

```python
from core.galaxy_module import GalaxyModel
from io.input_parser import InputParser

# Create galaxy model
galaxy = GalaxyModel()

# Read input parameters
parser = InputParser()
galaxy.model_params = parser.read_input('input_file.txt')

# Run calculations
galaxy.compute_sed()

# Write output
galaxy.write_output()
```

## Project Structure

```
src/python/
├── core/               # Core modules
│   ├── constants.py    # Physical constants
│   ├── galaxy_module.py # Main galaxy model
│   └── data_profiles.py # Data profile handling
├── io/                 # Input/Output modules
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

## Features

- Modern Python 3.8+ implementation
- Object-oriented design
- Comprehensive error handling
- Multiple input/output formats (standard, JSON, INI)
- Extensive logging
- Unit tests
- Type hints for better code maintainability

## Testing

Run tests with:
```bash
python -m pytest src/python/tests/
```

## Development

- Follow PEP 8 style guidelines
- Add type hints to all functions
- Write unit tests for new features
- Update documentation as needed

## Compatibility

This Python implementation aims to be compatible with the original Fortran version while providing:
- Improved code organization
- Better error handling
- More flexible I/O options
- Easier maintenance and extension

## License

Same as the original Starburst99 code.