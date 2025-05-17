# Starburst99 Data Conversion Tools

This directory contains tools for converting Starburst99/Galaxy data files to JSON format.

## Available Converters

### Fortran Converter (Text Files)

The Fortran converter is optimized for text-based data files and follows the same code style as the main application.

**Files:**
- `galaxy_dat2json.f90` - Module implementing conversion functions
- `convert_to_json.f90` - Main program with command-line interface
- `convert_all_to_json.sh` - Shell script to convert all data files

**Usage:**
```bash
# Run from main directory
make convert-json

# Or use directly
cd tools/converters
./convert_to_json irfeatures ../../auxil/irfeatures.dat output.json
./convert_to_json tracks ../../tracks/Z0020v00.txt output.json
./convert_to_json modfiles ../../tracks/modc001.dat output.json
```

### Python Converter (Binary & Text Files)

The Python converter offers more robust parsing capabilities, including handling binary files and malformed data.

**Files:**
- `convert_data_to_json.py` - Python script with advanced parsing capabilities
- `convert_all_data.sh` - Shell script to run the Python converter on all files

**Usage:**
```bash
# Run from main directory
make convert-json-py

# Or use directly
cd tools/converters
./convert_data_to_json.py ../../tracks/modc001.dat -o output.json
./convert_data_to_json.py ../../tracks/ -o ../../json_data
```

## Choosing a Converter

- **Fortran Converter**: Use for standard text-based data files when maintaining the same codebase style is important.
- **Python Converter**: Use for handling binary files or when dealing with problematic data formats.

## Creating New Converters

To add support for a new data format:

1. **For Fortran**: Add a new subroutine to `galaxy_dat2json.f90`
2. **For Python**: Add a new handler function to `convert_data_to_json.py`

Then update the file type detection logic in the conversion scripts to recognize the new format.

## Output Format

All converters generate JSON in the `json_data/` directory at the root of the project. The JSON format includes:

- Original data values
- Metadata about the source file
- Conversion information
- Error details (if conversion failed)

See `json_data/README.md` for more details on the output format and structure.