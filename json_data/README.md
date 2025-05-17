# Starburst99 JSON Data Files

This directory contains JSON versions of all data files used by the Starburst99 (formerly Galaxy) stellar population synthesis code. These JSON files are created automatically from the original `.dat` and `.txt` files by the `convert_to_json` utility included with the code.

## Purpose

The JSON format provides several advantages over the original fixed-format data files:

1. **Interoperability**: JSON is widely supported across programming languages and platforms
2. **Structure**: Clear hierarchical organization of data with named fields
3. **Readability**: Self-documenting format with explicit data types
4. **Web Integration**: Direct use in web applications and APIs
5. **Accessibility**: Easier to parse and use with modern data analysis tools

## File Categories

The JSON files in this directory are organized by their original data type:

### 1. Stellar Evolution Tracks

Files named like `Z0020v00.json` contain stellar evolution tracks. Each file includes:

```json
{
  "name": "Z0020v00",
  "columns": 24,
  "rows": 400,
  "tracks": [
    {
      "index": 1,
      "parameters": [1.77212388049011E+04, 119.975994, 6.226097, 4.769865, ...]
    },
    ...
  ]
}
```

### 2. Model Atmosphere Data

Files named like `modc001.json` and `mode001.json` contain model atmosphere data. Each file includes:

```json
{
  "data": [
    [parameter1, parameter2, ...],
    [parameter1, parameter2, ...],
    ...
  ]
}
```

### 3. Auxiliary Data

Files like `irfeatures.json` contain specialized data used for various calculations:

```json
{
  "wavelengths": [3000.00, 3200.00, ...],
  "extinctions": [0.00, 0.50, ...],
  "data": [
    [3.00, 2.43, 1.95, ...],
    ...
  ]
}
```

## Usage

These JSON files can be used directly with any programming language that supports JSON parsing. For example:

### Python

```python
import json

# Load track data
with open('Z0020v00.json', 'r') as f:
    track_data = json.load(f)

# Access parameters
for track in track_data['tracks']:
    age = track['parameters'][0]  # Age in years
    mass = track['parameters'][1]  # Mass in solar masses
    # Process data...
```

### JavaScript

```javascript
// Load track data
fetch('Z0020v00.json')
  .then(response => response.json())
  .then(data => {
    // Access parameters
    data.tracks.forEach(track => {
      const age = track.parameters[0];  // Age in years
      const mass = track.parameters[1]; // Mass in solar masses
      // Process data...
    });
  });
```

## Regenerating JSON Files

If the original data files are updated, you can regenerate all JSON files using one of two methods:

### 1. Fortran Converter (Text Files Only)

```bash
# From the main directory:
make convert-json
```

Or run the conversion script directly:

```bash
./tools/converters/convert_all_to_json.sh
```

### 2. Python Converter (Both Text and Binary Files)

```bash
# From the main directory:
make convert-json-py
```

Or run the Python conversion script directly:

```bash
./tools/converters/convert_all_data.sh
```

The Python converter provides better handling for binary files and more robust error handling for malformed data, but requires Python 3 with the NumPy package installed.

## File Format Details

For full details on the data structure and parameters contained in each file type, please refer to the original Starburst99 documentation and the `tools/converters/galaxy_dat2json.f90` module or `tools/converters/convert_data_to_json.py` script, which contain detailed comments on each data format.