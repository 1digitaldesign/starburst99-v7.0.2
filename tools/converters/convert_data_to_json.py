#!/usr/bin/env python3
"""
Starburst99/Galaxy Data to JSON Converter

This script converts various data file formats used by the Starburst99/Galaxy
code to JSON format. It handles binary and text files with robust error handling.
"""

import os
import sys
import json
import argparse
import struct
import numpy as np
from pathlib import Path
import logging
import re
from typing import List, Dict, Any, Union, Tuple, Optional

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger('data2json')

# Constants
LINE_LENGTH = 1024
MAX_HEADER_LINES = 10
OUTPUT_DIR = "json_data"


def is_binary_file(file_path: str) -> bool:
    """Check if a file is binary by reading the first chunk of data."""
    try:
        with open(file_path, 'r', encoding='utf-8') as file:
            file.read(1024)
        return False
    except UnicodeDecodeError:
        return True


def detect_file_type(file_path: str) -> str:
    """Detect the type of the data file based on its name and content."""
    file_name = os.path.basename(file_path)
    
    # Check file extension and name patterns
    if "irfeatures.dat" in file_path:
        return "irfeatures"
    elif file_name.startswith("Z") and file_name.endswith(".txt"):
        return "tracks"
    elif file_name.startswith("mod") and file_name.endswith(".dat"):
        return "modfiles"
    elif file_name.endswith(".dat"):
        # General .dat file - need to check content
        if is_binary_file(file_path):
            return "binary"
        else:
            return "generic_text"
    elif file_name.endswith(".txt"):
        return "generic_text"
    else:
        return "unknown"


def read_text_file(file_path: str) -> List[str]:
    """Read a text file and return its lines."""
    try:
        with open(file_path, 'r', encoding='utf-8') as file:
            return [line.strip() for line in file.readlines()]
    except UnicodeDecodeError:
        logger.warning(f"Error reading {file_path} as text. It may be binary.")
        return []


def read_binary_file(file_path: str) -> bytes:
    """Read a binary file and return its content."""
    with open(file_path, 'rb') as file:
        return file.read()


def convert_irfeatures(file_path: str) -> Dict[str, Any]:
    """Convert irfeatures.dat file to JSON format."""
    lines = read_text_file(file_path)
    if not lines:
        raise ValueError(f"Failed to read {file_path}")
    
    # Parse wavelengths (first line)
    try:
        wavelengths = [float(x) for x in lines[0].split()]
    except (ValueError, IndexError):
        logger.warning(f"Error parsing wavelengths in {file_path}")
        wavelengths = []
    
    # Parse extinctions (second line)
    try:
        extinctions = [float(x) for x in lines[1].split()]
    except (ValueError, IndexError):
        logger.warning(f"Error parsing extinctions in {file_path}")
        extinctions = []
    
    # Parse data matrix
    data_matrix = []
    for i in range(2, len(lines)):
        try:
            row = [float(x) for x in lines[i].split()]
            if row:  # Skip empty lines
                data_matrix.append(row)
        except ValueError:
            logger.warning(f"Error parsing data line {i+1} in {file_path}")
    
    return {
        "wavelengths": wavelengths,
        "extinctions": extinctions,
        "data": data_matrix
    }


def convert_tracks(file_path: str) -> Dict[str, Any]:
    """Convert track files to JSON format."""
    lines = read_text_file(file_path)
    if not lines:
        raise ValueError(f"Failed to read {file_path}")
    
    # Extract track name from first line
    track_name = lines[0].strip()
    
    # Extract dimensions from second line
    dimensions = lines[1].strip().split()
    try:
        num_cols = int(dimensions[0]) if len(dimensions) > 0 else 0
        num_rows = int(dimensions[1]) if len(dimensions) > 1 else 0
    except (ValueError, IndexError):
        logger.warning(f"Error parsing dimensions in {file_path}")
        num_cols, num_rows = 0, 0
    
    # Skip header lines (usually 3 more lines)
    data_start = 2
    for i in range(2, min(MAX_HEADER_LINES, len(lines))):
        if re.match(r'^\s*\d+\s+', lines[i]):
            data_start = i
            break
        # Otherwise it's a header line
    
    # Parse track data
    tracks = []
    for i in range(data_start, len(lines)):
        line = lines[i].strip()
        if not line:
            continue
            
        # Try to parse line with track data
        try:
            values = line.split()
            if len(values) > 1:
                index = int(values[0])
                parameters = [float(x) for x in values[1:]]
                tracks.append({
                    "index": index,
                    "parameters": parameters
                })
        except (ValueError, IndexError):
            logger.warning(f"Error parsing track data on line {i+1} in {file_path}")
    
    return {
        "name": track_name,
        "columns": num_cols,
        "rows": num_rows,
        "tracks": tracks
    }


def convert_generic_text(file_path: str) -> Dict[str, Any]:
    """Convert a generic text data file to JSON format."""
    lines = read_text_file(file_path)
    if not lines:
        raise ValueError(f"Failed to read {file_path}")
    
    # Try to determine if it's a header+data format or just data
    header = []
    data = []
    data_start = 0
    
    # Scan for potential headers
    for i, line in enumerate(lines[:min(MAX_HEADER_LINES, len(lines))]):
        # Check if line looks like data (has numbers)
        if re.match(r'^\s*[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?\s', line):
            data_start = i
            break
        else:
            header.append(line)
    
    # Parse data
    for i in range(data_start, len(lines)):
        line = lines[i].strip()
        if not line:
            continue
            
        # Try to convert all values to floats
        try:
            row = [float(x) for x in line.split()]
            data.append(row)
        except ValueError:
            # If conversion fails, treat as string data
            data.append(line)
    
    result = {"data": data}
    
    # Add header if present
    if header:
        result["header"] = header
    
    return result


def convert_binary_data(file_path: str) -> Dict[str, Any]:
    """Attempt to convert binary data file to JSON."""
    # Read the binary data
    binary_data = read_binary_file(file_path)
    if not binary_data:
        raise ValueError(f"Failed to read binary file {file_path}")
    
    # Try to interpret as Fortran binary format with record markers
    # Common Fortran binary format has 4 or 8-byte record markers
    try:
        # Try to determine if it's floating-point data with either single or double precision
        # This is experimental and may need adjustment for specific file formats
        
        # First attempt: Try as Fortran unformatted file with 4-byte record markers
        data = []
        pos = 0
        
        while pos < len(binary_data):
            # Try to read record length (4 bytes)
            if pos + 4 > len(binary_data):
                break
                
            record_length = struct.unpack('i', binary_data[pos:pos+4])[0]
            pos += 4
            
            # Sanity check for record length
            if record_length <= 0 or record_length > 10000:
                raise ValueError("Invalid record length")
                
            # Try to read as many 4-byte floats as will fit in the record
            num_floats = record_length // 4
            float_data = []
            
            for _ in range(num_floats):
                if pos + 4 > len(binary_data):
                    break
                value = struct.unpack('f', binary_data[pos:pos+4])[0]
                float_data.append(value)
                pos += 4
                
            # Skip the ending record marker
            pos += 4
            
            if float_data:
                data.append(float_data)
                
        if not data:
            raise ValueError("No valid data found")
            
        return {
            "format": "binary",
            "data": data,
            "_converted": "Experimental binary conversion - verify data integrity"
        }
            
    except Exception as e:
        # If structured interpretation fails, include a hex dump for troubleshooting
        logger.warning(f"Binary interpretation failed: {str(e)}")
        
        # Include a sample of the binary data as hex
        hex_sample = binary_data[:min(100, len(binary_data))].hex()
        
        return {
            "format": "binary",
            "conversion_error": "Unable to interpret binary format",
            "file_size_bytes": len(binary_data),
            "hex_sample": hex_sample,
            "_note": "Binary file could not be automatically converted. Manual conversion required."
        }


def convert_data_file(file_path: str, output_path: str) -> bool:
    """Convert a data file to JSON based on its type."""
    try:
        file_type = detect_file_type(file_path)
        logger.info(f"Converting {file_path} (detected type: {file_type})")
        
        if file_type == "irfeatures":
            data = convert_irfeatures(file_path)
        elif file_type == "tracks":
            data = convert_tracks(file_path)
        elif file_type == "modfiles":
            if is_binary_file(file_path):
                data = convert_binary_data(file_path)
            else:
                data = convert_generic_text(file_path)
        elif file_type == "binary":
            data = convert_binary_data(file_path)
        elif file_type == "generic_text":
            data = convert_generic_text(file_path)
        else:
            logger.warning(f"Unknown file type for {file_path}")
            data = {
                "error": "Unknown file type",
                "file_path": file_path
            }
        
        # Add metadata to all conversions
        data["_metadata"] = {
            "source_file": os.path.basename(file_path),
            "detected_type": file_type,
            "converted_by": "convert_data_to_json.py"
        }
        
        # Write JSON output
        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2)
            
        logger.info(f"Successfully converted to {output_path}")
        return True
        
    except Exception as e:
        logger.error(f"Error converting {file_path}: {str(e)}")
        
        # Create a JSON file with error information
        with open(output_path, 'w', encoding='utf-8') as f:
            error_data = {
                "error": str(e),
                "file_path": file_path,
                "_metadata": {
                    "source_file": os.path.basename(file_path),
                    "conversion_failed": True,
                    "error_type": type(e).__name__
                }
            }
            json.dump(error_data, f, indent=2)
            
        return False


def convert_directory(input_dir: str, output_dir: str, file_patterns: List[str] = None) -> Tuple[int, int]:
    """Convert all data files in a directory that match the patterns."""
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Default patterns if none provided
    if not file_patterns:
        file_patterns = ["*.dat", "*.txt"]
        
    # Find all files that match the patterns
    all_files = []
    for pattern in file_patterns:
        all_files.extend(list(Path(input_dir).glob(pattern)))
    
    success_count = 0
    failure_count = 0
    
    for file_path in all_files:
        input_path = str(file_path)
        base_name = os.path.basename(input_path)
        output_path = os.path.join(output_dir, f"{os.path.splitext(base_name)[0]}.json")
        
        if convert_data_file(input_path, output_path):
            success_count += 1
        else:
            failure_count += 1
    
    return success_count, failure_count


def main():
    parser = argparse.ArgumentParser(
        description='Convert Starburst99/Galaxy data files to JSON format'
    )
    
    # Define command-line arguments
    parser.add_argument(
        'path',
        help='File or directory to convert'
    )
    parser.add_argument(
        '-o', '--output',
        help='Output file (for single file) or directory (for directory conversion)',
        default=OUTPUT_DIR
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose logging'
    )
    parser.add_argument(
        '-p', '--pattern',
        action='append',
        help='File pattern to match (can be specified multiple times, e.g., -p "*.dat" -p "*.txt")'
    )
    
    args = parser.parse_args()
    
    # Configure logging level
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    # Handle path
    input_path = args.path
    
    if os.path.isfile(input_path):
        # Single file conversion
        output_path = args.output
        if os.path.isdir(output_path):
            # If output is a directory, place the file there
            base_name = os.path.basename(input_path)
            output_path = os.path.join(args.output, f"{os.path.splitext(base_name)[0]}.json")
        
        # Ensure output directory exists
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        
        # Convert the file
        success = convert_data_file(input_path, output_path)
        sys.exit(0 if success else 1)
        
    elif os.path.isdir(input_path):
        # Directory conversion
        output_dir = args.output
        success_count, failure_count = convert_directory(
            input_path, 
            output_dir, 
            args.pattern
        )
        
        logger.info(f"Conversion complete. Success: {success_count}, Failures: {failure_count}")
        
        if failure_count > 0:
            sys.exit(1)
        else:
            sys.exit(0)
    else:
        logger.error(f"Input path does not exist: {input_path}")
        sys.exit(1)


if __name__ == "__main__":
    main()