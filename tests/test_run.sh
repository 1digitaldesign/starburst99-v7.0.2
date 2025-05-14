#!/bin/bash
#==============================================================================
# Starburst99 Test Script
#==============================================================================
# This script runs a standard test of the Starburst99 code
#
# Usage: ./test_run.sh [input_file]
#   - Default: Uses standard.input1 in the current directory
#   - Specify a different input file as an argument
#==============================================================================

# Get directory of this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Project root directory
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
# Output directory
OUTPUT_DIR="$SCRIPT_DIR/output"
# Default input file
INPUT_FILE="$SCRIPT_DIR/standard.input1"

# Parse command-line arguments
if [ $# -ge 1 ]; then
  INPUT_FILE="$1"
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

# Set up symbolic links to data directories
echo "Setting up test environment..."
if [ ! -e tracks ]; then
  ln -sf "$ROOT_DIR/data/tracks" tracks
fi

if [ ! -e lejeune ]; then
  ln -sf "$ROOT_DIR/data/lejeune" lejeune
fi

if [ ! -e auxil ]; then
  ln -sf "$ROOT_DIR/data/auxil" auxil
fi

# Copy input file to current directory
cp "$INPUT_FILE" test.input1

# Create symbolic link for Fortran unit
if [ -e fort.1 ]; then
  rm fort.1
fi
ln -sf test.input1 fort.1

# Run the code
echo "Running Starburst99 test case..."
"$ROOT_DIR/bin/galaxy" > test.log 2>&1
RESULT=$?

if [ $RESULT -ne 0 ]; then
  echo "ERROR: Test failed with exit code $RESULT"
  echo "See test.log for details"
  exit $RESULT
fi

# Save output files with test prefix
echo "Saving test outputs..."
"$ROOT_DIR/scripts/save_output" test 1 >> test.log

# Verify outputs
SPECTRUM_FILE="test.spectrum1"
if [ -e "$SPECTRUM_FILE" ]; then
  SPECTRUM_SIZE=$(wc -l < "$SPECTRUM_FILE")
  echo "Test completed successfully."
  echo "Generated $SPECTRUM_SIZE lines of spectral data."
  echo "Output files are in: $OUTPUT_DIR"
else
  echo "ERROR: Failed to generate spectral output."
  exit 1
fi