#!/bin/bash
#==============================================================================
# Convert all data files to JSON format using Python
#==============================================================================
# This script uses the Python converter to handle both text and binary files

# Set up variables
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
TOOLS_DIR="${SCRIPT_DIR}/tools/converters"
JSON_DIR="${SCRIPT_DIR}/json_data"
PYTHON_SCRIPT="${TOOLS_DIR}/convert_data_to_json.py"
PYTHON="python3"

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

# Function to handle errors
handle_error() {
    echo -e "${RED}Error: $1${NC}" >&2
    exit 1
}

# Create JSON directory if it doesn't exist
mkdir -p "${JSON_DIR}"

# Check if Python is available
if ! command -v ${PYTHON} &> /dev/null; then
    handle_error "Python 3 not found. Please install Python 3 to use this converter."
fi

# Check if the Python script exists
if [ ! -f "${PYTHON_SCRIPT}" ]; then
    handle_error "Python converter script not found at ${PYTHON_SCRIPT}"
fi

# Make the Python script executable
chmod +x "${PYTHON_SCRIPT}"

# Function to convert a directory
convert_directory() {
    local dir="$1"
    local name="$2"
    
    echo -e "${YELLOW}Converting ${name} files from ${dir}...${NC}"
    "${PYTHON}" "${PYTHON_SCRIPT}" "${dir}" -o "${JSON_DIR}" -v || handle_error "Failed to convert ${name} files"
    echo
}

# Convert all data directories
echo -e "${GREEN}Starting conversion of all data files to JSON format...${NC}"
echo "Output directory: ${JSON_DIR}"
echo

# Convert track files
convert_directory "${SCRIPT_DIR}/data/tracks" "tracks"

# Convert auxil files
convert_directory "${SCRIPT_DIR}/data/auxil" "auxiliary"

# Convert lejeune files
convert_directory "${SCRIPT_DIR}/data/lejeune" "lejeune"

# Count the number of converted files
JSON_COUNT=$(find "${JSON_DIR}" -name "*.json" | wc -l)
echo -e "${GREEN}Conversion complete!${NC}"
echo -e "Successfully converted ${YELLOW}${JSON_COUNT}${NC} files to JSON format"
echo -e "JSON files are in: ${YELLOW}${JSON_DIR}${NC}"

# Check for error files (files that have "error" field in the JSON)
ERROR_COUNT=$(grep -l '"error":' "${JSON_DIR}"/*.json | wc -l)
if [ "${ERROR_COUNT}" -gt 0 ]; then
    echo -e "${RED}Warning: ${ERROR_COUNT} files had conversion errors.${NC}"
    echo "Check the JSON files with 'error' field for details."
fi