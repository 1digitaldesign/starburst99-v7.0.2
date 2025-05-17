#!/bin/bash
#==============================================================================
# Convert all data files to JSON format
#==============================================================================
# This script converts all GALAXY/Starburst99 data files to JSON format

# Set up variables
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
TOOLS_DIR="${SCRIPT_DIR}/tools/converters"
JSON_DIR="${SCRIPT_DIR}/json_data"
EXECUTABLE="${TOOLS_DIR}/convert_to_json"

# Create JSON directory if it doesn't exist
mkdir -p "${JSON_DIR}"

# Function to handle errors
handle_error() {
    echo "Error: $1" >&2
    exit 1
}

# Check if executable exists
if [ ! -x "${EXECUTABLE}" ]; then
    echo "Compiling converter program..."
    gfortran -o convert_to_json galaxy_dat2json.f90 convert_to_json.f90 || handle_error "Failed to compile converter"
    chmod +x convert_to_json
fi

# Convert irfeatures.dat
echo "Converting irfeatures.dat..."
"${EXECUTABLE}" irfeatures "${SCRIPT_DIR}/data/auxil/irfeatures.dat" "${JSON_DIR}/irfeatures.json" || handle_error "Failed to convert irfeatures.dat"

# Convert all track files
echo "Converting track files..."
for track_file in "${SCRIPT_DIR}"/data/tracks/Z*.txt; do
    base_name=$(basename "${track_file}" .txt)
    echo "  - ${base_name}"
    "${EXECUTABLE}" tracks "${track_file}" "${JSON_DIR}/${base_name}.json" || handle_error "Failed to convert ${track_file}"
done

# Convert all mod files
echo "Converting model files..."
for mod_file in "${SCRIPT_DIR}"/data/tracks/mod*.dat; do
    base_name=$(basename "${mod_file}" .dat)
    echo "  - ${base_name}"
    "${EXECUTABLE}" modfiles "${mod_file}" "${JSON_DIR}/${base_name}.json" || handle_error "Failed to convert ${mod_file}"
done

# Convert other data files in auxil
echo "Converting auxil data files..."
for data_file in "${SCRIPT_DIR}"/data/auxil/*.dat; do
    # Skip irfeatures.dat since we've already converted it
    if [[ "${data_file}" != *"irfeatures.dat" ]]; then
        base_name=$(basename "${data_file}" .dat)
        echo "  - ${base_name}"
        "${EXECUTABLE}" modfiles "${data_file}" "${JSON_DIR}/${base_name}.json" || handle_error "Failed to convert ${data_file}"
    fi
done

echo "All conversions complete! JSON files are in ${JSON_DIR}"
echo "Successfully converted $(find "${JSON_DIR}" -name "*.json" | wc -l) files"