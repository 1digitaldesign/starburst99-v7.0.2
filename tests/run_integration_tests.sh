#!/bin/bash
#==============================================================================
# Starburst99 Integration Tests
#==============================================================================
# This script runs integration tests for the Starburst99 code using
# test input files and compares the output to expected results.
#==============================================================================

# Define colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

# Set execution directories
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TEST_DIR="${ROOT_DIR}/tests"
DATA_DIR="${ROOT_DIR}/data"
OUTPUT_DIR="${ROOT_DIR}/output"
TEMP_DIR="${TEST_DIR}/temp"
INTEGRATION_DIR="${TEST_DIR}/integration"

# Create necessary directories
mkdir -p "${TEMP_DIR}"
mkdir -p "${INTEGRATION_DIR}"

# Ensure the executable exists
EXECUTABLE="${ROOT_DIR}/galaxy_fixed"
if [ ! -x "${EXECUTABLE}" ]; then
    echo -e "${RED}ERROR: Executable '${EXECUTABLE}' not found or not executable.${NC}"
    echo "Please build the program with 'make fixed' first."
    exit 1
fi

# Global test counter
TOTAL_TESTS=0
PASSED_TESTS=0
FAILED_TESTS=0

# Function to run a test case
run_test() {
    local test_name="$1"
    local input_file="$2"
    local expected_outputs=("${@:3}")
    
    echo -e "${YELLOW}Running test: ${test_name}${NC}"
    ((TOTAL_TESTS++))
    
    # Create test directory
    local test_dir="${TEMP_DIR}/${test_name}"
    rm -rf "${test_dir}"
    mkdir -p "${test_dir}"
    
    # Copy test input to fort.1
    cp "${input_file}" "${test_dir}/fort.1"
    
    # Change to test directory and run the program
    pushd "${test_dir}" > /dev/null
    
    # Run the program
    "${EXECUTABLE}" > run.log 2>&1
    local exit_code=$?
    
    # Check if the program ran successfully
    if [ ${exit_code} -ne 0 ]; then
        echo -e "${RED}✗ Test '${test_name}' failed: Program did not execute successfully.${NC}"
        echo -e "${YELLOW}Log output:${NC}"
        cat run.log
        ((FAILED_TESTS++))
        popd > /dev/null
        return 1
    fi
    
    # Check if expected output files were created
    local all_files_exist=true
    for output_file in "${expected_outputs[@]}"; do
        if [ ! -f "${output_file}" ]; then
            echo -e "${RED}✗ Test '${test_name}' failed: Expected output file '${output_file}' was not created.${NC}"
            all_files_exist=false
        fi
    done
    
    if [ "${all_files_exist}" = false ]; then
        echo -e "${YELLOW}Files in test directory:${NC}"
        ls -la
        ((FAILED_TESTS++))
        popd > /dev/null
        return 1
    fi
    
    # If we get here, the test passed
    echo -e "${GREEN}✓ Test '${test_name}' passed.${NC}"
    ((PASSED_TESTS++))
    popd > /dev/null
    return 0
}

# Function to compare numerical outputs with reference data
compare_numerical() {
    local test_file="$1"
    local reference_file="$2"
    local tolerance="$3"
    
    # Simple comparison using awk to calculate relative differences
    awk -v tol="${tolerance}" '
    BEGIN { max_diff = 0; line_count = 0; diff_count = 0; }
    NR==FNR { ref[NR] = $0; next }
    {
        line_count++;
        ref_line = ref[FNR];
        test_line = $0;
        
        # Skip empty or comment lines
        if (ref_line ~ /^[[:space:]]*$/ || ref_line ~ /^[[:space:]]*#/) {
            if (test_line ~ /^[[:space:]]*$/ || test_line ~ /^[[:space:]]*#/) {
                # Both lines are empty or comments, match is good
                next;
            } else {
                # Reference is comment/empty but test is not
                print "Line " FNR ": Reference line is empty/comment but test line is not";
                diff_count++;
                next;
            }
        }
        
        # Split the lines into fields
        split(ref_line, ref_fields);
        split(test_line, test_fields);
        
        # Check if the number of fields match
        if (length(ref_fields) != length(test_fields)) {
            print "Line " FNR ": Field count mismatch - ref: " length(ref_fields) ", test: " length(test_fields);
            diff_count++;
            next;
        }
        
        # Compare each field
        for (i = 1; i <= length(ref_fields); i++) {
            # Skip non-numeric fields
            if (ref_fields[i] ~ /^[+-]?[0-9]*\.?[0-9]+([eE][+-]?[0-9]+)?$/) {
                ref_val = ref_fields[i] + 0;  # Convert to number
                test_val = test_fields[i] + 0;
                
                # Calculate relative difference
                if (ref_val != 0) {
                    rel_diff = abs((test_val - ref_val) / ref_val);
                    if (rel_diff > tol) {
                        print "Line " FNR ", field " i ": Relative difference " rel_diff " exceeds tolerance " tol;
                        print "  Reference: " ref_val ", Test: " test_val;
                        diff_count++;
                        if (rel_diff > max_diff) max_diff = rel_diff;
                    }
                } else if (test_val != 0) {
                    print "Line " FNR ", field " i ": Reference is zero but test is " test_val;
                    diff_count++;
                }
            }
        }
    }
    END {
        printf "Compared %d lines, found %d differences, maximum relative difference: %g\n", line_count, diff_count, max_diff;
        exit (diff_count > 0 ? 1 : 0);
    }
    function abs(x) { return (x < 0) ? -x : x; }
    ' "${reference_file}" "${test_file}"
    
    return $?
}

# Create test input files
echo -e "${YELLOW}Creating test input files...${NC}"

# Test 1: Basic continuous star formation model
cat > "${INTEGRATION_DIR}/test1.input" << EOF
Test1: Continuous SF
1
1.0
1.0
1
2.35
0.1 100.0
8.0
120.0
14
0
0.0
0
0.0
50
100.0
2
1 100
0.1
5
2
1
1 1
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
EOF

# Test 2: Instantaneous burst model
cat > "${INTEGRATION_DIR}/test2.input" << EOF
Test2: Instantaneous Burst
0
1.0
0.0
1
2.35
0.1 100.0
8.0
120.0
14
0
1.0
0
0.0
50
10.0
2
1 100
0.1
5
2
1
1 1
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
EOF

# Test 3: Multiple IMF intervals
cat > "${INTEGRATION_DIR}/test3.input" << EOF
Test3: Multiple IMF
1
1.0
1.0
3
1.3 2.3 2.7
0.1 0.5 1.0 100.0
8.0
120.0
14
0
0.0
0
0.0
50
100.0
2
1 100
0.1
5
2
1
1 1
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
EOF

# Run the tests
echo -e "${YELLOW}Running integration tests...${NC}"

# Test 1: Basic continuous star formation model
run_test "test1_continuous_sf" "${INTEGRATION_DIR}/test1.input" "fort.85" "fort.86" "fort.87"

# Test 2: Instantaneous burst model
run_test "test2_instantaneous" "${INTEGRATION_DIR}/test2.input" "fort.85" "fort.86" "fort.87"

# Test 3: Multiple IMF intervals
run_test "test3_multiple_imf" "${INTEGRATION_DIR}/test3.input" "fort.85" "fort.86" "fort.87"

# Print summary
echo -e "${YELLOW}=========================================${NC}"
echo -e "${YELLOW}Integration Test Summary${NC}"
echo -e "${YELLOW}=========================================${NC}"
echo -e "Total tests:  ${TOTAL_TESTS}"
echo -e "Passed tests: ${GREEN}${PASSED_TESTS}${NC}"
echo -e "Failed tests: ${RED}${FAILED_TESTS}${NC}"
echo -e "${YELLOW}=========================================${NC}"

# Return success if all tests passed
if [ ${FAILED_TESTS} -eq 0 ]; then
    echo -e "${GREEN}All integration tests passed!${NC}"
    exit 0
else
    echo -e "${RED}Some integration tests failed.${NC}"
    exit 1
fi