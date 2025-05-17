#!/bin/bash
#==============================================================================
# Starburst99 Test Suite
#==============================================================================
# This script runs a comprehensive suite of tests for the Starburst99 code
#
# Usage: ./test_suite.sh [options]
#   Options:
#     --quick    Run only basic tests
#     --full     Run all tests including long-running ones
#     --help     Display this help message
#==============================================================================

# Set up environment
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
TEST_OUTPUT_DIR="$SCRIPT_DIR/output"
TEST_DATA_DIR="$SCRIPT_DIR/test_data"
TEST_CASES_DIR="$SCRIPT_DIR/test_cases"
REFERENCE_DIR="$SCRIPT_DIR/reference_outputs"

# Set terminal colors
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Default test mode
TEST_MODE="quick"

# Parse command-line arguments
for arg in "$@"; do
  case "$arg" in
    --quick)
      TEST_MODE="quick"
      ;;
    --full)
      TEST_MODE="full"
      ;;
    --help)
      echo "Usage: ./test_suite.sh [options]"
      echo "Options:"
      echo "  --quick    Run only basic tests (default)"
      echo "  --full     Run all tests including long-running ones"
      echo "  --help     Display this help message"
      exit 0
      ;;
    *)
      echo "Unknown option: $arg"
      echo "Use --help for usage information"
      exit 1
      ;;
  esac
done

# Create necessary directories
mkdir -p "$TEST_OUTPUT_DIR" "$TEST_DATA_DIR" "$TEST_CASES_DIR" "$REFERENCE_DIR"

# Clean previous test outputs
clean_outputs() {
  echo -e "${BLUE}Cleaning previous test outputs...${NC}"
  rm -f "$TEST_OUTPUT_DIR"/*
}

# Function to check if galaxy executable exists
check_executable() {
  if [ ! -e "$ROOT_DIR/galaxy" ] && [ ! -e "$ROOT_DIR/bin/galaxy" ]; then
    echo -e "${RED}ERROR: Galaxy executable not found.${NC}"
    echo "Build the code first with: cd $ROOT_DIR && make"
    exit 1
  fi
}

# Function to create test case input file
create_test_input() {
  local name="$1"
  local sf_mode="$2"  # Star formation mode
  local metallicity="$3"
  local imf_exp="$4"  # IMF exponent
  
  # Ensure the test cases directory exists
  mkdir -p "$TEST_CASES_DIR"
  
  local output_file="$TEST_CASES_DIR/${name}.input"
  
  cat > "$output_file" << EOF
$name     Model designation
$sf_mode     Star formation mode: >0 continuous, -1 instantaneous
1.0     Total stellar mass (10^6 Msun)
1.0     Star formation rate (Msun/yr)
1     Number of intervals for the IMF
$imf_exp     Exponent for the IMF
0.1 100.     Mass boundaries for the IMF (Msun)
8.0     SNII cutoff mass (Msun)
120.0     Black hole cutoff mass (Msun)
$metallicity     Metallicity + tracks identifier
0     Wind model
0.0    Initial time (10^6 yr)
0     Time scale: 0=linear, 1=logarithmic
1.0     Time step (10^6 yr) if jtime=0
10     Number of steps if jtime=1
100.0     Last grid point (10^6 yr)
0     Atmosphere model
0     Metalicity of high-resolution models
1     UV line spectrum (1=solar, 2=LMC/SMC)
0 0     RSG parameters for individual features
1 1 1 1 1     Output files 1-5
1 1 1 1 1     Output files 6-10
1 1 1 1 1     Output files 11-15
EOF
  
  echo "Created test input: $output_file"
  echo "$output_file"
}

# Function to run a single test case
run_test_case() {
  local test_name="$1"
  local input_file="$2"
  local output_prefix="$test_name"
  
  echo -e "${BLUE}Running test: ${YELLOW}$test_name${NC}"
  
  # Prepare test environment
  cd "$TEST_OUTPUT_DIR"
  
  # Create symlinks to data directories
  if [ ! -e tracks ]; then ln -sf "$ROOT_DIR/data/tracks" tracks; fi
  if [ ! -e lejeune ]; then ln -sf "$ROOT_DIR/data/lejeune" lejeune; fi
  if [ ! -e auxil ]; then ln -sf "$ROOT_DIR/data/auxil" auxil; fi
  
  # Copy input file to fort.1
  cp "$input_file" fort.1
  
  # Run the code - check in both root and bin directories
  echo "  Executing galaxy with input: $(basename "$input_file")"
  if [ -x "$ROOT_DIR/galaxy" ]; then
    "$ROOT_DIR/galaxy" > "${test_name}.log" 2>&1
  elif [ -x "$ROOT_DIR/bin/galaxy" ]; then
    "$ROOT_DIR/bin/galaxy" > "${test_name}.log" 2>&1
  else
    echo "  ERROR: Cannot find executable galaxy in $ROOT_DIR or $ROOT_DIR/bin"
    return 1
  fi
  local result=$?
  
  if [ $result -ne 0 ]; then
    echo -e "${RED}  ✗ Test failed with exit code $result${NC}"
    echo "  See ${test_name}.log for details"
    return 1
  fi
  
  # Save outputs
  "$ROOT_DIR/scripts/save_output" "$output_prefix" 1 >> "${test_name}.log"
  
  # Verify outputs
  if [ -e "${output_prefix}.spectrum1" ]; then
    echo -e "${GREEN}  ✓ Test passed${NC}"
    return 0
  else
    echo -e "${RED}  ✗ Test failed: No spectrum output generated${NC}"
    return 1
  fi
}

# Function to verify output against reference
verify_output() {
  local test_name="$1"
  local output_file="${TEST_OUTPUT_DIR}/${test_name}.spectrum1"
  local reference_file="${REFERENCE_DIR}/${test_name}.spectrum1"
  
  # If reference doesn't exist, create it
  if [ ! -e "$reference_file" ]; then
    echo "  Creating reference output for future tests"
    cp "$output_file" "$reference_file"
    return 0
  fi
  
  # Compare outputs (first line is header, skip it)
  echo "  Comparing with reference output"
  local diff_result=$(diff <(tail -n +2 "$output_file") <(tail -n +2 "$reference_file"))
  
  if [ -z "$diff_result" ]; then
    echo -e "${GREEN}  ✓ Output matches reference${NC}"
    return 0
  else
    echo -e "${RED}  ✗ Output differs from reference${NC}"
    echo "  Differences found: $(echo "$diff_result" | wc -l) lines"
    return 1
  fi
}

# Function to run a set of tests
run_test_group() {
  local group_name="$1"
  local total_tests=0
  local passed_tests=0
  
  echo -e "\n${BLUE}===== Running $group_name tests =====${NC}\n"
  
  # Run tests specific to this group
  case "$group_name" in
    "Basic")
      # Test different star formation modes
      test_input=$(create_test_input "test_instantaneous" "-1" "5" "2.35")
      run_test_case "test_instantaneous" "$test_input"
      [ $? -eq 0 ] && ((passed_tests++))
      ((total_tests++))
      
      test_input=$(create_test_input "test_continuous" "1" "5" "2.35")
      run_test_case "test_continuous" "$test_input"
      [ $? -eq 0 ] && ((passed_tests++))
      ((total_tests++))
      ;;
      
    "IMF")
      # Test different IMF slopes
      test_input=$(create_test_input "test_imf_salpeter" "-1" "5" "2.35")
      run_test_case "test_imf_salpeter" "$test_input"
      [ $? -eq 0 ] && ((passed_tests++))
      ((total_tests++))
      
      test_input=$(create_test_input "test_imf_top_heavy" "-1" "5" "1.5")
      run_test_case "test_imf_top_heavy" "$test_input"
      [ $? -eq 0 ] && ((passed_tests++))
      ((total_tests++))
      
      test_input=$(create_test_input "test_imf_bottom_heavy" "-1" "5" "3.0")
      run_test_case "test_imf_bottom_heavy" "$test_input"
      [ $? -eq 0 ] && ((passed_tests++))
      ((total_tests++))
      ;;
      
    "Metallicity")
      # Test different metallicities
      # 1=0.001, 2=0.004, 3=0.008, 4=0.02, 5=0.04 Z_sun (Geneva tracks)
      for z in 1 2 3 4 5; do
        test_input=$(create_test_input "test_z$z" "-1" "$z" "2.35")
        run_test_case "test_z$z" "$test_input"
        [ $? -eq 0 ] && ((passed_tests++))
        ((total_tests++))
      done
      ;;
  esac
  
  echo -e "\n${BLUE}===== $group_name tests: $passed_tests/$total_tests passed =====${NC}\n"
  return $((total_tests - passed_tests))
}

# Main function to run test suite
run_tests() {
  local total_failures=0
  
  # Check if executable exists
  check_executable
  
  # Clean previous outputs
  clean_outputs
  
  # Always run basic tests
  run_test_group "Basic"
  total_failures=$((total_failures + $?))
  
  # Run IMF tests
  run_test_group "IMF"
  total_failures=$((total_failures + $?))
  
  if [ "$TEST_MODE" = "full" ]; then
    # Run additional tests for full mode
    run_test_group "Metallicity"
    total_failures=$((total_failures + $?))
  fi
  
  # Summary
  echo -e "\n${BLUE}===== Test Summary =====${NC}"
  if [ $total_failures -eq 0 ]; then
    echo -e "${GREEN}All tests passed!${NC}"
    return 0
  else
    echo -e "${RED}$total_failures test(s) failed.${NC}"
    return 1
  fi
}

# Execute the test suite
run_tests
exit $?