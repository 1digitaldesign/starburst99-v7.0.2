#!/bin/bash
#==============================================================================
# Run All Starburst99 Tests
#==============================================================================
# This script runs all tests for the Starburst99 code, including:
# - Unit tests
# - Functional tests
# - Numerical tests
#
# Usage: ./run_all_tests.sh [options]
#   Options:
#     --quick    Run only essential tests (faster)
#     --full     Run all tests (more thorough)
#==============================================================================

# Set up environment
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

# Set terminal colors
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Build unit test program
build_unit_tests() {
  echo -e "${BLUE}Building unit tests...${NC}"
  
  # Compile unit tests
  cd "$ROOT_DIR"
  gfortran -std=f2018 -o "$SCRIPT_DIR/unit_test_runner" \
    "$ROOT_DIR/src/galaxy_module.f90" \
    "$ROOT_DIR/src/galaxy_module_error.f90" \
    "$SCRIPT_DIR/unit_tests.f90"
    
  local result=$?
  if [ $result -ne 0 ]; then
    echo -e "${RED}Failed to build unit tests${NC}"
    return 1
  fi
  
  echo -e "${GREEN}Unit tests built successfully${NC}"
  return 0
}

# Run unit tests
run_unit_tests() {
  echo -e "\n${BLUE}Running unit tests...${NC}"
  
  # Run unit test program
  cd "$SCRIPT_DIR"
  ./unit_test_runner
  
  local result=$?
  if [ $result -ne 0 ]; then
    echo -e "${RED}Unit tests failed${NC}"
    return 1
  fi
  
  echo -e "${GREEN}Unit tests passed${NC}"
  return 0
}

# Run functional tests
run_functional_tests() {
  echo -e "\n${BLUE}Running functional tests...${NC}"
  
  # Make test files executable
  chmod +x "$SCRIPT_DIR/test_suite.sh"
  
  # Run test suite
  "$SCRIPT_DIR/test_suite.sh" "$@"
  
  local result=$?
  if [ $result -ne 0 ]; then
    echo -e "${RED}Functional tests failed${NC}"
    return 1
  fi
  
  echo -e "${GREEN}Functional tests passed${NC}"
  return 0
}

# Run numerical tests
run_numerical_tests() {
  echo -e "\n${BLUE}Running numerical tests...${NC}"
  
  # Make test script executable
  chmod +x "$SCRIPT_DIR/numerical_tests.sh"
  
  # Run numerical tests
  "$SCRIPT_DIR/numerical_tests.sh"
  
  local result=$?
  if [ $result -ne 0 ]; then
    echo -e "${RED}Numerical tests failed${NC}"
    return 1
  fi
  
  echo -e "${GREEN}Numerical tests passed${NC}"
  return 0
}

# Main function
main() {
  local total_failures=0
  local test_args=("$@")
  
  # Build the main code if it doesn't exist
  if [ ! -e "$ROOT_DIR/galaxy" ]; then
    echo -e "${BLUE}Building main program...${NC}"
    cd "$ROOT_DIR"
    make
    
    if [ $? -ne 0 ]; then
      echo -e "${RED}Failed to build main program${NC}"
      exit 1
    fi
  fi
  
  # Build and run unit tests
  build_unit_tests
  if [ $? -eq 0 ]; then
    run_unit_tests
    total_failures=$((total_failures + $?))
  else
    total_failures=$((total_failures + 1))
  fi
  
  # Run functional tests
  run_functional_tests "${test_args[@]}"
  total_failures=$((total_failures + $?))
  
  # Run numerical tests if --full specified
  if [[ " ${test_args[*]} " == *" --full "* ]]; then
    run_numerical_tests
    total_failures=$((total_failures + $?))
  else
    echo -e "\n${YELLOW}Skipping numerical tests (use --full to run them)${NC}"
  fi
  
  # Summary
  echo -e "\n${BLUE}===== Test Summary =====${NC}"
  if [ $total_failures -eq 0 ]; then
    echo -e "${GREEN}All tests passed!${NC}"
    exit 0
  else
    echo -e "${RED}$total_failures test group(s) failed.${NC}"
    exit 1
  fi
}

# Run tests
main "$@"