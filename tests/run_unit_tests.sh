#!/bin/bash
#==============================================================================
# Run Starburst99 Unit Tests
#==============================================================================
# This script compiles and runs the unit tests for the Starburst99 core modules.
#
# Usage: ./run_unit_tests.sh
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

# Build and run unit tests
build_and_run_unit_tests() {
  echo -e "${BLUE}Building and running unit tests...${NC}"
  
  # Create test output directory if it doesn't exist
  mkdir -p "$SCRIPT_DIR/output"
  
  # Compile unit tests with extra validation
  cd "$ROOT_DIR"
  gfortran -std=f2018 -Wall -Wextra -fimplicit-none -fcheck=all -o "$SCRIPT_DIR/unit_test_runner" \
    "$ROOT_DIR/src/galaxy_module.f90" \
    "$ROOT_DIR/src/galaxy_module_error.f90" \
    "$SCRIPT_DIR/unit_tests.f90"
    
  local result=$?
  if [ $result -ne 0 ]; then
    echo -e "${RED}Failed to build unit tests${NC}"
    return 1
  fi
  
  # Run unit tests and capture output
  "$SCRIPT_DIR/unit_test_runner" > "$SCRIPT_DIR/output/unit_test_results.log" 2>&1
  result=$?
  
  # Display the test results
  cat "$SCRIPT_DIR/output/unit_test_results.log"
  
  if [ $result -ne 0 ]; then
    echo -e "${RED}Unit tests failed with exit code $result${NC}"
    return 1
  fi
  
  echo -e "${GREEN}All unit tests passed!${NC}"
  return 0
}

# Main function
main() {
  build_and_run_unit_tests
  exit $?
}

# Run tests
main