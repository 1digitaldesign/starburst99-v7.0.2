#!/bin/bash
#==============================================================================
# Starburst99 Build Script
#==============================================================================
# This script builds the Starburst99 population synthesis code and sets up
# the runtime environment.
#
# Usage: ./build.sh [clean|debug|release]
#   - No arguments: Standard build
#   - clean: Remove build artifacts
#   - debug: Build with debug options
#   - release: Build optimized release version
#==============================================================================

# Default build type
BUILD_TYPE="standard"

# Parse command-line arguments
if [ $# -ge 1 ]; then
  BUILD_TYPE="$1"
fi

# Set terminal colors
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}Starburst99 Build System${NC}"
echo "=============================="

# Execute the appropriate build command
case "$BUILD_TYPE" in
  "clean")
    echo -e "${YELLOW}Cleaning build artifacts...${NC}"
    make clean
    ;;
    
  "debug")
    echo -e "${YELLOW}Building with debug options...${NC}"
    make debug
    ;;
    
  "release")
    echo -e "${YELLOW}Building optimized release version...${NC}"
    make release
    ;;
    
  *)
    echo -e "${YELLOW}Building standard version...${NC}"
    make
    ;;
esac

# Setup the runtime environment
echo -e "\n${GREEN}Setting up runtime environment...${NC}"
make setup

echo -e "\n${GREEN}Build process complete${NC}"
echo "To run the code, use: make run"
echo "Or cd to the output directory: cd output && ./go_galaxy"