###############################################################################
# Makefile for the Starburst99 Population Synthesis Code
# - Modernized for Fortran 2018 standards
###############################################################################

#==============================================================================
# Project configuration
#==============================================================================
# Version information
VERSION = 7.1.0
# Executable name (maintains backward compatibility)
EXE = galaxy
# Directory structure
SRC_DIR = src
BIN_DIR = bin
SCRIPTS_DIR = scripts
TOOLS_DIR = tools
DATA_DIR = data
OUTPUT_DIR = output
DOCS_DIR = docs

#==============================================================================
# Source file structure
#==============================================================================
# Core module files (in order of dependency)
MODULE_SRCS = $(SRC_DIR)/galaxy_module.f90 $(SRC_DIR)/galaxy_module_error.f90 $(SRC_DIR)/galaxy_interface.f90 $(SRC_DIR)/data_profiles.f90 $(SRC_DIR)/data_profiles_init.f90

# Main program file (full version and fixed version)
MAIN_SRC = $(SRC_DIR)/starburst_main_fixed.f90

# Minimal version for testing
MINIMAL_SRC = $(SRC_DIR)/starburst_minimal.f90

# Default sources for standard compilation
SOURCES = $(MODULE_SRCS) $(MAIN_SRC)

# Sources for minimal test version
MINIMAL_SOURCES = $(MODULE_SRCS) $(MINIMAL_SRC)

# Object files derived from sources (keep in src directory for simplicity)
OBJECTS = $(SOURCES:.f90=.o)

# Object files for minimal version
MINIMAL_OBJECTS = $(MINIMAL_SOURCES:.f90=.o)

# Module files created by compilation
MODULES = $(MODULE_SRCS:.f90=.mod)

# JSON conversion program sources
JSON_SRCS = $(TOOLS_DIR)/converters/galaxy_dat2json.f90 $(TOOLS_DIR)/converters/convert_to_json.f90
JSON_OBJECTS = $(patsubst %.f90,%.o,$(JSON_SRCS))
JSON_EXE = $(TOOLS_DIR)/converters/convert_to_json

#==============================================================================
# Compiler settings
#==============================================================================
# Compiler selection
FC = gfortran

# Common compiler flags
COMMON_FLAGS = -std=f2018 -fimplicit-none -fall-intrinsics

# Production build flags
FFLAGS = -c -O2 $(COMMON_FLAGS) -fcoarray=single

# Debug build flags
DBGFLAGS = -c -g -Wall -Wextra $(COMMON_FLAGS) -fcheck=all -fbacktrace -pedantic -Wconversion

# Linker flags
LFLAGS = 

# External libraries
LIBS = 

#==============================================================================
# Build rules
#==============================================================================
# Default target builds the executable
all: $(EXE)

# Main build target
$(EXE): $(OBJECTS) 
	@echo "Linking $(EXE) (version $(VERSION))"
	$(FC) $(LFLAGS) $(OBJECTS) -o $(EXE) $(LIBS)
	@ln -sf ../$(EXE) $(BIN_DIR)/$(EXE)
	@echo "Build complete."

# Debug build with additional checks
debug: FFLAGS = $(DBGFLAGS)
debug: clean all
	@echo "Debug build complete with extra validation checks enabled."

# Release build with optimization
release: clean all
	@echo "Production build complete."
	@echo "Version: $(VERSION)"
	
# Build minimal version for testing
minimal: $(SRC_DIR)/galaxy_module.o $(SRC_DIR)/galaxy_module_error.o $(SRC_DIR)/starburst_minimal.o
	@echo "Building minimal version for testing..."
	$(FC) $(LFLAGS) $(MINIMAL_OBJECTS) -o $(EXE) $(LIBS)
	@ln -sf ../$(EXE) $(BIN_DIR)/$(EXE)
	@echo "Minimal version build complete."

# Build fixed version
fixed: $(SRC_DIR)/galaxy_module.o $(SRC_DIR)/galaxy_module_error.o $(SRC_DIR)/data_profiles.o $(SRC_DIR)/data_profiles_init.o src/starburst_main_fixed.o
	@echo "Building fixed version..."
	$(FC) $(LFLAGS) $^ -o $(EXE)_fixed $(LIBS)
	@ln -sf ../$(EXE)_fixed $(BIN_DIR)/$(EXE)_fixed
	@echo "Fixed version build complete."

#==============================================================================
# Module dependencies
#==============================================================================
# Module dependencies - explicit declaration to ensure correct build order
$(SRC_DIR)/galaxy_module_error.o: $(SRC_DIR)/galaxy_module.o
$(SRC_DIR)/galaxy_interface.o: $(SRC_DIR)/galaxy_module.o
$(SRC_DIR)/data_profiles.o: $(SRC_DIR)/galaxy_module.o
$(SRC_DIR)/data_profiles_init.o: $(SRC_DIR)/data_profiles.o
$(SRC_DIR)/starburst_main_fixed.o: $(SRC_DIR)/galaxy_module.o $(SRC_DIR)/galaxy_module_error.o $(SRC_DIR)/galaxy_interface.o $(SRC_DIR)/data_profiles.o

#==============================================================================
# Pattern rules
#==============================================================================
# Build rule for Fortran 2018 files
%.o: %.f90
	@echo "Compiling $<"
	$(FC) $(FFLAGS) -o $@ $<

#==============================================================================
# Utility targets
#==============================================================================
# Clean up intermediate files and executable
clean:
	@echo "Cleaning build artifacts..."
	rm -f $(SRC_DIR)/*.o $(SRC_DIR)/*.mod $(SRC_DIR)/*.smod $(EXE) $(BIN_DIR)/$(EXE)
	rm -f $(TOOLS_DIR)/converters/*.o $(JSON_EXE)
	@echo "Clean complete."

# Install the executable to user's bin directory
install: $(EXE) $(JSON_EXE)
	@echo "Installing $(EXE) and $(JSON_EXE) to ~/bin"
	mkdir -p ~/bin
	cp $(EXE) $(JSON_EXE) ~/bin/
	@echo "Installation complete."

# Setup directory structure and links
setup:
	@echo "Setting up runtime environment..."
	mkdir -p $(BIN_DIR) $(OUTPUT_DIR)
	ln -sf ../$(SCRIPTS_DIR)/go_galaxy $(OUTPUT_DIR)/go_galaxy
	@echo "Setup complete."

# Build the JSON converter program
json: $(JSON_EXE)

$(JSON_EXE): $(JSON_OBJECTS)
	@echo "Building JSON converter program"
	$(FC) $(LFLAGS) $(JSON_OBJECTS) -o $(JSON_EXE) $(LIBS)
	@echo "JSON converter built."

# Convert data files to JSON format (Fortran version)
convert-json: $(JSON_EXE)
	@echo "Converting data files to JSON format (Fortran converter)..."
	bash $(TOOLS_DIR)/converters/convert_all_to_json.sh
	@echo "Conversion complete."

# Convert data files to JSON format (Python version - handles binary files)
convert-json-py:
	@echo "Converting data files to JSON format (Python converter)..."
	bash $(TOOLS_DIR)/converters/convert_all_data.sh
	@echo "Conversion complete."

# Run the code
run: $(EXE) setup
	@echo "Running Starburst99..."
	cd $(OUTPUT_DIR) && ./go_galaxy
	@echo "Execution complete."

# Run all tests (requires minimal build for functional tests)
test: minimal test-unit
	@echo "Running Starburst99 tests..."
	cd tests && ./run_all_tests.sh --quick
	@echo "Tests complete."

# Run quick tests only (requires minimal build for functional tests)
test-quick: minimal test-unit
	@echo "Running quick Starburst99 tests..."
	cd tests && ./run_all_tests.sh --quick
	@echo "Tests complete."

# Run all tests including time-consuming ones (requires minimal build for functional tests)
test-full: minimal test-unit
	@echo "Running full Starburst99 test suite..."
	cd tests && ./run_all_tests.sh --full
	@echo "Tests complete."

# Run unit tests only
test-unit:
	@echo "Building and running unit tests..."
	./tests/run_unit_tests.sh
	@echo "Unit tests complete."

# Run single test
test-single: $(EXE)
	@echo "Running Starburst99 single test..."
	cd tests && ./test_run.sh
	@echo "Single test complete."

# Display help information
help:
	@echo "Starburst99 Makefile Help"
	@echo "-----------------------"
	@echo "Available targets:"
	@echo "  all           - Build the executable (default target)"
	@echo "  debug         - Build with debug flags and extra validation"
	@echo "  release       - Build optimized version for production use"
	@echo "  minimal       - Build minimal version for testing"
	@echo "  clean         - Remove all build artifacts"
	@echo "  setup         - Setup runtime environment"
	@echo "  install       - Install executable to ~/bin"
	@echo "  json          - Build the JSON converter program"
	@echo "  convert-json  - Convert all data files to JSON format (Fortran)"
	@echo "  convert-json-py - Convert all data files to JSON format (Python)"
	@echo "  run           - Run the Starburst99 code"
	@echo "  test          - Run all tests with default options"
	@echo "  test-quick    - Run only quick tests"
	@echo "  test-full     - Run all tests including time-consuming ones"
	@echo "  test-unit     - Run only unit tests"
	@echo "  test-single   - Run a single test case"
	@echo "  help          - Display this help message"

# Mark targets that don't produce files
.PHONY: all clean debug release minimal install setup help json convert-json convert-json-py run test test-quick test-full test-unit test-single