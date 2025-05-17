# Galaxy/Starburst99 Test Suite

This directory contains a comprehensive test suite for the Galaxy/Starburst99 code. The tests are designed to validate the functionality and correctness of the code from multiple perspectives:

## Test Types

1. **Unit Tests** - Test individual functions and modules
   - `unit_tests.f90` - Tests basic utility functions, module initialization, error handling, etc.
   - `run_unit_tests.sh` - Script to compile and run the unit tests

2. **Astronomical Tests** - Validate the physics calculations
   - `astro_tests.f90` - Tests blackbody radiation, stellar luminosity, IMF, etc.
   - `numerical_tests.sh` - Tests numerical precision and stability

3. **Integration Tests** - Test the system as a whole
   - `run_integration_tests.sh` - Creates test inputs and runs the complete program
   - `test_suite.sh` - Comprehensive functional testing

## Running Tests

The test suite can be run using the following commands:

```bash
# Run all tests
./run_all_tests.sh

# Run only quick tests
./run_all_tests.sh --quick

# Run comprehensive tests
./run_all_tests.sh --full
```

## Test Directories

- `output/` - Contains test output files
- `reference_outputs/` - Contains expected output files for comparison
- `test_cases/` - Contains specific test case definitions
- `test_data/` - Contains input data for tests

## Adding New Tests

To add new tests:

1. For unit tests, add a new test subroutine to `unit_tests.f90`
2. For astronomical tests, add a new test subroutine to `astro_tests.f90`
3. For integration tests, create a new test input file and add it to `run_integration_tests.sh`

## Test Coverage

The current test suite covers:

- Basic utilities (exp10, linear_interp, etc.)
- Module initialization and cleanup
- Error handling and reporting
- File I/O operations
- Blackbody radiation calculations
- Stellar luminosity calculations
- Initial Mass Function implementations
- End-to-end program execution with various inputs