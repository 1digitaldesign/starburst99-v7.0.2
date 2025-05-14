#!/usr/bin/env python3
"""Run all Starburst99 Python tests."""

import sys
import subprocess
import os

def run_tests():
    """Run pytest on the src/python directory."""
    
    # Get the root directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    root_dir = os.path.dirname(script_dir)
    python_dir = os.path.join(root_dir, 'src', 'python')
    
    # Run pytest with coverage
    cmd = [
        sys.executable, '-m', 'pytest',
        python_dir,
        '-v',
        '--cov=src.python',
        '--cov-report=term-missing',
        '--cov-report=html'
    ]
    
    print("Running Python tests with coverage...")
    print(f"Command: {' '.join(cmd)}")
    
    result = subprocess.run(cmd, cwd=root_dir)
    
    if result.returncode == 0:
        print("\n✅ All tests passed!")
    else:
        print("\n❌ Tests failed!")
        
    return result.returncode

if __name__ == '__main__':
    sys.exit(run_tests())