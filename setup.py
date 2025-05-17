"""Setup script for Starburst99 Python package"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="starburst99",
    version="7.0.2",
    author="Starburst99 Development Team",
    author_email="starburst99@example.com",
    description="A stellar population synthesis code",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/starburst99/starburst99",
    package_dir={"": "src/python"},
    packages=find_packages(where="src/python"),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "starburst99=starburst_main:main",
        ],
    },
    include_package_data=True,
    package_data={
        "starburst99": [
            "data/**/*",
            "config/*.yaml",
            "config/*.json",
        ],
    },
)