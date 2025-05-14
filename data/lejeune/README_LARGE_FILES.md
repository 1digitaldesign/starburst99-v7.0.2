# Large Data Files

Several large data files have been excluded from this repository due to GitHub file size limitations.

## Missing Files

The following files are not included in the repository:
- `allstarscont_m05.txt` (59.61 MB)
- `allstarscont_m10.txt` (58.62 MB)
- `allstarscont_p00.txt` (58.90 MB)
- `allstarscont_p03.txt` (58.20 MB)
- `allstarsflux_m05.txt` (59.23 MB)
- `allstarsflux_m10.txt` (59.37 MB)
- `allstarsflux_p00.txt` (59.66 MB)
- `allstarsflux_p03.txt` (58.94 MB)
- `allstarswave.txt` (Large file)

## How to Obtain These Files

These files can be obtained from the official Starburst99 website:
[http://www.stsci.edu/science/starburst99/](http://www.stsci.edu/science/starburst99/)

Or by contacting the authors directly.

## Installation

After obtaining these files, place them in this directory (`data/lejeune/`) to use the full functionality of the code.

## JSON Versions

The JSON versions of these files have also been excluded from the repository. They will be automatically generated if you run:

```bash
make convert-json
```

or

```bash
make convert-json-py
```

After placing the original data files in this directory.