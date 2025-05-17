# Large JSON Data Files

Several large JSON data files have been excluded from this repository due to GitHub file size limitations.

## Missing JSON Files

The following files are not included in the repository:
- `allstarscont_m05.json` (106.07 MB)
- `allstarscont_m10.json` (104.27 MB)
- `allstarscont_p00.json` (104.92 MB)
- `allstarscont_p03.json` (103.77 MB)
- `allstarsflux_m05.json` (103.54 MB)
- `allstarsflux_m10.json` (103.75 MB)
- `allstarsflux_p00.json` (104.14 MB)
- `allstarsflux_p03.json` (102.95 MB)
- `allstarswave.json` (Large file)

## How to Generate These Files

These JSON files can be automatically generated from the original data files. First, obtain the original data files as described in `data/lejeune/README_LARGE_FILES.md`, then run:

```bash
make convert-json
```

or

```bash
make convert-json-py
```

to convert all data files to JSON format.

## Purpose

These JSON files are used by the modern version of the Starburst99 code for more efficient data loading and processing. They are a transformed version of the original data files with the same scientific content but in a machine-readable format.