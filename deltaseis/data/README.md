# DeltaSEIS Data Files

This directory contains example and test data files used in tutorials and examples.

## Data Files

Due to their large size, binary data files are not tracked in the Git repository. This includes:
- Seismic data files (`.sgy`, `.seg-y`, `.seg2`, `.segd`, `.semd`, `.tdms`)
- Image files (`.png`, `.jpg`, `.jpeg`, `.gif`, `.bmp`, `.tiff`)
- Video files (`.avi`, `.mp4`, `.mov`)
- Processed data files (`.pkl`, `.pickle`)
- Grid data files (`.asc`)
- PDF documentation (`.pdf`)

## Getting the Data Files

### Option 1: Download from Release Assets
Data files will be provided as release assets. Check the [Releases page](https://github.com/Deltares-research/DeltaSEIS/releases) for downloadable data packages.

### Option 2: Generate Your Own
Many of the example data files can be generated using the tools provided in DeltaSEIS. Refer to the tutorials for instructions.

### Option 3: Use Your Own Data
DeltaSEIS is designed to work with your own seismic data files. Simply place your data files in this directory or specify the path when loading data.

## Data Files Previously Included

The following files were previously tracked in Git but have been removed to reduce repository size:
- `sfs.sgy` (~14 MB) - Seismic survey data
- `bergen0111a.sgy` (~13.5 MB) - Bergen survey data
- `bergen0111a_TUTORIAL.sgy` (~10.8 MB) - Tutorial seismic data
- `sfs_SH_survey.sgy` (~6.9 MB) - SH wave survey
- `bergen0111a_2.png` (~3.6 MB) - Survey visualization
- `bergen0111a_1.png` (~3.5 MB) - Survey visualization
- `bergen0111a_3.png` (~3.0 MB) - Survey visualization
- `graafgeluid.png` (~2.7 MB) - Visualization
- `bathy_clip_bergen0111a.asc` (~3.6 MB) - Bathymetry data
- `bathy_clip_bergen0111a_interpolator.pkl` (~1.3 MB) - Interpolator data
- `interferentie_beverhol.png` (~1.4 MB) - Visualization
- Video files in `examples/` directory

## Note for Contributors

Please do not commit large binary files to the repository. If you need to share example data:
1. Keep example data files small (< 100 KB if possible)
2. For larger datasets, provide instructions on how to generate or download them
3. Consider using Git LFS for essential binary files if absolutely necessary
