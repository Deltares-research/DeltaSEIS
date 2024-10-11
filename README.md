# deltaseis
[![License: MIT](https://img.shields.io/pypi/l/imod)](https://choosealicense.com/licenses/mit)
[![Lifecycle: experimental](https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![Build: status](https://gitlab.com/deltares/tgg-projects/subsurface-toolbox/pysst/badges/main/pipeline.svg)](https://gitlab.com/deltares/tgg-projects/subsurface-toolbox/pysst/-/pipelines)
[![Coverage](https://gitlab.com/deltares/tgg-projects/subsurface-toolbox/pysst/badges/main/coverage.svg)](https://gitlab.com/deltares/tgg-projects/subsurface-toolbox/pysst/-/pipelines)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)
[![Formatting: ruff](https://camo.githubusercontent.com/238a858d2190f028e2acde6cf05c6f71b67b3ad439e2788a518ae40f8fc2d0e3/68747470733a2f2f696d672e736869656c64732e696f2f62616467652f6c696e7465722d727566662d726564)](https://github.com/charliermarsh/ruff)

The Geophysical Modelling and Processing Toolbox (deltaseis) package is designed to handle all common formats of geophysical data (for now only seismic data). It provides selection, editing, processing, analysis, and export methods that can be applied generically to the loaded data. It is designed to connect with other Deltares developments such as [iMod](https://gitlab.com/deltares/imod) and [DataFusionTools](https://bitbucket.org/DeltaresGEO/datafusiontools/src/master/).

Both the **Seismic** and **Segy_editor** classes make heavey use of the [segyio](https://segyio.readthedocs.io/en/latest/index.html) and [obspy](https://docs.obspy.org/) packages as well as numpy ans scipy for data handling and signal processing.

## Installation process
The installation uses package manager pixi, for installation options see https://pixi.sh/latest/

To install pixi on windows, in powershell type:
```
winget install prefix-dev.pixi
```
Now clone deltaseis to your local drive using:
```
git clone https://gitlab.com/deltares/tgg-projects/subsurface-toolbox/deltaseis.git
```
Then navigate into that folder with:
```
cd deltaseis
```
To create a conda enviroment and install deltaseis in it type:
```
pixi run start
```

## Update deltaseis
To update deltaseis with the latest version from gitlab, open a shell in the deltaseis folder and:
```
git pull
```
And the same as with installation type:
```
pixi run start
```

## Tutorial
Activate the deltaseis_env in a command prompt and type 
```
pixi run notebook
```
In the browser that opens up, navigate to deltaseis/tutorials and click on segy_editing.ipynb to start a tutorial on how to edit seismic data files called seg-y using the class **Segy_editor**
## Usage
The same tutorial folder there is a .py script that has the same commands as the notebook, that you can use as a template for using the Segy_editor.
In the below example we load a segy file for which we like to change the record length to 55 ms and write a copy that reflects that change:
```
from deltaseis import Seg_editor
seismic_edit = Segy_editor(path/to/segy_infile)
seismic_edit.set_record_lenght(55)
seismic_edit.write(path/to/segy_outfile)
```
The Seismic class can be called similarly but has a simple 2D data array as input and performs more advanced processing to the data and not just edits. 
In the below example applies a time-squared gain and a bandpass filter between 80 and 5000 Hz to the data:
```
seismic = Seismic(data_array)
seismic.time_squared_gain()
seismic.bandpass_filter(lowcut=80, highcut=5000)
```
## Supported geophysics formats
- Seismic files
    - Post-stack seismic data (.seg-y)
    - Pre-stack seismic data (seg2, segd, dat)
    - Synthetic data from SPECFEM (.semd)
    - Distributed Acoustic Sensing (DAS) files (Silixa .tdms)


## Features
Todo

## Roadmap
Todo

## Contributing

You can contribute by testing, raising issues and making pull requests. Some general guidelines:

- Use new branches for developing new features or bugfixes. Use prefixes such as feature/ bugfix/ experimental/ to indicate the type of branch
- Add unit tests (and test data) for new methods and functions. We use pytest.
- Add Numpy-style docstrings
- Use Black formatting with default line lenght (88 characters)
- Update requirement.txt en environment.yml files if required

## License
MIT license
