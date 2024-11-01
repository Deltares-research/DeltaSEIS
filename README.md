# deltaseis
[![License: MIT](https://img.shields.io/pypi/l/imod)](https://choosealicense.com/licenses/mit)
[![Lifecycle: experimental](https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![Formatting: ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/charliermarsh/ruff)

The Deltares Seismic (DeltaSEIS) package is designed to handle all common formats of seismic data. This includes conventional seismic data, DAS fiber optic and simulated data. It provides selection, editing, processing, analysis, and export methods that can be applied generically to the loaded data. It is designed to connect with other Deltares developments such as [iMod](https://gitlab.com/deltares/imod) and [DataFusionTools](https://bitbucket.org/DeltaresGEO/datafusiontools/src/master/).

Both the **Seismic** and **Segy_editor** classes make heavey use of the [segyio](https://segyio.readthedocs.io/en/latest/index.html) and [obspy](https://docs.obspy.org/) packages as well as numpy ans scipy for data handling and signal processing.

## Installation process
The installation uses package manager pixi, for installation options see https://pixi.sh/latest/

To install pixi on windows, in powershell type:
```
winget install prefix-dev.pixi
```
Now clone deltaseis to your local drive using:
```
git clone https://github.com/Deltares-research/deltaseis.git
```
Then navigate into that folder with:
```
cd deltaseis
```
To create a conda enviroment and install deltaseis in it type:
```
pixi run install
```

## Update deltaseis
To update deltaseis with the latest version from gitlab, open a shell in the deltaseis folder and:
```
git pull
```
And the same as with installation type:
```
pixi run install
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
