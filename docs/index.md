---
layout: default
title: DeltaSEIS Documentation
---

# DeltaSEIS

Doing cool stuff with seismic data.

## Overview

DeltaSEIS is a Python package for seismic data processing and analysis. It provides tools for:

- SEG-Y file editing and manipulation
- Seismic data processing and visualization  
- Quality control and analysis
- Marine seismic pulse analysis

## Installation

```bash
pip install deltaseis
```

## Quick Start

```python
from deltaseis import Seismic, Segy_edit

# Load and process seismic data
segy_file = "your_data.sgy"
editor = Segy_edit(segy_file)
editor.plot()
```

## Documentation

For detailed documentation, see the [README.md](https://github.com/Deltares-research/DeltaSEIS/blob/main/README.md) file.

## Examples

Check out the [examples directory](https://github.com/Deltares-research/DeltaSEIS/tree/main/examples) for usage examples.

## License

This project is licensed under the MIT License.