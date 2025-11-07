# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 12:53:28 2025

@author: vermaas
"""

from pathlib import Path
from importlib import resources
import deltaseis
import numpy as np
import os
from matplotlib import pyplot as plt

segy_folder = Path(r'D:\Projects\DIS_Offshore\data\xstar')

# Ensure output subdirectory 'deltaseis' exists
output_folder = segy_folder / 'deltaseis'
output_folder.mkdir(parents=True, exist_ok=True)

# Output file suffix
output_suffix = "trace_av1_full"

grid_path = resources.files('deltaseis.data') / "bathy_clip_bergen0111a.asc"


#%% loop
segy_paths = [f for f in segy_folder.iterdir() if f.suffix.lower() in {'.sgy', '.seg', '.segy'}]

for i, segy_path in enumerate(segy_paths):
    print(f"{i+1}/{len(segy_paths)}: {segy_path.name}")

    s = deltaseis.Segy_edit(segy_path) 
    s.xstar_split('full')
    #s.set_record_length(30)

    # set the input crs WGS 84 and transform to ETRS89 UTM31N
    s.set_crs(4326)
    s.transform_coordinates(25831)

    # data processing
    seis = deltaseis.Seismic(np.array(s.trace_data).T,20,0.4)
    # seis.time_power_gain()
    seis.trace_averaging(1)
        
    s.trace_data = seis.data.T
    #s.trace_data = (32767/np.max(data))*data #normalize to fit in int16 format
    
    s.renumber_shotpoints()
    s.write(output_folder / f"{segy_path.stem}_{output_suffix}.sgy")
