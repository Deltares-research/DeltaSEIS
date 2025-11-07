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
output_suffix = "trace_av1_full_heave_corrected_time_power_gain3_decon_update"

grid_path = resources.files('deltaseis.data') / "bathy_clip_bergen0111a.asc"


#%% loop
segy_paths = [f for f in segy_folder.iterdir() if f.suffix.lower() in {'.sgy', '.seg', '.segy'}]

for i, segy_path in enumerate(segy_paths):
    print(f"{i+1}/{len(segy_paths)}: {segy_path.name}")

    s = deltaseis.Segy_edit(segy_path) 
    s.xstar_split('full')

    # set the input crs WGS 84 and transform to ETRS89 UTM31N
    s.set_crs(4326)
    s.transform_coordinates(25831)

    #read grid and extract values coinciding with seg-y coordinates
    s.read_grid(grid_path, 4326, horizon_name='bathy')
    s.get_seabed_pick(10, 100, 9, 3, truncate=10)
    s.plot(save_plot=True, clip=0.98)

    # filter heave from seabed pick, calculate the difference with the original and apply as vertical corrections to the data
    s.filter_horizon_savgol('seabed_pick', 'seabed_pick_savgol', 501, 4)
    s.calculate_difference_horizon('seabed_pick_savgol', 'seabed_pick', difference_horizon_name = 'heave')
    s.vertical_trace_corrections(s.heave)
    s.plot(save_plot=True, clip=0.98)

    # data processing
    trace_array = np.array(s.trace_data)
    print(f"Raw trace_data shape: {trace_array.shape}")
    # Seismic expects (samples × traces), so check if transpose is needed
    if trace_array.shape[0] == 12404 and trace_array.shape[1] == 3256:
        # Data is (traces × samples), need to transpose
        seis = deltaseis.Seismic(trace_array.T, fs=50_000, dx=0.4)
    else:
        # Data is already (samples × traces)
        seis = deltaseis.Seismic(trace_array, fs=50_000, dx=0.4)
    print(f"Seismic data shape: {seis.data.shape} (samples × traces)")
    seis.time_power_gain(3)
    seis.trace_averaging(1)
    seis.signature_deconvolution(trace_number=6415,
                                 start_time_ms=30.5,
                                 end_time_ms=31.7,
                                 method='wiener',
                                 epsilon=0.01,
                                 prewhiten=True,
                                 prewhiten_percent=1.0)     

    #seis.bandpass_filter(lowcut=1000, highcut=20000)  # Adjust to your data   
    s.trace_data = seis.data.T
    
    s.renumber_shotpoints()
    s.write(output_folder / f"{segy_path.stem}_{output_suffix}.sgy")

