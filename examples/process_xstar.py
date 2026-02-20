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

segy_folder = Path(r'C:\Projects\DIS_Offshore\data\AquaPulse')

# Ensure output subdirectory 'deltaseis' exists
output_folder = segy_folder / 'deltaseis'
output_folder.mkdir(parents=True, exist_ok=True)

# Output file suffix
output_suffix = "full_heave_gain_bpf_trav_mute"
color_map = "Greys"

# Read bathymetry grid
grid_path = resources.files('deltaseis.data') / "bathy_clip_bergen0111a.asc"


# Process each SEG-Y file in the folder
segy_paths = [f for f in segy_folder.iterdir() if f.suffix.lower() in {'.sgy', '.seg', '.segy'}]
segy_paths = [segy_paths[3]]

# Or select specific files
# segy_paths = [segy_folder / 'B01_BX03.003_echo_selection_envelope2.sgy',]

for i, segy_path in enumerate(segy_paths):
    print(f"{i+1}/{len(segy_paths)}: {segy_path.name}")

    s = deltaseis.Segy_edit(segy_path) 
    #s.xstar_split('envelope')
    s.set_record_length(70)

    # set the input crs WGS 84 and transform to ETRS89 UTM31N
    #s.set_crs(4326)
    #s.transform_coordinates(25831)
    s.set_crs(25831)
    
    #read grid and extract values coinciding with seg-y coordinates
    s.read_grid(grid_path, 4326, horizon_name='bathy')
    s.get_seabed_pick(10, 100, 9, 3, truncate=10)
    s.plot(save_plot=True, clip=0.1, cmap=color_map, show_horizons=False)

    # filter heave from seabed pick, calculate the difference with the original and apply as vertical corrections to the data
    s.filter_horizon_savgol('seabed_pick', 'seabed_pick_savgol', 501, 4)
    s.calculate_difference_horizon('seabed_pick_savgol', 'seabed_pick', difference_horizon_name = 'heave')
    s.vertical_trace_corrections(s.heave)
    s.plot(save_plot=True, clip=0.1, cmap=color_map, show_horizons=False)

    # data processing
    seis = deltaseis.Seismic(np.array(s.trace_data).T, fs=50_000, dx=0.4)
    seis.time_power_gain(1.2)

    s.trace_data = seis.data.T
    s.plot(save_plot=True, clip=0.1, cmap=color_map, show_horizons=False)

    # seis.signature_deconvolution(trace_number=6414,
    #                              start_time_ms=30.5,
    #                              end_time_ms=31.7,
    #                              method='wiener',
    #                              epsilon=0.5,
    #                              prewhiten=True,
    #                              prewhiten_percent=1.0)     

    
    seis.bandpass_filter(lowcut=2700, highcut=8000)
    s.trace_data = seis.data.T
    s.plot(save_plot=True, clip=4.0, cmap=color_map, show_horizons=False)

    seis.time_frequency_denoise(threshold=0.12)
    s.trace_data = seis.data.T
    s.plot(save_plot=True, clip=0.1, cmap=color_map, show_horizons=False)

    seis.adaptive_noise_filter(reference_trace=213, 
                               reference_start_ms=60, 
                               reference_end_ms=67,
                               filter_length=128, 
                               max_lag_ms=2.0, 
                               mu=0.02)
    s.trace_data = seis.data.T
    s.plot(save_plot=True, clip=0.1, cmap=color_map, show_horizons=False)
    
    seis.trace_averaging(1)
    s.trace_data = seis.data.T
    # s.plot(save_plot=True, clip=0.1, cmap=color_map, show_horizons=False)
    
    # seis.top_mute(s.seabed_pick_savgol, shift_ms=0)
    # s.trace_data = seis.data.T
    # s.plot(save_plot=True, clip=0.1, cmap=color_map, show_horizons=False)   
    
    s.renumber_shotpoints()
    s.write(output_folder / f"{segy_path.stem}_{output_suffix}.sgy")

