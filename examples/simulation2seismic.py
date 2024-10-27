# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 17:22:24 2024

Convert simulated data files to industry standard seismic files. 
No example model data distributed with package due to size restrictions.

1) seismic data simulated in SPECFEM2D (.semd files) are read 
2) for model stability the semd files are strongly oversampled so they have to be downsampled.
3) the resulting data is exported to industry standard (.seg-y) files for further processing as deliverable.
4) if snapshot images of the wavefield were output during the simulation, these will be made into a movie (.avi file)

@author: nieboer
"""
from deltaseis import read_semd, resample, export_sgy
from deltaseis.tools.moviemaker import create_video_from_images
from importlib import resources
from pathlib import Path 
import matplotlib.pyplot as plt
from re import findall

# required input parameters
component = "BXX"           # choose receiver component BXX for horizontal, BXZ for vertical
extension = ".semd"         # extension of simulation files (most commonly ".semd")
new_sampling_rate = 40000          # Hertz, choose new sampling rate
receiver_spacing = 0.5                    # meter, set receiver spacing
shotpoint_interval = 0.0    # meter, set shotpoint interval if relevant

# set folder that contains model runs (no model data distrubted with package due to size restrictions)
data_folder = Path(resources.files('deltaseis.data'))
#data_folder = Path(r'\\wsl.localhost\Ubuntu\home\rlnd\projects\DiggerDAS') # my specfem project folder in WSL2

sub_directories = sorted([d for d in data_folder.rglob('*') if d.name == 'OUTPUT_FILES' and d.is_dir()])
#sub_directories =[root_directory/'hole4/OUTPUT_FILES']  #use this line to select specific subdirectories

for directory in sub_directories:

    print(directory)

    # list of simulated data filepaths
    file_list = sorted([file for file in directory.iterdir() if file.is_file() and file.name.endswith(extension)])
    
    out_folder = Path(file_list[0]).parents[2]
    shot_name = Path(file_list[0]).parents[1].name
    
    # (1) read simulated data and convert to deltaseis format
    data, sampling_rate, t_start = read_semd(file_list, verbose=True)

    # (2) resample data (honor your Nyquist)
    data_resampled = resample(data, sampling_rate, new_sampling_rate)

    # (3) export deltaseis format to seismic data files (seg-y files)
    shot_number = int(findall(r'\d+', shot_name)[0])
    out_sgy =(out_folder / f"{shot_name}_{component}").with_suffix(".sgy")
    export_sgy(data_resampled, new_sampling_rate, receiver_spacing, shot_number, shotpoint_interval, out_sgy)

    # (4) create video from the snapshots created by specfem
    try:
        create_video_from_images(directory, directory.parent/f"movie_{directory}_{component}.avi", fps=5, cut_percentage=60)  
    except FileNotFoundError:
        print(f"No snapshot images found to be made into a video in folder: {directory}")



    
    
