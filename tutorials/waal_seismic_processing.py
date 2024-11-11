
"""
This script processes seismic data from the river Waal using the DeltaSEIS library.
The script performs the following steps for each SEG-Y file:
1. Reads the SEG-Y file and sets the coordinate reference system (CRS) and endian format.
2. Applies vertical corrections to the trace data based on recording delays.
3. Processes the trace data, including applying time power gain and trace averaging.
4. Despikes the navigation data using smoothing and median filtering.
5. Picks the seabed horizon and applies vertical corrections based on bathymetry data.
6. Filters the seabed horizon and calculates corrections for differences due to changes in water level or sensor depth.
7. Generates quality control (QC) plots and writes the processed data to a new SEG-Y file.
Dependencies:
- deltaseis
- matplotlib
- pathlib
- numpy
- xarray
Inputs:
- List of SEG-Y files to be processed.
- Bathymetry raster file.
Outputs:
- Processed SEG-Y files with "_PRC" suffix.
- QC plots saved as images.
"""
from deltaseis import Segy_edit, Seismic
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import xarray as xr

segy_files = [Path(r"D:\Projects\PES waal\segy\segy_full_renamed\S2\S2_20210623_3_N.sgy"), 
              Path(r"D:\Projects\PES waal\segy\segy_full_renamed\S2\S2_20210623_3_ns.sgy"),
              Path(r"D:\Projects\PES waal\segy\segy_full_renamed\S2\S2_20210623_3_sn.sgy"),
              Path(r"D:\Projects\PES waal\segy\segy_full_renamed\S2\S2_20210623_3_S.sgy")]

raster = Path(r"D:\Projects\PES waal\qgis\20210623PES waal S2_Area 3_RD_NAP_25cm.nc")
bathy = xr.open_dataset(raster)['__xarray_dataarray_variable__']

for i, segy_file in enumerate(segy_files):

    print(f"{i+1}/{len(segy_files)}: processing {segy_file.stem}")
    edit = Segy_edit(segy_file)
    edit.set_crs(28992)
    edit.set_endian('little')

    #vertical corrections due to diffent recording delays
    manual_recording_delay = 4 #ms
    shift = edit.recording_delay - manual_recording_delay
    edit.vertical_trace_corrections(shift)
    edit.recording_delay = np.full(len(edit.recording_delay), manual_recording_delay)
                
    #trace data processing
    dx_mean = edit.factor*edit.shot_point_interval.mean()
    fs = edit.sampling_rate
    data = np.array(edit.trace_data).T

    print("Start data processing")

    seis = Seismic(data, fs, dx_mean)
    seis.time_power_gain(2)
    seis.trace_averaging(1)
    seis.convert_to_trace_data(edit.data_sample_format)
    edit.trace_data = seis.data
    edit.gain_type = 2 #t2 gain

    #despike navigation
    edit.fix_navigation_median(maximum_dist=1, save_figures=True)

    print("Start seabed pick")
    ## import bathymetry and verical correction with seabed pick
    edit.get_seabed_pick(10, 100, 9, 3, truncate = 10)
    
    print("Start Horizon processing")
    edit.filter_horizon_savgol('seabed_pick', 'seabed_pick_savgol', 51, 4)
    edit.add_horizon("seabed_multiple", 2*edit.seabed_pick_savgol + edit.recording_delay)
    edit.read_grid(raster, epsg_grid=28992, horizon_name='bathy', velocity=1489)

    # set the seismic datum 2 m (= 2.66 ms TWTT) above NAP
    edit.bathy = edit.bathy + (2/0.75)

    # calculate and apply correction for differences due to changes in water level and/or sensor depth and apply long wavelength filter
    edit.calculate_difference_horizon('bathy', 'seabed_pick', difference_horizon_name ='difference_stage_raw')
    edit.filter_horizon_savgol('difference_stage_raw', 'difference_stage', 40001, 1)
    edit.vertical_trace_corrections(edit.difference_stage)

    print("Start QC plot and write to SEG-Y")

    # save QC plot and write to SEG-Y
    edit.plot(save_plot=True, clip=0.9)
    processed_file = segy_file.with_stem(f"{segy_file.stem}_PRC")
    edit.write(processed_file)