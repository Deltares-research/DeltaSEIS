'''
"""
This script processes SEGY files to select seismic data traces that intersect with a user-defined cross section through the river Waal.
The script performs the following steps:
1. Reads SEGY files from specified paths.
2. Loads crosspoints from a shapefile which contains the trace numbers of the crosspoints.
3. For each SEGY file, calculates the distance of traces from the crosspoints.
4. Selects traces within a specified distance from the crosspoints.
5. Writes the selected traces to new SEGY files.
Modules used:
- deltaseis.Segy_edit: For editing SEGY files.
- matplotlib.pyplot: For plotting (not used in this script but imported).
- pathlib.Path: For handling file paths.
- numpy: For numerical operations.
- geopandas: For reading shapefiles.
Variables:
- segy_files: List of paths to SEGY files.
- crosspoints: Path to the shapefile containing crosspoints.
- crosspoints_gdf: GeoDataFrame containing the crosspoints data.
- trace_numbers: Array of trace numbers from the crosspoints.
- distance: Distance in meters within which traces will be selected.
"""

'''	
#%%
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np

from deltaseis import Segy_edit

segy_files = [Path(r"D:\Projects\PES waal\segy\S2\S2_20210623_3_N_PRC.sgy"), 
              Path(r"D:\Projects\PES waal\segy\S2\S2_20210623_3_ns_PRC.sgy"),
              Path(r"D:\Projects\PES waal\segy\S2\S2_20210623_3_S_PRC.sgy"),
              Path(r"D:\Projects\PES waal\segy\S2\S2_20210623_3_sn_PRC.sgy")]

# get the crosspoints of the profile perpendicular to the segy trackes and the associated trace numbers
crosspoints = Path(r"D:\Projects\PES waal\qgis\PES_waal_survey2_area3_Trackpoints_crossings.shp")
crosspoints_gdf = gpd.read_file(crosspoints)
trace_numbers = crosspoints_gdf['Trace numb'].values.astype(int)

# traces within  'distance' meters from the crosspoints will be select and written to a new segy file
distance = 50

for i, segy_file in enumerate(segy_files):

    print(f"{i+1}/{len(segy_files)}: processing {segy_file.stem}")
    edit = Segy_edit(segy_file)
    edit.set_crs(28992)

    # set the recording delay as the difference between the seismic datum and crs datum (NAP) in this case 2 m
    edit.recording_delay = np.full(len(edit.recording_delay), 2)

    # apply conversion to depth using a constant velocity of v=15000 m/s (1.5 m/ms) by updating the sampling interval
    # this is for the export of depth converted seisnc data. Note that edit.sampling_interval is assumed in ms, not conventional microseconds.
    velocity = 1500 # m/s
    edit.sampling_interval = int(np.round(edit.sampling_interval * (0.5 * -velocity*1e-3), 0))   #the sampling interval is now in meters
    print(f"Depth conversion applied using constant velocity = {velocity} m/s (depth sampling interval: {edit.sampling_interval} mm)")

    distance_along_line = np.cumsum(edit.shot_point_interval)
    distance_from_cross = edit.line_distance - edit.line_distance[trace_numbers[i]]
    indices = np.where(np.abs(distance_from_cross) <= distance / edit.factor)[0]

    edit.select_traces(indices)

    selection_file = segy_file.with_stem(f"{segy_file.stem}_DEPTH_SEL")
    edit.write(selection_file)