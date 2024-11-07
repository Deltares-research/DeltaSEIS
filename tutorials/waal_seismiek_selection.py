'''
- calculate distance
- trace selection based on 4 crosspoints (set manually by plotting the cross section with the qgis exports of the tracks)
- find closest point, calculate distance of all points from THAT specific point. Select traces based on distance < x e.g 50 on both sides
- conversion of all the raw data to make a clean repo would be an excellent idea (raw and processed)
'''	
#%%
from deltaseis import Segy_edit
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import geopandas as gpd

segy_files = [Path(r"D:\Projects\PES waal\segy\segy_full_renamed\S2\S2_20210623_3_N_PRC.sgy"), 
              Path(r"D:\Projects\PES waal\segy\segy_full_renamed\S2\S2_20210623_3_ns_PRC.sgy"),
              Path(r"D:\Projects\PES waal\segy\segy_full_renamed\S2\S2_20210623_3_S_PRC.sgy"),
              Path(r"D:\Projects\PES waal\segy\segy_full_renamed\S2\S2_20210623_3_sn_PRC.sgy")]

# get the crosspoints of the profile perpendicular to the segy trackes and the associated trace numbers
crosspoints = Path(r"D:\Projects\PES waal\qgis\PES_waal_survey2_area3_Trackpoints_crossings.shp")                       
crosspoints_gdf = gpd.read_file(crosspoints)
trace_numbers = crosspoints_gdf['Trace numb'].values.astype(int)
print(trace_numbers)

#%% traces within  'diastance' meters from the crosspoints will be select and written to a new segy file
distance = 50

for i, segy_file in enumerate(segy_files):

    print(f"{i+1}/{len(segy_files)}: processing {segy_file.stem}")
    edit = Segy_edit(segy_file)
    edit.set_crs(28992)


    distance_along_line = np.cumsum(edit.shot_point_interval)
    distance_from_cross = edit.line_distance - edit.line_distance[trace_numbers[i]]
    indices = np.where(np.abs(distance_from_cross) <= distance / edit.factor)[0]

    edit.select_traces(indices)

    selection_file = segy_file.with_stem(f"{segy_file.stem}_SEL")
    edit.write(selection_file)