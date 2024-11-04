#%%
from deltaseis import Segy_edit, Seismic
from pathlib import Path
import numpy as np
import xarray as xr

#segy_file = Path(r"D:\Projects\PES waal\segy\PES data Surv 2 Area 3\S2_A3_N_23.06.2021\merged.sgy")
segy_file = Path(r"D:\Projects\PES waal\segy\PES data Surv 2 Area 3\S2_A3_ns_23.06.2021\merged.sgy")
#segy_file = Path(r"D:\Projects\PES waal\segy\PES data Surv 2 Area 3\S2_A3_sn_23.06.2021\merged.sgy")
#segy_file = Path(r"D:\Projects\PES waal\segy\PES data Surv 2 Area 3\S2_A3_S_23.06.2021\merged.sgy")

bathy_path = Path(r"D:\Projects\PES waal\segy\PES data Surv 2 Area 3\Bathymetry\Bathymetry.nc")

edit = Segy_edit(segy_file)
edit.set_crs(28992)
data = np.array(edit.trace_data).T

manual_recording_delay = 4 #ms

shift = edit.recording_delay - manual_recording_delay
edit.vertical_trace_corrections(shift)
edit.recording_delay = manual_recording_delay
            
#%% data processing, this part needs some profiling and optimization
shot_point_interval_mean = edit.factor*edit.shot_point_interval.mean()

seis = Seismic(data, edit.sampling_rate, shot_point_interval_mean)
seis.time_squared_gain()
#seis.bandpass_filter(200, 20000) 

##?? where is the trace averager Tommer!?
## data derived decon would be cool here
## kirchhoff migration would be cool here
## multiple suppression would be cool here

## end data processing

edit.gain_type = 2 #t2 gain
edit.trace_data = seis.data.astype('int32')

## import bathymetry and verical correction with seabed pick
bath = xr.open_dataset(bathy_path)

# set correct crs? what is this again fro pes data?
#crosspoint_x = kaas
#%%crosspoint_y = kaas2


#edit.trace_selection(selection)
edit.write(segy_file.parent / 'processed.sgy')


'''
NEXT STEPS:
- find the bottleneck for processing (bandpass filter is probably the slowest and not really needed)
- seabedpick and and correct with bath (using a netcdf file, meaning the get_bathy function should be updated)
- get
- trace selection based on 4 crosspoints (set manually by plotting the cross section with the qgis exports of the tracks)
 find closest point, calculate distance of all points from THAT specific point. Select traces based on distance < x e.g 50 on both sides
'''	