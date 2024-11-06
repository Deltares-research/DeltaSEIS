#%%
from deltaseis import Segy_edit, Seismic
from pathlib import Path
import numpy as np
import xarray as xr

#segy_file = Path(r"D:\Projects\PES waal\segy\PES data Surv 2 Area 3\S2_A3_N_23.06.2021\merged.sgy")
segy_file = Path(r"D:\Projects\PES waal\segy\PES data Surv 2 Area 3\S2_A3_ns_23.06.2021\Site 3_20210623_090246_RAW_LF.sgy")
#segy_file = Path(r"D:\Projects\PES waal\segy\PES data Surv 2 Area 3\S2_A3_sn_23.06.2021\merged.sgy")
#segy_file = Path(r"D:\Projects\PES waal\segy\PES data Surv 2 Area 3\S2_A3_S_23.06.2021\merged.sgy")

bathy_path = Path(r"D:\Projects\PES waal\qgis\20210623PES waal S2_Area 3_RD_NAP_25cm.nc")

edit = Segy_edit(segy_file)
edit.set_crs(28992)

manual_recording_delay = 4 #ms

shift = edit.recording_delay - manual_recording_delay
edit.vertical_trace_corrections(shift)
edit.recording_delay = np.full(len(edit.recording_delay), manual_recording_delay).astype('int32')
            
#%% data processing, this part needs some profiling and optimization
shot_point_interval_mean = edit.factor*edit.shot_point_interval.mean()

seis = Seismic(data, edit.sampling_rate, shot_point_interval_mean)
#seis.time_squared_gain()
seis.agc_gain()
#seis.trace_averaging(3)
#seis.bandpass_filter(200, 20000) #wajoo duurt veel te lang jonge

## data derived decon would be cool here
## kirchhoff migration would be cool here
## multiple suppression would be cool here



## import bathymetry and verical correction with seabed pick
bathy = xr.open_dataset(bathy_path)

# set correct crs? what is this again fro pes data?
#crosspoint_x = kaas
#%%crosspoint_y = kaas2


#edit.trace_selection(selection)
edit.write(segy_file.parent / 'processed.sgy')


'''
NEXT STEPS:
- find the bottleneck for processing (bandpass filter is probably the slowest and not really needed)
- seabedpick and and correct with bathy (using a netcdf file, meaning the get_bathy function should be updated)
- trace selection based on 4 crosspoints (set manually by plotting the cross section with the qgis exports of the tracks)
- find closest point, calculate distance of all points from THAT specific point. Select traces based on distance < x e.g 50 on both sides
- conversion of all the raw data to make a clean repo would be an excellent idea (raw and processed)
'''	