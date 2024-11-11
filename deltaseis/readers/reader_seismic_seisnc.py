import xarray as xr
import numpy as np

seisnc_file = r"D:\Projects\PES waal\segy\segy_full_renamed\S2\S2_20210623_3_N_PRC_DEPTH_SEL.seisnc"
ds = xr.open_dataset(seisnc_file)


source_x = ds.source_x.values
source_y = ds.source_y.values

distances = np.sqrt(np.diff(source_x)**2 + np.diff(source_y)**2)
cumulative_distances = np.insert(np.cumsum(distances), 0, 0)

ds['distance'] = xr.DataArray(cumulative_distances, dims='shot')

print(ds.offset)
