

#%%
from importlib import resources
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from deltaseis import Segy_edit
from deltaseis.tools.navigation import compute_offset_points

# refer to data path
#data_folder = Path(resources.files('deltaseis.data'))                                       #part of data set for example (and testing)
data_folder = Path(r"D:\Projects\Universiteit Utrecht - Seismic Field School\deltaseis")     #full data seet for production

segy_file = data_folder / "sfs_p_shots_63-87.segy"

nav_file = data_folder / "sfs_GNSS_geophones.txt"
far_shots_file = data_folder / "sfs_GNSS_SHOTS_FAR_OFFSET.txt"

# set the instance for the bergen segy file
sfs = Segy_edit(segy_file)

# set the manual order for the segy file to get it correctly sorted
shot_order = (
            list(range(16, 24)) + 
            list(range(8, 16)) + 
            list(range(0, 8)) + 
            list(range(47, 23, -1)) +
            list(range(48, 72)) +
            list(range(95, 71, -1))
)

number_of_shots = len(np.unique(sfs.ffid)) #number of shots in file
manual_order = [shot + i * len(shot_order) for i in range(number_of_shots) for shot in shot_order] #reorder for all shots in file

# apply the sorting
sfs.sort_traces(manual_order=manual_order)
sfs.trace_sequence.sort()

# add coordinates and elevation to the trace headers
sfs.set_crs(28992)          # set to RD
sfs.set_crs_vertical(7415)  # set to NAP

sfs.set_scalar(-1000)            # set to millimeter precision
sfs.set_scalar_vertical(-1000)   # set to millimeter precision

# get the coordinates from the navigation file
coordinates = np.genfromtxt(nav_file, delimiter=',', skip_header=1, usecols=range(4))
receiver_x = coordinates[:, 1]
receiver_y = coordinates[:, 2]
receiver_z = coordinates[:, 3]
source_z = receiver_z[::2]

# compute 2 m offsets for the shot positions
source_x, source_y = compute_offset_points(receiver_x, receiver_y, crossline_distance=2, inline_distance=0)
source_x = source_x[::2]
source_y = source_y[::2]

# now add the coordinates for all traces
number_of_lines = 2
ffid_start_of_line = sfs.ffid[0]
receiver_x_header = np.tile(receiver_x, number_of_lines*number_of_shots)
receiver_y_header = np.tile(receiver_y, number_of_lines*number_of_shots)
receiver_z_header = np.tile(receiver_z, number_of_lines*number_of_shots)

# set the ffids so that the increment by one (this can be off due to discarded shots)
sfs.ffid = np.repeat(np.arange(sfs.ffid[0], sfs.ffid[0] + len(np.unique(sfs.ffid))), sfs.number_of_channels)

# add shot locations, 
source_x_header = np.take(source_x, sfs.ffid - ffid_start_of_line)
source_y_header = np.take(source_y, sfs.ffid - ffid_start_of_line)
source_z_header = np.take(source_z, sfs.ffid - ffid_start_of_line)
source_z_header = np.take(source_z, sfs.ffid - ffid_start_of_line)

sfs.add_coordinates(source_x_header, source_y_header, receiver_x_header, receiver_y_header, source_z_header, receiver_z_header)

#%%trace header math
sfs.ffid = sfs.ffid + 1000 # set first to 1000 to distinguish from original shotpoints (industry standard)
offsets = np.sqrt((source_x_header - receiver_x_header)**2 + (source_y_header - receiver_y_header)**2) # set offset to 2 m

cdpx = (source_x_header + receiver_x_header) / 2
cdpy = (source_y_header + receiver_y_header) / 2

#cdp number calculation (unique for this survey setup)
cdps = np.tile(np.arange(1, 49), sfs.number_of_shots * 2) + np.repeat(np.arange(sfs.number_of_shots), 96)

#%%QC plot
plt.figure(figsize=(8, 10))

# Plot receiver and shot points
plt.plot(receiver_x, receiver_y, 'k.', label='receiver points', markersize=2)
plt.plot(source_x, source_y, 'r.', label='shot points', markersize=2)
plt.plot(cdpx, cdpy, 'b.', label='cdp points', markersize=2)
plt.grid()
plt.legend()
plt.title(f'Receiver and Shot Points (EPSG:{sfs.crs.to_epsg()})')
plt.xlabel('X Coordinate')
plt.ylabel('Y Coordinate')
plt.axis('equal')
plt.savefig(segy_file.with_suffix('.png'), dpi=600)

# write as segy file header format
sfs.offsets = np.round(offsets / sfs.factor).astype(np.int32)
sfs.cdpx = np.round(cdpx / sfs.factor).astype(np.int32)
sfs.cdpy = np.round(cdpy / sfs.factor).astype(np.int32)
sfs.cdps = np.round(cdps).astype(np.int32)

#%% set selection for traces to write

# P hammer survey
select_lf_vertical_component = [i for i, channel in enumerate(sfs.channel_numbers) if 1 <= channel <= 48]
sfs.select_traces(select_lf_vertical_component)
segy_outfile = segy_file.with_stem(segy_file.stem + "_V4_5HZ_SORTED")
sfs.write(segy_outfile)

select_vertical_component =   [i for i, channel in enumerate(sfs.channel_numbers) if 49 <= channel <= 96]
sfs.select_traces(select_vertical_component)
segy_outfile = segy_file.with_stem(segy_file.stem + "_V10HZ_SORTED")
sfs.write(segy_outfile)