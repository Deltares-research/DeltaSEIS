'''
#first: clean up data folder on D and set proper linking scheme
#second: find order required

#(1) all .dat files opened with geo2view and saved as a single seg-y file: sfs2024.sgy 
#(2) load with deltaseis Segyeditor
(3) sort with manual order (where did you save that) - on teams

(5) export as different segy
(6) add coordinates, base in gnss data (make add_coordinates method in Segy_edit)
- channel postions x, y and z, offset by 10? cm for '2nd' line, which was which
- shot potition x, y and Z, lateral offset by 2m, was there an inline offset as well? (or alligned with rec positions?)
    - if so compute offset points has to be adapted to take into account the inline offset
- figure out which bytes these should be in the segy file and properly handle both scalars (set to -1000 for millimeter precision)
(4) select to get different cathegories (P-reflection and rayleighwaves, SH reflection and love waves?, diffraction (from P-wave traces but based on shotpoints in logsheet))
(7) optional apply an f-k filter (to split p-waves from rayleigh waves and SH-waves from love waves)
(8) ask edwin to help with the surface wave processing and refraction seismic.
(9) due to the near co-located data sets, the data can serve as a use case for multi-physics inverse physics informed neural IMPINN

'''

#%%
from importlib import resources
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from deltaseis import Segy_edit
from deltaseis.tools.navigation import compute_offset_points

#%% refer to data path
data_folder = Path(resources.files('deltaseis.data'))                                           #part of data set for example (and testing)
#data_folder = Path(r"D:\Projects\Universiteit Utrecht - Seismic Field School\deltaseis")       #full data seet for production

segy_file = data_folder / "sfs.sgy"
segy_outfile = segy_file.with_stem(segy_file.stem + "_SORTED")

nav_file = data_folder / "sfs_GNSS_geophones.txt"
far_shots_file = data_folder / "sfs_GNSS_SHOTS_FAR_OFFSET.txt"


#%% set the instance for the bergen segy file
sfs = Segy_edit(segy_file)
sfs.set_endian('little')


#%% set the manual order for the segy file to get it correctly sorted
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

#%% add coordinates and elevation to the trace headers

sfs.set_crs(28992)          # set to RD
sfs.set_crs_vertical(7415)  # set to NAP

sfs.set_output_scalar(-1000)            # set to millimeter precision
sfs.set_vertical_output_scalar(-1000)   # set to millimeter precision

# get the coordinates from the navigation file
coordinates = np.genfromtxt(nav_file, delimiter=',', skip_header=1, usecols=range(4))
rx = coordinates[:, 1]
ry = coordinates[:, 2]
rz = coordinates[:, 3]
z = rz

# compute 2 m offsets for the shot positions
x_offset, y_offset = compute_offset_points(rx, ry, crossline_distance=2, inline_distance=0)
#OFFSETS ARE NOT CORRECT YET, THERE IS AN INLINE EFFECT AND ONE POINT MISSES

# Quick QC plot
plt.figure(figsize=(7, 8))

# Plot receiver and shot points
plt.plot(rx, ry, 'k.', label='receiver points')
plt.plot(x_offset, y_offset, 'r.', label='shot points')
plt.grid()
plt.legend()
plt.title(f'Receiver and Shot Points (EPSG:{sfs.crs.to_epsg()})')
plt.xlabel('X Coordinate')
plt.ylabel('Y Coordinate')
plt.axis('equal')
plt.show()

#%% now add the coordinates to the segy trace headers
number_of_lines = 2
ffid_start_of_line = 6
rx_header = np.tile(rx, number_of_lines*number_of_shots)
ry_header = np.tile(ry, number_of_lines*number_of_shots)
rz_header = np.tile(rz, number_of_lines*number_of_shots)

x_header = np.take(x_offset, sfs.ffid - ffid_start_of_line)
y_header = np.take(y_offset, sfs.ffid - ffid_start_of_line)
z_header = np.take(z, sfs.ffid - ffid_start_of_line)

sfs.add_coordinates(x_header, y_header, z_header, rx_header, ry_header, rz_header)


sfs.set_channels(channel_list)
sfs.set_offsets(offsets_list)
sfs.set_cdps(cdp, cdp_x, cdp_y)


#
#%%
#sfs.select


#%% other header math might be required, e.g reset shotpoints? calc offsets, cdps (look in radexpro manual)

sfs.write(segy_outfile)



# detail: use compute offset point to get 10 cm offset for 2nd line
# sfs select traces to split different lines
# coordinate scalar should be made easier (in/out etc, just use set_ ?)
# also make a more regular shot interval (exact staright line, exact 2m spacing, exact 2 m offset to compare)