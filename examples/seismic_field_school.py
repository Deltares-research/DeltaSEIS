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
import deltaseis
from deltaseis import Segy_edit
from deltaseis.tools.navigation import compute_offset_points

#%% refer to deltaseis data path
data_folder = Path(resources.files('deltaseis.data'))
segy_file = data_folder / "sfs.sgy"
segy_outfile = segy_file.with_stem(segy_file.stem + "_SORTED")

nav_file = data_folder / "sfs_GNSS_geophones.txt"
far_shots_file = data_folder / "sfs_GNSS_SHOTS_FAR_OFFSET.txt"


#%% set the instance for the bergen segy file
sfs = Segy_edit(segy_file)

sfs.set_endian('little')
sfs.set_crs(28992)

#%% set the manual order for the segy file to get it correctly sorted
shot_order = (
            list(range(16, 24)) + 
            list(range(8, 16)) + 
            list(range(0, 8)) + 
            list(range(47, 23, -1)) +
            list(range(48, 72)) +
            list(range(95, 71, -1))
)

number_of_shots = sfs.ffid.max()

manual_order = []
for i in range(number_of_shots):
    manual_order.extend([shot + i * len(shot_order) for shot in shot_order])

#%% methods to we written
sfs.sort_traces(manual_order=manual_order)

#%%

data = np.genfromtxt(nav_file, delimiter=',', skip_header=1, filling_values=np.nan)
#CONTINUE HERE:
#load x, y and z coordiantes in the apporpriate segy headers
elevation_scalar_byte = 69
coordinate_scalar_byte = 71

source_x_byte = 73
source_y_byte = 77
source_z_byte = 57

receiver_x_byte = 81
receiver_y_byte = 85
receiver_z_byte = 53



#%% Get the offset points

# Assuming original_points is the data loaded from nav_file
original_points = np.loadtxt(nav_file, delimiter=',', skiprows=1)

#offset_points = compute_offset_points(original_points, offset_distance)
#sfs.add_coordinates(nav_file)

#sfs.set_input_scalar(-1000) # set to millimeter precision
#%%
#sfs.select


#%%

sfs.write(segy_outfile)