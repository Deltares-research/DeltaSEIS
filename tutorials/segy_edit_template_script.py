#%%

from pathlib import Path
import numpy as np
from deltaseis import Segy_edit

segy_file = Path(r"P:\430-tgg-data\Tutorial\bergen0111a.sgy")

#set the path of the grid, from which you would like to extract the horizon
grid_path = Path(r"P:\430-tgg-data\Tutorial\bathy19msl_dis_buff12.asc")

#set the segy file path
segy_outfile = segy_file.with_stem(segy_file.stem + "_TUTORIAL")

# set the instance for the bergen segy file
bergen = Segy_edit(segy_file)

bergen.set_record_length(60)
bergen.set_endian('little')
bergen.fix_navigation_median()
bergen.set_input_scalar(-100)

bergen.transform_coordinates(4326, 32631)

bergen.read_grid(grid_path, 32631, 4326, horizon_name='bathy')
bergen.get_seabed_pick(10, 100, 9, 3, truncate=10) #see if you can find default for part of the input parameters to simplify
bergen.plot(save_plot=True, clip=0.1)

#filter heave from seabed pick, calculate the difference with the original and apply as vertical corrections to the data
bergen.filter_horizon_savgol('seabed_pick', 'seabed_pick_savgol', 501, 4)
bergen.calculate_difference_horizon('seabed_pick_savgol', 'seabed_pick', difference_horizon_name = 'heave')
bergen.vertical_trace_corrections(bergen.heave)
bergen.plot(save_plot=True, clip=0.1)

#filter the tide grid so that it becomes a smooth curve for tide and sensor depth correction
bergen.calculate_difference_horizon('bathy', 'seabed_pick', difference_horizon_name ='tide_raw')
bergen.filter_horizon_savgol('tide_raw', 'tide', 2001, 2)
bergen.vertical_trace_corrections(bergen.tide)
bergen.plot(save_plot=True, clip=0.1)

#miscelaneous tace header operations
bergen.set_trace_number_in_field_record()
bergen.renumber_shotpoints()
bergen.copy_source_coordinates_to_group()
bergen.set_output_scalar(-1000)

#write to segy
bergen.write(segy_outfile)


#print all atrribute of 'bergen'
# print(dir(bergen))



