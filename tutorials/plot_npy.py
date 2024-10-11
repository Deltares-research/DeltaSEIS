

import numpy as np
import matplotlib.pyplot as plt

#numpy_file = r"D:\Projects\PES waal\segy\segy_full_renamed\compare marios\S1_1_N_RAW_LF_loc_stack.npy"
#numpy_file = r"D:\Projects\PES waal\segy\segy_full_renamed\compare marios\S1_1_N_RAW_LF_stack_agg.npy"
#numpy_file = r"D:\Projects\PES waal\segy\segy_full_renamed\compare marios\S1_1_N_RAW_LF_stack_agg_trim.npy"
numpy_file = r"D:\Projects\PES waal\segy\segy_full_renamed\compare marios\S1_20210531_1_N.npy"

vmin = -0.01
vmax = 0.01

outfile = numpy_file.replace(".npy", ".png")

data = np.load(numpy_file)
data = data / np.max(data)
print(data.shape)

print(outfile)

plt.imshow(data.T, vmin=vmin, vmax=vmax, cmap="RdBu", aspect='auto')

plt.savefig(outfile)
plt.show()
