#%%
from deltaseis import Segy_edit, Seismic
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import matplotlib.gridspec as gridspec
import numpy as np


file_intact = Path(r"\\wsl.localhost\Ubuntu\home\rlnd\projects\DiggerDAS\base3.sgy")
file_hole   = Path(r"\\wsl.localhost\Ubuntu\home\rlnd\projects\DiggerDAS\hole3.sgy")

outfile = file_intact.parent/"difference_plot"

#use Segy editor to retrieve data and sampling rate
intact = Segy_edit(file_intact)
data_intact = np.array(intact.trace_data).T
fs_intact = intact.sampling_rate

hole   = Segy_edit(file_hole)
data_hole = np.array(hole.trace_data).T
fs_hole = hole.sampling_rate

#because no trace coordinates were defined during the simulation, set this manually
dx = 0.5

#instantiate seismic objects for preprocessing
seis_intact = Seismic(data_intact, fs_intact, dx)
seis_hole = Seismic(data_hole, fs_hole, dx)

#calculate agc
seis_intact_agc = seis_intact.time_squared_gain(np.array(seis_intact.data))
seis_hole_agc = seis_hole.time_squared_gain(np.array(seis_hole.data))

# Compute the difference plot
seis_diff = seis_intact_agc - seis_hole_agc

#%% Create a figure with 3 subplots in a row

import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import matplotlib.gridspec as gridspec

# Conversion factor for inches to cm
inch = 2.54

# Create a gridspec with 4 columns; the fourth column is for the colorbar
fig = plt.figure(figsize=(14.6/inch, 7/inch))
gs = gridspec.GridSpec(1, 4, width_ratios=[1, 1, 1, 0.15], wspace=0.2)  # Adjust wspace for desired spacing

vmin = -0.01
vmax = 0.01

extent = [0, seis_intact_agc.shape[1]*dx, seis_intact_agc.shape[0]/fs_intact, 0]

# Plot the first image (data_intact)
ax1 = fig.add_subplot(gs[0])
im1 = ax1.imshow(seis_intact_agc, cmap='viridis', aspect='auto', extent=extent, vmin=vmin, vmax=vmax, interpolation='kaiser')
ax1.set_title('Data Intact', fontsize=10)
ax1.set_ylabel("Two-way travel time (s)", fontsize=7)
ax1.set_xlabel("Receiver position (m)", fontsize=7)
ax1.tick_params(axis='both', which='major', labelsize=7)

# Plot the second image (data_hole)
ax2 = fig.add_subplot(gs[1])
im2 = ax2.imshow(seis_hole_agc, cmap='viridis', aspect='auto', extent=extent, vmin=vmin, vmax=vmax, interpolation='kaiser')
ax2.set_title('Data Hole', fontsize=10)
ax2.set_xlabel("Receiver position (m)", fontsize=7)
ax2.tick_params(axis='both', which='major', labelsize=7)
ax2.set_yticklabels([])  # Hide y-axis tick labels for the second subplot

# Plot the difference image (data_diff)
ax3 = fig.add_subplot(gs[2])
im3 = ax3.imshow(seis_diff, cmap='viridis', aspect='auto', extent=extent, vmin=vmin, vmax=vmax, interpolation='kaiser')
ax3.set_title('Difference', fontsize=10)
ax3.set_xlabel("Receiver position (m)", fontsize=7)
ax3.tick_params(axis='both', which='major', labelsize=7)
ax3.set_yticklabels([])  # Hide y-axis tick labels for the third subplot

# Create a colorbar in the fourth column of the grid
cbar_ax = fig.add_subplot(gs[3])
cbar = plt.colorbar(im3, cax=cbar_ax)
cbar.ax.tick_params(labelsize=7)  # Set colorbar tick label size
cbar.set_label('Amplitude', fontsize=7)

# Use ScalarFormatter to display colorbar labels in scientific notation or shorter format
formatter = ScalarFormatter(useMathText=True)
formatter.set_powerlimits((-2, 2))  # Adjust the power limits for when scientific notation is used
cbar.formatter = formatter
cbar.update_ticks()  # Update the colorbar ticks to apply the formatter

# Adjust the font size and alignment of the exponent
cbar.ax.yaxis.get_offset_text().set_fontsize(7)  # Match the font size of the exponent to the ticks
cbar.ax.yaxis.get_offset_text().set_position((1, 0))  # Adjust the position of the exponent (move slightly closer to the ticks)

# Manually adjust the layout to provide extra padding
fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.15)

# Show the plot
plt.show()

# Save the figure with high resolution
fig.savefig(outfile, dpi=600)
