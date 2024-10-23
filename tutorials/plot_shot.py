#%%


from deltaseis import Segy_edit, Seismic
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import numpy as np


file_hole   = Path(r"\\wsl.localhost\Ubuntu\home\rlnd\projects\DiggerDAS\segy\hole4_BXX.sgy")

outfile = file_hole.parent/"shot_plot_BXX"

#use Segy editor to retrieve data and sampling rate
hole   = Segy_edit(file_hole)

#cut the record lenght to the relevant interval (ms)
hole.set_record_length(20)

data_hole = np.array(hole.trace_data).T
fs_hole = hole.sampling_rate

#because no trace coordinates were defined during the simulation, set this manually
dx = 0.5

#instantiate seismic sobjects for preprocessing
seis_hole = Seismic(data_hole, fs_hole, dx)

#calculate agc
seis_hole_agc = seis_hole.time_squared_gain(np.array(seis_hole.data))

#%% Create a figure with 3 subplots in a row

# Conversion factor for inches to cm
inch = 2.54

# Create figure and subplots (1 row, 2 columns)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14.6/inch, 12/inch), sharey=True)

# Define color limits and extent
vmin, vmax = -1e-4, 1e-4
extent = [0, seis_hole_agc.shape[1] * dx, seis_hole_agc.shape[0] / fs_hole, 0]

# Plot the data on ax1
im1 = ax1.imshow(seis_hole_agc, cmap='viridis', aspect='auto', extent=extent, vmin=vmin, vmax=vmax)
ax1.set(title='Graafgeluid data', xlabel='Receiver position (m)', ylabel='Two-way travel time (s)')
ax1.tick_params(axis='both', labelsize=7)

# Plot additional data on ax2 (replace with your actual data)
im2 = ax2.imshow(seis_hole_agc, cmap='viridis', aspect='auto', extent=extent, vmin=vmin, vmax=vmax)  # Example with different colormap
ax2.set(title='Interpretatie', xlabel='Receiver position (m)')
ax2.tick_params(axis='both', labelsize=7)

# Create a colorbar for ax1
cbar = plt.colorbar(im2, ax=ax2, fraction=0.05, pad=0.04)
cbar.ax.tick_params(labelsize=8)
cbar.set_label('Amplitude', fontsize=8)

# Format colorbar labels
formatter = ScalarFormatter(useMathText=True)
formatter.set_powerlimits((-2, 2))
cbar.formatter = formatter
cbar.update_ticks()
cbar.ax.yaxis.get_offset_text().set_fontsize(8)
cbar.ax.yaxis.get_offset_text().set_position((3, 0))

# Adjust layout to prevent overlap
plt.tight_layout()

# Show and save the plot
plt.show()
fig.savefig(outfile, dpi=600)


