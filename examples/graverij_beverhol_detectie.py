'''
Door het vergelijken van een meetsimulatie van (1) een intacte dijk en (2) een dijk met een bevergraafgang
kan worden onderzocht of aan de hand van seismische interferometrie met glasvezel (DAS) gedetecteerd kan
worden wanneer er een nieuwe graafgang ontstaat.

De simulaties zijn doorgerekend met SPECFEM2D en zijn reproduceerbaar via:
https://github.com/Deltares-research/specfem-examples/tree/main/DiggerDAS
(hier is ook een animatie van het golfveld door de dijk te vinden: hol_interferentie.avi)

De resultaten zijn geconverteerd naar seismische data met simulation2seismic.py van DeltaSEIS: 
https://github.com/Deltares-research/DeltaSEIS)

De onderstaande code reproduceert het figuur dat bij de interferometrische aanpak hoort
in Deltares rapportnummer 11210320-010-BGS-0001*

(1) de gesimuleerde seismische data voor een intacte dijk en een dijk met hol wordt ingeladen
(2) de afname van amplitude met afstand door geometrische verspreiding wordt gecorribeerd (T-square gain)
(3) het verschil tussen de intacte-dijk-data en de dijk-met-hol dat wordt berekend
(4) tot slot wordt een figuur met 3 subplots geproduceert met intact, hol en verschil, deze hebben
gelijke assen en amplitudes (door vmin en vmax aan te passen kan de nadruk verlegd worden op zwakkere of sterkere signalen)

*de annotatie en interpretatie in het rapport is apart ingetekend

'''
from deltaseis import Segy_edit, Seismic
from pathlib import Path
from importlib import resources
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import matplotlib.gridspec as gridspec
import numpy as np

#set color limits
vmin = -0.002
vmax = 0.002

# deltaseis data folder
data_folder = Path(resources.files("deltaseis.data"))

# paths to sgy data generated with simulation2seismic
file_intact = data_folder / "geen_hol_data.sgy"     # dijk zonder hol
file_hole =   data_folder / "hol_data.sgy"     # dijk met hol

# set output filename
outfile = file_intact.parent/"interferentie_beverhol"

# use Segy editor to retrieve data and sampling rate
intact = Segy_edit(file_intact)
data_intact = np.array(intact.trace_data).T
sampling_rate_intact = intact.sampling_rate

hole   = Segy_edit(file_hole)
data_hole = np.array(hole.trace_data).T
sampling_rate_hole = hole.sampling_rate

# because no trace coordinates were defined during the simulation, set this manually
receiver_spacing = 0.5

# instantiate seismic objects for preprocessing
seis_intact = Seismic(data_intact, sampling_rate_intact, receiver_spacing)
seis_hole = Seismic(data_hole, sampling_rate_hole, receiver_spacing)

# calculate agc
seis_intact_T2 = seis_intact.time_squared_gain(np.array(seis_intact.data))
seis_hole_T2 = seis_hole.time_squared_gain(np.array(seis_hole.data))

# compute the difference plot
seis_diff = seis_intact_T2 - seis_hole_T2

#%% Create a figure with 3 subplots in a row

# conversion factor for inches to cm
inch = 2.54

# create a gridspec with 4 columns; the fourth column is for the colorbar
fig = plt.figure(figsize=(14.6/inch, 7/inch))
gs = gridspec.GridSpec(1, 4, width_ratios=[1, 1, 1, 0.15], wspace=0.2)  # Adjust wspace for desired spacing

extent = [0, seis_intact_T2.shape[1]*receiver_spacing, seis_intact_T2.shape[0]/sampling_rate_intact, 0]

# plot the first image (data_intact)
ax1 = fig.add_subplot(gs[0])
im1 = ax1.imshow(seis_intact_T2, cmap='viridis', aspect='auto', extent=extent, vmin=vmin, vmax=vmax, interpolation='kaiser')
ax1.set_title('Data Intact', fontsize=10)
ax1.set_ylabel("Two-way travel time (s)", fontsize=7)
ax1.set_xlabel("Receiver position (m)", fontsize=7)
ax1.tick_params(axis='both', which='major', labelsize=7)

# plot the second image (data_hole)
ax2 = fig.add_subplot(gs[1])
im2 = ax2.imshow(seis_hole_T2, cmap='viridis', aspect='auto', extent=extent, vmin=vmin, vmax=vmax, interpolation='kaiser')
ax2.set_title('Data Hole', fontsize=10)
ax2.set_xlabel("Receiver position (m)", fontsize=7)
ax2.tick_params(axis='both', which='major', labelsize=7)
ax2.set_yticklabels([])  # Hide y-axis tick labels for the second subplot

# plot the difference image (data_diff)
ax3 = fig.add_subplot(gs[2])
im3 = ax3.imshow(seis_diff, cmap='viridis', aspect='auto', extent=extent, vmin=vmin, vmax=vmax, interpolation='kaiser')
ax3.set_title('Difference', fontsize=10)
ax3.set_xlabel("Receiver position (m)", fontsize=7)
ax3.tick_params(axis='both', which='major', labelsize=7)
ax3.set_yticklabels([])  # Hide y-axis tick labels for the third subplot

# create a colorbar in the fourth column of the grid
cbar_ax = fig.add_subplot(gs[3])
cbar = plt.colorbar(im3, cax=cbar_ax)
cbar.ax.tick_params(labelsize=7)  # Set colorbar tick label size
cbar.set_label('Amplitude', fontsize=7)

# use ScalarFormatter to display colorbar labels in scientific notation or shorter format
formatter = ScalarFormatter(useMathText=True)
formatter.set_powerlimits((-2, 2))  # Adjust the power limits for when scientific notation is used
cbar.formatter = formatter
cbar.update_ticks()  # Update the colorbar ticks to apply the formatter

# adjust the font size and alignment of the exponent
cbar.ax.yaxis.get_offset_text().set_fontsize(7)  # Match the font size of the exponent to the ticks
cbar.ax.yaxis.get_offset_text().set_position((1, 0))  # Adjust the position of the exponent (move slightly closer to the ticks)

# manually adjust the layout to provide extra padding
fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.15)

# show the plot
plt.show()

# save the figure with high resolution
fig.savefig(outfile, dpi=600)