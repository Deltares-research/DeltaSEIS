'''
Door het vergelijken van een meetsimulatie van (1) een intacte dijk en (2) een dijk met een bevergraafgang
kan worden onderzocht of aan de hand van seismische interferometrie met glasvezel (DAS) gedetecteerd kan
worden wanneer er een nieuwe graafgang ontstaat.

De simulaties zijn doorgerekend met SPECFEM2D en zijn reproduceerbaar via:
https://github.com/Deltares-research/specfem-examples/tree/main/DiggerDAS
(hier is ook een animatie van het graafgeluid vanuit de dijk te vinden: bever_geluid.avi)

De resultaten zijn geconverteerd naar seismische data met simulation2seismic.py van DeltaSEIS: 
https://github.com/Deltares-research/DeltaSEIS)

De onderstaande code reproduceert het figuur dat bij de graafgeluid aanpak hoort
in Deltares rapportnummer 11210320-010-BGS-0001*

(1) de gesimuleerde seismische data graafgeluid uit een hol in de dijk wordt geladen
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
import numpy as np

#set color limits
vmin = -1e-4
vmax = 1e-4

# deltaseis data folder
data_folder = Path(resources.files('deltaseis.data'))

# paths to sgy data generated with simulation2seismic
file_hole = data_folder / 'bever_geluid_data.sgy'

# set output filename
outfile = file_hole.parent/"graafgeluid"

#use Segy editor to retrieve data and sampling rate
hole   = Segy_edit(file_hole)

#cut the record lenght to the relevant interval (ms)
hole.set_record_length(20)

data_hole = np.array(hole.trace_data).T
sampling_rate_hole = hole.sampling_rate

#because no trace coordinates were defined during the simulation, set this manually
dx = 0.5

#instantiate seismic sobjects for preprocessing
seis_hole = Seismic(data_hole, sampling_rate_hole, dx)

#calculate agc
seis_hole_agc = seis_hole.time_squared_gain(np.array(seis_hole.data))

#%% Create a figure with 3 subplots in a row

# Conversion factor for inches to cm
inch = 2.54

# Create figure and subplots (1 row, 2 columns)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14.6/inch, 12/inch), sharey=True)

# Define extent
extent = [0, seis_hole_agc.shape[1] * dx, seis_hole_agc.shape[0] / sampling_rate_hole, 0]

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


