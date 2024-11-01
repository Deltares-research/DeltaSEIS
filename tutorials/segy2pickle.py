'''
This script show how to export the trace data (nothing else) of segy files to a pickles file. 
The exported pickles contain matrices with data containing the amplitude data for all traces.
This data can be used to generate input pictures for AI deep learning algorithms.

TODO: this script should become part of export_seismic_pickle

'''

from deltaseis import Segy_edit, Seismic
import pickle
import numpy as np
from pathlib import Path

folder = Path(r'D:\Projects\PES waal\segy\segy_full_renamed\S2')        
segy_files = txt_files = folder.glob('*.segy')
segy_files = [Path(path) for path in segy_files]

#%%

def export_pickle(infile):

    outfile = infile.with_suffix(".npy")
    print(outfile)
    s = Segy_edit(infile)
    print(np.array(s.trace_data).shape)

    seismic = Seismic(s.trace_data, s.sampling_rate, s.spi.mean())
    print(s.sampling_rate)
    print(s.spi)
    print(s.spi.mean())

    with open(outfile, 'wb') as f:
        np.save(f, np.array(seismic.data))


for segy_file in segy_files:
    export_pickle(segy_file)
# %%

