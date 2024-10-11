'''
This script show how to export the trace data (nothing else) of segy files to a pickles file. 
The exported pickles contain matrices with data containing the amplitude data for all traces.
This data can be used to generate input pictures for AI deep learning algorithms.

TODO: this script should become part of export_seismic_pickle

'''

from deltaseis import Segy_edit
import pickle
import numpy as np
from pathlib import Path

folder = Path(r'D:\Projects\PES waal\segy\segy_full_renamed')        
segy_files = txt_files = folder.glob('*.sgy')
segy_files = [Path(path) for path in segy_files]

#%%

def export_pickle(infile):

    outfile = infile.with_suffix(".npy")
    print(outfile)
    seismic = Segy_edit(infile)
    print(np.array(seismic.trace_data).shape)

    with open(outfile, 'wb') as f:
        np.save(f, np.array(seismic.trace_data))


for segy_file in segy_files:
    export_pickle(segy_file)
# %%

