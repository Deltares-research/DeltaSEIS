# -*- coding: utf-8 -*-
"""
The scripts reads pre-stack seismic data e.g. as recorded by geodes
The user can select .sg2 or .dat file via a pop-up file dialog. 

The scripts reads the data into a deltaseis-compatable format for further use.

For now it doesn't. I want to use this to edit seg2, so i don't do the step towards
gmpt compatable, but i stay in a copy of sg2 and output to segy with fully populated headers

edit_sg2 does not want to read sg2/dat files as far as i know, as a workaround i now
used geo2view (http://geo2x.com/wp/seismic-viewer-geo2view/) to open de sg2 and 
save as segy from the file menu.

"""
def header_geometry(ntraces, dx, shot_number, spi):

    
    import numpy as np
        
    geometry_dict = {}
    
    geometry_dict['zeros']          = [0]*ntraces
    geometry_dict['spn']            = [shot_number]*ntraces
    geometry_dict['ffid']           = [spn + 1000 for spn in geometry_dict['spn']]
    geometry_dict['source_number']  = [shot_number]*ntraces
    geometry_dict['channels']       = [channel + 1 for channel in range(ntraces)] 
    geometry_dict['x_rec']          = [spi + (i-1)*dx for i in geometry_dict['channels']]
    geometry_dict['y_rec']          = geometry_dict['zeros']
    geometry_dict['x_src']          = [(shot_number - 1)*spi]*ntraces
    geometry_dict['y_src']          = geometry_dict['zeros']
    geometry_dict['offset']         = [np.sqrt((y_rec - y_src)**2 + (x_rec - x_src)**2) for y_rec, y_src, x_rec, x_src in zip(geometry_dict['y_rec'],  geometry_dict['y_src'],  geometry_dict['x_rec'],  geometry_dict['x_src'])]
    geometry_dict['x_cdp']          = [0.5*(x_rec - x_src) + x_src for x_rec, x_src in zip(geometry_dict['x_rec'],geometry_dict['x_src'])]
    geometry_dict['y_cdp']          = geometry_dict['zeros']
    geometry_dict['cdp']            = [int(x_cdp/(0.5*dx)) for x_cdp in geometry_dict['x_cdp']]
    
    c0 = int( geometry_dict['x_rec'][0]/dx) #trick to let the first CDP number start at 1
    
    geometry_dict['cdp']        = [cdp - c0 + 1 for cdp in geometry_dict['cdp']]  #set first cdp number to 1   
    
    return geometry_dict


def edit_sg2(sg2_file, dx, shot_number, spi):
    import segyio
      
    
    with segyio.open(sg2_file, mode= 'r', endian='big', ignore_geometry=True) as src:
        spec = segyio.tools.metadata(src)
        
        ntraces = spec.tracecount
        g = header_geometry(ntraces, dx, shot_number, spi)
        
        sg2_file_edited = sg2_file.replace(".sgy", "_edit.sgy")
   
        with segyio.create(sg2_file_edited, spec) as dst:
            dst.text[0] = src.text[0]
            dst.bin = src.bin
            dst.header[:] = src.header[:]
            dst.trace = src.trace
            
            for i in range(ntraces):
                dst.header[i] = src.header[i]
                dst.trace[i] = src.trace[i]
                #populate trace headers
                dst.header[i][9]      = g['ffid'][i]                  #FFID - original field record number
                dst.header[i][17]     = g['spn'][i]                   #energy source point number
                dst.header[i][21]     = g['cdp'][i]                   #cdp ensemble number
                dst.header[i][71]     = src.header[i][71]             #scaler to all coordinates
                dst.header[i][73]     = int(100*g['x_src'][i])        #source x coordinate
                dst.header[i][77]     = int(100*g['y_src'][i])        #source y coordinate
                dst.header[i][81]     = int(100*g['x_rec'][i])        #group x coordinate (receiver)
                dst.header[i][85]     = int(100*g['y_rec'][i])        #group y coordinate (receiver)
                dst.header[i][109]    = src.header[i][109]            #recording delay
                dst.header[i][115]    = len(spec.samples)             #number of samples in this trace
                dst.header[i][117]    = src.header[i][117]            #sample interval in ms for this trace
                dst.header[i][181]    = int(100*g['x_cdp'][i])        #CDP X
                dst.header[i][185]    = int(100*g['y_cdp'][i])        #CDP Y
                
 

#%%------------------------------apply to files--------------------------
from pathlib import Path
import tkinter.filedialog, tkinter.simpledialog
import tkinter as tk
#root = tk.Tk().withdraw()
               
dx = 2 #m
spi = 2 #m
 
 
#sg2_files = tk.filedialog.askopenfilenames(filetypes = (("seg-y files",("*.sgy","*.segy")), ("all files","*.*")))
sg2_files = [r"C:\Users\nieboer\OneDrive - Stichting Deltares\Desktop\6.segy"]

for i, sg2_file in enumerate(sg2_files):
    print("{}/{}: {}".format(i+1, len(sg2_files), sg2_file))
    shot_number = int(Path(sg2_file).stem)
    edit_sg2(sg2_file, dx, shot_number, spi)
    











# %%
