# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 12:53:28 2025

@author: vermaas
"""

from pathlib import Path
import deltaseis
import numpy as np
import os
from matplotlib import pyplot as plt

indir = Path(r'D:\Projects\DIS_Offshore\marsdiep2015_xstar')
outdir = Path(r'D:\Projects\DIS_Offshore\marsdiep2015_xstar\output')

#%% loop
segy_files = os.listdir(indir)
for i, segy_file in enumerate(segy_files):
    print(f"{i+1}/{len(segy_files)}: {segy_file}")

    f = Path(indir,segy_file)
    s = deltaseis.Segy_edit(f) 
    s.xstar_split()
    s.set_record_length(30)
    s.set_crs(4326)

    #s.set_scalar(-100)   
    s.transform_coordinates(25831) #to RD (to WGS84 UTM31N = 32631)
    s.renumber_shotpoints()
    
    data = np.array(s.trace_data).T #transpose zodat tijd op 'y' (rijen) trace op 'x' (kolommen)
    seis = deltaseis.Seismic(data,12.5,0.5)
    # seis.time_power_gain()
    seis.trace_averaging()
    data2 = seis.data.T
    
    s.trace_data = (32767/np.max(data2))*data2 #normalize to fit in int16 format
    
    # s.plot()
    # plt.savefig(Path(outdir,'pngs',f"{f.stem}.png"))
    # plt.close()
    
    s.set_endian('big')
    
    s.write(Path(outdir,f"{f.stem}_edit{f.suffix}y"))
