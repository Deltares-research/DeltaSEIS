'''
Read pre-stack from legacy data shot north of the Dutch Wadden Islands
'''
#%%
from pathlib import Path
import numpy as np
from deltaseis import Segy_edit
import matplotlib.pyplot as plt


segy_file = Path(r"D:\else9854_3_1777.segy")
number_of_channels = 12
channel = 1


#initiate object
prestack = Segy_edit(segy_file)
prestack.set_endian('little')

#add channel number to (near) trace gather output file
segy_outfile = segy_file.with_stem(segy_file.stem + f"_NTG{channel}")

#extract single channel gather
prestack.extract_near_trace_gather(number_of_channels=number_of_channels, channel=channel)
#write to new seg-y
prestack.write(segy_outfile)


