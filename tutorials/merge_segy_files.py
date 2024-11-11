#%%
from deltaseis import merge_segys
from pathlib import Path
import matplotlib.pyplot as plt

folder = Path(r"D:\Projects\PES waal\raw data\PES data Surv 2 Area 2\S2_A2_sn_22.06.2021")


folders = [f for f in folder.iterdir() if f.suffix in ('.sgy', '.segy')]

s = merge_segys(folders, make_plot=False, record_length=14)

outfile = folder / 'merged.sgy'
s.write(outfile)

