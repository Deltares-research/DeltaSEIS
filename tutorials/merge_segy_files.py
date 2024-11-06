#%%
from deltaseis import merge_segys
from pathlib import Path
import matplotlib.pyplot as plt

#folder = Path(r"D:\Projects\PES waal\segy\PES data Surv 2 Area 3\S2_A3_N_23.06.2021")
#folder = Path(r"D:\Projects\PES waal\segy\PES data Surv 2 Area 3\S2_A3_ns_23.06.2021")
#folder = Path(r"D:\Projects\PES waal\segy\PES data Surv 2 Area 3\S2_A3_sn_23.06.2021")
folder = Path(r"D:\Projects\PES waal\segy\PES data Surv 2 Area 3\S2_A3_S_23.06.2021")


folders = [f for f in folder.iterdir() if f.suffix in ('.sgy', '.segy')]

s = merge_segys(folders, make_plot=False, record_length=14)

outfile = folder / 'merged.sgy'
s.write(outfile)

