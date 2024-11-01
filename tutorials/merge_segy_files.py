#%%
from deltaseis import merge_segys
from pathlib import Path

folder = Path(r'D:\Projects\PES waal\segy\PES data Surv 2 Area 3\S2_A3_ns_23.06.2021')
folders = [f for f in folder.iterdir() if f.suffix in ('.sgy', '.segy')]

#%%



s = merge_segys(folders, record_length=14)

#%%

outfile = Path("D:\Projects\PES waal\segy\PES data Surv 2 Area 3\S2_A3_ns_23.06.2021\merged.sgy")
s.write(outfile)
# %%
