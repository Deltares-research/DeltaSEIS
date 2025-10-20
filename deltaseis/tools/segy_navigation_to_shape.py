# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 10:42:19 2019

Creates point and/or line track shapefiles using SEGY files from specified folder as input. 
Start- and end_filename between which the line name is found is important to set as well as epsg for the shapefile prj file

BUG: if anaconda modules are not installed correctly, the PRJ file holding the coordinate system will be blank. Always check this.
ESSENTIAL QC: check if the prj file is correctly populated.

goto: https://spatialreference.org/ search your coordinate system (e.g. by EPGS code)
after creation of the shapefiles, the prj should hold a text identical to the WKT in the grey box half way the page (vital QC check). You can also copy this to the blank PRJ file.

todo: put into GUI (make tab?)
todo: remove outliers

2024: this needs refactor - make simpler, add entrypoint where function is run (segy load parallely?, but in the end one shapefile) but also loadable as module 

@author: r.nieboer
"""

#Import modules
import time
import tkinter as tk

# import re
import tkinter.filedialog
import tkinter.simpledialog
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import segyio
from shapely.geometry import LineString, Point

root = tk.Tk().withdraw()

#input files with dialog
segy_files= tk.filedialog.askopenfilenames(filetypes = (("SEG-Y files","*.sgy;*.segy"),("SEG files","*.seg"),("all files","*.*")))

#or manual input list of files
#segy input folder
# segy_files=["D:/Projects/DIS/interpretations/zeeland2019/Seismiek_S2_S3_test/test/beton_2001_lijn01.sgy"]

#%%
outfile="HVS" #outfile basename for shapefiles
trackpoints=False #output track point shape file
tracklines=True #output track line shape file
point_type="Trace number" #e.g. shotpoint or CDP depending on which point_byte is used
point_byte=1 #check textual header. Generally: shotpoint=17, CDP=21
easting_byte=73 #check textual header. Genearlly:  for m3ultichannel,CDP_X=181 (or 81) for SBP src_x=rec_x so source x=73 is mostly used 
northing_byte=77 #check textual header. Gnerally: for multichannel, CDP=185 (or 85) for SBP src_y=rec_y so source y=77 is mostly used
scalar_byte=71
pos_in_arcsec=False #true if postions are in arc seconds, for easting/northings use False
replace_shotpoint=True #if non unique shotsource available in any of the bytes, this adds a number sequence from 1 to #tracecount

scalar_corr=1 #multiplier if the coordinate scalar in the segy traceheader is incorrect, default should be 1

# goto: https://spatialreference.org/ search your coordinate system (e.g. by EPGS code)
# after creation of the shapefiles, the prj should hold a text identical to the WKT in the grey box half way the page (vital QC check)

#coordinate and vertical reference systems, commonly used:

# WGS 84 geographical 4326
# ED 50 geographical 4230 
# WGS 84 / UTM 31N 32631 
# WGS 84 / UTM 28N 32628
# ETRS89 / UTM 31N 25831
# ED50 UTM 31N: 23031
# RD new: 28992

epsg_in='28992'
epsg_out='28992'

vdatum='NAP'

# start_filename="beton_2001_lijn"
# start_filename="maasvl_1997_lijn"
# end_filename="_PRC"

#end of user inputs
#%%

folder = Path(segy_files[0]).parent

crs_in='epsg:{}'.format(epsg_in)
crs_out='epsg:{}'.format(epsg_out)

#initiate timer
t0=time.time()

def opensegy(file, endian,count):
    global df1, sp0, spn
    from os.path import basename
    
    with segyio.open(file, ignore_geometry=True, endian=endian) as f:
    #with segyio.open(file, mode="r+", iline = 189, xline = 193, endian='big') as f:
            
            if replace_shotpoint==True:
                sp=np.arange(1,len(f.trace)+1)
            else: 
                sp= f.attributes(point_byte)[:]
                
            if f.attributes(scalar_byte)[:].any()==0:
                scalar=scalar_corr
            else:
                scalar=f.attributes(scalar_byte)[:]
                
            if scalar[10]<0:
                scalar=-1/scalar[10]
            
            
            sp=sp.tolist() #convert to native python interger, required to create shapefile
            x = f.attributes(easting_byte)[:]*scalar
            y = f.attributes(northing_byte)[:]*scalar
            
            
            if pos_in_arcsec==True:
                x=x/3600; y=y/3600
                
            #qc step, print min/max values
            # print("xmin/xmax = {}/{}) and ymin/ymax = {}/{}   ||| {}".format(min(x),max(x),min(y),max(y), file))
            
            #Extract linename  
            filename = basename(file)
            line=[filename]*len(x) #use to insert filename as a whole
            
            sp0=sp[0]
            spn=sp[-1]
            
            #store results in dataframe
            df1=pd.DataFrame([line,sp,x,y])
         
            
            print("{}/{}: {} || first/last {}: {}/{}".format(count,len(segy_files),file,point_type,sp0,spn))
            
            
            return df1,sp0,spn

#get filenames and loop over all files
count=0
df=pd.DataFrame()
first_sp=[]; last_sp=[]
for file in segy_files:
    count=count+1
    try: opensegy(file, endian='big',count=count) 
    except: opensegy(file, endian='little', count=count)
    df=pd.concat([df,df1],axis=1)
    first_sp.append(sp0)
    last_sp.append(spn)
        

df=df.T
df.columns= columns=["Line name", point_type,"Easting","Northing"]

df.Easting=df.Easting.apply(float)
df.Northing=df.Northing.apply(float)

df=df[df.Northing>0]
df=df[df.Easting>0]


t1=int(time.time()-t0)
print("Data extracted to dataframe in "+ str(t1) + " seconds")

#%%
#output point shapefile

if Path(folder / 'GIS').is_dir():
    print("Writing files to existing GIS folder...")

(folder / 'GIS').mkdir(parents=True, exist_ok=True)

geometry=[Point(xy) for xy in zip(df.Easting, df.Northing)]
df = gpd.GeoDataFrame(df, crs=crs_in, geometry=geometry).reset_index(drop=True)
df=df.to_crs(crs_out)

if trackpoints is True:
    df.to_file(folder / 'GIS' / f'{outfile}_Trackpoints.shp', driver='ESRI Shapefile')

t2=int(time.time()-t1-t0)

if trackpoints is True:
    print("Trackpoint shapefile created in "+ str(t2) + " seconds")

#%%
#output line shapefile
if tracklines is True:
    dfl= df.groupby('Line name', as_index=False).agg({'geometry': lambda x: LineString(x.tolist())})

    dfl = gpd.GeoDataFrame(dfl, crs=crs_out, geometry='geometry')
    try:
        dfl['First {}'.format(point_type)]=first_sp
        dfl['Last  {}'.format(point_type)]=last_sp
    except: print("First and last {} cannot be added to Trackline shapefile, possibly due to line name problems (check splitter)".format(point_type))

    outfile_shape = '{}_Tracklines_{}.shp'.format(outfile, epsg_out)
    dfl.to_file(folder / 'GIS' / outfile_shape, driver='ESRI Shapefile')

    t3=int(time.time()-t2-t0)

    print("Trackline shapefile created in "+ str(t3) + " seconds")

