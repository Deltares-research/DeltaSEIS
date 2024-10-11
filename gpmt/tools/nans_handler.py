# -*- coding: utf-8 -*-
"""
Created on December 2020 @author: r.nieboer: 
    
"""
def prep_for_filter(horizon):
    
    

#%%
def padding_symmetric(vector, padding_length, truncate_data=0):
    '''
    Pads values to a vector, so that the trend of a line is followed ('linearly' extrapolated),
    this avoids edge effects in later filtering with scipy medfilt, that per default
    applies zero-padding, which creates edge effects. Optionally,the input data can be truncated at the start 
    and the end of the line if anomoulous values occur there (peaks at start and/or end of line)

    Parameters
    ----------
    vector : list with e.g. x or y coordinates
       
    padding_length : int with number of values to be padded on both sides of the vector
    this can be any value for instance equal to the length of the planned filter window
    
    truncate_data: int. Removes points from the start and the end of line if peaks occur in one and/or the other

    Returns
    -------
    vector_padded : the padded vector which has as lenght len(vector) + 2*padding_length

    '''
    import numpy as np
    
    if truncate_data > 0:
        vector = vector[truncate_data: -truncate_data]
    elif truncate_data == 0:
        vector = vector

    vector_padded = np.pad(vector, (padding_length + truncate_data), mode='symmetric', reflect_type='odd')
    
    return vector_padded


#%%
def get_fit(x, deg): 
    
    '''
    Return a polygonial fit of a chosen degree through a 1D vector
    and print the residual of the least squares fit
    
    Parameters
    ----------
    x numpy list with values to be fitted
    deg: polynomial degree to be used (residual resulting from choice will be printed)
    
    Returns
    -------
    
    f(k): curve fitting the data with k being the index of the original list
    '''
    
    import numpy as np
    
    k = np.arange(len(x))

    p, res, _, _, _= np.polyfit(k,x,deg, full=True)
    
    f = np.poly1d(p)
    
    print("Residual of the least square fit with degree {} is {}".format(deg,res))
    
    return f(k)

#%%
def make_polyfit(x, degree_polyfit, truncate_data = 0):
    
    '''
    Replace the original navigation data with a polynomial fit of that data.
    This is a last resort to be able to get some sensible line navigation,
    but it does not represent true navigation

    Parameters
    ----------
    x : numpy list to be fit\n
    degree: degree of polynomial fit that will replace the navigation data \n
    truncate_data: int. default = 0
    Removes points from the start and the end of line if peaks occur in one and/or the other
    
    Returns
    -------
    x1 : numpy list holding the polyfitted curve \n
   
    '''
    
    #padd the data symmetrically
    x_padded = padding_symmetric(x, degree_polyfit, truncate_data)
 
    #find a polynomially fitted curve through the data
    x_polyfit_padded = get_fit(x_padded, degree_polyfit)

    #truncate the fitted curve to original length of the data
    x1 = x_polyfit_padded[degree_polyfit:-degree_polyfit]
    
    return x1
   

#%%
def get_bathy(x,y,spl, offset=0):
    
    if transform_coordinates==True:
        from pyproj import Transformer
        transformer = Transformer.from_crs(epsg_seismic, epsg_bathy)
        
        if invert_xy==True:
            y, x = transformer.transform(x, y)
        else:
            x, y = transformer.transform(x, y)
            
    depths = spl(x,y,grid=False)

    depths = -np.round(depths,2)

    
    #truncate data to remove edge effects and 'extrapolate' by padding symmetrically
    depths = padding_symmetric(depths, padding_length=0, truncate_data=20)
    
    #if there is a misfit between the seismic and bathy seabed, bathy can be shifted by 'offset' to match
    if offset > 0:
        depths = padding_symmetric(depths, padding_length=abs(offset), truncate_data=0)[2*offset:]
    
    elif offset < 0: 
        depths = padding_symmetric(depths, padding_length=abs(offset), truncate_data=0)[:2*offset]
        
    return depths

#%%

def longest_valid_section(x):
    
        
    # find where null
    w = np.where(np.isnan(x))[0]
    w = np.append(w, len(x)-1)
    w = np.insert(w, 0, 0)
    
    try:
    # diff to find length of stretch and argmax to find where largest stretch
        v = np.diff(w).argmax()
    
    except ValueError:
        return print("No nans encountered in input vector"), x
    
    z = x[w[v]+1:w[v+1]]   
        
    # return longest stretch
    return z
    
    

#%% 
def find_shift(a,b):
    '''
    Find the shift between 2 similar signals using cross-correlation

    Parameters
    ----------
    a : list
        first signal
    b : list
        second signal

    Returns
    -------
    lag: int
    shift required in a to match up with b
    (or the other way around, whatever)

    '''
    import numpy as np
    from scipy.signal import correlate, correlation_lags
    
    #find longest sections without nans
    
    a = longest_valid_section(a.tolist())
    b = longest_valid_section(b.tolist())
    
    a = np.divide(a, 0.75) # convert bathy to m
    a = a - np.mean(a) #remove dc component
    b = b - np.mean(b) #remove dc component

    #calculate cross correlation and lags
    corr = correlate(a,b, mode='full')
    lags = correlation_lags(len(a),len(b), mode='full')
    
    # plt.plot(lags, corr)
    
    #find location of maximum value for correlation
    lag_idx = np.argmax(corr)
    # 
    #find corresponding lag
    lag = lags[lag_idx]
        
    return lag

#%%
def vertical_correct_segy(filepath, shifts, endian='big'):
    
    outpath=filepath.replace(".sgy", "_tideSD.sgy")

    with segyio.open(filepath, ignore_geometry=True, endian=endian) as src:
        spec = segyio.tools.metadata(src)
        sampling_interval = src.bin[3217]
        
               
        with segyio.create(outpath, spec) as dst:
            dst.text[0] = src.text[0]
            dst.bin = src.bin
            dst.header = src.header
            
            for i, trace in enumerate(src.trace):
               shift_samples=int(round(float(shifts[i]/(1e-3*sampling_interval)))) #in millisecond
                
               if shift_samples>0:
                   trace = np.pad(trace,(shift_samples,0),'constant')[:-shift_samples]
               
               elif shift_samples < 0:                   
                   trace = np.pad(trace, (0, abs(shift_samples)),'constant')[-shift_samples:]
               
               else: trace = trace
               
               dst.trace[i]=trace

#%%  import grid and calculate spline, or loads earlier saved spline
if bathydepth==True:
            spl = import_bathy(asciigrid)

#%%

#set and make output dirs if they don't exist

bottompath=os.path.dirname(filepaths[0])+"/"+"QC"
try: os.mkdir(bottompath)
except FileExistsError: print("Writing bottom to existing {} output folder".format(bottompath))

if output_bathy_segy_diff==True:
    if bathydepth==False: 
        output_bathy_segy_diff=False
        print("Bathydepth is set to False, tideSD cannot be calculated")

#loop over all input segy files0
count=0
tide_all=pd.DataFrame()
for i, filepath in enumerate(filepaths): 
    
    filename=os.path.basename(filepaths[i])
    print("{} of {}, {}".format(count+1,len(filepaths), filepath))
    count+=1
    
    try:
        results=seabed_pick(filepath, peak_threshold,endian='big',manyear=manyear,manday=manday,div=div)
    except RuntimeError:
        results=seabed_pick(filepath, peak_threshold,endian='little', manyear=manyear,manday=manday,div=div)
    
    if filterpick==True:
                results['seabed(ms)']=savgol_filter(results['seabed(ms)'],window_length,polyorder)
    results['seabed(ms)']=results['seabed(ms)'].round(2)
    
    if bathydepth==True: #write to function in v8, instead of reapeating same steps
        
        results["bathy(m)"]=get_bathy(results.lon,results.lat,spl)
        results["bathy(m)"]=np.where(results["bathy(m)"]<mindepth,np.nan,results["bathy(m)"]) #bath < mindepth to nan 
        results["bathy(m)"].interpolate(limit=limit, inplace=True) #interpolate the smaller gabs
        
               
        #remove interpolation edge effects by setting bathy values before and after a gab to nan and then interpolate and/or extrapolate if set by user
        nanedges=np.diff(np.where(np.isnan(results["bathy(m)"]),0,1)) #1 where valid bathy goes to gab and -1 where gab goes to valid bathy
      
        for idx in np.argwhere(nanedges==-1).flatten(): #remove values before gabs
            b=idx-taper
            if idx-taper<0: b=0
            results.loc[b:idx,'bathy(m)']=np.nan
                
        for idx in np.argwhere(nanedges==1).flatten(): #remove values after gabs
            a=idx+taper
            if idx+taper>len(results): a=len(results)
            results.loc[idx:a,'bathy(m)']=np.nan
            
        if shiftbathy == True:
        
            offset = find_shift(results["bathy(m)"],results["seabed(ms)"])
            print("offset = {} ".format(offset))
        
        else:
            offset = 0
                   
        
        results["bathy(m)"]=get_bathy(results.lon,results.lat,spl, offset=offset)
        results["bathy(m)"]=np.where(results["bathy(m)"]<mindepth,np.nan,results["bathy(m)"]) #bath < mindepth to nan 
        results["bathy(m)"].interpolate(limit=limit, inplace=True) #interpolate the smaller gabs
        
        #remove interpolation edge effects by setting bathy values before and after a gab to nan and then interpolate and/or extrapolate if set by user
        nanedges=np.diff(np.where(np.isnan(results["bathy(m)"]),0,1)) #1 where valid bathy goes to gab and -1 where gab goes to valid bathy
      
        for idx in np.argwhere(nanedges==-1).flatten(): #remove values before gabs
            b=idx-taper
            if idx-taper<0: b=0
            results.loc[b:idx,'bathy(m)']=np.nan
                
        for idx in np.argwhere(nanedges==1).flatten(): #remove values after gabs
            a=idx+taper
            if idx+taper>len(results): a=len(results)
            results.loc[idx:a,'bathy(m)']=np.nan
            
        invalid=results['bathy(m)'].isna().sum() #count nan values
  
        
        if invalid>100:
            print("WARNING: no bathymetry found for {} pings ({}%)".format(invalid, int(100*invalid/len(results))))
   
        results["bathy(ms)"]=(results["bathy(m)"]/0.75).round(2)
        
        
         
    if output_bathy_segy_diff==True:
        
        st=results["bathy(m)"].first_valid_index()
        end=results["bathy(m)"].last_valid_index()
   
        results["bathy(m)"].interpolate(inplace=True, limit_area=('inside'))
        nanedges_trim=nanedges[st:end]
        
        startnans=np.argwhere(nanedges_trim==-1).flatten()
        endnans=np.argwhere(nanedges_trim==1).flatten()
        
        for i,val in enumerate(startnans):
            try:
                print ("WARNING: bathymetry interpolated in pings range {}-{}".format(startnans[i],endnans[i]))
            except IndexError:
                print ("ERROR, no bathymetry coverage along this line")
                
        results["bathy(m)"].loc[results["bathy(m)"].first_valid_index():results["bathy(m)"].last_valid_index()]
        results["bathy(ms)"]=(results["bathy(m)"]/0.75).round(2)
        
        results["diff(ms)"]=results['bathy(ms)']-results['seabed(ms)']
        
        #set difference grid to nan where bathy has gabs
        for i,val in enumerate(startnans):
            results["diff(ms)"].values[startnans[i]-trim:endnans[i]+trim] = np.nan 
        
     
        try:
            if poly_diff ==True:
                    results["tideSD(ms)"] = np.nan
                    results["diff(ms)"].interpolate(method='pchip', order=6,inplace=True, limit_area=('inside'))
                    idx = np.isfinite(results["diff(ms)"])
                    results["tideSD(ms)"][idx] = make_polyfit(results["diff(ms)"][idx],polyfit_order)
            else:
                diff_savgol=savgol_filter(results["diff(ms)"].loc[st:end],diff_window_length,diff_polyorder)
                results["tideSD(ms)"]=np.pad(diff_savgol,(st,len(results)-end-1), 'constant', constant_values=(np.nan,np.nan))
                
        except:
            results["tideSD(ms)"]=[0]*len(results) #set all bathy to zero, if no coverage
            if len(results)<diff_window_length:
                print("Error: diff_window_length ({}) should be set smaller than number of pings ({}), tidesd set to zero".format(diff_window_length,len(results)))
            
        #construct datetime-tide file every x minutes as defined by t_int 
        tide=pd.DataFrame()
        tide['tideSD(m)']=(-results['tideSD(ms)']*0.75).round(2)
        tide['datetime']=pd.to_datetime(results['date']+' '+results['time'],dayfirst=True,errors='coerce')
        tide=tide.set_index('datetime')
        tide=tide.resample(t_int).median()
        tide['date']=tide.index.strftime("%d/%m/%Y")
        tide['time']=tide.index.strftime("%H:%M:%S")
        tide=tide[["date","time","tideSD(m)"]]
        
        tide.to_csv(bottompath+"/"+os.path.splitext(filename)[0]+"_tideSD.txt", header=True, index=False, sep=" ",quoting=csv.QUOTE_NONE, escapechar=' ')   
        tide_all =  pd.concat([tide_all, tide])
   
        plt.ioff() #don't show but save directly
         
        if bathydepth==False:results.plot(x='ping', y='seabed(ms)', grid=True)
        elif output_bathy_segy_diff==False:results.plot(x='ping', y=['bathy(ms)','seabed(ms)'], grid=True)
        else:results.plot(x='ping', y=['bathy(ms)','seabed(ms)','diff(ms)','tideSD(ms)'], grid=True)
        plt.gca().invert_yaxis()
        plt.savefig(bottompath+"/"+os.path.splitext(filename)[0]+"_QCplot.png")
        plt.grid() 
        
        if qc_plots==True:
            plt.show()
        else: plt.close()
         
    if bathydepth==True: 
        bottom=results[["ping","seabed(ms)","bathy(ms)"]] 
        if amplitude==True: 
            bottom=results[["ping","seabed(ms)","bathy(ms)","amplitude"]]
    if bathydepth==False: 
        bottom=results[["ping","seabed(ms)"]] 
        if amplitude==True: 
            bottom=results[["ping","seabed(ms)","amplitude"]] 
    
    #write results to bottom
    bottom.to_csv(bottompath+"/"+os.path.splitext(filename)[0]+"_BOTTOM.txt",header=True, index=False, sep=" ")   
    
    if correct_segy==True:
        shifts = results["diff(ms)"]
        try:
            vertical_correct_segy(filepath, shifts, endian='big')
        except RuntimeError:
            vertical_correct_segy(filepath, shifts, endian='little')
            

if output_bathy_segy_diff==True:
    tide_all.sort_index(inplace=True) #sort chronologically
    
    j=0
    while os.path.exists(bottompath+"/tideSD_all%s.txt" % j): j+=1  
    tide_all.to_csv(bottompath+"/tideSD_all{}.txt".format(j),header=True, index=False,sep=" ")


    
    
    
print('Done, total runtime: {} seconds.'.format(int(time.time()-start)))