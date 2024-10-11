# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import medfilt


def docstring():
    """
    Edit the navigation meta data of seg-y files, e.g. remove outliers due to GPS wander
    
    Parameters
    ----------
    segy_files : select seg-y files to corrected via the file dialog, alternatively
    files can be specified as a list, but this is deactivated per default ("#' commented)
    Extentions include: .sgy, .seg, segy
    
    remove_junk (boolean): some seg-y files have junk null bytes at the end of the file
    so that the python package segyio cannot load it. If this is the case can be expected in seisee. 
    If remove_junk is set to True, these null bytes are removed, but it will increase the processing time.
    If set to False (default), the program wil suggest to set to true if the seg-y does not load. 
    
    show_figures (boolean): if set to True, this will show figures with the results of before and after correction.
    Can be set to False when processing large batches of files but this is not recommended.
    
    write_to_segy (boolean): if set to true (default) the program will output the seg-y with the corrected navigation data.
    Can be set to False while experimenting with different 'outlier_mode' options and filter/polyfit settings
    
    outlier_mode (string): different techniques can be used to get rid of the outliers in the navigation data, which one
    is most useful depends on the input data and might be subject to trial-and-error
    
        - 'median': creates a curve applying median filter with a long window (window_outlier_filter) through the input data points (using scipy.medfilt) and discards
          points in the input data that are a certain distance (maximum_dist) away from this curve. The empty points are then linearly interpolated.
        
        - 'polyfit': creates a curve by making a polynomial fit (of degree 'degree_polyfit') of the input data points (using numpy.polyfit) and discards
          points in the input data that are a certain distance (maximum_dist) away from this curve. The empty points are then linearly interpolated.
        
        - 'replace': replaces the input data with the curve of a polynomial fit (of degree 'degreee_polyfit') of the input data points (using numpy.polyfit).
          This is a last resort to gain a seg-file with loadable navigation, though it is in no sense the true navigation, therefore this setting is not recommended.
            
        - 'none': no outliers will be discarded, this can be used if only a 'smooth' of the navigational data is required
        
    smooth (boolean): if set to True, the input data is replaced by a smoothed version of it. 
    A median filter with window 'window_smooth_filter' (using scipy.medfilt) is applied
    
    window_outlier_filter (interger): (long) filter window used for the 'median' outlier_mode, value should be odd
    
    degree_polyfit (integer): degree of the polynomial to be fitted through the data for the 'polyfit' or 'replace' outlier_mode options
    
    maximum_dist (integer): distance to discard outliers in the outlier_mode options. the distance can have different unit depending on what
    the trace headers of the orignal seg-y file was in so it can be meters, but also arcseconds (or anything else).
    
    truncate_data (integer): advanced option to remove points at the start and end of the data in case spikes occur there. this
    makes sure that the symmetric padding that is applied before any filtering/fitting follows the general line trend.  
    
    endian (string): state the byte order of the seg-y input file (inspect in seisee)
    
    x_byte (int): location of the x variable in the seg-y trace header (default = 73)
    
    y_bute (int): location of the y variable in the seg-y trace header (default = 77)
    
    scalar_byte (int) = location of the coordinate scalar in the seg-y trace header (default = 71)
    
    Program Returns
    -------
    seg-y file (.sgy) with no junk bytes at the end of file and corrected navigation added to the trace headers
    figures showing the results of selected the navigation correction mode ('outlier_mode' and 'smooth')
    the corrected x and y coordinates as list for further analysis or plotting
    
    
    @author: Roeland Nieboer @Deltares
    """  
     
            
    #%%   

def plot_xy_changed(x,y,x1,y1,title='',xlabel='original data',x1label='changed data'):
        fig, axs = plt.subplots(1, 3, figsize=(18, 6))
        fig.canvas.manager.set_window_title(title)
        fig.suptitle(title, fontsize=16)

        axs[0].set_title("x coordinates")
        axs[0].plot(x, label=xlabel)
        axs[0].plot(x1, label=x1label) 

        axs[0].legend()
        axs[0].grid()

        axs[1].set_title("y coordinates")    
        axs[1].plot(y)
        axs[1].plot(y1)

        axs[1].grid()

        axs[2].set_title("x versus y (equal axis)")    
        axs[2].plot(x, y)
        axs[2].plot(x1, y1)

        axs[2].grid()

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
        
    k = np.arange(len(x))

    p, res, _, _, _= np.polyfit(k,x,deg, full=True)
    
    f = np.poly1d(p)
    
    print("Residual of the least square fit with degree {} is {}".format(deg,res))
    
    return f(k)

#%%

def remove_outliers_median(x, y , order, max_dist, truncate_data, show_figures, save_path=None):
    
    '''
    Discard and interpolate coordinate points that are more than max_dist meters away
    from a strongly median filtered version of the same points

    Parameters
    ----------
    x : numpy list of faulty x coordinates \n
    y : numpy list of faulty y coordinates \n
    order: n-th order median filter that will be applied to the data. Distance of the original
    data to the filtered output will determine if points will be discarded as outliers
    truncate_data: int. Removes points from the start and the end of line if peaks occur in one and/or the other
    
    Returns
    -------
    x1 : numpy list of corrected x coordinates \n
    y1 : numpy list of corrected y coordinates \n
   

    '''
    
    #padd the data symmetrically
    x_padded = padding_symmetric(x, order, truncate_data)
    y_padded = padding_symmetric(y, order, truncate_data)    
    
    #median filter the padded data
    x_filtered_padded = medfilt(x_padded, order)
    y_filtered_padded = medfilt(y_padded, order)
    
    #truncate the filtered data to original length
    x_filtered = x_filtered_padded[order:-order]
    y_filtered = y_filtered_padded[order:-order]    

    #calculate the distance between the original data and filtered data
    x_diff = abs(x-x_filtered)
    y_diff = abs(y-y_filtered)  
    
    #dataframe to use pandas functionalities
    df = pd.DataFrame(data = {'x': x, 'y': y, 'x_diff': x_diff, 'y_diff': y_diff})  
     
    #discard outliers and interpolate linearly   
    df.loc[df['x_diff'] > max_dist, ['x']] = np.nan           
    df.loc[df['y_diff'] > max_dist, ['y']] = np.nan

    #keep first and last point, even if they were discard, to still be able to interpolate
    if np.isnan(df['x'].iloc[0]) == True:
        df.loc[0, 'x'] = x[0] 
        print('WARNING: first point of x was discarded as outlier, it was reinserted but consider updating order_outlier_filter for better results')
    if np.isnan(df['x'].iloc[-1]) == True: 
        df.loc[df.index[-1], 'x'] = x[-1]
        print('WARNING: last point of x was discarded as outlier, it was reinserted but consider updating order_outlier_filter for better results')
 
    if np.isnan(df['y'].iloc[0]) == True:
        df.loc[df.index[0], 'y'] = y[0]
        print('WARNING: first point of y was discarded as outlier, it was reinserted but consider updating order_outlier_filter for better results')
    if np.isnan(df['y'].iloc[-1]) == True:
        df.loc[df.index[-1], 'y'] = y[-1]
        print('WARNING: last point of y was discarded as outlier, it was reinserted but consider updating order_outlier_filter for better results')

    #interpolate nans    
    df[['x','y']] = df[['x', 'y']].interpolate()              
    
    #results from pandas to list
    x1 = df['x'].tolist() 
    y1 = df['y'].tolist() 
    
    
    if show_figures == True:
        plot_xy_changed(x,y,x1,y1,title='Remove outliers using a median filtered curve',x1label=f"corrected data (max_dist={max_dist})")        

    if save_path is not None:
        plt.savefig(save_path)
  
    return x1, y1 

#%%   
def remove_outliers_polyfit(x, y, degree_polyfit, max_dist, truncate_data, show_figures):
    
    '''
    Discard and interpolate coordinate points that are more than max_dist meters away
    from a curve fitted in a least squares sense to the original points

    Parameters
    ----------
    x : numpy list of faulty x coordinates \n
    y : numpy list of faulty y coordinates \n
    degree: degree of polynomial fitting curve. Distance of the original
    data to the fitted curve will determine if points will be discarded as outliers \n
    truncate_data: int. Removes points from the start and the end of line if peaks occur in one and/or the other
    
    Returns
    -------
    x1 : numpy list of corrected x coordinates \n
    y1 : numpy list of corrected y coordinates \n
   
    '''    
    #padd the data symmetrically
    x_padded = padding_symmetric(x, degree_polyfit, truncate_data)
    y_padded = padding_symmetric(y, degree_polyfit, truncate_data) 
        
    #find a polynomially fitted curve through the data
    x_polyfit_padded = get_fit(x_padded, degree_polyfit)
    y_polyfit_padded = get_fit(y_padded, degree_polyfit)
    
    #truncate the fitted curve to original length of the data
    x_polyfit = x_polyfit_padded[degree_polyfit:-degree_polyfit]
    y_polyfit = y_polyfit_padded[degree_polyfit:-degree_polyfit]    

    #calculate the distance between the original data and fitted curve
    x_diff = abs(x-x_polyfit)
    y_diff = abs(y-y_polyfit)  
    
    #dataframe to use pandas functionalities
    df = pd.DataFrame(data = {'x': x, 'y': y, 'x_diff': x_diff, 'y_diff': y_diff})  
     
    #discard outliers and interpolate linearly   
    df.loc[df['x_diff'] > max_dist, ['x']] = np.nan           
    df.loc[df['y_diff'] > max_dist, ['y']] = np.nan

    #keep first and last point, even if they were discard, to still be able to interpolate
    if np.isnan(df['x'].iloc[0]) == True:
        df['x'].iloc[0] = x[0]
        print('WARNING: first point of x was discarded as outlier, it was reinserted but consider updating degree_polyfit for better results')
    if np.isnan(df['x'].iloc[-1]) == True:
        df['x'].iloc[-1] = x[-1]
        print('WARNING: last point of x was discarded as outlier, it was reinserted but consider updating degree_polyfit for better results')
  
    if np.isnan(df['y'].iloc[0]) == True:
        df['y'].iloc[0] = y[0]  
        print('WARNING: first point of y was discarded as outlier, it was reinserted but consider updating degree_polyfit for better results')
    if np.isnan(df['y'].iloc[-1]) == True:
        df['y'].iloc[-1] = y[-1] 
        print('WARNING: last point of y was discarded as outlier, it was reinserted but consider updating degree_polyfit for better results')
        
    #interpolate nans   
    df[['x','y']] = df[['x', 'y']].interpolate()  
         
    
    #results from pandas to list
    x1 = df['x'].tolist() 
    y1 = df['y'].tolist() 
    
    if show_figures == True:
        plot_xy_changed(x,y,x1,y1,title='Remove outliers using a polyfit curve',x1label=f"corrected data (max_dist={max_dist})") 
    
    return x1, y1 

#%%
def smoothXY(x, y, order, truncate_data, show_figures):
    
    '''
    Smooth the coordinate point along a survey line 
    
    Parameters
    ----------
    x : numpy list of faulty x coordinates \n
    y : numpy list of faulty y coordinates \n
    order: n-th order median filter that will be applied to the data \n
    truncate_data: int. Removes points from the start and the end of line if peaks occur in one and/or the other
    
    Returns
    -------
    x1 : numpy list of corrected x coordinates \n
    y1 : numpy list of corrected y coordinates \n
    

    '''
    
    #padd the data symmetrically
    x_padded = padding_symmetric(x, order, truncate_data)
    y_padded = padding_symmetric(y, order, truncate_data)    
    
    #median filter the padded data
    x_filtered_padded = medfilt(x_padded, order)
    y_filtered_padded = medfilt(y_padded, order)
    
    #truncate the filtered data to original length
    x_filtered = x_filtered_padded[order:-order]
    y_filtered = y_filtered_padded[order:-order]   

    x1 = x_filtered.tolist()
    y1 = y_filtered.tolist()
    
    if show_figures == True:
        plot_xy_changed(x,y,x1,y1,title='Smooth navigation data using a median filter',x1label=f"filtered data (window={order})") 
       
    return x1, y1

#%%

def replace_with_polyfit(x, y, degree_polyfit, truncate_data, show_figures):
    
    '''
    Replace the original navigation data with a polynomial fit of that data.
    This is a last resort to be able to get some sensible line navigation,
    but it does not represent true navigation

    Parameters
    ----------
    x : numpy list of faulty x coordinates \n
    y : numpy list of faulty y coordinates \n
    degree: degree of polynomial fit that will replace the navigation data \n
    truncate_data: int. Removes points from the start and the end of line if peaks occur in one and/or the other
    
    Returns
    -------
    x1 : numpy list of corrected x coordinates \n
    y1 : numpy list of corrected y coordinates \n
   
    '''
        
    #padd the data symmetrically
    x_padded = padding_symmetric(x, degree_polyfit, truncate_data)
    y_padded = padding_symmetric(y, degree_polyfit, truncate_data) 
        
    #find a polynomially fitted curve through the data
    x_polyfit_padded = get_fit(x_padded, degree_polyfit)
    y_polyfit_padded = get_fit(y_padded, degree_polyfit)
    
    #truncate the fitted curve to original length of the data
    x1 = x_polyfit_padded[degree_polyfit:-degree_polyfit]
    y1 = y_polyfit_padded[degree_polyfit:-degree_polyfit]    

    
    if show_figures == True:
        plot_xy_changed(x,y,x1,y1,title='Replace navigation data a polyfit curve (last resort)',x1label=f"polynomial fit (deg={degree_polyfit})")
    
    return x1, y1 

def distance_from_xy(x,y):
    '''
    Calculate distance array from arrays with x an y coordinates

    Parameters
    ----------
    x : numpy array of (real) x coordinates
    y : numpy array of (real) y coordinates
    
    Returns
    -------
    array with distances between points and cumulative distance

    '''
    dx = np.diff(x).astype('float64')
    dy = np.diff(y).astype('float64')

    d = (dx**2 + dy**2)**0.5

    return d

def remove_start_end(x,y,max_distance,max_truncate_length):
    '''
    Function used to remove errorneous coordinates at start or end of line

    Parameters
    ----------
    x : numpy array
        real x coordinates
    y : numpy array
        real y coordinates
    max_distance : int or float
        maximum distance between two traces in real units (same as coordinate system)
    max_truncate_length : int
        maximum number of traces at start and end to find bad coordinates in 
        (outliers within line should be removed with e.g. fix_navigation_median() function)

    Returns
    -------
    array with distances between points and cumulative distance

    '''
    d = distance_from_xy(x,y)
    i_start = [i for i in range(max_truncate_length) if d[i]>max_distance]
    dflip = np.flip(d)
    i_end = [i for i in range(max_truncate_length) if dflip[i]>max_distance]

    if i_start: 
        n_bad_traces_start = i_start[0] + 1
    else:
        n_bad_traces_start = 0
    if i_end: 
        n_bad_traces_end = i_end[0] + 1
    else:
        n_bad_traces_end = 0

    return n_bad_traces_start, n_bad_traces_end