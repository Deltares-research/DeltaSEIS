# -*- coding: utf-8 -*-
"""
The scripts reads data from the LUNA dss system into a deltaseis standard format. 
Sampling rate (fs) and gauge length (dx) are infered from the file.

"""

def read_semd(file_list, verbose=True):
    """
    Read SPECFEM acoustic-elastic modeling output files and return in a
    deltaseis-compatable format

    Parameters
    ----------
    file_list : list
        A list of filepaths of individual semd files, each containing a trace
    verbose : boolean, optional
        Shows which semd file is being processed. The default is True.

    Returns
    -------
    data : 2d array (pandas dataframe)
        Resulting seismic prestack data (shot)
    fs : float
        Sampling rate (Hertz) of the input data
    t_start : TYPE
        Starting time of the data. This is not always zero but often
        negative (which would be acausal). If this is the case the
        script throws a warning

    """
    
    import pandas as pd

    
    
    trace_list = []

    for i, file in enumerate(file_list):
        
        #print filename in progress
        if verbose==True:
            print("{}/{}: {}".format(i+1, len(file_list),file))
        
        #read the trace for each semd file
        trace = pd.read_csv(file, usecols=[1], header=None, delim_whitespace=True)
        
        #store traces in list
        trace_list.append(trace)
        
        #retrieve the time vector and calcualte fs using last semd file
        if file==file_list[-1]:
            
            time_vector = pd.read_csv(file, usecols=[0], header=None, delim_whitespace=True)
            
            num_samples = len(time_vector)                  #number of samples in each trace
            dt = time_vector[0][1] - time_vector[0][0]      #sampling interval in seconds
            fs = 1/dt                                       #sampling rate in hz
            
            t_start = time_vector[0][0]                     # start time in semd files in seconds
            
    #concatenate traces in trace_list to generate data (shot record)
    data = pd.concat(trace_list, axis=1)
    data = data.values
    
    #print relevant variables for QC
    print("\nSampling interval = {} microseconds, sampling rate = {} Hertz".format(round(dt*1e6, 2), int(fs)))
    print("Number of samples = {}, Record length = {} seconds".format(num_samples, round(num_samples*dt, 4)))


    if t_start != 0:
         print("\nWARNING time vector starts at {} seconds".format(time_vector[0][1]))   
    
 
    return data, fs, t_start
