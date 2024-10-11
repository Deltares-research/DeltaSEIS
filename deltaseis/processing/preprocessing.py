# -*- coding: utf-8 -*-

"""
The scripts reads a deltaseis-compatable seismic format, resamples using a
user defined new sampling interval and return it again in a deltaseis format.

The scripts concatenates the data into a deltaseis-compatable format for further use. 

TODO: add more obspy preprocessing steps and make this

"""


def resample(data, fs, fs_new):
    """
    Resampling seismic data

    Parameters
    ----------
    data : deltaseis 2D dataframe
        Any seismic data in the deltaseis format (pre-stack shots or post-stack sections)
    fs: float
        Sampling rate of the input seismic data
    fs_new : float
        Desired new sampling rate. Should be in the order of 10x the high frequency in data
  
    Returns
    -------
    data_resampled : deltaseis 2D dataframe
        Downsampled seismic data in deltaseis format

    """
    
    
    #converting to obspy stream object    
    st = deltaseis_to_obspy(data, fs)
    
    num_samples_input = data.shape[0]
    
    #resample by factor
    st.resample(fs_new)
    
    #convert back to deltaseis-compatable
    data_resampled = obspy_to_deltaseis(st)
    num_samples_output = data_resampled.shape[0]
    
    factor = int(num_samples_input / num_samples_output)
    
    print("\nThe data was downsampled by factor = {}".format(factor))
    
    #print relevant variables for QC
    print("\nSampling interval = {} microseconds, sampling rate = {} Hertz".format(round(1e6/fs_new, 2), int(fs_new)))
    print("Number of samples = {}, Record length = {} seconds".format(num_samples_output, round((num_samples_output-1)/fs_new, 4)))
    
    return data_resampled
    
    
def deltaseis_to_obspy(data, fs):
    """
    Converts a gmpt seismic format to obspy for specialized pre-processing

    Parameters
    ----------
    data : deltaseis 2D dataframe
        Any seismic data in the deltaseis format (pre-stack shots or post-stack sections)
    fs : float
        Sampling rate of the seismic data

    Returns
    -------
    st : object
        obspy stream object
        
    """
    
    from obspy import Trace, Stream, UTCDateTime
    import numpy as np
    
    size = np.shape(data)
    st = Stream()

    for i in range(0, size[1]):

        trace = Trace(header={'station': 'seismic', 'channel': 'Z'+str(i)})
        trace.data = data[:, i]

        # --Default values
        trace.stats.starttime = UTCDateTime(2011, 11, 11, 11, 11, 0)
        trace.stats.npts = len(trace.data)
        trace.stats.sampling_rate = fs
        trace.stats.station = str(i)
        trace.stats.channel = 'Z'+str(i)
        trace.stats.distance = i
        trace.stats.network = 'seismic'

        st.append(trace)

    return st


def obspy_to_deltaseis(st):
    """
    Converts an obspy stream object into the deltaseis-compatable format

    Parameters
    ----------
    st : object
        obspy stream object

    Returns
    -------
    data : deltaseis 2D dataframe
        Seismic data in the deltaseis format.

    """

    import numpy as np
    
    size = np.shape(st)
    data = np.zeros((size[1], size[0]))

    for ii in range(0, size[0]):

        data[:, ii] = st[ii].data

    return data