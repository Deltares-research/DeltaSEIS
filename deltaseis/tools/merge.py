# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 16:17:17 2024

variety of functions not fitting in another object


@author: vermaas
"""

from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt

from deltaseis import Segy_edit


def merge_segys(filelist, make_plot=True, record_length=None):
    """
    Returns one Segy_edit() object of a list of segy files. The segy's are merged based on their order in the list.

    TODO: automatically detect order by checking start/end coordinates of segy's (mind erroneous coordinates, often first/last SP's...)

    Parameters
    ----------
    filelist : list of strings or Paths
        full paths to segy files (.sgy, .seg, .segy) to merge, in correct order
    make_plot : boolean
        flag to plot coordinates of merged segy files using rainbow color, to allow quick visual check if order is correct (default=True)

    """
    
    # empty variables
    x = []
    y = []
    z = []
    cdps = []
    cdpx = []
    cdpy = []
    ffid = []
    channel_numbers = []
    groupx = []
    groupy = []
    groupz = []
    data = []
    fold = []
    offsets = []
    recording_delay = []

    if make_plot:
        colors = plt.cm.jet(np.linspace(0, 1, len(filelist)))
        n = 0

    # loop over segy files and append variables
    for file in filelist:
        print(f"reading file {file}")
        s = Segy_edit(file)

        if record_length is not None:
            s.set_record_length(record_length)
            s.spec.samples = s.spec.samples[:len(s.trace_data[0])]

        ffid = np.append(ffid, s.ffid)
        channel_numbers = np.append(channel_numbers, s.channel_numbers)    
        x = np.append(x, s.x)
        y = np.append(y, s.y)
        z = np.append(z, s.z)
        groupx = np.append(groupx, s.groupx)
        groupy = np.append(groupy, s.groupy)
        groupz = np.append(groupz, s.groupz)
        offsets = np.append(offsets, s.offsets)
        cdps = np.append(cdps, s.cdps)
        cdpx = np.append(cdpx, s.cdpx)
        cdpy = np.append(cdpy, s.cdpy)
        data = data + s.trace_data
        fold = np.append(fold, s.fold)
        recording_delay = np.append(recording_delay, s.recording_delay)
       
        
                        
        if make_plot:
            plt.plot(s.x, s.y, color=colors[n])
            n = n + 1
    
    if make_plot:
        plt.axis("equal")
        plt.grid()

    # overwrite data in last segy and update e.g. shotpointnr
    s.ffid = ffid.astype("int32")
    s.channel_numbers = channel_numbers.astype("int32")
    s.x = x.astype("int32")
    s.y = y.astype("int32")
    s.z = z.astype("int32")
    s.groupx = groupx.astype("int32")
    s.groupy = groupy.astype("int32")
    s.groupz = groupz.astype("int32")
    s.offsets = offsets.astype("int32")
    s.cdps = cdps.astype("int32")
    s.cdpx = cdpx.astype("int32")
    s.cdpy = cdpy.astype("int32")
    s.trace_data = data
    s.fold = fold.astype("int32")
    s.trace_number = len(s.x)
    s.spec.tracecount = len(s.x)
    s.indices = np.arange(s.trace_number)
    s.trace_sequence = np.arange(s.trace_number)
    s.recording_delay = recording_delay.astype("int32")
    s.renumber_shotpoints(0)

    return s
