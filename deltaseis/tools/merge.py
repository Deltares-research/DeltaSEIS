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


def merge_segys(filelist, make_plot=True):
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
    groupx = []
    groupy = []
    data = []

    if make_plot:
        colors = plt.cm.jet(np.linspace(0, 1, len(filelist)))
        n = 0

    # loop over segy files and append variables
    for file in filelist:
        print(f"reading file {file}")
        s = Segy_edit(file)
        x = np.append(x, s.x)
        y = np.append(y, s.y)
        groupx = np.append(groupx, s.groupx)
        groupy = np.append(groupy, s.groupy)
        data = data + s.trace_data

        if make_plot:
            plt.plot(s.x, s.y, color=colors[n])
            n = n + 1
    plt.axis("equal")
    plt.grid()

    # overwrite data in last segy and update e.g. shotpointnr
    s.x = x.astype("int32")
    s.y = y.astype("int32")
    s.groupx = groupx.astype("int32")
    s.groupy = groupy.astype("int32")
    s.trace_data = data

    s.trace_number = len(x)
    s.spec.tracecount = len(x)
    s.renumber_shotpoints(0)

    return s
