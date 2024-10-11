# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 13:24:26 2023
Read a grid and extract values at coordinates

Example: find the bathymetric depth along the navigation line of seismic data

@author: nieboer
edited by vermaas
"""

import numpy as np
from geost.spatial import get_raster_values


def get_grid(grid_path, x, y, velocity):
    '''Extract the bathymetric depth at specific x,y locations, the epsg from the
    grid coordinate data should be the same as the x, y coordinates'''

    depths = get_raster_values(x,y,grid_path)

    depths = -np.round(depths, 2)
    two_way_times = depth_time_conversion(depths, velocity)

    return depths, two_way_times


def depth_time_conversion(depths, velocity):
    '''Convert a depth grid to time using a constant velocity (m/s)'''
    two_way_time = depths/(0.5*velocity*1e-3)

    return two_way_time
