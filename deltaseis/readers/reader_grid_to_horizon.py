# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 13:24:26 2023
Read a grid and extract values at coordinates

Example: find the bathymetric depth along the navigation line of seismic data
TODO: there has to be a better solution for the x, y projected or y x geographic order problem

@author: nieboer
edited by vermaas
"""

import numpy as np
from geost.spatial import get_raster_values
from scipy.spatial import cKDTree


def get_grid(grid_path, x, y, velocity, grid_format="raster"):
    """Extract depths and two-way times at specified x,y locations.

    Parameters
    ----------
    grid_path : str or pathlib.Path
        Path to a raster grid or an XYZ point file.
    x, y : array-like
        Coordinates in the grid's coordinate reference system.
    velocity : float
        P-wave velocity in m/s for the depth-to-time conversion.
    grid_format : {"raster", "xyz"}, optional
        Use raster sampling for raster grids, or nearest-point lookup for
        headerless, whitespace-separated XYZ files. Defaults to "raster".
    """
    if grid_format == "raster":
        try:
            depths = get_raster_values(x, y, grid_path)
        except Exception:
            depths = get_raster_values(y, x, grid_path)
    elif grid_format == "xyz":
        points = np.loadtxt(grid_path, dtype=float, ndmin=2)
        if points.shape[1] != 3:
            raise ValueError(
                f"Expected an XYZ file with exactly 3 columns; found {points.shape[1]}."
            )

        valid_points = np.isfinite(points).all(axis=1)
        if not valid_points.any():
            raise ValueError("XYZ file contains no rows with finite X, Y, and Z values.")

        tree = cKDTree(points[valid_points, :2])
        query_coordinates = np.column_stack((x, y))
        depths = np.full(query_coordinates.shape[0], np.nan)
        valid_queries = np.isfinite(query_coordinates).all(axis=1)
        depths[valid_queries] = points[valid_points, 2][
            tree.query(query_coordinates[valid_queries])[1]
        ]
    else:
        raise ValueError(
            f"Unsupported grid_format {grid_format!r}; expected 'raster' or 'xyz'."
        )

    depths = -np.round(depths, 2)
    two_way_times = depth_time_conversion(depths, velocity)

    return depths, two_way_times


def depth_time_conversion(depths, velocity):
    '''Convert a depth grid to time using a constant velocity (m/s)'''
    two_way_time = depths/(0.5*velocity*1e-3)

    return two_way_time
