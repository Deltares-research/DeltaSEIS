# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 16:03:50 2023

Segy_editor using segyio library. It performs common edits and writes these back
to segy. The main operations it performs are:
    - set record length, resample, set endianness and scalar
    - transform coordinates, fix coordinate spikes
    - vertical trace corrections, e.g. for bulk shift or tide and senosor depth corrections
    (add other functions)


lower priority
near_trace_gather: use taiwan_near_trace_gather.py
output different data format e.g. from 4 byte float to 16 byte int (narrowing)

@author: nieboer
"""
import math
import warnings
from itertools import cycle

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import segyio
from obspy.signal.trigger import classic_sta_lta
from pyproj import Transformer
from pyproj import CRS
from scipy.signal import find_peaks, hilbert, savgol_filter
from shapely.geometry import LineString, Point

from deltaseis import deltaseis_to_obspy, resample
from deltaseis.readers.parser_seismic_segy import parse_segy
from deltaseis.readers.reader_grid_to_horizon import get_grid
from deltaseis.tools.segy_navigation_despiker import (
    remove_outliers_median,
    remove_outliers_polyfit,
    replace_with_polyfit,
    smoothXY,
    remove_start_end, 
)


class Segy_edit:
    """Segy editor for both tracedata and header meta data"""

    def __init__(self, segy_filepath):
        self.segy_filepath, self.endian, self.data_sample_format = parse_segy(
            segy_filepath
        )

        self.src = segyio.open(
            self.segy_filepath, ignore_geometry=True, endian=self.endian
        )

        self.spec = segyio.tools.metadata(self.src)
        self.text_header = self.src.text[0]
        self.bin_header = self.src.bin
        self.trace_header = self.src.header
        self.trace_data = [
            tr.astype(np.dtype(self.data_sample_format)) for tr in self.src.trace
        ]

        self.trace_number = len(self.src.trace)
        self.shotpoint_numbers = self.src.attributes(17)[:]
        self.trace_samples = len(self.spec.samples)

        self.sampling_interval = self.src.bin[3217]
        self.sampling_rate = 1 / (self.sampling_interval * 1e-6)

        self.trace_number_in_field_record = self.src.bin[3213]
        self.record_length = self.trace_samples * 1e-3 * self.sampling_interval

        self.x = self.src.attributes(73)[:]
        self.y = self.src.attributes(77)[:]
        self.scalar = self.src.attributes(71)[10][0]
        self.out_scalar = self.scalar

        self.groupx = self.src.attributes(81)[:]
        self.groupy = self.src.attributes(85)[:]

        self.horizons = {}
        self.shift_samples = np.zeros(self.trace_number)
        self.plot_number = 0

        self.crs = None
        print('Warning: CRS is not defined, please define with set_crs() if needed')

    def fix_multiple_scalars(self):
        """Fixes the scalar and coordinates in the segy file when multiple scalars are present in one file"""

        for i in range(self.trace_number):
            if self.trace_header[i][71] != self.scalar:
                # correct the coordinates to match the new scalar
                factor = self.scalar / self.trace_header[i][71]
                self.x[i] = (factor * self.x[i]).astype(
                    np.dtype(self.data_sample_format)
                )
                self.y[i] = (factor * self.y[i]).astype(
                    np.dtype(self.data_sample_format)
                )

    def interpolate_trace_header(self, trace_header, start_index=0, end_index=0):
        """Interpolate trace header values, e.g. to fill in missing values or to replace erroneous values"""

        trace_header[start_index:end_index] = np.interp(
            np.arange(start_index, end_index),
            np.r_[0:start_index, end_index : len(trace_header)],
            np.r_[trace_header[0:start_index], trace_header[end_index:]],
        )
        return trace_header

    def set_record_length(self, milliseconds):
        """Shorten the record length to decrease file size, this number should
        be smaller than the current record length, else the output is equal to the input file
        """

        new_sample_length = int(milliseconds / (self.sampling_interval * 1e-3) + 1)
        self.spec.samples = self.spec.samples[:new_sample_length]
        self.trace_data = [tr[:new_sample_length] for tr in self.trace_data]
        self.record_length = self.sampling_interval * 1e-3 * new_sample_length
        self.trace_samples = new_sample_length

    def resample(self, new_sampling_rate):
        """Resamples the seismic data when writing to file"""

        # cast into a 2D array and transpose for resampling with obspy
        data = np.array(self.trace_data).T
        resampled_data = resample(data, self.sampling_rate, new_sampling_rate)

        # set spec.samples for segyio write
        new_sample_length = resampled_data.shape[0]
        self.spec.samples = self.spec.samples[:new_sample_length]

        # bring back trace data to original format that can be written by segyio
        self.trace_data = list(resampled_data.T)
        self.sampling_rate = new_sampling_rate
        self.sampling_interval = math.ceil(1 / (self.sampling_rate * 1e-6))

    def set_endian(self, endian):
        """
        Set the endianness for the output file to 'big' or 'little'
        """
        self.spec.endian = endian

    def set_input_scalar(self, manual_scalar):
        """
        Set the coordinate scalar if not correct in the segy (e.g. zero)
        """
        self.scalar = manual_scalar
        self.out_scalar = self.scalar

    def set_crs(self,epsg):
        """
        Change crs attribute based on epsg code

        Parameters
        ----------
            - epsg : int
                epsg code for coordinate reference system

        Commonly used projections for the Netherlands
        - WGS 84 geographical: 4326
        - ED 50 geographical: 4230
        - WGS 84 / UTM 31N: 32631 (North Sea)
        - ETRS89 / UTM 31N: 25831 (North Sea)
        - ED50 UTM 31N: 23031 (North Sea)
        - RD new: 28992 (on land)
        """

        self.crs = CRS(epsg)

    def set_output_scalar(self, manual_out_scalar):
        """
        Set the coordinate scale if not correct in the segy (e.g. zero)
        """
        self.out_scalar = manual_out_scalar

    def xy_to_real(self):
        """
        Returns real x and y coordinates

        Parameters
        ----------
            - no parameters required

        Commonly used projections for the Netherlands
        - WGS 84 geographical: 4326
        - ED 50 geographical: 4230
        - WGS 84 / UTM 31N: 32631 (North Sea)
        - ETRS89 / UTM 31N: 25831 (North Sea)
        - ED50 UTM 31N: 23031 (North Sea)
        - RD new: 28992 (on land)
        """

        if not self.crs:
            print("CRS is not defined, cannot perform function - first use set_crs() to define CRS")

        else:
            # correct coordinates with scalar, assuming 0 should be 1:
            if self.scalar == 0:
                factor = 1
                print('Warning: scalar is 0, interpreted as 1, make sure this is correct or change mannually')
            elif self.scalar > 0:
                factor = self.scalar
            elif self.scalar < 0:
                factor = -1 / self.scalar

            if self.crs.is_geographic:
                x = factor * self.y / 3600
                y = factor * self.x / 3600
            else:
                x = factor * self.x
                y = factor * self.y

            return x,y

    def transform_coordinates(
        self, epsg_out, manual_in_scalar=0, apply_new_coordinates=True
    ):
        """
        Transform the navigation coordinates from the segy file

        Commonly used projections for the Netherlands

        - WGS 84 geographical: 4326
        - ED 50 geographical: 4230
        - WGS 84 / UTM 31N: 32631 (North Sea)
        - ETRS89 / UTM 31N: 25831 (North Sea)
        - ED50 UTM 31N: 23031 (North Sea)
        - RD new: 28992 (on land)

        Last to didgets of projected systems refer to the UTM zone
        Refer to spatialreference.org for other epsg codes
        """

        if manual_in_scalar != 0:
            scalar = self.set_input_scalar(manual_in_scalar)
            print(f"In-scalar of segy changed to {manual_in_scalar}")

        # apply the transformations
        transformer = Transformer.from_crs(self.crs, epsg_out)

        x, y = self.xy_to_real()
        x_transformed, y_transformed = transformer.transform(x,y)

        # write back into seg-y header required format
        if apply_new_coordinates == True:
            if self.out_scalar > 0:
                out_factor = self.out_scalar
            elif self.out_scalar < 0:
                out_factor = -1 / self.out_scalar

            with warnings.catch_warnings():
                warnings.filterwarnings("error", category=RuntimeWarning)

                try:
                    if CRS(epsg_out).is_projected:
                        self.x = (x_transformed / out_factor).astype(np.int32)
                        self.y = (y_transformed / out_factor).astype(np.int32)
                    elif CRS(epsg_out).is_geographic:
                        self.y = (3600 * x_transformed / out_factor).astype(np.int32)
                        self.x = (3600 * y_transformed / out_factor).astype(np.int32)
                    self.crs = CRS(epsg_out)

                except:
                    raise ValueError(
                        "Coordinate transformation failed, check if the epsg codes are correct or set a different out_scalar"
                    )
        else: # apply_new_coordinates == False
            return x_transformed, y_transformed

    def calc_max_dist(self, maximum_dist):
        """Get the input maximum distance in meters or arcseconds
        in the format of the coordinates in the trace haders"""
        if self.scalar < 0:
            max_dist = int(abs(maximum_dist * self.scalar))
        elif self.scalar > 0:
            max_dist = int(abs(maximum_dist / self.scalar))

        return max_dist

    def fix_navigation_median(
        self,
        order=501,
        maximum_dist=30,
        truncate_data=0,
        show_figures=True,
        save_figures=False,
    ):
        """
        Discard and interpolate coordinate points that are more than max_dist meters away
        from a strongly median filtered version of the same points

        Parameters
        ----------
        order: n-th order median filter that will be applied to the data. Distance of the original
        data to the filtered output will determine if points will be discarded as outliers
        truncate_data: int. Removes points from the start and the end of line if peaks occur in one and/or the other
        show_figures: displays figures of original and the corrected coordinates, default=True
        save_figures: saves figures of original and the corrected coordinates, default=False

        Returns
        -------
        self.x : numpy list of corrected x coordinates \n
        self.y : numpy list of corrected y coordinates \n

        """

        max_dist = self.calc_max_dist(maximum_dist)

        nav_path = self.segy_filepath.with_stem(
            self.segy_filepath.stem + "_NAVIGATION"
        ).with_suffix(".png")

        if "_TEMPORARY" in str(nav_path):
            nav_string_path = nav_path.name.replace("_TEMPORARY", "")
            nav_path = nav_path.with_name(nav_string_path)

        if save_figures == True:
            self.x, self.y = remove_outliers_median(
                self.x,
                self.y,
                order,
                max_dist,
                truncate_data,
                show_figures,
                save_path=nav_path,
            )
        else:
            self.x, self.y = remove_outliers_median(
                self.x, self.y, order, max_dist, truncate_data, show_figures
            )

        self.x = np.array(self.x).astype(np.int32)
        self.y = np.array(self.y).astype(np.int32)

    def fix_navigation_polyfit(
        self, degree_polyfit=3, maximum_dist=30, truncate_data=0, show_figures=True
    ):
        """
        Discard and interpolate coordinate points that are more than max_dist meters away
        from a curve fitted in a least squares sense to the original points

        Parameters
        ----------
        degree: degree of polynomial fitting curve. Distance of the original
        data to the fitted curve will determine if points will be discarded as outliers \n
        truncate_data: int. Removes points from the start and the end of line if peaks \n
        occur in one and/or the other. Default = 0
        show_figures: displays figures of original and the corrected coordinates, default=True

        Returns
        -------
        self.x : numpy list of corrected x coordinates \n
        self.y : numpy list of corrected y coordinates \n

        """

        max_dist = self.calc_max_dist(maximum_dist)
        self.x, self.y = remove_outliers_polyfit(
            self.x, self.y, degree_polyfit, max_dist, truncate_data, show_figures
        )

        self.x = np.array(self.x).astype(np.int32)
        self.y = np.array(self.y).astype(np.int32)

    def fix_navigation_smooth(self, order=501, truncate_data=0, show_figures=True):
        """
        Smooth the coordinate point along a survey line

        Parameters
        ----------
        x : numpy list of faulty x coordinates \n
        y : numpy list of faulty y coordinates \n
        order: n-th order median filter that will be applied to the data \n
        truncate_data: int. Removes points from the start and the end of line if peaks occur in one and/or the other
        show_figures: displays figures of original and the corrected coordinates, default=True

        Returns
        -------
        self.x : numpy list of corrected x coordinates \n
        self.y : numpy list of corrected y coordinates \n


        """

        self.x, self.y = smoothXY(self.x, self.y, order, truncate_data, show_figures)

        self.x = np.array(self.x).astype(np.int32)
        self.y = np.array(self.y).astype(np.int32)

    def fix_navigation_replace(
        self, degree_polyfit=3, truncate_data=0, show_figures=True
    ):
        """
        Replace the original navigation data with a polynomial fit of that data.
        This is a last resort to be able to get some sensible line navigation,
        but it does not represent true navigation

        Parameters
        ----------
        x : numpy list of faulty x coordinates \n
        y : numpy list of faulty y coordinates \n
        degree: degree of polynomial fit that will replace the navigation data \n
        truncate_data: int. Removes points from the start and the end of line if peaks occur in one and/or the other
        show_figures: displays figures of original and the corrected coordinates, default=True

        Returns
        -------
        self.x : numpy list of replaced x coordinates \n
        self.y : numpy list of replaced y coordinates \n

        """

        self.x, self.y = replace_with_polyfit(
            self.x, self.y, degree_polyfit, truncate_data, show_figures
        )

        self.x = np.array(self.x).astype(np.int32)
        self.y = np.array(self.y).astype(np.int32)

    def fix_navigation_remove_startend(self,max_distance=200,max_truncate_length=20):
        """
        Removes traces at start and/or end of line with bad coordinates. 

        Parameters
        ----------
        max_distance : int or float
            maximum distance between two traces in real units (same as coordinate system)
        max_truncate_length : int
            maximum number of traces at start and end to find bad coordinates in 
            (outliers within line should be removed with e.g. fix_navigation_median() function)

        Returns
        -------
        Changes Segy_edit object without leading and trailing bad coordinates
        """
        x, y = self.xy_to_real()
        n_bad_traces_start, n_bad_traces_end = remove_start_end(x,y,max_distance,max_truncate_length)

        if (n_bad_traces_start>0) or (n_bad_traces_end>0):
            print(f"Removing {n_bad_traces_start} traces from start of line and/or {n_bad_traces_end} traces from end of line")

            # change coordinates, shotpoint and other general info
            self.x = self.x[n_bad_traces_start:self.trace_number-n_bad_traces_end]
            self.y = self.y[n_bad_traces_start:self.trace_number-n_bad_traces_end]
            self.trace_data = self.trace_data[n_bad_traces_start:self.trace_number-n_bad_traces_end]

            self.renumber_shotpoints(0)
            self.trace_number = len(self.trace_data)
            self.spec.tracecount = len(self.trace_data)
        else:
            print("No traces found to remove")

    def renumber_shotpoints(self, start_shotpoint=1001):
        """renumber the shotpoint numbers is byte 17 and 21 to be sequential
        you can define your starting shotpoint (default=1001)"""

        self.shotpoint_numbers = np.array(
            (np.arange(self.trace_number) + start_shotpoint), dtype="int32"
        )

    def set_trace_number_in_field_record(self, num=1):
        """set the trace number in field when to go e.g after extracting near-trace
        from multi-channel data or after decimating the new 3 trace x-star data.
        The default value is 1"""

        self.trace_number_in_field_record = num

    def copy_source_coordinates_to_group(self):
        """Copy the coordinates found in the source postion bytes over to the
        group/receiver position bytes"""

        self.groupx = self.x
        self.groupy = self.y

    def read_grid(
        self, grid_path, epsg_grid, horizon_name, velocity=1500
    ):
        """
        Retrieves depths and two-way traveltimes corresponding the coordinates
        of the seismic navigation line.

        Parameters
        ----------
        grid_path : pathlib.Path object
            path to the input ascii grid
        epsg_grid : int
            epsg code of the input ascii grid
        velocity : int, optional
            The p-wave velocity to be used to convert the horizon to two-way time.
            The default is 1500 m/s

        """

        # transform the seismic navigation coordinates to match that of the grid_file (fastest, because less data point)
        x_transformed, y_transformed = self.transform_coordinates(
            epsg_grid,
            manual_in_scalar=self.out_scalar,
            apply_new_coordinates=False,
        )

        # extrect horizon depth from the grid associated with each seismic navigation line point, also calculated horizon two-way time
        horizon_depths, horizon_two_way_times = get_grid(
            grid_path, x_transformed, y_transformed, velocity
        )

        setattr(self, horizon_name, horizon_two_way_times)

        self.add_to_horizons(horizon_name)

    def get_seabed_pick(
        self,
        nsta,
        nlta,
        threshold,
        peak_threshold,
        pick_peaks=True,
        start_window=0,
        end_window=0,
        truncate=0,
    ):
        """
        Pick the seabed in the seismic sections

        Parameters
        ----------
        trigger_threshold : float
            Trigger threshold (amplitude ratio)
        sta : float
            Short-term average window (in milliseconds)
        lta : float
            Long-term average window (in milliseconds)
            
        """

        start_window_sample = int(start_window / (self.sampling_interval * 1e-3))

        # define window to look for seabed pick
        if end_window == 0:
            end_window_sample = self.trace_samples
        else:
            end_window_sample = int(end_window / (self.sampling_interval * 1e-3))

        # cast data in obspy readible format
        data = np.array(self.trace_data)[:, start_window_sample:end_window_sample]

        # convert to obspy stream object
        st = deltaseis_to_obspy(data.T, self.sampling_rate)

        onsets = []
        peaks = []

        for i, tr in enumerate(st):
            sta_lta = classic_sta_lta(tr, nsta=nsta, nlta=nlta)
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=RuntimeWarning)
                rms = np.sqrt(np.mean(np.square(tr)))
            try:
                onset_samp = np.where(sta_lta > threshold)[0][0]
                peak_samps = find_peaks(tr[onset_samp:], peak_threshold * rms)[0]
                peak_samp = peak_samps[0] + onset_samp

                onsets.append(onset_samp)
                peaks = np.append(peaks, peak_samp)

            except IndexError:
                onsets.append(np.nan)
                peaks = np.append(peaks, np.nan)

        # add constant value start_window_sample to peaks
        peaks = peaks + start_window_sample

        # onset or peak of seabed pick
        if pick_peaks == True:
            seabed_pick = peaks * self.sampling_interval * 1e-3
        else:
            seabed_pick = np.array(onsets) * self.sampling_interval * 1e-3

        if truncate != 0:
            self.seabed_pick = self.truncate_non_nan_signal(seabed_pick, truncate)
        else:
            self.seabed_pick = seabed_pick

        self.add_to_horizons("seabed_pick")

    def truncate_non_nan_signal(self, horizon, truncate_samples=0):
        """truncate the non-nan part of a vector and set to nan to avoid edge effects
        e.g. from seabed picking"""

        # Find the index of the first non-NaN value and the last non-nan
        first_non_nan_idx = np.where(~np.isnan(horizon))[0][0]
        last_non_nan_idx = np.where(~np.isnan(horizon))[0][-1]

        horizon[first_non_nan_idx : first_non_nan_idx + truncate_samples] = np.nan
        horizon[last_non_nan_idx - truncate_samples : last_non_nan_idx] = np.nan

        return horizon

    def calculate_difference_horizon(
        self, horizon1_name, horizon2_name, difference_horizon_name
    ):
        """
        Calculate the difference between 2 horizons. E.g. to get the tide by calculating the difference
        between the seabed pick and a time converted bathymetry horizon. Or to determine (two-way time)
        thickness of units (isopachs).

        Parameters
        ----------
        horizon1_name : (string) first input horizon name
        horizon2_name : (string) second input horizon name
        difference_horizon_name : (string) resulting difference horizon name

        Returns
        -------
        The difference horizon name is stored as an attribute of the object

        """

        if hasattr(self, horizon1_name) and hasattr(self, horizon2_name):
            setattr(
                self,
                difference_horizon_name,
                getattr(self, horizon1_name) - getattr(self, horizon2_name),
            )

        else:
            print(
                f"{horizon1_name} and/or {horizon2_name} does not exist, {difference_horizon_name} not calculated"
            )

        self.add_to_horizons(difference_horizon_name)

    def handle_nans(self, horizon):
        """Interpolate nans values in a horizon, the prepare it for filtering, set nan at the start
        and end of horizon equal to the (next to) non-nan edge value"""

        # find nans in horizon
        horizon_no_nans = np.copy(horizon)
        nan_indices = np.isnan(horizon_no_nans)

        # Create an array of indices for non-NaN values
        non_nan_indices = np.arange(len(horizon_no_nans))

        # replace nans at the start and the end with the (next to) edge value
        # this avoids extrapolation edge effects when attempting to interpolate the nans

        first_non_nan_idx = np.where(~np.isnan(horizon_no_nans))[0][1]
        last_non_nan_idx = np.where(~np.isnan(horizon_no_nans))[0][-2]

        horizon_no_nans[: first_non_nan_idx - 1] = horizon_no_nans[first_non_nan_idx]
        horizon_no_nans[last_non_nan_idx + 1 :] = horizon_no_nans[last_non_nan_idx]

        # now interpolate the internal nans
        horizon_no_nans[nan_indices] = np.interp(
            non_nan_indices[nan_indices],
            non_nan_indices[~nan_indices],
            horizon_no_nans[~nan_indices],
        )

        return horizon_no_nans

    def filter_horizon_savgol(
        self, horizon_name, savgol_horizon_name, window_length, polyorder
    ):
        """Applies a smoothing filter to a horizon, this can be used to correct for static
        that are for instance caused by heave, but I might also filter seabed features like
        mega ripples, so use with care"""

        horizon = np.copy(getattr(self, horizon_name))
        horizon = self.handle_nans(horizon)

        savgol_horizon = savgol_filter(horizon, window_length, polyorder, mode="mirror")

        setattr(self, savgol_horizon_name, savgol_horizon)
        self.add_to_horizons(savgol_horizon_name)


    def add_to_horizons(self, horizon_name):
        """Checks if horizon exists and gives a warning when that it overwrites. If it
        does not exist, it is added to the horizons dictionary"""

        if horizon_name in self.horizons:
            print(f"Warning: overwriting existing {horizon_name} horizon")
            self.horizons[horizon_name] = getattr(self, horizon_name)
        else:
            self.horizons[horizon_name] = getattr(self, horizon_name)

    def vertical_trace_corrections(self, corrections):
        """
        Applies vertical corrections to traces by:
        Padding zeros at the start of the trace and cutting the same amount of values from the end, to shift down
        Cutting values at the start of the trace and padding the same amojunt of zeros to the end, to shift up

        Parameters
        ----------
        shifts : 1-D array
            Contains the shift in milliseconds required for each trace. The shifts vector
            should have the same length as the amount of traces. Shift can
            have different values for each trace, or the same for all traces (bulk shift)

        """

        shifts = np.copy(corrections)
        shifts[np.isnan(shifts)] = 0

        shift_samples = np.array(shifts / (1e-3 * self.sampling_interval)).astype(
            int
        )  # in millisecond

        trace_data = np.array(self.trace_data, dtype=self.data_sample_format)

        for i, tr in enumerate(trace_data):
            if shift_samples[i] > 0:
                trace_data[i] = np.pad(
                    trace_data[i], (shift_samples[i], 0), "constant"
                )[: -shift_samples[i]]
            elif shift_samples[i] < 0:
                trace_data[i] = np.pad(
                    trace_data[i], (0, abs(shift_samples[i])), "constant"
                )[-shift_samples[i] :]

        self.trace_data = trace_data.tolist()

    def plot(
        self,
        cmap="Greys",
        clip=0.8,
        fontsize=16,
        show_horizons=True,
        save_plot=False,
        show_time_axis=True,
    ):
        """Quick plot of seismic section, optionally plots imported or picked horizons on top if they
        exist (default=True), optionally saves the plot (default=False)"""

        if show_time_axis == True:
            extent = [
                0,
                self.trace_number,
                self.trace_samples * (self.sampling_interval * 1e-3),
                0,
            ]
            ylabel = "Two-way time (ms)"
        else:
            extent = [0, self.trace_number, self.trace_samples, 0]
            ylabel = "Sample number"

        plt.figure(figsize=(20, 8))
        plt.title(self.segy_filepath.name, fontsize=fontsize)
        vmin = clip * np.min(self.trace_data)
        vmax = clip * np.max(self.trace_data)
        if vmin<0: #if data is full wave form
            vmin = -clip * np.mean(np.abs(self.trace_data))
            vmax = clip * np.mean(np.abs(self.trace_data))
            cmap = 'bwr'
        im = plt.imshow(
            np.array(self.trace_data).T,
            cmap=cmap,
            aspect="auto",
            extent=extent,
            vmin=vmin,
            vmax=vmax,
        )

        # Increase font size for axis labels
        plt.xlabel("Trace number", fontsize=fontsize)
        plt.ylabel(ylabel=ylabel, fontsize=fontsize)
        # ax2.set_ylabel('Two-way time (ms)', fontsize=fontsize)
        # ax2.set_ylim(0, self.trace_samples*(self.sampling_interval*1e-3))

        # Increase font size for axis tick labels
        plt.xticks(fontsize=fontsize)
        plt.yticks(fontsize=fontsize)

        cbar = plt.colorbar(im)
        cbar.set_label("Amplitudes", fontsize=fontsize)
        cbar.ax.tick_params(labelsize=fontsize)

        # set colorbar values alwasys as exponents
        cbar.formatter.set_powerlimits((0, 0))

        # Increase the font size for the exponent text
        cbar.ax.yaxis.get_offset_text().set_fontsize(fontsize)

        # plot existing (difference) horizons and seabed pick
        if show_horizons == True:
            colors = cycle(
                ["red", "darkorange", "yellow", "gold", "limegreen", "springgreen"]
            )

            for horizon in self.horizons:
                if show_time_axis == True:
                    horizon_two_way_time = getattr(self, horizon)
                    plt.plot(
                        horizon_two_way_time,
                        color=next(colors),
                        linewidth=1,
                        label=horizon,
                    )
                else:
                    horizon_samples = getattr(self, horizon) / (
                        self.sampling_interval * 1e-3
                    )
                    plt.plot(
                        horizon_samples, color=next(colors), linewidth=1, label=horizon
                    )

            plt.legend(fontsize=fontsize, loc="lower right")

        plt.tight_layout()

        if save_plot == True:
            self.plot_number += 1
            plot_outfile = self.segy_filepath.with_stem(
                self.segy_filepath.stem + f"_{self.plot_number}"
            ).with_suffix(".png")

            if "_TEMPORARY" in str(plot_outfile):
                string_plot_outfile = plot_outfile.name.replace("_TEMPORARY", "")
                plot_outfile = plot_outfile.with_name(string_plot_outfile)

            plt.savefig(plot_outfile, dpi=300)

    def plot_nav(self):
        """
        Simple plot of (real) coordinates for a quick check

        Parameters
        ----------
        - none
        """
        x, y = self.xy_to_real()
        
        plt.plot(x,y)
        if self.crs.is_projected:
            plt.xlabel('x coordinate')
            plt.ylabel('y coordinate')
        else: 
            plt.xlabel('longitude')
            plt.ylabel('latitude')
        plt.grid()
        plt.axis('equal')

    def write(self, segy_outpath):
        # cast trace data to correct format
        self.trace_data = np.array(self.trace_data).astype(
            np.dtype(self.data_sample_format)
        )

        with segyio.create(segy_outpath, self.spec) as dst:
            dst.text[0] = self.text_header
            dst.bin = self.bin_header
            dst.header = self.trace_header

            dst.bin[3221] = len(self.spec.samples)
            dst.bin[3213] = self.trace_number_in_field_record
            dst.bin[3217] = self.sampling_interval

            for i in range(self.trace_number):
                dst.header[i][13] = self.trace_number_in_field_record
                dst.header[i][17] = self.shotpoint_numbers[i]
                dst.header[i][71] = self.out_scalar
                dst.header[i][73] = self.x[i]
                dst.header[i][77] = self.y[i]
                dst.header[i][81] = self.groupx[i]
                dst.header[i][85] = self.groupy[i]
                dst.trace[i] = self.trace_data[i]

            if "_TEMPORARY" in str(self.segy_filepath):
                self.src.close()
                if self.segy_filepath.exists():
                    self.segy_filepath.unlink()
            else:
                self.src.close()

    def xstar_split(self, mode="envelope"):
        """
        Produces a segy containing only the analytical (envelope), or the full signal in xstar data

        Parameters
        ----------
        mode : string
            'envelope' (default) or 'full'
        """

        if mode == "envelope":
            self.trace_data = self.trace_data[::3]
        elif mode == "full":
            new_trace_data = []
            for i in range(int(self.trace_number / 3)):
                im = self.trace_data[i * 3 + 1]
                re = self.trace_data[i * 3 + 2]

                # Perform Fourier transform on the real and imaginary signals
                fft_real = np.fft.fft(re)
                fft_imag = np.fft.fft(im)

                # Combine in the frequency domain
                combined_frequency = fft_real + 1j * fft_imag

                # Compute the inverse Fourier transform to get the time-domain signal
                combined_waveform = np.real(np.fft.ifft(combined_frequency))

                new_trace_data.append(combined_waveform)

            self.trace_data = new_trace_data

        # change coordinates, shotpoint and other general info
        self.x = self.x[::3]
        self.y = self.y[::3]
        self.renumber_shotpoints()
        self.trace_number_in_field_record = 1  # set number of traces to 1
        self.trace_number = int(self.trace_number / 3)
        self.spec.tracecount = int(self.spec.tracecount / 3)

    def full2envelope(self):
            """
            Transforms full waveform data to envelope data using Hilberts transformation

            """
            d = self.trace_data
            if np.all(np.stack(d)>0):
                print("""WARNING: ALL DATA VALUES ARE POSITIVE, this does not look like full wave form data, DATA IS NOT TRANSFORMED!!""")
            else:
                for i in range(len(d)):
                    d[i] = np.abs(hilbert(d[i]))
                self.trace_data = d

    def extract_near_trace_gather(self, number_of_channels, channel=1):
        """
        Produces the near trace gather from pre-stack seg-y data. 

        Parameters
        ----------
        number_of_channels: int
            The number of channels in the streamer or spread
        channel : int
            The number of the channel for which you would like to extract the gather 
            (default=1, which is the nearest trace)
        """

        self.trace_data = self.trace_data[channel -1::number_of_channels]

        # change coordinates, shotpoint and other general info
        self.x = self.x[channel - 1::number_of_channels - 1]
        self.y = self.y[channel - 1::number_of_channels - 1]
        self.renumber_shotpoints()
        self.trace_number_in_field_record = 1  # set number of traces to 1
        self.trace_number = int(self.trace_number / number_of_channels)
        self.spec.tracecount = int(self.spec.tracecount / number_of_channels)

    def remove_empty_traces(self):
        """
        Removes traces with only 0 values from line

        """
        traces_with_data = [
            i for i in range(self.trace_number) if np.any(self.trace_data[i])
        ]

        if len(traces_with_data) < self.trace_number:
            print(f"removing {self.trace_number - len(traces_with_data)} empty traces")

            # change coordinates, shotpoint and other general info
            self.x = self.x[traces_with_data]
            self.y = self.y[traces_with_data]

            self.trace_data = [self.trace_data[i] for i in traces_with_data]

            self.renumber_shotpoints(0)
            self.trace_number = len(traces_with_data)
            self.spec.tracecount = len(traces_with_data)
        else:
            print("no empty traces to remove")

    def to_geodataframe(self, geomtype="line", linename="", linename_prefix=""):
        """
        Returns a geodataframe containing the coordinates of the segy, either as a line
        or as points

        Parameters
        ----------
        geomtype : string
            'line' (default) or 'point'
            geometry type, returning one 'Linestring' object for a single segy, or 
            ping-coordinates as points (number of points depending on thinning factor)
        linename : string
            string to be used in the attribute 'linename'. if not given,
            the filename without extention will be used
        linename_prefix : string
            string added as prefix to the linename, e.g. 'IJmuidenGronden_',
            only used if linename is empty

        """

        x, y = self.xy_to_real()

        # create string for linename
        if not linename:
            linename = linename_prefix + self.segy_filepath.stem

        # create dataframe with linename and geometry
        if geomtype == "point":
            geom = gpd.points_from_xy(x, y)
            df = pd.DataFrame(
                {
                    "linename": [linename for i in range(len(geom))],
                    "geometry": geom,
                }
            )
        else:
            geom = LineString([Point(xyi) for xyi in zip(x, y)])
            df = pd.DataFrame({"linename": [linename], "geometry": [geom]})

        # create geopandas dataframe
        gdf = gpd.GeoDataFrame(df, geometry="geometry", crs=self.crs)

        return gdf
