import pandas as pd
import re
import warnings
# -*- coding: utf-8 -*-
"""
This script reads CSV files in the Digilent WaveForms Oscilloscope format (e.g., acq0010.csv).
It extracts the sample rate from the header (line starting with '#Sample rate: ...'),
then reads the time series data starting at row index 10 (after the header).

Returns the data as a pandas DataFrame, the sample rate (fs), and the start time.
"""
import pandas as pd
import re
import warnings

def read_waveforms_csv(file_list, verbose=True):
    """
    Read one or more Digilent WaveForms Oscilloscope CSV files and concatenate their time series.

    Parameters
    ----------
    file_list : list of str or Path
        List of CSV file paths (e.g., [acq0010.csv, acq0011.csv, ...])
    verbose : bool
        Print extracted sample rate and QC info

    Returns
    -------
    data : pandas.DataFrame
        DataFrame with columns ['Time (s)', 'Channel 1 (V)', ...] (one column per file)
    fs : float
        Sample rate in Hz (must be identical for all files)
    t_start : float
        Starting time in seconds (warning if not identical)
    """

    trace_list = []
    fs_list = []
    t_start_list = []

    for i, filepath in enumerate(file_list):
        # Read header lines
        header_lines = []
        with open(filepath, 'r') as f:
            for j in range(10):
                line = f.readline()
                header_lines.append(line)

        # Extract sample rate
        fs = None
        for line in header_lines:
            m = re.match(r"#Sample rate:\s*([\deE\+\.]+)Hz", line)
            if m:
                fs = float(m.group(1))
                break
        if fs is None:
            raise ValueError(f"Sample rate not found in header of file {filepath}.")
        fs_list.append(fs)

        # Read time series data (skip first 9 rows)
        data = pd.read_csv(filepath, skiprows=9)
        t_start = data.iloc[0, 0]
        t_start_list.append(t_start)

        # Only keep the voltage column
        trace_list.append(data.iloc[:, 1].reset_index(drop=True))

        if verbose:
            print(f"{i+1}/{len(file_list)}: {filepath}")
            print(f"  Sample rate: {fs} Hz, Start time: {t_start} s, Samples: {len(data)}")

    # Check sample rate consistency
    if len(set(fs_list)) > 1:
        raise ValueError(f"Sample rate (fs) differs between files: {fs_list}")
    fs = fs_list[0]

    # Check t_start consistency
    if len(set(t_start_list)) > 1:
        t_start_ms = [float(x) * 1000 for x in t_start_list]
        min_t = min(t_start_ms)
        max_t = max(t_start_ms)
        if abs(max_t - min_t) > 0.01:
            warnings.warn(f"Start time (t_start) differs between files: min={min_t:.6f} ms, max={max_t:.6f} ms")
    t_start = t_start_list[0]

    # Concatenate traces into DataFrame
    data = pd.concat(trace_list, axis=1)
    data.columns = [f"Channel_{i+1}" for i in range(len(trace_list))]

    # Add time column (from first file)
    time_data = pd.read_csv(file_list[0], skiprows=9).iloc[:, 0].reset_index(drop=True)
    data.insert(0, "Time (s)", time_data)

    return data, fs, t_start
