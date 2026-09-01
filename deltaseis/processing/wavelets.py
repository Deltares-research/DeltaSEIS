from pathlib import Path

import numpy as np

from deltaseis.base_seismic import Seismic
from deltaseis.tools.segy_editor import Segy_edit


def extract_wavelet_from_segy(segy_file, trace_number, start_time_ms,
                              end_time_ms, taper=0.2):
    """Extract a signature wavelet from a reference SEG-Y file."""
    print(f"Extracting signature wavelet from {Path(segy_file).name}")
    edit = Segy_edit(segy_file)
    data = np.array(edit.trace_data).T
    dx_mean = edit.factor * edit.shot_point_interval.mean()
    seismic = Seismic(data, edit.sampling_rate, dx_mean)
    wavelet = seismic.extract_wavelet(
        trace_number, start_time_ms, end_time_ms, taper=taper
    )
    return wavelet, seismic.fs
