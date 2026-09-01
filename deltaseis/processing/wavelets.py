import json
from pathlib import Path

import numpy as np

from deltaseis.base_seismic import Seismic
from deltaseis.tools.segy_editor import Segy_edit


def load_or_extract_wavelet(wavelet_file, segy_file, trace_number,
                            start_time_ms, end_time_ms, taper=0.2):
    """Load a cached wavelet or extract it from a reference SEG-Y file."""
    wavelet_file = Path(wavelet_file)
    metadata_file = wavelet_file.with_suffix(".json")

    if not wavelet_file.exists() or not metadata_file.exists():
        print(f"Extracting signature wavelet from {Path(segy_file).name}")
        edit = Segy_edit(segy_file)
        data = np.array(edit.trace_data).T
        dx_mean = edit.factor * edit.shot_point_interval.mean()
        seismic = Seismic(data, edit.sampling_rate, dx_mean)
        wavelet = seismic.extract_wavelet(
            trace_number, start_time_ms, end_time_ms, taper=taper
        )
        np.save(wavelet_file, wavelet)
        metadata_file.write_text(json.dumps({
            "file": Path(segy_file).name,
            "trace_number": trace_number,
            "start_time_ms": start_time_ms,
            "end_time_ms": end_time_ms,
            "fs": float(seismic.fs),
        }, indent=2))

    wavelet = np.load(wavelet_file)
    wavelet_fs = json.loads(metadata_file.read_text())["fs"]
    return wavelet, wavelet_fs
