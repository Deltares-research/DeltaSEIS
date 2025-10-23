from deltaseis import Segy_edit
from deltaseis.tools.qc import analyze_seismic_data
from pathlib import Path
import numpy as np


# Use raw string (r"") to avoid escape sequence issues
segy_folder = Path(r'D:\Projects\DIS_Offshore\qc\example_sparker_data_texel')

qc_folder = segy_folder / 'qc_reports'
qc_folder.mkdir(parents=True, exist_ok=True)

# Check if segy folder exists
if not segy_folder.exists():
    print(f"Warning: SEGY folder {segy_folder} does not exist!")
    print("Using local test data instead...")
    segy_folder = Path('deltaseis/data')

# List all SEGY files first
segy_files = list(Path(segy_folder).glob('*.sgy'))
print(f"Found {len(segy_files)} SEGY files in {segy_folder}")

for i, segy_file in enumerate(segy_files):
    print(f"\nFile {i+1}/{len(segy_files)}: {segy_file.name}")
    print(f"File size: {segy_file.stat().st_size / (1024*1024):.1f} MB")
    
    # only load segy files when qc pdf files do not exist
    if not any(qc_folder.glob(f"{segy_file.stem}_*.pdf")):
        print(f"Processing {segy_file.name}...")

        # Load SEGY file with progress indication
        print("Loading SEGY file...")
        segy = Segy_edit(segy_file)
        print(f"SEGY loaded: {segy.trace_number} traces, {segy.trace_samples} samples per trace")
        
        print("Converting to numpy array...")
        data = np.array(segy.trace_data).T
        print(f"Data shape: {data.shape} ({data.nbytes / (1024*1024):.1f} MB)")
        
        fs = segy.sampling_rate         # sampling frequency
        signal_window = (500, 1500)     # example signal window
        noise_window = (0, 400)         # example noise window
        
        # Use fast analysis with trace sampling
        analyze_seismic_data(data, fs, signal_window, noise_window,
                             output_pdf=qc_folder / f"{segy_file.stem}_qc_report.pdf",
                             max_traces=50)  # Limit to 50 traces for speed
    else:
        print(f"Skipping {segy_file.name} - QC report already exists")