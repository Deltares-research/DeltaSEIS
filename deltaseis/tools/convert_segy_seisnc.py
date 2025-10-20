import importlib.resources
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from tkinter import Tk, filedialog

import segysak
import xarray as xr

from deltaseis import config

#notes: 
'''
- perhaps a set of jsons if required for loading different types of data (e.g. x-star, pre-stack, tdms etc)
- perhaps the jsons also need a bin_header bytes dictonary
- perhaps ths jsons also need to include the text header somehow? or could be loaded using scrape, find a way to edit this easily (ascii not ebcdic)
'''

def convert_segy_to_seisnc(segy_file_path, config_path=None):

    if config_path is None:
        # Default to the bundled resource in the package
        resource_path = importlib.resources.files(config) / 'seismic_headers.json'
    else:
        # If a custom path is provided, use it
        resource_path = Path(config_path)
    
    # Load SEG-Y header fields from the JSON config file
    with resource_path.open('r') as f:
        fields = json.load(f)

    # Load SEG-Y file with xarray
    ds = xr.open_dataset(
         segy_file_path,
         dim_byte_fields = fields["dim_bytes"],
         extra_byte_fields = fields["trace_header_bytes"],
         segyio_kwargs = {'endian': 'little'}
     )
    
        
    # Scale coordinates with SEG-Y coordinate scalar
    ds.segysak.scale_coords(-100)
    
    # Write to SeisNC format
    output_file = Path(segy_file_path).with_suffix('.seisnc')
    
    ds.to_netcdf(output_file)
    print(f"Converted {segy_file_path} to {output_file}")

def main():
    # Hide the root window
    root = Tk()
    root.withdraw()
    # Open file dialog to select multiple SEGY files
    segy_file_paths = filedialog.askopenfilenames(title="Select SEGY files", filetypes=[("SEGY files", "*.sgy *.segy")])
 
    if segy_file_paths:
        # Use ThreadPoolExecutor to process files in parallel
        with ThreadPoolExecutor() as executor:
            executor.map(convert_segy_to_seisnc, segy_file_paths)

if __name__ == "__main__":
    main()