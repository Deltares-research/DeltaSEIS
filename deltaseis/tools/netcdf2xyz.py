import tkinter as tk
from multiprocessing import Pool
from tkinter import filedialog

import numpy as np
import xarray as xr


def netcdf_to_xyz(netcdf_file, output_file, variable_name, no_value=-9999):
    # Open the NetCDF file using xarray
    dataset = xr.open_dataset(netcdf_file)

    # Extract the variable data
    variable_data = dataset[variable_name].values

    # Replace NaNs with no_value
    variable_data = np.nan_to_num(variable_data, nan=no_value)

    # Extract the coordinates
    y = dataset['y'].values
    x = dataset['x'].values

    # Open the output file
    with open(output_file, 'w') as f:
        # Loop over the data and write to the file
        for i in range(len(x)):
            for j in range(len(y)):
                f.write(f"{y[j]} {x[i]} {variable_data[j, i]}\n")

    # Close the dataset
    dataset.close()


def process_file(file_path):
    output_file = file_path.replace('.nc', '.xyz')
    variable_name = '__xarray_dataarray_variable__'  # set to actual variable name
    netcdf_to_xyz(file_path, output_file, variable_name)

def main():
    root = tk.Tk()
    root.withdraw()  # Hide the root window
    file_paths = filedialog.askopenfilenames(
        title="Select NetCDF files",
        filetypes=[("NetCDF files", "*.nc")]
        )

    if file_paths:
        with Pool() as pool:
            pool.map(process_file, file_paths)

if __name__ == "__main__":
    main()
