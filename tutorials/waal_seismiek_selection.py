'''
NEXT STEPS:
- clean segy folder, naming convention
- check the data processing steps
- get shape file of segy coordinates
- find the crosspoints
- find the closest point
- calculate distance
- find the bottleneck for processing (bandpass filter is probably the slowest and not really needed)
- seabedpick and and correct with bathy (using a netcdf file, meaning the get_bathy function should be updated)
- trace selection based on 4 crosspoints (set manually by plotting the cross section with the qgis exports of the tracks)
- find closest point, calculate distance of all points from THAT specific point. Select traces based on distance < x e.g 50 on both sides
- conversion of all the raw data to make a clean repo would be an excellent idea (raw and processed)
'''	