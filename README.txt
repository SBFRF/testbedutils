this group of files creates netCDF's from the traditional .txt (grid) and .csv (transect) data.
written by Spicer Bak, 10/11/16 
to run from command line:
python surveyToNetCDF.py fname

where fname can be either a .txt (grid) or a .csv (transect) file  This script will 
also create a log called 'Bathy_LOG.log' placed in the "path" file location.  This script assumes a folder located 
parallel to the script that there is a folder "yamlFiles" that contains the following files (needed for netCDF metadata)

grid_Global.yml         # grid global yaml for metadata for netCDF files 
grid_variables.yml      # variables meta data for netCDF files 
transect_Global.yml     # global yaml for the transect netCDF files
transect_variables.yml  # transect meta data for the variables in the netCDF files  