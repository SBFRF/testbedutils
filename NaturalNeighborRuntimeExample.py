# this code makes natural neighbor grids
# cliff change this to where the FRF_natneighbor.py lives if not in same folder
sys.path.append('/home/spike/repos/frf_bathy_interp')
import FRF_natneighbor as nn
import getopt, sys, glob, makenc
import netCDF4 as nc
import sblib as sb


globalYaml = '/home/spike/CMTB/yaml_files/Survey/Transect_Survey_LARC_Global.yml' # location of yaml files for transect data
varYaml = 'yaml_files/Survey/Transect_Survey_variables.yml' # location of yaml files for transect data

grid_global = 'yaml_files/Survey/grid_Global.yml'  # location of yaml files for gridded data
grid_var = 'yaml_files/Survey/grid_variables.yml'  # location of yaml files for gridded data


filelist = glob.glob('/home/spike/repos/frf_bathy_interp/data')

for fname in filelist:
    try:
        os.remove(file[:-4] + '.nc')
        os.remove(file[:-4] + '_grid.nc')
    except OSError:
        continue
    # first make transect
    TransectDict = sb.load_FRF_Transect(fname)
    ofname = fname[:-3] + 'nc'
    TransectDict['time'] = nc.date2num(TransectDict['time'], 'seconds since 1970-01-01')
    makenc.makenc_FRFTransect(bathyDict=TransectDict, ofname=ofname, globalYaml=globalYaml, varYaml=varYaml)


gridList = glob.glob()
for fname in gridList:
    # second make grid
    gridDict = sb.importFRFgrid(fname)

    gridDict = nn.frf_grid_product(fname, dxdy=10)
    ofname = fname[:-4]+'_grid.txt'
    nn.write_grid(ofname=ofname, grid_dict=gridDict)
    # put plotting functing that saves file here

if __name__ == "__main__":
    opts, args = getopt.getopt(sys.argv[1:], "h", ["help"])
    #.
    #
    fname = args[0]