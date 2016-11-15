import sys, getopt, os
sys.path.append('/home/spike/CMTB/')
import makenc, glob
from matplotlib import pyplot as plt
import netCDF4 as nc
import numpy as np
import sblib as sb

def convertText2NetCDF(globPath):
    """
    This function searches the given path for both grids and transect files present at the FRF to make the data into
    netCDF files
    the below yaml's have to be in the same folder as the
    :param globPath: a path to find the yaml files and the data
    :return:
    """
    ## INPUTS
    yamlPath = 'yamlFiles/'
    gridGlobalYaml = yamlPath + 'grid_Global.yml'  # grid Global yaml
    gridVarYaml = yamlPath + 'grid_variables.yml'  # grid yaml location
    transectGlobalYaml = yamlPath + 'transect_Global.yml'
    transectVarYaml = yamlPath + 'transect_variables.yml'

    gridList = glob.glob(globPath + 'FRF_*latlon.txt') # searching for grid files
    filelist = glob.glob(globPath + 'FRF_*.csv') # searching for transect Files

    logFile = globPath + 'Bathy_LOG.log'

    errorFname, errors = [],[]
    # creating a list of transect files to look for
    print 'Converting %d Transect (csv) file to netCDF ' % len(filelist)
    for transectFname in filelist:
        try:
            Tofname = transectFname[:-3] + 'nc'
            print 'Making %s ' % Tofname
            # first make transect
            TransectDict = sb.load_FRF_Transect(transectFname)
            TransectDict['time'] = nc.date2num(TransectDict['time'], 'seconds since 1970-01-01')
            makenc.makenc_FRFTransect(bathyDict=TransectDict, ofname=Tofname, globalYaml=transectGlobalYaml, varYaml=transectVarYaml)
        except Exception, e:
            print e
            errors.append(e)
            errorFname.append(transectFname)

    print 'Converting %d Grid (txt) file to netCDF ' % len(gridList)
    for gridFname in gridList:
        try:
            ofname = gridFname.split('.')[0] + '.nc'
            print 'Making %s ' %ofname
            makenc.convert_FRFgrid(gridFname, ofname, gridGlobalYaml, gridVarYaml, plotFlag=False)
            # ncfile = nc.Dataset(ofname)
            #elevation_points = np.ma.array(ncfile['elevation'][0,:,:], mask=np.isnan(ncfile['elevation'][0,:,:]))
            # remove -999's
            # maskedElev = (ncfile['elevation'][0,:,:]==-999.0)
            # elevation_points = np.ma.array(ncfile['elevation'][0,:,:], mask=maskedElev)
            #
            # plt.pcolor(ncfile['FRF_Xshore'][:], ncfile['FRF_Yshore'][:], elevation_points)
            # plt.colorbar()
            # plt.title('FRF GRID %s'% ofname[:-3].split('/')[-1])
            # plt.savefig(ofname[:-4]+'.png')
            # plt.close()
        except Exception, e:
            print e
            errors.append(e)
            errorFname.append(gridFname)

    # log errors
    # log errors that were encountered during file creation
    f = open(logFile, 'w')  # opening file
    f.write('File, Error\n')  # writing headers
    for aa in range(0, len(errorFname)):  # looping through errors
        f.write('%s,\n %s\n----------------------------\n\n' %(errorFname, errors))
    f.close()

if __name__ == "__main__":
    opts, args = getopt.getopt(sys.argv[1:], "h", ["help"])
    #  Location of where to look for files to convert
    globPath = args[0]
    if globPath[-1] != '/':
        globPath = globPath + '/'
    convertText2NetCDF(globPath)