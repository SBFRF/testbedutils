import datetime as DT
import os
import netCDF4 as nc
import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import griddata
import makenc
from getdatatestbed.getDataFRF import getObs
from sblib import geoprocess as gp, sblib as sb
from sblib.anglesLib import geo2STWangle


def frf2ij(xfrf, yfrf, x0, y0, dx, dy, ni, nj):
    """
    Convert FRF coordinates to ij grid locations.

    Matthew P. Geheran
    01 December 2017

    Parameters
    ----------
    xfrf : float
        FRF x-coordinate to convert
    yfrf : float
        FRF y-coordinate to convert
    x0 : float
        Grid origin x-coordinate (FRF)
    y0 : float

        Grid origin y-coordinate (FRF)
    dx : float
        Grid resolution in x-direction.
    dy : float
        Grid resolution in y-direction.
    ni : int
        Number of grid cells in the i direction.
    nj : int
        Number of grid cells in the j direction.

    Returns
    -------
    i : int
        Grid location in the i direction
    j: int
        Grid location in the j direction
    """
    dx_is_single_value = isinstance(dx, (float, int, long))
    dy_is_single_value = isinstance(dy, (float, int, long))

    # This routine does not handle variable grid spacing.
    if not dx_is_single_value or not dy_is_single_value:
        raise NotImplementedError('This routine does not handle variable grid spacing')

    xFRFgrid = x0 - np.arange(ni - 1)*dx - 0.5*dx
    yFRFgrid = y0 - np.arange(nj - 1)*dy - 0.5*dy

    i = np.abs(xfrf - xFRFgrid).argmin()
    j = np.abs(yfrf - yFRFgrid).argmin()

    # Convert from python base 0 indexing to STWAVE base 1.
    i += 1
    j += 1

    # Assign -99999 to i_sensor and j_sensor if the locations are
    # outside the grid.
    x_is_outside = xfrf < xFRFgrid.min() or xfrf > xFRFgrid.max()
    y_is_outside = yfrf < yFRFgrid.min() or yfrf > yFRFgrid.max()
    if x_is_outside or y_is_outside:
        i = -99999
        j = -99999

    return i, j

def convertGridNodesFromStatePlane(icoords, jcoords):
    """
    this function converts nodes of a grid coordinate in state plane to FRF coordinates using FRFcoord function

    :param icoords: an array of the i coordinates of a grid (easting, northing)
    :param jcoords: an array of the j coordinates of a grid (easting, northing)

    :return: array of frf coordinates for I and J of the grid
    """

    out = gp.FRFcoord(icoords[0], icoords[1])
    outIfrf =  np.array((out['xFRF'], out['yFRF'])).T

    out = gp.FRFcoord(jcoords[0], jcoords[1])
    outJfrf = np.array((out['xFRF'], out['yFRF'])).T

    return outIfrf, outJfrf

def makeTimeMeanBackgroundBathy(dir_loc, dSTR_s=None, dSTR_e=None, scalecDict=None, splineDict=None, plot=None):
    """

    This function will create a time-averaged background surface.
    It takes in a background netcdf file and adds in every survey between the start and end dates.
    Each survey is converted to a grid using scaleCinterpolation.
    These grids are all stacked on top of each other and averaged.
    (note - the original background grid is only counted once;
        in areas where it is the only data point the other values are nan)
    This final grid is smoothed using the scale-C interpolation at the end then written to a netcdf file.

    :param dSTR_s: string that determines the start date of the times of the surveys you want to use to update the DEM
                    format is  dSTR_s = '2013-01-04T00:00:00Z'
                    no matter what you put here, it will always round it down to the beginning of the month
    :param dSTR_e: string that determines the end date of the times of the surveys you want to use to update the DEM
                    format is dSTR_e = '2014-12-22T23:59:59Z'
                    no matter what you put here, it will always round it up to the end of the month
    :param dir_loc: place where you want to save the .nc files that get written
                    the function will make the year directories inside of this location on its own.
    :param scalecDict: keys are:
                        x_smooth - x direction smoothing length for scalecInterp
                        y_smooth - y direction smoothing length for scalecInterp

                        if not specified it will default to:
                        x_smooth = 100
                        y_smooth = 200
    :param splineDict: keys are:
                        splinebctype
                            options are....
                            2 - second derivative goes to zero at boundary
                            1 - first derivative goes to zero at boundary
                            0 - value is zero at boundary
                            10 - force value and derivative(first?!?) to zero at boundary
                        lc - spline smoothing constraint value (integer <= 1)
                        dxm -  coarsening of the grid for spline (e.g., 2 means calculate with a dx that is 2x input dx)
                                can be tuple if you want to do dx and dy separately (dxm, dym), otherwise dxm is used for both
                        dxi - fining of the grid for spline (e.g., 0.1 means return spline on a grid that is 10x input dx)
                                as with dxm, can be a tuple if you want separate values for dxi and dyi
                        targetvar - this is the target variance used in the spline function.
                        wbysmooth - y-edge smoothing length scale
                        wbxsmooth - x-edge smoothing length scale

                        if not specified it will default to:
                        splinebctype = 10
                        lc = 4
                        dxm = 2
                        dxi = 1
                        targetvar = 0.45
                        wbysmooth = 300
                        wbxsmooth = 100
    :param plot: do I want to plot this or not? 1 for yes, 0 for no
    :return:  netCDF file of the time mean bathymetry
    """
    # TODO add directions as to where to import these or how to get them, where they should be located ....
    from bsplineFunctions import bspline_pertgrid
    from scaleCinterp_python.DEM_generator import DEM_generator

    #HARD CODED VARIABLES!!!
    filelist = ['http://134.164.129.55/thredds/dodsC/FRF/geomorphology/elevationTransects/survey/surveyTransects.ncml']
    # this is just the location of the ncml for the transects!!!!!

    nc_b_loc = 'C:\Users\dyoung8\Desktop\David Stuff\Projects\CSHORE\Bathy Interpolation\TestNCfiles'
    nc_b_name = 'backgroundDEMt0.nc'
    # these together are the location of the standard background bathymetry that we started from.

    # Yaml files for my .nc files!!!!!
    global_yaml = 'C:\Users\dyoung8\PycharmProjects\makebathyinterp\yamls\BATHY\FRFt0_global.yml'
    var_yaml = 'C:\Users\dyoung8\PycharmProjects\makebathyinterp\yamls\BATHY\FRFt0_TimeMean_var.yml'

    # CS-array url - I just use this to get the position, not for any other data
    cs_array_url = 'http://134.164.129.55/thredds/dodsC/FRF/oceanography/waves/8m-array/2017/FRF-ocean_waves_8m-array_201707.nc'
    # where do I want to save any QA/QC figures
    fig_loc = 'C:\Users\dyoung8\Desktop\David Stuff\Projects\CSHORE\Bathy Interpolation\Test Figures\TimeWA_figs'

    #check scalecDict
    if scalecDict is None:
        x_smooth = 100  # scale c interp x-direction smoothing
        y_smooth = 200  # scale c interp y-direction smoothing
    else:
        x_smooth = scalecDict['x_smooth']  # scale c interp x-direction smoothing
        y_smooth = scalecDict['y_smooth']  # scale c interp y-direction smoothing

    #check dSTR_s
    if dSTR_s is None:
        dSTR_s = '1970-01-01T00:00:00Z' # set it to before the first survey
    else:
        pass

    #check dSTR_e
    if dSTR_e is None:
        dSTR_e = DT.datetime.strftime(DT.datetime.now(), '%Y-%m-%dT%H:%M:%SZ')  # set it to right now
    else:
        pass

    # force the survey to start at the first of the month and end at the last of the month!!!!
    dSTR_s = dSTR_s[0:7] + '-01T00:00:00Z'
    if dSTR_e[5:7] == '12':
        dSTR_e = str(int(dSTR_e[0:4]) + 1) + '-01' + '-01T00:00:00Z'
    else:
        dSTR_e = dSTR_e[0:5] + str(int(dSTR_e[5:7]) + 1).zfill(2) + '-01T00:00:00Z'

    d_s = DT.datetime.strptime(dSTR_s, '%Y-%m-%dT%H:%M:%SZ')
    d_e = DT.datetime.strptime(dSTR_e, '%Y-%m-%dT%H:%M:%SZ')

    # ok, I just need to go through and find all surveys that fall in this date range
    bathy = nc.Dataset(filelist[0])
    # pull down all the times....
    times = nc.num2date(bathy.variables['time'][:], bathy.variables['time'].units, bathy.variables['time'].calendar)
    all_surveys = bathy.variables['surveyNumber'][:]

    # find some stuff here...
    mask = (times >= d_s) & (times < d_e)  # boolean true/false of time
    idx = np.where(mask)[0]

    # what surveys are in this range?
    surveys = np.unique(bathy.variables['surveyNumber'][idx])

    # get rid of any surveys with rounded middle time not in my range
    for tt in range(0, len(surveys)):
        ids = (all_surveys == surveys[tt])
        surv_times = times[ids]
        # pull out the mean time
        surv_timeM = surv_times[0] + (surv_times[-1] - surv_times[0]) / 2
        # round it to nearest 12 hours.
        surv_timeM = sb.roundtime(surv_timeM, roundTo=1 * 12 * 3600)

        # if the rounded time IS in the month, great
        if (surv_timeM >= d_s) and (surv_timeM < d_e):
            pass
        else:
            # if not set it to a fill value
            surveys[tt] == -1000
    # drop all the surveys that we decided are not going to use
    surveys = surveys[surveys >= 0]


    # pull the original background DEM
    old_bathy = nc.Dataset(os.path.join(nc_b_loc, nc_b_name))
    Zi = old_bathy.variables['elevation'][:]
    xFRFi_vec = old_bathy.variables['xFRF'][:]
    yFRFi_vec = old_bathy.variables['yFRF'][:]

    # read out the dx and dy of the background grid!!!
    # assume this is constant grid spacing!!!!!
    dx = abs(xFRFi_vec[1] - xFRFi_vec[0])
    dy = abs(yFRFi_vec[1] - yFRFi_vec[0])

    xFRFi, yFRFi = np.meshgrid(xFRFi_vec, yFRFi_vec)
    rows, cols = np.shape(xFRFi)

    # pre-allocate my netCDF dictionary variables here....
    elevation = np.nan*np.zeros((len(surveys)+1, rows, cols))
    weights = np.nan * np.zeros((len(surveys) + 1, rows, cols))
    xFRF = np.zeros(cols)
    yFRF = np.zeros(rows)

    # ok, now that I have the list of the surveys I am going to keep.....
    for tt in range(0, len(surveys)):

        # get the times of each survey
        ids = (all_surveys == surveys[tt])

        # pull out this NC stuf!!!!!!!!
        dataX, dataY, dataZ = [], [], []
        dataX = bathy['xFRF'][ids]
        dataY = bathy['yFRF'][ids]
        dataZ = bathy['elevation'][ids]
        profNum = bathy['profileNumber'][ids]
        survNum = bathy['surveyNumber'][ids]
        stimes = nc.num2date(bathy.variables['time'][ids], bathy.variables['time'].units,
                             bathy.variables['time'].calendar)
        # pull out the mean time
        stimeM = min(stimes) + (max(stimes) - min(stimes)) / 2
        # round it to nearest 12 hours.
        stimeM = sb.roundtime(stimeM, roundTo=1 * 12 * 3600)

        assert len(np.unique(survNum)) == 1, 'makeTimeMeanBackgroundBathy error: You have pulled down more than one survey number!'
        assert isinstance(dataZ, np.ndarray), 'makeTimeMeanBackgroundBathy error: Script only handles np.ndarrays for the transect data at this time!'

        # build my new bathymetry from the FRF transect files

        # what are my subgrid bounds?
        surveyDict = {}
        surveyDict['dataX'] = dataX
        surveyDict['dataY'] = dataY
        surveyDict['profNum'] = profNum

        gridDict = {}
        gridDict['dx'] = dx
        gridDict['dy'] = dy
        gridDict['xFRFi_vec'] = xFRFi_vec
        gridDict['yFRFi_vec'] = yFRFi_vec

        temp = subgridBounds2(surveyDict, gridDict, maxSpace=249)
        x0 = temp['x0']
        x1 = temp['x1']
        y0 = temp['y0']
        y1 = temp['y1']
        del temp

        # if you wound up throwing out this survey!!!
        if x0 is None:
            newZi = np.nan * np.zeros(np.shape(Zi))

        else:
            print np.unique(survNum)
            dict = {'x0': x0,  # gp.FRFcoord(x0, y0)['Lon'],  # -75.47218285,
                    'y0': y0,  # gp.FRFcoord(x0, y0)['Lat'],  #  36.17560399,
                    'x1': x1,  # gp.FRFcoord(x1, y1)['Lon'],  # -75.75004989,
                    'y1': y1,  # gp.FRFcoord(x1, y1)['Lat'],  #  36.19666112,
                    'lambdaX': dx,
                    # grid spacing in x  -  Here is where CMS would hand array of variable grid spacing
                    'lambdaY': dy,  # grid spacing in y
                    'msmoothx': x_smooth,  # smoothing length scale in x
                    'msmoothy': y_smooth,  # smoothing length scale in y
                    'msmootht': 1,  # smoothing length scale in Time
                    'filterName': 'hanning',
                    'nmseitol': 0.75,
                    'grid_coord_check': 'FRF',
                    'grid_filename': '',  # should be none if creating background Grid!  becomes best guess grid
                    'data_coord_check': 'FRF',
                    'xFRF_s': dataX,
                    'yFRF_s': dataY,
                    'Z_s': dataZ,
                    }

            out = DEM_generator(dict)

            # read some stuff from this dict like a boss
            Zn = out['Zi']
            xFRFn_vec = out['x_out']
            yFRFn_vec = out['y_out']
            MSEn = out['MSEi']
            targetvar = 0.45
            wb = 1 - np.divide(MSEn, targetvar + MSEn) # these are my weights from scale C

            try:
                x1 = np.where(xFRFi_vec == min(xFRFn_vec))[0][0]
                x2 = np.where(xFRFi_vec == max(xFRFn_vec))[0][0]
                y1 = np.where(yFRFi_vec == min(yFRFn_vec))[0][0]
                y2 = np.where(yFRFi_vec == max(yFRFn_vec))[0][0]

                newZi = np.nan * np.zeros(np.shape(Zi))
                newZi[y1:y2 + 1, x1:x2 + 1] = Zn

                new_wb = np.nan * np.zeros(np.shape(Zi))
                new_wb[y1:y2 + 1, x1:x2 + 1] = wb

            except IndexError:
                newZi = np.nan * np.zeros(np.shape(Zi))
                new_wb = np.nan * np.zeros(np.shape(wb))

            elevation[tt, :, :] = newZi
            weights[tt, :, :] = new_wb

            """
            # plot each newZi to see if it looks ok
            fig_name = 'backgroundDEM_' + str(surveys[tt]) + '.png'
            plt.pcolor(xFRFi_vec, yFRFi_vec, elevation[tt, :, :], cmap=plt.cm.jet, vmin=-13, vmax=5)
            cbar = plt.colorbar()
            cbar.set_label('(m)')
            plt.scatter(dataX, dataY, marker='o', c='k', s=1, alpha=0.25, label='Transects')
            plt.xlabel('xFRF (m)')
            plt.ylabel('yFRF (m)')
            plt.legend()
            plt.savefig(os.path.join(fig_loc, fig_name))
            plt.close()
            """

            """
            # plot where it is nan....
            nan_loc = np.isnan(elevation[tt, :, :])
            fig_name = 'backgroundDEM_' + str(surveys[tt]) + 'NaN_loc' + '.png'
            plt.pcolor(xFRFi_vec, yFRFi_vec, nan_loc, cmap=plt.cm.jet, vmin=0, vmax=1)
            cbar = plt.colorbar()
            plt.xlabel('xFRF (m)')
            plt.ylabel('yFRF (m)')
            plt.legend()
            plt.savefig(os.path.join(fig_loc, fig_name))
            plt.close()
            """

    # drop in my original bathymetry as the last index!
    elevation[-1, :, :] = Zi
    weights[-1, :, :] = np.ones((rows, cols))
    xFRF = xFRFi[0, :]
    yFRF = yFRFi[:, 1]

    cleaned_elevation = np.ma.masked_array(elevation, np.isnan(elevation))
    cleaned_weights = np.ma.masked_array(weights, np.isnan(weights))

    # do a nanmean on the elevation!!!!
    Z = np.ma.average(cleaned_elevation, axis=0, weights=cleaned_weights)

    # plot the mean to see if that is the problem?
    fig_name = 'backgroundDEM_' + 'TimeMean_NoScaleC' + '.png'
    plt.pcolor(xFRFi_vec, yFRFi_vec, Z[:, :], cmap=plt.cm.jet, vmin=-13, vmax=5)
    cbar = plt.colorbar()
    cbar.set_label('(m)')
    plt.xlabel('xFRF (m)')
    plt.ylabel('yFRF (m)')
    plt.legend()
    plt.savefig(os.path.join(fig_loc, fig_name))
    plt.close()


    # run this through the DEM_generator function to smooth it....
    xFRF_mesh, yFRF_mesh = np.meshgrid(xFRF, yFRF)
    # reshape them and my Z...

    dataX = np.reshape(xFRF_mesh, (np.shape(xFRF_mesh)[0] * np.shape(xFRF_mesh)[1], 1)).flatten()
    dataY = np.reshape(yFRF_mesh, (np.shape(yFRF_mesh)[0] * np.shape(yFRF_mesh)[1], 1)).flatten()
    dataZ = np.reshape(Z, (np.shape(Z)[0] * np.shape(Z)[1], 1)).flatten()

    dict = {'x0': max(xFRF),
            'y0': max(yFRF),
            'x1': min(xFRF),
            'y1': min(yFRF),
            'lambdaX': dx,
            # grid spacing in x  -  Here is where CMS would hand array of variable grid spacing
            'lambdaY': dy,  # grid spacing in y
            'msmoothx': int(x_smooth),  # smoothing length scale in x
            'msmoothy': int(2*y_smooth),  # smoothing length scale in y
            'msmootht': 1,  # smoothing length scale in Time
            'filterName': 'hanning',
            'nmseitol': 0.75,
            'grid_coord_check': 'FRF',
            'grid_filename': '',  # should be none if creating background Grid!  becomes best guess grid
            'data_coord_check': 'FRF',
            'xFRF_s': dataX,
            'yFRF_s': dataY,
            'Z_s': dataZ,
            }

    out2 = DEM_generator(dict)

    # read some stuff from this dict like a boss
    del Z
    del xFRF
    del yFRF
    Z = out2['Zi']
    MSEn = out2['MSEi']
    xFRF = out2['x_out']
    yFRF = out2['y_out']

    # do we want to spline the ends?
    if splineDict is None:
        pass
    else:
        # we do spline the ends....
        splinebctype = splineDict['splinebctype']
        lc = splineDict['lc']
        dxm = splineDict['dxm']
        dxi = splineDict['dxi']
        targetvar = splineDict['targetvar']

        # get the difference!!!!
        Zdiff = Z - Zi

        # spline time?
        wb = 1 - np.divide(MSEn, targetvar + MSEn)
        newZdiff = bspline_pertgrid(Zdiff, wb, splinebctype=splinebctype, lc=lc, dxm=dxm, dxi=dxi)
        newZ = Zi + newZdiff

        del Z
        Z = newZ

    # save this to an nc file?
    # write the nc_file for this month, like a boss, with greatness
    nc_dict = {}
    nc_dict['elevation'] = Z
    nc_dict['xFRF'] = xFRF
    nc_dict['yFRF'] = yFRF

    nc_name = 'backgroundDEMt0_TimeMean' + '.nc'
    makenc.makenc_t0BATHY(os.path.join(dir_loc, nc_name), nc_dict, globalYaml=global_yaml, varYaml=var_yaml)

    if plot is None:
        pass
    else:
        # plot the bathymetry before and after....
        # where is the cross shore array?
        test = nc.Dataset(cs_array_url)
        Lat = test['latitude'][:]
        Lon = test['longitude'][:]
        # convert to FRF
        temp = gp.FRFcoord(Lon, Lat)
        CSarray_X = temp['xFRF']
        CSarray_Y = temp['yFRF']

        # original data
        fig_name = 'backgroundDEM_orig' + '.png'
        plt.figure()
        plt.pcolor(xFRF, yFRF, Zi, cmap=plt.cm.jet, vmin=-13, vmax=5)
        cbar = plt.colorbar()
        cbar.set_label('(m)')
        plt.plot(CSarray_X, CSarray_Y, 'rX', label='8m-array')
        plt.xlabel('xFRF (m)')
        plt.ylabel('yFRF (m)')
        plt.legend()
        plt.savefig(os.path.join(fig_loc, fig_name))
        plt.close()

        # new time-mean data
        fig_name = 'backgroundDEM_TimeMean' + '.png'
        plt.figure()
        plt.pcolor(xFRF, yFRF, Z, cmap=plt.cm.jet, vmin=-13, vmax=5)
        cbar = plt.colorbar()
        cbar.set_label('(m)')
        plt.plot(CSarray_X, CSarray_Y, 'rX', label='8m-array')
        plt.xlabel('xFRF (m)')
        plt.ylabel('yFRF (m)')
        plt.legend()
        plt.savefig(os.path.join(fig_loc, fig_name))
        plt.close()

def createGridNodesinFRF(x0, y0, dx, dy, ni, nj):
    """
    This function assumes azimuth of the grid is the same as that of the FRF coordinate system
    code developed for CMS wave and
    :param x0: origin of x in FRF coords
    :param y0: origin of grid in FRF coords
    :param dx: Array of dx values
    :param dy: Array of dy values
    :param ni: number of cells in i
    :param nj: number of cells in j
    :return:  array of i coords, array of j coordinates
    """
    assert dx.shape[0] == ni-1, 'This function assumes that there are n-1 dx values'

    if np.mean(np.diff(dx)) != np.mean(dx):  # vairable spacing cell array
        icoord = np.zeros(ni)  # assume
        jcoord = np.zeros(nj)
        icoord[0] = x0
        jcoord[0] = y0
        for xx, dxx in enumerate(dx):
            icoord[xx+1] = icoord[xx] - dxx  # assumes offshore origin
        for yy, dyy in enumerate(dy):
            jcoord[yy+1] = jcoord[yy] - dyy
    else:
        raise NotImplementedError

    return icoord, jcoord

def makeBackgroundBathyAzimuth(origin, geo_ang, dx, dy, ni, nj, coord_system='FRF'):

    """
    This function makes the grid nodes using the origin and the azimuth

    :param origin: this is the origin of your new grid in the form (xFRF, yFRF), (Lat, Lon), (easting, northing)
    :param geo_ang: angle of the x-axis of your grid clockwise relative to true north
    :param dx: x-direction spacing between your grid nodes in m
    :param dy: y-direction spacing between your grid nodes in m
    :param ni: number of nodes in the x-direction
    :param nj: number of nodes in the y-direction
    :param coord_system: 'FRF', 'utm', 'stateplane', 'LAT/LON'

    :return: dictionary with keys containing
            2D arrays of x & y grid nodes in the coordinate system you specify (easting/northing, lat/lon)
            2D array of bottom elevation at those node locations from the background dem
    """

    assert len(origin) == 2, 'makeBackgroundBathy Error: invalid origin input.  origin input must be of form (xFRF, yFRF), (easting, northing), or (LAT, LON)'

    # first check the coord_system string to see if it matches!
    coord_list = ['FRF', 'stateplane', 'utm', 'Lat/Lon']
    import pandas as pd
    import string
    exclude = set(string.punctuation)
    columns = ['coord', 'user']
    df = pd.DataFrame(index=range(0, np.size(coord_list)), columns=columns)
    df['coord'] = coord_list
    df['user'] = coord_system
    df['coordToken'] = df.coord.apply(lambda x: ''.join(ch for ch in str(x) if ch not in exclude).strip().upper())
    df['coordToken'] = df.coordToken.apply(lambda x: ''.join(str(x).split()))
    df['userToken'] = df.user.apply(lambda x: ''.join(ch for ch in str(x) if ch not in exclude).strip().upper())
    df['userToken'] = df.userToken.apply(lambda x: ''.join(str(x).split()))
    userToken = np.unique(np.asarray(df['userToken']))[0]
    assert df['coordToken'].str.contains(userToken).any(), 'makeBackgroundBathy Error: invalid coord_system string.  Acceptable strings include %s' % coord_list

    # convert origin to stateplane if it isn't already....
    if userToken == 'FRF':
        temp = gp.FRF2ncsp(origin[0], origin[1])
        x0 = temp['StateplaneE']
        y0 = temp['StateplaneN']

    elif userToken == 'STATEPLANE':
        x0 = origin[0]
        y0 = origin[1]

    elif userToken == 'UTM':
        temp = gp.utm2ncsp(origin[0], origin[1], 18, 'S')
        x0 = temp['easting']
        y0 = temp['northing']

    elif userToken == 'LATLON':
        temp = gp.LatLon2ncsp(origin[1], origin[0])
        x0 = temp['StateplaneE']
        y0 = temp['StateplaneN']

    else:
        pass

    # convert my geographic coordinate angle to azimuth!!
    azi = geo2STWangle(geo_ang, zeroAngle=71.8)
    # azi = geo_ang

    # note: I just striaght up pulled this bit of code from CreateGridNodesInStatePlane

    # calculating change in alongshore coordinate for northing and easting
    # given the associated dx dy
    dE_j = dy * np.cos(np.deg2rad(azi + 90))
    dN_j = dy * np.sin(np.deg2rad(azi + 90))

    # calculating change in cross-shore coordinate for northing and easting
    dE_i = dx * np.cos(np.deg2rad(azi))
    dN_i = dx * np.sin(np.deg2rad(azi))

    easting = np.zeros((ni, nj))
    northing = np.zeros((ni, nj))

    for ii in range(0, ni):
        for jj in range(0, nj):
            easting[ii, jj] = x0 + ii * dE_i + jj * dE_j
            northing[ii, jj] = y0 + ii * dN_i + jj * dN_j


    #convert all my new points to utm!
    east_vec = easting.reshape((1, easting.shape[0] * easting.shape[1]))[0]
    north_vec = northing.reshape((1, northing.shape[0] * northing.shape[1]))[0]

    # convert them to UTM
    temp = gp.ncsp2utm(east_vec, north_vec)
    utmE = temp['utmE']
    utmN = temp['utmN']


    # pull out the piece of the DEM I need!
    # these are just some random times I made up because the getObs class requires it.  They have no effect on the
    # bathymetry that is pulled, so put whatever you want in here...
    d_s = DT.datetime.strptime('2015-06-20T12:00:00Z', '%Y-%m-%dT%H:%M:%SZ')
    d_e = DT.datetime.strptime('2015-06-20T12:00:00Z', '%Y-%m-%dT%H:%M:%SZ')
    frf_bathy = getObs(d_s, d_e)
    buffer = 20 # buffer around my grid in m to make sure I pull at least one point to the outside

    bathyDEM = frf_bathy.getBathyRegionalDEM(min(utmE) - buffer, max(utmE) + buffer, min(utmN) - buffer, max(utmN) + buffer)
    assert np.size(np.where(bathyDEM == -9999)) <= 1, 'makeBackgroundDEM Error:  Your domain contains areas with no background DEM data!'


    # interpolate the bottom elevation onto my new nodes!!!!
    utmEdem = bathyDEM['utmEasting'].reshape((1, bathyDEM['utmEasting'].shape[0] * bathyDEM['utmEasting'].shape[1]))[0]
    utmNdem = bathyDEM['utmNorthing'].reshape((1, bathyDEM['utmNorthing'].shape[0] * bathyDEM['utmNorthing'].shape[1]))[0]

    points = (utmEdem, utmNdem)
    values = bathyDEM['bottomElevation'].reshape((1, bathyDEM['bottomElevation'].shape[0] * bathyDEM['bottomElevation'].shape[1]))[0]


    # do the interpolation
    bottomElevation_vec = griddata(points, values, (utmE, utmN), method='linear')
    # reshape it back to 2D array!
    bottomElevation = bottomElevation_vec.reshape((easting.shape[0], easting.shape[1]))


    # now convert my stateplane grid back into the coordinates specified!!!!
    if userToken == 'FRF':
        temp = gp.ncsp2FRF(east_vec, north_vec)
        x_vec = temp['xFRF']
        y_vec = temp['yFRF']

    elif userToken == 'STATEPLANE':
        x_vec = east_vec
        y_vec = north_vec

    elif userToken == 'UTM':
        x_vec = utmE
        y_vec = utmN

    elif userToken == 'LATLON':
        temp = gp.ncsp2LatLon(east_vec, north_vec)
        x_vec = temp['lon']
        y_vec = temp['lat']

    else:
        pass


    # reshape them back
    x = x_vec.reshape((easting.shape[0], easting.shape[1]))
    y = y_vec.reshape((easting.shape[0], easting.shape[1]))

    # return the grid in the coordinate system of the origin
    out = {}
    out['bottomElevation'] = bottomElevation

    if userToken == 'FRF':
        out['xFRF'] = x
        out['yFRF'] = y

    elif userToken == 'STATEPLANE':
        out['easting'] = x
        out['northing'] = y

    elif userToken == 'UTM':
        out['utmEasting'] = x
        out['utmNorthing'] = y

    elif userToken == 'LATLON':
        out['longitude'] = x
        out['latitude'] = y

    else:
        pass

    return out

def makeBackgroundBathyCorners(LLHC, URHC, dx, dy, coord_system='FRF'):

    """
    This function makes grid nodes using the corners of the grid using different coordinate systems

    :param LLHC: tuple: lower left hand corner of the desired domain (xFRF, yFRF) (easting, northing) or (Lat, Lon)
    :param URHC: tuple: upper right hand corner of the desired domain (xFRF, yFRF) (easting, northing) or (Lat, Lon)
    :param dx: x-direction grid spacing in m - lat/lon corners get converted to utm!!!
    :param dy: y-direction grid spacing in m - lat/lon corners get converted to utm!!!
    :param coord_system: string containing the coordinate system for your corners ('FRF' 'utm', 'stateplane', or 'LAT/LON')
    :return: dictionary containing 2D arrays of:
            xFRF (or easting or longitude)
            yFRF (or northing or Latitude)
            bottomElevation at those points interpolated from background DEM onto desired grid
    """

    # first check the coord_system string to see if it matches!
    coord_list = ['FRF', 'LAT/LON', 'utm', 'stateplane']
    import pandas as pd
    import string
    exclude = set(string.punctuation)
    columns = ['coord', 'user']
    df = pd.DataFrame(index=range(0, np.size(coord_list)), columns=columns)
    df['coord'] = coord_list
    df['user'] = coord_system
    df['coordToken'] = df.coord.apply(lambda x: ''.join(ch for ch in str(x) if ch not in exclude).strip().upper())
    df['coordToken'] = df.coordToken.apply(lambda x: ''.join(str(x).split()))
    df['userToken'] = df.user.apply(lambda x: ''.join(ch for ch in str(x) if ch not in exclude).strip().upper())
    df['userToken'] = df.userToken.apply(lambda x: ''.join(str(x).split()))
    userToken = np.unique(np.asarray(df['userToken']))[0]
    assert df['coordToken'].str.contains(userToken).any(), 'makeBackgroundBathy Error: invalid coord_system string.  Acceptable strings include %s' % coord_list

    # second, check the format of the corner inputs
    LLHC = np.asarray(LLHC)
    URHC = np.asarray(URHC)
    assert len(LLHC) == len(URHC) == 2, 'makeBackgroundBathy Error: invalid corner input.  corner inputs must be of form (xFRF, yFRF) (easting, northing) or (LAT, LON)'


    # make my new grid first!!!
    x_pts = [LLHC[0], URHC[0]]
    y_pts = [LLHC[1], URHC[1]]

    if userToken == 'LATLON':
        # if corners are in LAT/LON then we convert directly to UTM and work from that
        temp = gp.LatLon2utm(x_pts, y_pts)
        x_vec = np.arange(min(temp['utmE']), max(temp['utmE']), dx)
        y_vec = np.arange(min(temp['utmN']), max(temp['utmN']), dy)

    else:
        x_vec = np.arange(x_pts[0], x_pts[1], dx)
        y_vec = np.arange(y_pts[0], y_pts[1], dy)


    xv, yv = np.meshgrid(x_vec, y_vec)
    # reshape my points
    xv_vec = xv.reshape((1, xv.shape[0] * xv.shape[1]))[0]
    yv_vec = yv.reshape((1, yv.shape[0] * yv.shape[1]))[0]


    # convert all my points to UTM
    utmE = np.zeros(len(xv_vec))
    utmN = np.zeros(len(yv_vec))
    if userToken == 'FRF':
        for ii in range(0, len(xv_vec)):
            # note: I didn't use FRFcoord fxn because I don't want the code to "guess" what coordinate system I am in.
            temp = gp.FRF2ncsp(xv_vec[ii], yv_vec[ii])
            spE = temp['StateplaneE']
            spN = temp['StateplaneN']
            temp2 = gp.ncsp2utm(spE, spN)
            utmE[ii] = temp2['utmE']
            utmN[ii] = temp2['utmN']

    elif userToken == 'LATLON':
        utmE = xv_vec
        utmN = yv_vec

    elif userToken == 'UTM':
        utmE = xv_vec
        utmN = yv_vec

    elif userToken == 'STATEPLANE':
        temp = gp.ncsp2utm(xv_vec, yv_vec)
        utmE = temp['utmE']
        utmN = temp['utmN']

    else:
        pass


    # these are just some random times I made up because the getObs class requires it.  They have no effect on the
    # bathymetry that is pulled, so put whatever you want in here...
    d_s = DT.datetime.strptime('2015-06-20T12:00:00Z', '%Y-%m-%dT%H:%M:%SZ')
    d_e = DT.datetime.strptime('2015-06-20T12:00:00Z', '%Y-%m-%dT%H:%M:%SZ')
    frf_bathy = getObs(d_s, d_e)
    buffer = 20  # buffer in m(grid spacing is 10m, so this will make sure you always
    # have at least one node to each side

    bathyDEM = frf_bathy.getBathyRegionalDEM(min(utmE) - buffer, max(utmE) + buffer, min(utmN) - buffer, max(utmN) + buffer)

    # the getBathyDEM function will check to see if you are too close to the bounds.
    # All you have to do now is check to see if any piece of this sub-DEM has fill values instead of data!

    assert np.size(np.where(bathyDEM == -9999)) <= 1, 'makeBackgroundDEM Error:  Your domain contains areas with no background DEM data!'

    """
    # check to see if this actually worked...
    import matplotlib.pyplot as plt
    fig_name = 'DEMsubgrid.png'
    fig_loc = 'C:\Users\RDCHLDLY\Desktop\David Stuff\Projects\CSHORE\Bathy Interpolation\Test Figures'
    plt.contourf(bathyDEM['utmEasting'], bathyDEM['utmNorthing'], bathyDEM['bottomElevation'])
    plt.axis('equal')
    plt.savefig(os.path.join(fig_loc, fig_name))
    plt.close()
    """

    # reshape my DEM into a list of points
    utmEdem = bathyDEM['utmEasting'].reshape((1, bathyDEM['utmEasting'].shape[0] * bathyDEM['utmEasting'].shape[1]))[0]
    utmNdem = bathyDEM['utmNorthing'].reshape((1, bathyDEM['utmNorthing'].shape[0] * bathyDEM['utmNorthing'].shape[1]))[0]

    points = (utmEdem, utmNdem)
    values = bathyDEM['bottomElevation'].reshape((1, bathyDEM['bottomElevation'].shape[0] * bathyDEM['bottomElevation'].shape[1]))[0]
    # do the interpolation
    bottomElevation_vec = griddata(points, values, (utmE, utmN), method='linear')
    # reshape it back to 2D array!
    bottomElevation = bottomElevation_vec.reshape((xv.shape[0], xv.shape[1]))

    # so now I have xv, yv, bottomElevation on a rectangular grid in my new coordinate system. I think
    """
    # check to see if this actually worked...
    import matplotlib.pyplot as plt
    fig_name = 'newGrid.png'
    fig_loc = 'C:\Users\RDCHLDLY\Desktop\David Stuff\Projects\CSHORE\Bathy Interpolation\Test Figures'
    plt.contourf(xv, yv, bottomElevation)
    plt.axis('equal')
    plt.xlabel('xFRF')
    plt.ylabel('yFRF')
    plt.savefig(os.path.join(fig_loc, fig_name))
    plt.close()
    """

    # now return my stuff to the user....
    out = {}
    if userToken == 'FRF':
        out['xFRF'] = xv
        out['yFRF'] = yv
        out['bottomElevation'] = bottomElevation

    elif userToken == 'STATEPLANE':
        out['easting'] = xv
        out['northing'] = yv
        out['bottomElevation'] = bottomElevation

    elif userToken == 'UTM':
        out['utmEasting'] = xv
        out['utmNorthing'] = yv
        out['bottomElevation'] = bottomElevation

    elif userToken == 'LATLON':
        # if it is lat lon I have to convert all my points back from UTM!!!!
        temp = gp.utm2LatLon(xv_vec, yv_vec, 18, 'S')

        lat_vec = temp['lat']
        lon_vec = temp['lon']

        out['latitude'] = lat_vec.reshape((xv.shape[0], xv.shape[1]))
        out['longitude'] = lon_vec.reshape((xv.shape[0], xv.shape[1]))
        out['bottomElevation'] = bottomElevation

    else:
        pass

    return out

def CreateGridNodesInStatePlane(x0, y0, azi, dx, dy, ni, nj):
    """
    this function takes in a sim file and creates tuples of grid locations
    in state plane, can further be converted to lat/lon
    stateplane sp3200

    :param x0: integer/float describing origin in x (easting)
    :param y0: integer/float describing origin in y (northing)
    :param azi: grid azimuth defining rotation of grid
    :param dx: can be integer/float or numpy array/list describing cell width in x direction (i)
    :param dy: can be integer/float or numpy array/list describing cell with in y direction (j)
    :param ni: integer/float describing number of cells in i
    :param nj: integer/float describing  number of cells in j

    :return: tuples of i/j coords, jStatePlane in stateplane sp3200
    """
    # calculating change in alongshore coordinate for northing and easting
    # given the associated dx dy
    dE_j = dy * np.cos(np.deg2rad(azi + 90))
    dN_j = dy * np.sin(np.deg2rad(azi + 90))
    # calculating change in cross-shore coordinate for northing and easting
    dE_i = dx * np.cos(np.deg2rad(azi))
    dN_i = dx * np.sin(np.deg2rad(azi))
    # create Easting & Northing coordinates for
    # cross shore location (in grid space)
    try:  # this works for when dE_i is not an array .... ie regularly spaced grid nodes
        easting_i = np.linspace(x0, x0 + ni * dE_i, num=ni, endpoint=True)
        northing_i = np.linspace(y0, y0 + ni * dN_i, num=ni, endpoint=True)
        # create Northing and Easting coords for Along-shore location
        easting_j = np.linspace(x0, x0 + nj * dE_j, num=nj, endpoint=True)
        northing_j = np.linspace(y0, y0 + nj * dN_j, num=nj, endpoint=True)
    except ValueError:  # for instances when grid nodes are irregularly spaced
        easting_i, northing_i = [x0], [y0]  # seeding the origin for the first value in the coordinates
        easting_j, northing_j = [x0], [y0]
        for ii in range(0, ni-1):  # first doing the i coordinate (start at origin add the change in e/n for each grid cell
            easting_i.append(easting_i[ii] + dE_i[ii])
            northing_i.append(northing_i[ii] + dN_i[ii])
        for jj in range(0, nj-1):
            easting_j.append(easting_j[jj] + dE_j[jj])
            northing_j.append(northing_j[jj] + dN_j[jj])
        # convert lists to arrays
        easting_i = np.array(easting_i)
        easting_j = np.array(easting_j)
        northing_i = np.array(northing_i)
        northing_j = np.array(northing_j)
        assert easting_j.shape[0] == nj, 'len of cstateplane sp3200oordinates are not the same as the number of cells'
        assert northing_i.shape[0] == ni, 'len of coordinates are not the same as the number of cells'

    icoords = np.array([easting_i, northing_i])
    jcoords = np.array([easting_j, northing_j])

    return icoords, jcoords