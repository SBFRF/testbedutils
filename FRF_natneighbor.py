# -*- coding: utf-8 -*-
"""
This projects is created to test out the creation of a tool for 
the natural neighbor interpolation algorithm  3 tools are created, based in linux
the frf_grid_product creates the grid using the pyNGL library, then it is written
the voroni creates a tesselation of the surface for each file

assumes
"""

try:
    import pytess as pt
except ImportError:
    print 'no pytess module found, pip installing now'
    import os
    os.system('pip install pytess')
    import pytess as pt
import numpy as np
import csv
from collections import Counter
import glob
import datetime as DT
from matplotlib import pyplot as plt
try:
    from PyNGL import Ngl as ngl
except ImportError:
    raise ImportError("download the PyNGL library from https://www.pyngl.ucar.edu/Download/"\
     "then place in site-packages of python distro. only works with linux and macOS")

def frf_grid_product(fname_in, dxdy=5, header=0, **kwargs):
    """
    This function will create a gridded product using the FRF coordinates.  It 
    is designed to currently flow into the standard work flow of the FRF geomorphology
    workunits work flow currently in place.  Designed to replace the kriging algorithm
    currently used to create the gridded product.  
    the function relies on the PyNGL library created by the National Center for Atmospheric
    Resarch, Supported by NSF.  the library is available here: https://www.pyngl.ucar.edu/Download/
    once the package is downloaded, place the PyNGL folder into .../X/lib/site-packages/ folder
    where .../X is your python distribution folder.  eg /home/user/anaconda2/lib/python2.7/site-packages
    
    This script was created by Spicer Bak at USACE CHL FRF 
    contact: spicer.bak@usace.army.mil
    
    ::INPUT::
    :param fname_in: is a csv file formatted as standard frf transect data
    :param dxdy:  This is the output grid size in meters, default is 5 meters
            if dxdy = False, then unique x and y grid cells can be defined 
            using kwargs of dx and dy 
    :param header: this is the length of file header to be skipped upon csv import
    
    **kwarg
       dx = grid x spacing if not equal to dy
       dy = grid y spacing if not equal to dx
       plot = 1 turns a display plot on
       xmax = FRF max Xshore value to grid (w/o kwarg default 950)
       xmin = FRF min Xshore value to grid (w/o kwarg default 50)
       ymax = FRF max Yshore value to grid (w/o kwarg default 1100)
       xmin = FRF min Yshore value to grid (w/o kwarg default -100)
    
    :return:  data dictionary with x coords, y coords, z grid, and a tuple for x and y for voroni
        polygon creation
    """
    miniSurveyThresh = 22 # number of lines that are needed to grid a full survey
    # counting number of fields to ensure
    fileFields = fname_in.replace('/','_').split('_')
    field_count = Counter(fileFields)

    if field_count['FRF'] >= 2:
        Full = True  # its a full survey
    else:
        Full = False
    try:
        dx = kwargs['dx']
        dy = kwargs['dy']
    except (NameError, KeyError):
        dx = dy = dxdy    
    # import file with given path above
    raw_x = []
    raw_y = []
    lnNum, raw_z = [], []
    frf=[]
    with open(fname_in, 'rb') as f:
        reader = csv.reader(f)
        for row in reader:
            frf.append(row[0])
            lnNum.append(row[1])  # line number (named from y in frf coord)
            raw_x.append(row[7])  # x in string format
            raw_y.append(row[8])  # y in string format
            raw_z.append(row[9])  # z in string format
    # initializing values to make strictly numbers, imported as strings
    num_x = np.zeros(len(raw_x) - header)
    num_y = np.zeros(len(raw_y) - header)
    num_z = np.zeros(len(raw_z) - header)
    tup = []
    # making strings into floating point numbers
    for ii in range(header, len(raw_x)):
        num_x[ii-1] = float(raw_x[ii])
        num_y[ii-1] = float(raw_y[ii])
        num_z[ii-1] = float(raw_z[ii])
        tup.append((float(raw_x[ii]), float(raw_y[ii])))  #, float(raw_z[ii])))
    # count line numbers in file
    lineCount = len(Counter(lnNum).keys())


# assigning coordinate bounds for grid
    try:
        xmax = kwargs['xmax']
    except (NameError, KeyError):
        if lineCount > miniSurveyThresh and Full == True:
            xmax = 950
        else:
            xmax = np.max(num_x)
    try:
        xmin = kwargs['xmin']
    except (NameError, KeyError):
        if lineCount > miniSurveyThresh and Full == True:
            xmin = 50
        else:
            xmin = np.min(num_x)
    try:
        ymax = kwargs['ymax']
    except (KeyError, NameError):
        if lineCount > miniSurveyThresh and Full == True:
            ymax = 1100
        else:
            ymax = np.max(num_y)
    try:
        ymin = kwargs['ymin']
    except (KeyError, NameError):
        if lineCount > miniSurveyThresh and Full == True:
            ymin = -100
        else:
            ymin = np.min(num_y)

        #
## do some data checks to ensure proper formatting being done

    assert frf[0] == 'FRF', 'input file may not be in the appropriate format, check file'
    assert frf[-1] == 'FRF', 'input file may not be in the appropriate format, check file'
    ##
    xgrid = np.arange(xmin, xmax, dx) #np.min(num_x), np.max(num_x), dx) # array of y coords
    ygrid = np.arange(ymin, ymax, dy) #np.min(num_y), np.max(num_y), dy)  # array of x coords
    grid_data = ngl.natgrid(num_x, num_y, num_z, xgrid, ygrid)
    try:
        if kwargs['plot'] == 1 or kwargs['plot'] == True:
            plt.figure(figsize=(5, 5))
            plt.contourf(xgrid, ygrid, grid_data.T)
            cbar = plt.colorbar()
            cbar.set_label('elevation [m]')
            plt.xlabel('Cross Shore Distance [m]')
            plt.ylabel('Alongshore Distance [m]')
            plt.title('Gridded Bathymetry at the FRF')
            plt.savefig(fname_in+'plot.png')
            plt.close()
            print 'saved figure here: %s' % (fname_in[:-4] + 'plot_%s.png' %dx)
    except (NameError, KeyError):
        pass
    # packaging dict to return 
    product = {'xgrid': xgrid,
               'ygrid': ygrid,
               'zgrid' : grid_data,
               'vor_tup': tup}
    return product

def write_grid(ofname, grid_dict):
    """
    This function takes the gridded product created by frf_grid_prod and writes
    the data to a csv file, it requires a dictionary in with the grid, and a list
    for the x coordinates (xgrid) and a list for the y coordinates (ygrid)
    
    """
    xcoord = grid_dict['xgrid']
    ycoord = grid_dict['ygrid']
    grid = grid_dict['grid'].T  # this grid produced is the
                #transpose of the mesh grids produced below
    xx, yy = np.meshgrid(xcoord, ycoord)
    f= open(ofname,'w')
    for iy in range(0, np.size(grid, axis=0)):
        for ix in range(0, np.size(grid, axis=1)):
            f.write("%f, %f, %f\n" % (xx[iy, ix], yy[iy, ix], grid[iy, ix]))
    f.close()

def creat_vorpoly(tup,ofname):
    """
    This function created voroni polygons based on a scattered data set
    it relys on the pytess library which can be pip installed eg pip install pytess
    tup is a list of tuples of the scatter data [(x_i,y_i)...,(x_i+1, y_i+1)] where
    i is the point number see help pt.voroni
    """     
    vorpolys = pt.voronoi(tup)
        
    from matplotlib import pyplot as plt
    for ii in range(0, np.shape(vorpolys)[0]):
        if vorpolys[ii][0] != None:
            #plt.plot(vorpolys[ii][0], 'rd')
            for pp in range(0, np.shape(vorpolys[ii][1])[0]):
                if pp == np.shape(vorpolys[ii][1])[0]-1:
                    xs = [vorpolys[ii][1][pp][0], vorpolys[ii][1][0][0]]
                    ys = [vorpolys[ii][1][pp][1], vorpolys[ii][1][0][1]]
                else:
                    xs = [vorpolys[ii][1][pp][0], vorpolys[ii][1][pp+1][0]]
                    ys = [vorpolys[ii][1][pp][1], vorpolys[ii][1][pp+1][1]]
                #print pp, ii
                plt.plot(xs, ys, 'k-')
        plt.figure(figsize=(8,8), edgecolor='red')
        plt.ylim([200, 1500])
        plt.xlim([-200, 1100])
        plt.title('example Veroni Tesselation of FRF survey data')
        plt.xlabel('xshore loc')
        plt.savefig(ofname)
        plt.close()
    
# ### actual code
# # path to file, fed to tool
# #fname =  'data/FRF_20140930_1093_FRF_NAVD88_LARC_GPS_UTC_v20151127.csv'
# flist = glob.glob('data/*LARC*.csv')
# #flist = [fname]
# grid_spacing = [50,25,10, 5,1]
# for i in range(0, len(flist)):
#     print 'doing %s ' %flist[i]
#     for x in range(0, len(grid_spacing)):
#
#         fname = flist[i]
#         dxdy = grid_spacing[x]
#         ofname = fname[:-4]+'_grid%sm.xyz' %dxdy
#         print '\n-\n-doing grid spacing of: %s and file %s' %(dxdy,ofname)
#         DT.datetime.now()
#         grid_dict= frf_grid_product(fname, dxdy=dxdy, header=0, plot=1)
#         DT.datetime.now()
#         print 'writing file: %s' %ofname
#         write_grid(ofname, grid_dict)
#         #if x < 2:
#            # print 'starting voroni poly plotting'
#             #creat_vorpoly(grid_dict['vor_tup'], ofname+'vor_%s.png' %dxdy)
#         #DT.datetime.now()
filelist = glob.glob(files)
for fname in filelist:
    gridDict = nn.frf_grid_product(fname, dxdy=10)
    ofname = fname[:-4] + '_grid.txt'
    nn.write_grid(ofname=ofname, grid_dict=gridDict)
    # put plotting functing that saves file here