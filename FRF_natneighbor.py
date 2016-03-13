# -*- coding: utf-8 -*-
"""
This projects is created to test out the creation of a tool for 
the natural neighbor interpolation algorithm  3 tools are created, based in linux
the frf_grid_product creates the grid using the pyNGL library, then it is written
the voroni creates a tesselation of the surface for each file

assumes
"""

import pytess as pt
import numpy as np
import csv
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
    
    :return:
    
    """
    xmax = 950
    xmin = 50
    ymax = 1100
    ymin = -100
    if dxdy == False:
        dx = kwargs['dx']
        dy = kwargs['dy']
    else:
        dx = dy = dxdy    
    # import file with given path above
    #lines = np.genfromtxt(fname, delimiter=',', dtype=None)
    raw_x = []
    raw_y = []
    raw_z = []
    frf=[]
    with open(fname, 'rb') as f:
        reader = csv.reader(f)
        try:        
            for row in reader:
                frf.append(row[0])
                raw_x.append(row[7])  # x in string format
                raw_y.append(row[8])  # y in string format
                raw_z.append(row[9])  # z in string format
            odd = 0
        except IndexError:  # this is for the subset data sets, they have a different format
            for row in reader:
                raw_x.append(row[0])
                raw_y.append(row[1])
                raw_z.append(row[2])
            odd = 1     
    # initializing values to make strictly nubmers, imported as strings
    num_x = np.zeros(len(raw_x) - header)
    num_y = np.zeros(len(raw_y) - header)
    num_z = np.zeros(len(raw_z) - header)
    tup=[]
    # making strings into floating point numbers
    for ii in range(header,len(raw_x)):
        num_x[ii-1] = float(raw_x[ii])
        num_y[ii-1] = float(raw_y[ii])
        num_z[ii-1] = float(raw_z[ii])
        tup.append((float(raw_x[ii]),float(raw_y[ii])))  #, float(raw_z[ii])))

#
## do some data checks to ensure proper formatting being done
# these are wrong below
    assert frf[0] == 'FRF', 'input file may not be in the appropriate format, check file'
    assert frf[-1] == 'FRF', 'input file may not be in the appropriate format, check file'
    ##
# Remove for tool
    if odd == 1:
        xmin = np.min(num_x)
        xmax = np.max(num_x)
        ymin = np.min(num_y)
        ymax = np.max(num_y)
    xgrid = np.arange(xmin, xmax, dx) #np.min(num_x), np.max(num_x), dx) # array of y coords
    ygrid = np.arange(ymin, ymax, dy) #np.min(num_y), np.max(num_y), dy)  # array of x coords
    grid_data = ngl.natgrid(num_x, num_y, num_z, xgrid, ygrid)
    try:
        if kwargs['plot']==1:
            plt.figure(figsize=(5,5))            
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
    product = { 'xgrid': xgrid,
                'ygrid': ygrid,
                'grid' : grid_data,
                'vor_tup': tup,
              }
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
            f.write("%f, %f, %f\n" % (xx[iy,ix], yy[iy,ix], grid[iy,ix]))
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
    
### actual code
# path to file, fed to tool
#fname =  'data/FRF_20140930_1093_FRF_NAVD88_LARC_GPS_UTC_v20151127.csv'
flist = glob.glob('data/*LARC*.csv')
#flist = [fname]
grid_spacing = [50,25,10, 5,1]
for i in range(0, len(flist)):
    print 'doing %s ' %flist[i]    
    for x in range(0, len(grid_spacing)):
        
        fname = flist[i]
        dxdy = grid_spacing[x]
        ofname = fname[:-4]+'_grid%sm.xyz' %dxdy
        print '\n-\n-doing grid spacing of: %s and file %s' %(dxdy,ofname)
        DT.datetime.now()
        grid_dict= frf_grid_product(fname, dxdy=dxdy, header=0, plot=1)
        DT.datetime.now()        
        print 'writing file: %s' %ofname
        write_grid(ofname, grid_dict)
        #if x < 2:
           # print 'starting voroni poly plotting'
            #creat_vorpoly(grid_dict['vor_tup'], ofname+'vor_%s.png' %dxdy)
        #DT.datetime.now()
print 'program finished'