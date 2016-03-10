# -*- coding: utf-8 -*-
"""
This projects is created to test out the creation of a tool for the natural neighbor interpolation algorithm

assumes
"""
import scipy.interpolate as interp
import pytess as pt
import numpy as np
import csv
import arcpy


# path to file, fed to tool
fname =  'data/FRF_20151014_subset.csv'
header = 1  # length of header lines, place to start finding data


# import file with given path above
#lines = np.genfromtxt(fname, delimiter=',', dtype=None)
raw_x = []
raw_y = []
raw_z = []
with open(fname, 'rb') as f:
    reader = csv.reader(f)
    for row in reader:
        raw_x.append(row[0])
        raw_y.append(row[1])
        raw_z.append(row[2])
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
#assert lines[0][0] == 'x', 'File does not follow established format, check file'
#assert lines[0][1] == 'y', 'File does not follow established format, check file'
#assert lines[0][2] == 'z', 'File does not follow established format, check file'
#assert size (raw_x) == same as size of reader  # not more lines , check make sure file imported proper
#
##
#



# make grid mesh
dx = 5 # output grid resolution in x
dy = 5 # output grid resolution in y
xgrid = np.arange(np.min(num_x), np.max(num_x), dx)
ygrid = np.arange(min(num_y), max(num_y), dy)
xx, yy = np.meshgrid(xgrid, ygrid)

ff = interp.interp2d(xx, yy, num_z, kind = 'cubic')


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
plt.figure(figsize=(5,5), edgecolor='red')
plt.ylim([200, 1500])
plt.xlim([-200, 1100])
plt.title('example Veroni Tesselation of FRF survey data')
plt.xlabel('xshore loc')






#for tt, line in enumerate(vorpolys)
#http://gis.stackexchange.com/questions/16122/creating-shapefile-from-lat-long-values-using-arcpy
#http://help.arcgis.com/en/arcgisdesktop/10.0/help/index.html#//00170000002p000000