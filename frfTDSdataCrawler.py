"""Module associated with collecting and querying gauge locations from FRF/CHL TDS server"""
import csv, pickle, urllib.request, warnings, sys, time
import numpy as np
import xml.etree.ElementTree as ET
from netCDF4 import Dataset, num2date
import progressbar
import datetime as DT


def query(startDate, endDate, type, sensor=None, inputName='database', outputName=None):
    """Will querey gauge location lookup table to provide locations of data within.

    Args:
        startDate: datetime object that is to be inclusively searched (>=) for finding locations
        endDate: datetime object that is to be inclusively searched (<=) for finding locations
        type: type of data to be searched (eg currents, waves, etc). defines which folder to search through on oceanography.
        sensor: gauge string name to look through
        inputName: pickle file name that has data stored in it (default='database.p')
        outputName: will save query as csv file if file name is given. program will append csv suffix. if None
            program will not export.  To export, must be text file extension (default=None)

    Returns:
        dictionary with keys ['DateStart', 'DateEnd',
         'Type': type of data
         'Sensor': sensor name string
         'Lat': lat of gauge location
         'Lon': lon of gague location
         'Url': OPeNDAP querey location

    """
    print('Querying')

    with open(inputName + '.p', 'rb') as outfile:
        database = pickle.load(outfile)
    # convert datetime objects to numerical year
    assert isinstance(startDate, DT.datetime)
    assert isinstance(endDate, DT.datetime)
    # 1st query = Date
    dateStartList = np.array(database['DateStart']).astype(DT.datetime)
    dateEndList = np.array(database['DateEnd']).astype(DT.datetime)
    I1 = np.logical_and(dateEndList >= startDate, dateEndList <= endDate)
    I2 = np.logical_and(dateStartList >= startDate, dateStartList <= endDate)
    I3 = np.logical_and(dateStartList <= startDate, dateEndList >= endDate)
    I = np.logical_or(I1, I2)
    I = np.logical_or(I, I3)

    # 2nd query = Type
    typeList = np.array(database['Type'])
    I = np.logical_and(typeList == type, I)
    allList = [dateStartList[I], dateEndList[I], typeList[I]]

    sensorList = np.array(database['Sensor'])
    if sensor is not None:
        # 3rd query = Sensor
        I = np.logical_and(sensorList == sensor, I)
        # allList = [dateStartList[I], dateEndList[I], typeList[I], sensorList[I]]

    latList = np.array(database['Lat'])
    # # 4th query = lat
    # I_ = np.logical_and(latList >= lat[0], latList <= lat[1])
    # I = np.logical_and(I, I_)
    #
    lonList = np.array(database['Lon'])
    # # 5th query = lon
    # I_ = np.logical_and(lonList >= lon[0], lonList <= lon[1])
    # I = np.logical_and(I, I_)

    urlList = np.array(database['Url'])[I]
    allList.append(urlList)
    allList = [dateStartList[I], dateEndList[I], typeList[I], sensorList[I],
               latList[I], lonList[I], urlList]

    queryData = dict()

    i = 0
    for key in database:
        queryData[key] = allList[i]
        i += 1

    # data are already in datetime
    # convert DateStart and DateEnd back to datetime objects
    # for var in ['DateStart', 'DateEnd']:
        #queryData[var] = np.array([DT.datetime.strptime(str(date), '%Y%m') for date in queryData[var]])
        
    # Stop and save data as necessary
    if outputName is not None:
        print('Saving query')

        with open(outputName + '.csv', 'w') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=queryData.keys())
            writer.writeheader()
            for i in range(len(urlList)):
                row = dict()
                for key in queryData:
                    row[key] = queryData[key][i]
                writer.writerow(row)

    return queryData


def buildLookupTable(server, outputName='database'):
    """main code structure"""
    t = DT.datetime.now()
    datatype = ['waves', 'currents']
    print('Get URLs')
    urlList = getUrls(server, datatype)
    print('building database')
    database, errorbase = buildDatabase(urlList)
    print('sorting database')
    database = sortDatabase(database)
    print('collecting Lon/Lat info')
    database = collectLatLon(database)

    print('saving data')
    saveBinary(outputName, database)

    saveCsv(outputName, database)
    
    showErrors(errorbase)
    print("process took {:.1f} minutes".format((DT.datetime.now()-t).total_seconds()/60))

def getUrls(server, datatype):

    urlList = [None] * 999999
    istart = 0
    
    for datatype_ in datatype:
        urlList, istart = getUrlsEachType(server, datatype_, urlList, istart)
    
    urlList[istart:] = []
    
    return urlList
    
def getUrlsEachType(server, datatype, urlList, istart):
    
    if server == 'http://134.164.129.55':
        urlMain = (server + '/thredds/catalog/FRF/oceanography/{}/catalog.xml'
            .format(datatype))
    elif server == 'https://chldata.erdc.dren.mil':
        urlMain = (server + '/thredds/catalog/frf/oceanography/{}/catalog.xml'
            .format(datatype))
    else:
        print('    Unknown server')
        quit()
    
    tree = ET.parse(urllib.request.urlopen(urlMain))
    root = tree.getroot()
    
    i = istart
    bar = progressbar.ProgressBar(maxval=len(root[-1]),
                    widgets=[progressbar.Bar('.', '[', ']'), ' ', progressbar.Percentage()])
    bar.start()
    for ii, child in enumerate(root[-1]):
        bar.update(ii+1)
        if child.tag[-10:] != 'catalogRef':
            continue
        urlChild = urlMain[:-11] + child.attrib[child.keys()[0]]
        tree = ET.parse(urllib.request.urlopen(urlChild))
        root = tree.getroot()
        for gchild in root[-1]:
            if gchild.tag[-10:] != 'catalogRef':
                continue
            urlGchild = urlChild[:-11] + gchild.attrib[gchild.keys()[0]]
            t, itMax = 0, 10
            while t < itMax:
                try:
                    tree = ET.parse(urllib.request.urlopen(urlGchild))
                    break
                except:
                    time.sleep(10)
                    t += 1
                    continue
            root = tree.getroot()
            for ggchild in root[-1]:
                if ggchild.tag[-7:] != 'dataset':
                    continue
                urlList[i] = '{}/thredds/dodsC/{}'.format(
                    server, ggchild.attrib['urlPath'])
                # print('Found {}'.format('_'.join(
                #     (urlList[i].split('/')[-1].split('_')[2:]))))
                i += 1
    bar.finish()
    return urlList, i
    
def buildDatabase(urlList):

    database = dict()

    headers = ['DateStart', 'DateEnd', 'Type', 'Sensor', 'Lat', 'Lon', 'Url']

    for header in headers:
        database[header] = [None] * len(urlList)
    
    errorbase = dict()
    
    headers = ['OpeningError', 'LatLonError']
    
    for header in headers:
        errorbase[header] = [None] * 999
    
    i = 0
    j = 0
    k = 0
    bar = progressbar.ProgressBar(maxval=len(urlList),
                    widgets=[progressbar.Bar('.', '[', ']'), ' ', progressbar.Percentage()])
    bar.start()
    for ii, url in enumerate(urlList):
        for attempt in range(10):
            try:
                # print('Parsing {}'.format('_'.join(
                #             (url.split('/')[-1].split('_')[2:]))))
                bar.update(ii+1)
                try:
                    rootgrp = Dataset(url)
                except OSError:
                    print('    Opening error')
                    print('    Url = ' + url)
                    errorbase['OpeningError'][j] = url
                    j += 1
                    continue

                varList = list(rootgrp.variables)
                
                if 'latitude' in varList:
                    lat = rootgrp['latitude'][:]
                    lon = rootgrp['longitude'][:]
                elif 'lat' in varList:
                    lat = rootgrp['lat'][:]
                    lon = rootgrp['lon'][:]
                elif 'lidarLatitude' in varList:
                    lat = rootgrp['lidarLatitude'][:]
                    lon = rootgrp['lidarLongitude'][:]
                else:
                    print('    Lat/lon error')
                    errorbase['LatLonError'][k] = url
                    k += 1
                    continue

                if (type(lat) is np.ma.core.MaskedConstant or type(lon) is
                np.ma.core.MaskedConstant):
                    print('    Lat/lon error')
                    errorbase['LatLonError'][k] = url
                    k += 1
                    continue

                DateStart = num2date(rootgrp['time'][0], 'seconds since 1970-01-01')
                DateEnd = num2date(rootgrp['time'][-1], 'seconds since 1970-01-01')
        
                rootgrp.close()

                # date = int(url.split('_')[-1][:-3])
                sensor = url.split('_')[2]

                database['DateStart'][i] = DateStart
                database['DateEnd'][i] = DateEnd
                database['Type'][i] = url.split('_')[1]
                database['Sensor'][i] = sensor
                database['Lat'][i] = round(float(lat), 4)
                database['Lon'][i] = round(float(lon), 4)
                database['Url'][i] = url
                i += 1
            except:
                print('\n\r retrying {}'.format('_'.join((url.split('/')[-1].split('_')[2:]))))
                continue
            else:
                break
    bar.finish()
    for key in database:
        database[key][i:] = []
        
    errorbase['OpeningError'][j:] = []
    errorbase['LatLonError'][k:] = []

    return database, errorbase

def sortDatabase(database):

    # get sorting indices - sort by Sensor, then by Date
    ind = np.lexsort((database['DateStart'], database['Sensor'],
        database['Type']))[::-1]
    bar = progressbar.ProgressBar(maxval=len(database),
                    widgets=[progressbar.Bar('.', '[', ']'), ' ', progressbar.Percentage()])
    bar.start()
    # sort each column using sorting indices
    for ii, header in enumerate(database):
        database[header] = list(np.array(database[header])[ind])
        bar.update(ii+1)
    bar.finish()
    return database

def collectLatLon(database):

    dbNew = dict()

    headers = ['DateStart', 'DateEnd', 'Type', 'Sensor', 'Lat', 'Lon', 'Url']

    for header in headers:
        dbNew[header] = [None] * 999

    dbNew['DateEnd'][0] = database['DateEnd'][0]
    j = 0
    
    if server == 'http://134.164.129.55':
        frfTag = 'FRF'
    else:
        frfTag = 'frf'

    bar = progressbar.ProgressBar(maxval=len(database['Url']),
                    widgets=[progressbar.Bar('.', '[', ']'), ' ', progressbar.Percentage()])
    bar.start()
    for i in range(len(database['Url'])):
        bar.update(i+1)
        if (i == len(database['Url']) - 1 or database['Lat'][i] !=
        database['Lat'][i + 1] or database['Lon'][i]
        != database['Lon'][i + 1]):
            dbNew['DateStart'][j] = database['DateStart'][i]
            dbNew['Type'][j] = database['Type'][i]
            dbNew['Sensor'][j] = database['Sensor'][i]
            dbNew['Lat'][j] = database['Lat'][i]
            dbNew['Lon'][j] = database['Lon'][i]
            
            url = ('{}/thredds/catalog/{}/oceanography/waves/{}/catalog.xml'
                .format(server, frfTag, database['Sensor'][i]))
            t, itMax = 0, 10
            while t < itMax:
                try:
                    tree = ET.parse(urllib.request.urlopen(url))
                    break
                except:
                    time.sleep(10)
                    t += 1
                    continue
            root = tree.getroot()
            foundNcml = False
            for child in root[-1]:
                if child.tag[-7:] == 'dataset':
                    dbNew['Url'][j] = '{}/thredds/dodsC/{}'.format(server,
                        child.attrib['urlPath'])
                    foundNcml = True
                    break
      
            if not foundNcml:
                dbNew['Url'][j] = database['Url'][i]
            
            if i != len(database['Url']) - 1:
                j += 1
                dbNew['DateEnd'][j] = database['DateEnd'][i + 1]
            else:
                break
    bar.finish()
    for key in dbNew:
        dbNew[key][j + 1:] = []

    return dbNew

def saveBinary(outputName, database):

    print('    Saving binary')

    with open(outputName + '.p', 'wb') as outfile:
        pickle.dump(database, outfile)

def saveCsv(outputName, database):

    print('    Saving csv')

    with open(outputName + '.csv', 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=database.keys())
        writer.writeheader()
        for i in range(len(database['Url'])):
            row = dict()
            for header in database:
                row[header] = database[header][i]
            writer.writerow(row)

def showErrors(errorbase):
    print('The following could not be opened:')
    for url in errorbase['OpeningError']:
        print(url)
    
    print('Lat/lon could not be found in the following:')
    for url in errorbase['LatLonError']:
        print(url)

if __name__ == '__main__':
    assert len(sys.argv) > 1, 'needs input'
    assert sys.argv[1].lower()  in ['chl', 'frf'], "input argument must be in ['chl', 'frf']"
    if sys.argv[-1].lower() == 'chl':
        server = 'https://chldata.erdc.dren.mil'
    else: # sys.argv[-1].lower() == 'frf':
        server = 'http://134.164.129.55'

    # if len(sys.argv) > 2:
    #     dtype = sys.argv[2:]

    buildLookupTable(server)
