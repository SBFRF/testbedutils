import numpy as np
import pyproj

def FRF2ncsp(xFRF, yFRF):
    """
    this function makes NC stateplane out of X and Y FRF coordinates,
    based on kent Hathaway's code, bill birkmeir's calculations .
         #  15 Dec 2014
    #  Kent Hathaway.
    #  Translated from Matlab to python 2015-11-30 - Spicer Bak
    #
    #  Uses new fit (angles and scales) Bill Birkemeier determined in Nov 2014
    #
    #  This version will determine the input based on values, outputs FRF, lat/lon,
    #  and state plane coordinates.  Uses NAD83-2011.
    #
    #  IO:
    #  p1 = FRF X (m), or Longitude (deg + or -), or state plane Easting (m)
    #  yFRF = FRF Y (m), or Latitude (deg), or state plane Northing (m)
    #
    #  X = FRF cross-shore (m)
    #  Y = FRF longshore (m)
    #  ALat = latitude (decimal degrees)
    #  ALon = longitude (decimal degrees, positive, or W)
    #  spN = state plane northing (m)
    #  spE = state plane easting (m)

    NAD83-86	2014
    Origin Latitude          36.1775975
    Origin Longitude         75.7496860
    m/degLat             110963.357
    m/degLon              89953.364
    GridAngle (deg)          18.1465
    Angle FRF to Lat/Lon     71.8535
    Angle FRF to State Grid  69.9747
    FRF Origin Northing  274093.1562
    Easting              901951.6805

    #  Debugging values
    xFRF=566.93;  yFRF=515.11;  % south rail at 1860
    ALat = 36.1836000
    ALon = 75.7454804
    p2= 36.18359977;  p1=-75.74548109;
    SP:  p1 = 902307.92; 	p2 = 274771.22;
    :param xFRF:
    :param yFRF:
    :return:
    """
    r2d = 180.0 / np.pi;

    Eom = 901951.6805;  # % E Origin State Plane
    Nom = 274093.1562;  # % N Origin State Plane
    spAngle = (90 - 69.974707831) / r2d
    X = xFRF
    Y = yFRF

    R = np.sqrt(X ** 2 + Y ** 2)
    Ang1 = np.arctan2(X, Y)  # % CW from Y
    #  to state plane
    Ang2 = Ang1 - spAngle
    AspN = R * np.cos(Ang2)
    AspE = R * np.sin(Ang2)
    spN = AspN + Nom
    spE = AspE + Eom
    out = {'xFRF': xFRF,  'yFRF': yFRF, 'StateplaneE': spE, 'StateplaneN': spN}
    return out

def ncsp2FRF(p1, p2):
    """
    this function converts nc StatePlane (3200 fips) to FRF coordinates
    based on kent Hathaways Code
        #
    #  15 Dec 2014
    #  Kent Hathaway.
    #  Translated from Matlab to python 2015-11-30 - Spicer Bak
    #
    #  Uses new fit (angles and scales) Bill Birkemeier determined in Nov 2014
    #
    #  This version will determine the input based on values, outputs FRF, lat/lon,
    #  and state plane coordinates.  Uses NAD83-2011.
    #
    #  IO:
    #  p1 = FRF X (m), or Longitude (deg + or -), or state plane Easting (m)
    #  p2 = FRF Y (m), or Latitude (deg), or state plane Northing (m)
    #
    #  X = FRF cross-shore (m)
    #  Y = FRF longshore (m)
    #  ALat = latitude (decimal degrees)
    #  ALon = longitude (decimal degrees, positive, or W)
    #  spN = state plane northing (m)
    #  spE = state plane easting (m)

    NAD83-86	2014
    Origin Latitude          36.1775975
    Origin Longitude         75.7496860
    m/degLat             110963.357
    m/degLon              89953.364
    GridAngle (deg)          18.1465
    Angle FRF to Lat/Lon     71.8535
    Angle FRF to State Grid  69.9747
    FRF Origin Northing  274093.1562
    Easting              901951.6805

    #  Debugging values
    p1=566.93;  p2=515.11;  % south rail at 1860
    ALat = 36.1836000
    ALon = 75.7454804
    p2= 36.18359977;  p1=-75.74548109;
    SP:  p1 = 902307.92; 	p2 = 274771.22;

    :param spE:
    :param spN:
    :return:
    """
    r2d = 180.0 / np.pi;

    Eom = 901951.6805;  # % E Origin State Plane
    Nom = 274093.1562;  # % N Origin State Plane
    spAngle = (90 - 69.974707831) / r2d

    spE = p1
    spN = p2  # designating stateplane vars
    # to FRF coords
    spLengE = p1 - Eom
    spLengN = p2 - Nom
    R = np.sqrt(spLengE ** 2 + spLengN ** 2)
    Ang1 = np.arctan2(spLengE, spLengN)
    Ang2 = Ang1 + spAngle
    # to FRF
    X = R * np.sin(Ang2)
    Y = R * np.cos(Ang2)
    # to Lat Lon
    ans = {'xFRF': X,
           'yFRF': Y,
           'StateplaneE': spE,
           'StateplaneN': spN}
    return ans

def ncsp2LatLon(spE, spN):
    """
    This function uses pyproj to convert


    test points taken from conversions made in USACE SMS modeling system

    nc stateplane  meters NAD83
    spE1 = 901926.2 m
    spN1 = 273871.0 m
    Lon1 = -75.75004989
    Lat1 =  36.17560399

    spE2 = 9025563.9 m
    spN2 = 276229.5 m
    lon2 = -75.47218285
    lat2 =  36.19666112
    :param spE: easting - assumed north carolina state plane Meters
    :param spN: northing - assumed north carolina state plane meters
    :return: dictionary with original coords and output of latitude and longitude.
    """
    EPSG = 3358  # taken from spatialreference.org/ref/epsg/3358
    # NC stateplane NAD83
    spNC = pyproj.Proj(init="epsg:%s" %EPSG)
    LL = pyproj.Proj(init='epsg:4269')  # epsg for NAD83 projection
    lon, lat = pyproj.transform(spNC, LL, spE, spN)
    ans = {'lon': lon, 'lat': lat, 'StateplaneE': spE, 'StateplaneN': spN}
    return ans

def LatLon2ncsp(lon, lat):
    """
      This function uses pyproj to convert longitude and latitude to stateplane

      test points taken from conversions made in USACE SMS modeling system

        nc stateplane  meters NAD83
        spE1 = 901926.2 m
        spN1 = 273871.0 m
        Lon1 = -75.75004989
        Lat1 =  36.17560399

        spE2 = 9025563.9 m
        spN2 = 276229.5 m
        lon2 = -75.47218285
        lat2 =  36.19666112
      :param lon: geographic longitude (NAD83)  decimal degrees
      :param lat: geographic longitude (NAD83)  decimal degrees
    :return:  output dictionary with original coords and output of NC stateplane FIPS 3200
    """
    EPSG = 3358  # taken from spatialreference.org/ref/epsg/3358
    # NC stateplane NAD83
    spNC = pyproj.Proj(init="epsg:%s" %EPSG)
    LL = pyproj.Proj(init='epsg:4269')  # epsg for NAD83 projection
    spE, spN = pyproj.transform(LL, spNC, lon, lat)
    ans = {'lon': lon, 'lat': lat, 'StateplaneE': spE, 'StateplaneN': spN}
    return ans

def FRFcoord(p1, p2):
    """
    updated FRF coord in python, using kent's original code as guide but converting pyproj for all
    conversions between state plane and Lat Lon,  Then all conversions between stateplane and
    FRF coordinates are done using kents original geometry.

    #
    #  15 Dec 2014
    #  Kent Hathaway.
    #  Translated from Matlab to python 2015-11-30 - Spicer Bak
    #
    #  Uses new fit (angles and scales) Bill Birkemeier determined in Nov 2014
    #
    #  This version will determine the input based on values, outputs FRF, lat/lon,
    #  and state plane coordinates.  Uses NAD83-2011.
    #
    #  IO:
    #  p1 = FRF X (m), or Longitude (deg + or -), or state plane Easting (m)
    #  p2 = FRF Y (m), or Latitude (deg), or state plane Northing (m)
    #
    #  X = FRF cross-shore (m)
    #  Y = FRF longshore (m)
    #  ALat = latitude (decimal degrees)
    #  ALon = longitude (decimal degrees, positive, or W)
    #  spN = state plane northing (m)
    #  spE = state plane easting (m)

    NAD83-86	2014
    Origin Latitude          36.1775975
    Origin Longitude         75.7496860
    m/degLat             110963.357
    m/degLon              89953.364
    GridAngle (deg)          18.1465
    Angle FRF to Lat/Lon     71.8535
    Angle FRF to State Grid  69.9747
    FRF Origin Northing  274093.1562
    Easting              901951.6805

    #  Debugging values
    p1=566.93;  p2=515.11;  % south rail at 1860
    ALat = 36.1836000
    ALon = 75.7454804
    p2= 36.18359977;  p1=-75.74548109;
    SP:  p1 = 902307.92; 	p2 = 274771.22;
    """

    # Determine Data type
    if np.floor(abs(p1)) == 75 and np.floor(p2) == 36:  # lat/lon input
        sp = LatLon2ncsp(p1, p2)
        frf = ncsp2FRF(sp['StateplaneE'], sp['StateplaneN'])
        coordsOut = {'xFRF': frf['xFRF'], 'yFRF': frf['yFRF'], 'StateplaneE':sp['StateplaneE'],
                     'StateplaneN': sp['StateplaneN'], 'Lat': p2, 'Lon':p1}

    elif (p1 > 800000) and p2 > 200000:  # state plane input
        frf = ncsp2FRF(p1, p2)     # convert state plane to FRF
        ll = ncsp2LatLon(p1, p2)  # convert state plane to Lat Lon
        coordsOut = {'xFRF': frf['xFRF'], 'yFRF': frf['yFRF'], 'StateplaneE': p1,
                     'StateplaneN': p2, 'Lat': ll['lat'], 'Lon':ll['lon']}

    elif (p1 > -10000 and p1 < 10000) and (p2 > -10000 and p2 < 10000):  # FRF input
        # this is FRF in
        sp = FRF2ncsp(p1, p2)
        ll = ncsp2LatLon(sp['StateplaneE'], sp['StateplaneN'])
        coordsOut = {'xFRF': p1, 'yFRF': p2, 'StateplaneE': sp['StateplaneE'],
                     'StateplaneN': sp['StateplaneN'], 'Lat': ll['lat'], 'Lon':ll['lon']}

    else:
        print '<<ERROR>> Cound not determine input type, returning NaNs'
        coordsOut = {'xFRF': float('NaN'), 'yFRF': float('NaN'), 'StateplaneE': float('NaN'),
             'StateplaneN': float('NaN'), 'Lat': float('NaN'), 'Lon':float('NaN')}
    return coordsOut