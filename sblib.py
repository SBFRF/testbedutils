# -*- coding: utf-8 -*-
"""
Author: spicer bak, phd
Contact: spicer.bak@usace.army.mil
Association: USACE CHL Field Research Facility


my own library with useful functions and tips from which i found to be helpful
on various general functions that don't fit into specific other codes
this folder needs to be added to sys.path to use
documented 12/2/17
"""
import numpy as np
import datetime as DT
import csv, warnings
import math
import netCDF4 as nc


def makegif(flist, ofname, size=None, dt=0.5):
    """
    This function uses imageio to create gifs from a list of images

    kwargs for mimwrite http://imageio.readthedocs.org/en/latest/format_gif.html#gif

    :param flist: a sorted list of files to be made into gifs (including path)
    :param ofname: output gif filename (including path)
    :param size: size of pictures (default not resized)
    :param loop: number of loops to do, 0 is default and infinite
    :return:
        will write a gif to ofname location
    """
    # images = [Image.open(fn) for fn in flist]
    #
    # for im in images:
    #     im.thumbnail(size, Image.ANTIALIAS)
    # images2gif.writeGif(ofname, images, duration=dt, nq=15)
    import imageio
    images = []
    if size != None:
        for im in images:
            im.thumbnail(size, Image.ANTIALIAS)
    for filename in flist:
        images.append(imageio.imread(filename))
    imageio.mimwrite(ofname, images, duration=dt)

def find_nearest(array, value):
    """
    Function looks for value in array and returns the closest array value
    (to 'value') and index of that value
    :param array: array to to find things in
    :param value: value to search against
    :return returns a number from the array value, that is closest to value input
    """
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx

def findUniqueFromTuple(a, axis=1):
    """
    This function finds the unique values of a multi dimensional tuple (quickly)
    :param a: an array of multidimensional size
    :param axis: this is the axis that it looks for unique values using (default is horizontal)
    :return: array of unique tuples
        warning, values are not sorted
    """
    if axis == 0:
        a = a.T

    b = np.ascontiguousarray(a).view(np.dtype((np.void, a.dtype.itemsize * a.shape[1])))
    _, idx = np.unique(b, return_index=True)

    unique_a = a[idx]
    return unique_a

def FRFcoord(p1, p2):
    """
    place older
    """
    raise ImportError('please use the function in sblib.geoprocess')

def findbtw(data, lwth, upth, type=0):
    """
    This function finds both values and indicies of a list values between two values
    :TODO: probably could be improved by using boolean compares


    :param upth: = upper level threshold
    :param lwth: = lower level threshold
    :param data: = list (or numpy array?)
    :param type:
        0 = non inclusive  ie. lwth < list <  upth
        1 = low incluisve  ie. lwth <=list <  upth
        2 = high inclusive ie. lwth < list <= upth
        3 = all inclusive  ie  lwth <=list <= upth
    :return
        indicies returns idices that meet established criteria
        values   returns associated values from da
    """
    indices = []
    vals = []
    shp = np.shape(data)
    if len(shp) == 2:
        for i, in enumerate(data):
            for j, elem in enumerate(range):
                if type == 0:
                    if elem < upth and elem > lwth:
                        indices.append((i, j))
                        vals.append(elem)
                elif type == 1:
                    if elem < upth and elem >= lwth:
                        indices.append((i, j))
                        vals.append(elem)
                elif type == 2:
                    if elem <= upth and elem > lwth:
                        indices.append((i, j))
                        vals.append(elem)
                elif type == 3:
                    if elem <= upth and elem >= lwth:
                        indices.append((i, j))
                        vals.append(elem)
    if len(shp) == 1:
        for j, elem in enumerate(data):
            if type == 0:
                if elem < upth and elem > lwth:
                    indices.append((j))
                    vals.append(elem)
            elif type == 1:
                if elem < upth and elem >= lwth:
                    indices.append((j))
                    vals.append(elem)
            elif type == 2:
                if elem <= upth and elem > lwth:
                    indices.append((j))
                    vals.append(elem)
            elif type == 3:
                if elem <= upth and elem >= lwth:
                    indices.append((j))
                    vals.append(elem)

    return indices, vals

class Bunch(object):
    """
    allows user to access dictionary data from 'object'
    instead of object['key']
    do x = sblib.Bunch(object)
    x.key

    do x = Bunch(object)
    x.key
    :param object: input dictionary
    :return class: with objects being the same as the keys from dictionary
    """

    def __init__(self, aDict):
        self.__dict__.update(aDict)

def myround(x, base=5):
    """
    This function will round any value to a multiple of the base,

    :param x: values to be rounded:
    :param base: this is the value by which x is rounded to a multiple of
        ie base = 10  x = [4, 8, 2, 12]  returns [0,10,0,10]
    :return: np.array of floating point numbers rounded to a multiple of base
    """
    x = np.array(x, dtype=float)
    return base * np.round(x/base)

def roundtime(timeIn=None, roundTo=60):
    """"
    Round a datetime object to any time lapse in seconds

    :param dt : datetime.datetime object, default now.
    :param roundTo : Closest number of seconds to round to, default 1 minute.

       Author: Thierry Husson 2012 - Use it as you want but don't blame me.

       modified by SB to include lists of datetime dataList,
       returned as a list if it came in as a list, if it came in as a datetime object
       it is returned as such
    :return
        list/array of times rounded to 'round to value'
    """
    # making dt a list

    if np.size(timeIn) > 1:
        dtlist = timeIn
    elif np.size(timeIn) == 1:
        if type(timeIn) == np.ndarray and timeIn.shape == (1,):
            dtlist = timeIn.tolist()
        else:
            dtlist = [timeIn]
    elif np.size(timeIn) == None:
        dtlist = [DT.datetime.now()]
        # checking to make datetime
        # if type(dt[0] != DT.datetime):
        # trying to convert epoch time to datetime object  if need be but doen'st solve problem currently working on

    # looping through list
    for ii, dt in enumerate(dtlist):
        seconds = (dt - dt.min).seconds
        # // is a floor division, not a comment on following line:
        rounding = (seconds + roundTo / 2) // roundTo * roundTo
        dtlist[ii] = dt + DT.timedelta(0, rounding - seconds, -dt.microsecond)
    # return data as it came in
    if len(dtlist) == 1:
        if type(timeIn) == np.ndarray:
            dtlist = np.array(dtlist)  # if it came in as an array, return it as an array
        else:
            dtlist = dtlist[0]  # if it came in as a single element, return it as such
    return dtlist

def cart2pol(x, y):
    """
        this translates from cartesian coords to polar coordinates (radians)

    :param x: x componant
    :param y: y componant
    :return:
        r radial componant
        theta angular compoanat (returned in radian)
    """
    r = np.sqrt(x ** 2 + y ** 2)
    theta = np.arctan2(y, x)
    return r, theta

def pol2cart(r, theta):
    """
    this translates from polar coords (radians) to polar coordinates
    assumed radian input for theta

    :param r: speed, magnatude
    :param theta:  direction (in radians)
    :return:
        x - componant
        y - componant
    """
    if (np.max(theta) > 2 * np.pi).any():
        print 'Warning polar2cart assumes radian direction in, angles found above 2pi'
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y

def angle_correct(angle_in, rad=0):
    """
    this function takes angles in that are both positve and negative
    and corrects them to posivitve only

    :param angle_in:
    :param rad: radian =0 input angles are in degrees
                radian =1 input anglesa are in radian
    :return:
        array of corrected angles in
    """
    angle_in = np.array(angle_in)
    try:
        assert (angle_in == 0).all() is not True, 'All of the Angles are 0, cannot correct'
    except AssertionError:
        return angle_in
    if rad == 0:
        if (angle_in == 0).all():
            warnings.warn('WARNING - Correcting angles of Zero')
        elif (np.abs(angle_in) < 2 * np.pi).all():
            warnings.warn(' WARNING angles are all < 2Pi , ensure that angles are in degrees not radians')

        shape = np.shape(angle_in)
        if len(shape) == 0:
            posmask = angle_in >= 360
            negmask = angle_in < 0
            while negmask.any() or posmask.any():
                if negmask.any() == True:
                    angle_in += 360
                elif posmask.any() == True:
                    angle_in -= 360
                posmask = angle_in >= 360
                negmask = angle_in < 0
        if len(shape) == 1:
            posmask = angle_in >= 360
            negmask = angle_in < 0
            while negmask.any() or posmask.any():
                if negmask.any():  # filter negs out
                    idxneg = np.where(negmask)
                    angle_in[idxneg] += 360
                if posmask.any():  # filter overly positives out
                    idxpos = np.where(posmask)
                    angle_in[idxpos] -= 360
                posmask = angle_in >= 360
                negmask = angle_in < 0
        elif len(shape) == 2:
            for ii in range(0, np.size(angle_in, axis=0)):
                angle_in_2 = np.zeros((np.size(angle_in[ii, :])))  # initializing
                angle_in_2 = angle_in[ii, :]  # taking small chunk 1D array
                posmask = angle_in_2 >= 360  # seeing what's over 360
                negmask = angle_in_2 < 0  # seeing what's under 0
                while negmask.any() or posmask.any():
                    if negmask.any():  # filter negs out
                        idxneg = np.where(negmask)  # finding ids of where 
                        if np.size(angle_in_2) == 1 and negmask == True:  # if there's only 1 instance
                            angle_in_2 += 360
                        else:
                            angle_in_2[idxneg] += 360
                    if posmask.any():  # filter overly positives out
                        idxpos = np.where(posmask)
                        if np.size(angle_in_2) == 1 and posmask == True:
                            angle_in_2 -= 360
                        else:
                            angle_in_2[idxpos] -= 360
                    posmask = angle_in_2 >= 360
                    negmask = angle_in_2 < 0
                angle_in[ii, :] = angle_in_2

        elif len(shape) == 3:
            for yy in range(0, np.size(angle_in, axis=1)):
                angle_in_3 = np.zeros(np.size(angle_in, axis=1))
                angle_in_3 = angle_in[:, yy, :]
                for ii in range(0, np.size(angle_in, axis=0)):
                    angle_in_2 = np.zeros((np.size(angle_in_3[ii, :])))  # initializing
                    angle_in_2 = angle_in_3[ii, :]  # taking small chunk 1D array
                    posmask = angle_in_2 >= 360  # seeing what's over 360
                    negmask = angle_in_2 < 0  # seeing what's under 0
                    while negmask.any() or posmask.any():
                        if negmask.any():  # filter negs out
                            idxneg = np.where(negmask)  # finding ids of where 
                            if np.size(angle_in_2) == 1 and negmask == True:  # if there's only 1 instance
                                angle_in_2 += 360
                            else:
                                angle_in_2[idxneg] += 360
                        if posmask.any():  # filter overly positives out
                            idxpos = np.where(posmask)
                            if np.size(angle_in_2) == 1 and posmask == True:
                                angle_in_2 -= 360
                            else:
                                angle_in_2[idxpos] -= 360
                        posmask = angle_in_2 >= 360
                        negmask = angle_in_2 < 0
                    angle_in_3[ii, :] = angle_in_2
                angle_in[:, yy, :] = angle_in_3
    else:
        print '<<ERROR>> this function only takes angles in as degrees right now'
        raise
    assert (angle_in < 360).all() and (angle_in >= 0).all(), 'The angle correction function didn''t work properly'
    return angle_in

def statsBryant(observations, models):
    """
    This function does Non-Directional Statsistics
    These statistics are from the Bryant Wave stats CHETN - I - 91
    :param observations: array of observational data
    :param models:  array of model data
    :return: dictiornay
        :key  'bias' average of residuals
        :key 'RMSEdemeaned': RMSEdemeaned,
        :key 'RMSE': R oot mean square error
        :key  'RMSEnorm': normalized root Mean square error also Percent RMSE
        :key  'scatterIndex': ScatterIndex
        :key  'symSlope': symetrical slope
        :key  'corr': R^2 --- coefficient of determination
        :key  'PscoreWilmont': performance score developed by Wilmont
        :key  'PscoreIMEDS':  see Hanson 2007 for description
        :key  'residuals': model - observations
        :
    """
    obsNaNs = np.argwhere(np.isnan(observations))
    modNaNs = np.argwhere(np.isnan(models))
    if len(obsNaNs) > 0:
        print 'warning found nans in bryant stats'
        observations = np.delete(observations, obsNaNs)
        models = np.delete(models, obsNaNs)  # removing corresponding model data, that cnanot be compared
        modNaNs = np.argwhere(np.isnan(models))
        if len(modNaNs) > 0:
            observations = np.delete(observations, modNaNs)
            models = np.delete(models, modNaNs)  # removing c
    elif len(modNaNs) > 0:
        print 'warning found nans in bryant stats'
        models = np.delete(models, modNaNs)
        observations = np.delete(observations, modNaNs)
        obsNaNs = np.argwhere(np.isnan(observations))
        if len(obsNaNs) > 0:
            observations = np.delete(observations, obsNaNs)
            models = np.delete(models, obsNaNs)  # removing cor
    assert len(observations) == len(models), 'these data must be the same length'

    residuals = models - observations
    bias = np.nansum(residuals) / len(observations)

    ## RMSE's
    # demeaned RMSE
    RMSEdemeaned = np.sqrt(np.sum((residuals - bias) ** 2) / (len(observations) - 1))
    # regular RMSE
    RMSE = np.sqrt(np.sum(residuals ** 2) / len(observations))
    # normalized RMSE or percentage
    RMSEnorm = np.sqrt(np.sum(residuals ** 2) / np.sum(observations ** 2))
    # scatter index - a normalize measure of error often times presented as %
    ScatterIndex = RMSE / np.mean(observations)
    # symetric Slope
    symr = np.sqrt((models ** 2).sum() / (models ** 2).sum())
    # coefficient of determination
    r2 = np.sum((observations - observations.mean()) * (models - models.mean())) \
         / (np.sqrt(  ((observations - observations.mean()) ** 2).sum()) *
            np.sqrt(  ((models       - models.mean()      ) ** 2).sum()))
    # SSres = (residuals ** 2).sum()  ## from wiki
    # SStot = ((observations - observations.mean()) ** 2).sum()
    # r2 = 1 - SSres/SStot
    # wilmont 1985
    topW = np.abs(models - (observations).sum())
    botW = (np.abs(models - (observations).mean()) + np.abs(np.nansum(observations - (observations).mean())))
    Wilmont = 1 - topW / botW

    xRMS = np.sqrt((observations).sum() ** 2 / len(observations))
    pBias = 1 - np.abs(bias) / xRMS
    IMEDS = (xRMS + pBias) / 2
    stats = {'bias': bias,
             'RMSEdemeaned': RMSEdemeaned,
             'RMSE': RMSE,
             'RMSEnorm': RMSEnorm,
             'scatterIndex': ScatterIndex,
             'symSlope': symr,
             'corr': r2,
             'PscoreWilmont': Wilmont,
             'PscoreIMEDS': IMEDS,
             'residuals': residuals,
             'meta': 'please see Bryant, et al.(2016). Evaluation Statistics computed for the WIS ERDC/CHL CHETN-I-91'}

    return stats

def weightedAvg(toBeAveraged, weights, avgAxis=0, inverseWeight=False):
    """
    This function does a weighted average on a multidimensional array

    :param toBeAveraged: values to be averaged (array)
    :param weights: values to be used as weights (does not have to be normalized)
    :param avgAxis: axis over which to average
    :param inverseWeight: if true will invert the weights, high values weighted higher with marked FALSE
            to weight high values lower, inverseWeight must be True
    :return: an array of weighted average
    """
    if inverseWeight == True:
        weights = 1/weights
    assert toBeAveraged.shape == weights.shape, 'data and weights need to be the same shapes to be averaged'
    averagedData = np.sum(weights * toBeAveraged, axis=avgAxis) / (weights).sum(axis=avgAxis)
    return averagedData

def timeMatch(obs_time, obs_data, model_time, model_data):
    """
    This is the time match function from the IMEDs lite version created by ASA
    This has been removed from the IMEDS package to simplify use.
    This method returns the matching model data to the closest obs point.

    similar to time match imeds

    Time Matching is done by creating a threshold by taking the median of the difference of each time
       then taking the minimum of the difference between the two input times divided by 2.
       a small, arbitrary (as far as i know) factor is then subtracted from that minimum to remove the possiblity
       of matching a time that is exactly half of the sampling interval.

    :param obs_time: observation times, in
    :param obs_data: matching observation data, any shape
    :param model_time:  modeling time
    :param model_data:  modeling data (any shape)

    :return: time, (as float)
             obs_data_s -  data as input
             model_data_s  - data as input
    """

    # try to convert it from datetime to epochtime
    # this will fail if it already is in epochtime, so it wont do anything.
    dt_check = False
    try:
        timeunits = 'seconds since 1970-01-01 00:00:00'
        obs_time_n = nc.date2num(obs_time, timeunits)
        del obs_time
        obs_time = obs_time_n
        del obs_time_n
    except:
        pass
    try:
        timeunits = 'seconds since 1970-01-01 00:00:00'
        model_time_n = nc.date2num(model_time, timeunits)
        del model_time
        model_time = model_time_n
        del model_time_n
        dt_check = True
    except:
        pass

    assert type(obs_time[0]) != DT.datetime, 'time in must be numeric, try epoch!'
    assert type(model_time[0]) != DT.datetime, 'time in must be numeric, try epoch!'

    time = np.array([])
    obs_data_s = np.array([])
    model_data_s = np.array([])
    # 43 seconds here makes it 43 seconds less than 1/2 of smallest increment
    threshold = min(np.median(np.diff(obs_time)) / 2.0 - 43,
                    np.median(np.diff(model_time)) / 2.0 - 43)

    # Loop through model records
    for data, record in zip(model_data, model_time):
        in1 = np.where(obs_time <= record)[0]
        in2 = np.where(obs_time >= record)[0]

        if in1.size == 0 or in2.size == 0:
            continue

        if in1[-1] == in2[0]:  # Times match up
            indx = in2[0]
        else:
            d1 = record - obs_time[in1[-1]]
            d2 = obs_time[in2[0]] - record
            if min(d1, d2) > threshold:
                continue
            elif (d1 <= d2):
                indx = in1[-1]
            elif (d2 < d1):
                indx = in2[0]

        if (np.isnan(obs_data[indx]).all() or np.isnan(data).all()):
            continue

        time = np.append(time, record)
        obs_data_s = np.append(obs_data_s, obs_data[indx])
        model_data_s = np.append(model_data_s, data)

    # check to see if these are empty, and if so, interpolate
    if np.size(time) == 0:
        # this means it is empty!!!  interpolate obs onto model time
        obs_data_s = np.interp(model_time, obs_time, obs_data)
        model_data_s = model_data
        time = model_time
    else:
        pass

    # if I was handed a datetime, convert back to datetime
    if dt_check:
        timeunits = 'seconds since 1970-01-01 00:00:00'
        calendar = 'gregorian'
        time_n = nc.num2date(time, timeunits, calendar)
        del time
        time = time_n
        del time_n
    else:
        pass


    return time, obs_data_s, model_data_s

def waveStat(spec, dirbins, frqbins, lowFreq=0.05, highFreq=0.5):
    """     
    this function will calculate the mean direction from a full spectrum
    only calculates on one 2D spectrum at a time 
    defaults to 0.05 hz to 0.5 hz frequency for the statistics
    Input:
        %     spec  Frequency-direction spectra (2D)       shape(record,frqbin,dirbin)
        %  frqbins  Frequency vector (not assumed constant)
        %  dirbins  Direction vector (assumed constant)
        %
    Outputs (MKS, Hz, degrees, degrees CW from true north):
        %   Hmo   Significant wave height
        %    Tp   Period of the peak energy in the frequency spectra, (1/Fp).  AKA Tpd, not to be
        %           confused with parabolic fit to spectral period
        %    Tm02   Mean spectral period (Tm0,2, from moments 0 & 2), sqrt(m0/m2)
        %    Tm01   Average period, frequency sprectra weighted, from first moment (Tm0,1)
        %    Dp   Peak direction at the peak frequency
        %   Dmp   Mean direction at the peak frequency
        %    Dm   Mean wave direction  
        %  sprdF  Freq-spec spread (m0*m4 - m2^2)/(m0*m4)  (one definition)
        %  sprdD  Directional spread (m0*m4 - m2^2)/(m0*m4)  (one definition, Kuik 1988, buoys), 
        %         total sea-swell
        %         sprdD = r2d * sqrt(2.0 * (1.0 - sqrt(Xcomp.^2 + Ycomp^2)));
        %         where  Xcomp = sum(sin(Drad) .* Ds .* dwdir) ./ sum(Ds .* dwdir); 
        %                Ycomp = sum(cos(Drad) .* Ds .* dwdir) ./ sum(Ds .* dwdir);
        % sprdDhp  half-power direction width in direction spectra at peak freq (not currently incorporated)
        %  Tm10 - Mean Absolute wave Period from -1 moment

            return order [ Hm0, Tp, TmSecondMoment, Tm01,  Dp, Dm, Dmp, sprdF, sprdD, stats], Tm10
        Code Translated by Spicer Bak from: fd2BulkStats.m written by Kent Hathaway
            
    """
    warnings.warn('This function is depricated, The development should be moved to sb.waveLib version!!!!')
    assert type(frqbins) in [np.ndarray, np.ma.MaskedArray], 'the input frqeuency bins must be a numpy array'
    assert type(dirbins) in [np.ndarray, np.ma.MaskedArray], 'the input DIRECTION bins must be a numpy array'
    try:
        assert np.array(spec).ndim == 3, 'Spectra must be a 3 dimensional array'
    except AssertionError:
        spec = np.expand_dims(spec, axis=0)
    try:
        assert (spec != 0).all() is not True, 'Spectra must have energy to calculate statistics, all values are 0'
    except AssertionError:
        return 0
    # finding delta freqeucny (may change as in CDIP spectra)
    frq = np.array(np.zeros(len(frqbins) + 1))  # initializing frqbin bucket
    frq[0] = frqbins[0]
    frq[1:] = frqbins

    df = np.diff(frq, n=1)  # change in frequancy banding
    dd = np.abs(np.median(np.diff(dirbins)))  # dirbins[2] - dirbins[1]  # assume constant directional bin size
    # finding delta degrees
    # frequency spec
    fspec = np.sum(spec, axis=2) * dd  # fd spectra - sum across the frequcny bins to leave 1 x n-frqbins
    # doing moments over 0.05 to 0.33 Hz (3-20s waves) (mainly for m4 sake)
    [idx, vals] = findbtw(frqbins, lowFreq, highFreq, type=3)

    m0 = np.sum(fspec * df, axis=1)  # 0th momment
    m1 = np.sum(fspec[:, idx] * df[idx] * frqbins[idx], axis=1)  # 1st moment
    m2 = np.sum(fspec[:, idx] * df[idx] * frqbins[idx] ** 2, axis=1)  # 2nd moment
    m3 = np.sum(fspec[:, idx] * df[idx] * frqbins[idx] ** 3, axis=1)  # 3rd moment
    m4 = np.sum(fspec[:, idx] * df[idx] * frqbins[idx] ** 4, axis=1)  # 4th moment
    m11 = np.sum(fspec[:, idx] * df[idx] * frqbins[idx] ** -1, axis=1)  # negitive one moment

    # sigwave height
    Hm0 = 4 * np.sqrt(m0)
    # period stuff
    ipf = fspec.argmax(axis=1)  # indix of max frequency
    Tp = 1 / frqbins[ipf]  # peak period
    Tm02 = np.sqrt(m0 / m2)  # mean period
    Tm01 = m0 / m1  # average period - cmparible to TS Tm
    Tm10 = m11 / m0
    # directional stuff
    Ds = np.sum(spec * np.tile(df, (len(dirbins), 1)).T, axis=1)  # directional spectra (directional Spred)
    Dsp = []
    for ii in range(0, len(ipf)):
        Dsp.append(spec[ii, ipf[ii], :])  # direction spread at peak-f (ipf)
    Dsp = np.array(Dsp)
    idp = Dsp.argmax(axis=1)  # index of direction at peak frquency
    Dp = dirbins[idp]  # peak direction

    Drad = np.deg2rad(dirbins)  # making a radian degree bin
    # mean wave direction (e.g. Kuik 1988, used by USACE WIS)
    Xcomp = np.sum(np.cos(Drad) * Ds, axis=1)  # removed denominator as it canceles out in calculation
    Ycomp = np.sum(np.sin(Drad) * Ds, axis=1)  # removed denominator as it canceles out in calculation
    Dm = np.rad2deg(np.arctan2(Ycomp, Xcomp))
    Dm = angle_correct(Dm, rad=False)  # fixing directions above or below 360
    # Vector Dm (Hesser)
    sint = np.sin(Drad)  # sine of dirbins
    cost = np.cos(Drad)  # cosine of dirbins

    sint2 = np.tile(sint, [len(frqbins), 1])  # 2d diretion size of the spectra
    cost2 = np.tile(cost, [len(frqbins), 1])
    delsq = np.tile(df, [len(dirbins), 1]).T

    for tt in range(0, spec.shape[0]):
        xsum[tt] = sum(np.sum(cost2 * delsq * spec[tt, :, :], axis=1))  # summing along Xcomponant directions, then
        ysum[tt] = sum(np.sum(sint2 * delsq * spec[tt, :, :], axis=1))  # y componant

    vavgdir = np.rad2deg(np.arctan2(ysum, xsum))
    vavgdir = angle_correct(vavgdir)
    # assert vavgdir == Dm, 'Dm is calculated wrong ... at least once'
    # Mean direction at the peak frequency
    Dmp = np.rad2deg(np.arctan2(np.sum(np.sin(Drad) * Dsp, axis=1),
                                np.sum(np.cos(Drad) * Dsp, axis=1)))  # converting back to degrees
    Dmp = angle_correct(Dmp)
    # f-spec spread
    sprdF = (m0 * m4 - m2 ** 2) / (m0 * m4)

    # fd-spec spread
    sprdD = np.rad2deg(np.sqrt(2.0 * (1.0 - np.sqrt(Xcomp ** 2 + Ycomp ** 2))))

    ##### Exceprt from Kent's code for spreading - not sure how to handle
    # fd-spec spread, do a linear interp to get closer to half-power
    # % from the delta-deg increments
    # hp = np.max(Dsp)/2;

    ##### Exceprt from Kent's code for spreading - not sure how to handle
    #        % fd-spec spread, do a linear interp to get closer to half-power
    # % from the delta-deg increments
    # hp = np.max(Dsp)/2;
    # fd-spec spread, do a linear interp to get closer to half-power from the
    # delta-deg increments
    ##### Exceprt from Kent's code for spreading - not sure how to handle
    #        % fd-spec spread, do a linear interp to get closer to half-power
    # % from the delta-deg increments
    # hp = max(Dsp)/2;
    # ihp=find(Dsp > hp);
    #
    #  % Left (d1) and right (d2) interps: Y=Dir, X=E   
    # d1=interp1([Dsp(ihp(1)-1) Dsp(ihp(1)+1)], [dwdir(ihp(1)-1) dwdir(ihp(1)+1)], hp);
    # d2=interp1([Dsp(ihp(end)-1) Dsp(ihp(end)+1)], [dwdir(ihp(end)-1) dwdir(ihp(end)+1)], hp);
    # sprdDhp = d2 - d1;

    # wrapping up data into dictionary
    meta = 'Tp - peak period, Tm - mean period, Tave - average period, comparable to Time series mean period, Dp - peak direction, Dm - mean direction, Dmp - mean direction at peak frequency, vavgdir - Vector Averaged Mean Direction,sprdF - frequency spread, sprdD - directional spread'
    stats = {'Hm0': Hm0,
             'Tp': Tp,
             'Tm': Tm02,
             'Tave': Tm01,
             'Dp': Dp,
             'Dm': Dm,
             'Dmp': Dmp,
             'VecAvgMeanDir': vavgdir,
             'sprdF': sprdF,
             'sprdD': sprdD,
             'meta': meta
             }
    # print meta
    return Hm0, Tp, Tm02, Tm01, Dp, Dm, Dmp, vavgdir, sprdF, sprdD, stats, Tm10

def geo2STWangle(geo_angle_in, zeroAngle=70., METin=1, fixanglesout=0):
    """
    This rotates an angle (angle_in) from geographic Meterological convention 0= True North
    and puts it to an STWAVE angle 0 is onshore
    variable pierang is the angle from shore to 90 degrees (shore coordinates) in geographic convention
    ie the pier at Duck NC is at an angle of 70 degrees TN (stateplane) and this is assumed to be shore perpendicular

    :param geo_angle_in:  an array or list of angles to be rotated from MET convention of angle from
    :param zeroAngle:  the angle of the pier, from this the azimuth is calculated (MET CONVENTION)
    :param METin:  = 1 if the input angle is in MET convention (angle from)
    :param fixanglesout: if set to 1, will correct out angles to +/-180

    :return:
        angle_out corrected angle back out, into math space
    """
    # assert len(np.shape(geo_angle_in)) <= 1, 'function geo2STWangle not tested in more than 1 dimension'
    azimuth = 270 - zeroAngle  # the control of the zero for rotation of the grid in TN coords
    geo_angle_in = np.array(geo_angle_in, dtype=float)  # making sure floating calcs are used
    if METin == 1:  # flip the from/to convention
        ocean_angle_in = angle_correct(geo_angle_in + 180)  # to 'ocean' from 'MET' convention
    else:
        ocean_angle_in = geo_angle_in
    rotate = angle_correct(90 - azimuth)  # moving geo
    STWangle = angle_correct(rotate - ocean_angle_in)  # rotation of angles to grid convention
    #  putting into +/- 180 degrees
    if fixanglesout == 1:
        flip = np.argwhere(STWangle > 180)  # indicies that need to be flipped
        STWangle[flip] -= 360
    return STWangle

def STWangle2geo(STWangle, pierang=70, METout=1):
    """
    This is the complementary function to geo2STWangle,  It takes STWAVE angles (local coordinate system with a towards
     definition and + CCW)
    and translates them into geospatial grid angles (with a MET -from- convention and a CW+ convention)

    :rtype: 1D array of angles in geographic convention both met or ocean convention
    :param gridangle: an array or list of angles to be rotated
    :param pierang:  the (MET CONVENTION)
    :param METout:  if left 1, this creates output into a MET conveinton with the definition in the from
    :return:
        angle_out array of angles returned back to geographic convention (true north, clockwise positive)
    """
    assert len(np.shape(STWangle)) <= 3, 'STWangle2geo has not been tested in greater than 3dimensions'
    azimuth = 270 - pierang  # rotation of the Grid in local coordinate
    rotate = angle_correct(90 - azimuth)  # putting azimuth into ocean (towards) convention
    angle_out = rotate - STWangle
    if METout == 1:
        angle_out += 180
    angle_out = angle_correct(angle_out)  # correcting to < +360
    return angle_out

def whatIsYesterday(now=DT.date.today(), string=1, days=1):
    """
    this function finds what yesterday's date string is in the format
    of yyy-mm-dd
    :params:
        now:: the date to start counting backwards from

        string:: (1) output is in stiring format (default)
                 (2) output is in datetime format
        days:: how many days to count backwards from
           default = 1 
    :return
        the date of yesterday ( "days" previous to input 'now')
    """

    yesterday = now - DT.timedelta(days)
    if string == 1:
        yesterday = DT.date.strftime(yesterday, '%Y-%m-%d')
    return yesterday

def createDateList(start, end, delta):
    """
    creates a generator of dates
    :param start: date to start (date time object ... i think)
    :param end: date to end (datetime object .. i think)
    :param:  the amount to change between dates in the list ( time delta object)
    :return
        generator for datelists with dates separated by delta

    """
    curr = start
    while curr <= end:
        yield curr
        curr += delta

def importXYZ(fname_in):
    '''
    This function imports a file comma seperated and returns a dictionary with keys x, y, z
    the file to be imported must be x y z order

    :param fname_in: file name (incl path for xyz file)
    :return: dictionary of x, y, z values

    '''
    import csv
    # import file with given path above

    raw_x, raw_y = [], []
    ptCt, raw_z = [], []

    with open(fname_in, 'rb') as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) == 1:
                ptCt = row[0]
            elif len(row) == 3:
                raw_x.append(row[0])  # x in string format
                raw_y.append(row[1])  # y in string format
                raw_z.append(row[2])  # z in string format
            elif len(row) == 5:
                raw_z.append(row[2])  # z in string format
                raw_x.append(row[3])  # x in string format
                raw_y.append(row[4])  # y in string format

    # initializing values to make strictly numbers, imported as strings
    num_x = np.zeros(len(raw_x))
    num_y = np.zeros(len(raw_y))
    num_z = np.zeros(len(raw_z))
    # making strings into floating point numbers
    for ii in range(0, len(raw_x)):
        num_x[ii] = float(raw_x[ii])
        num_y[ii] = float(raw_y[ii])
        num_z[ii] = float(raw_z[ii])
    # gridding


    # format dictionary for output
    out = {'x': num_x,
           'y': num_y,
           'z': num_z,
           }
    return out

def write_frfGridFileText(FRFgrid, ofname, LatLon=False):
    """This function will write the cBathy grid (output from get data get_cBathyFromNc) to
    a text file grid of the same standard and format to that of the standard grid

    :param FRFgrid: this is the output from getdata.get_cBathyFromNc
        :key 'xm'; x coordinates in meters
        :key 'ym': y coordinates in meters
        :key 'depth': depths
    :param ofname: file output name
    :param LatLon: want to write the lat/lon too? boolean T/F (Default value = False)

    """
    import geoprocess as gp
    datestring = ofname[ofname.rfind('/') + 1:-4]
    version_prefix = ofname.split('/')[2]
    ym = FRFgrid['xFRF']
    xm = FRFgrid['yFRF']
    depth = FRFgrid['depth']

    if LatLon == True:
        LatGrid, LonGrid = np.meshgrid(xm, ym)
        grid = depth[0, :, :]  # depths here are negitive values

        f = open(ofname, 'w')
        for iy in range(0, np.size(grid, axis=0)):
            for ix in range(0, np.size(grid, axis=1)):
                pack = gp.FRFcoord(xm[ix], ym[iy])
                # print 'x coord:%d  Y coord: %d ' %(cBathy_grid['xm'][ix], cBathy_grid['ym'][iy])
                # print 'Lon: %f  Lat: %f' %(pack['Lon'], pack['Lat'])
                LatGrid[iy, ix] = pack['Lat']
                LonGrid[iy, ix] = -pack['Lon']

                # for ix in range(0, np.size(grid, axis=0)):
                #     for iy in range(0, np.size(grid, axis=1)):
                # these indicies were flipped for fitting purposes - the last one had backwards shape
                f.write("%f, %f, %f\n" % (
                -pack['Lon'], pack['Lat'], grid[iy, ix].T))  # (yy[ix, iy], xx[ix, iy], grid[ix, iy]))
        f.close()

    else:
        xcoord = xm
        ycoord = ym

        grid = depth.T  # made neg to fit into work flow
        # transpose of the mesh grids produced below
        yy, xx = np.meshgrid(xcoord, ycoord)  # putting x and y coords into grid
        f = open(ofname, 'w')
        for iy in range(0, np.size(grid, axis=0)):
            for ix in range(0, np.size(grid, axis=1)):
                # these indicies were flipped for fitting purposes - the last one had backwards shape
                f.write("%f, %f, %f\n" % (xx[ix, iy], yy[ix, iy], grid[iy, ix]))
        f.close()
    # chopping off the prefix of the grid name to conform with parent code
    outname = datestring + '.txt'
    return outname

def imedsObsModelTimeMatch(obs_time, obs_data, model_time, model_data, checkMask=True):
    """
    This method returns the matching model data to the closest obs point.
        it does this by using the smallest dt (in obs_time or model_time) and matching the other within that threshold
        then returns associated data


    :param obs_time: epoch (numeric) time
    :param obs_data: data associated with model_time
    :param model_time: epoch (numeric) time
    :param model_data: data associated with model_time
    :param checkMask:   this will not return masked data/time/ or NaN's with value marked True
    :return:
        time, time matched time
        obs_data_s, data associated with newly matched time
        model_data_s, data associated with newly matched time
    """
    assert np.array(obs_time).shape[0] == np.array(obs_data).shape[
        0], 'observational data must have same length as time '
    assert np.array(model_data).shape[0] == np.array(model_data).shape[0], ' model data must have same length as time'
    time = np.array([])
    obs_data_s = np.array([])
    model_data_s = np.array([])

    if len(model_time) == 1:
        threshold = (np.median(np.diff(obs_time)) / 2.0)
    else:
        threshold = min(np.median(np.diff(obs_time)) / 2.0 - 43,  # subtract 43 seconds to offset
                        np.median(np.diff(model_time)) / 2.0 - 43)

    # Loop through model records
    rc, rcc = 0, 0
    for data, record in zip(model_data, model_time):
        rcc += 1  # record counter for start of loop
        in1 = np.where(obs_time <= record)[0]
        in2 = np.where(obs_time >= record)[0]

        if in1.size == 0 or in2.size == 0:
            continue

        if in1[-1] == in2[0]:  # Times match up
            indx = in2[0]
        else:
            d1 = record - obs_time[in1[-1]]
            d2 = obs_time[in2[0]] - record
            if min(d1, d2) > threshold:
                continue
            elif (d1 <= d2):
                indx = in1[-1]
            elif (d2 < d1):
                indx = in2[0]
        # don't append data that is masked or NaN
        if checkMask == True:
            if (np.isnan(obs_data[indx]) or np.isnan(data)):
                continue
            elif type(obs_data[indx]) == np.ma.core.MaskedConstant and obs_data[indx].mask:
                continue
            elif type(data) == np.ma.core.MaskedConstant and data.mask:
                continue

        rc += 1  # record counter, how many loops
        time = np.append(time, record)
        obs_data_s = np.append(obs_data_s, obs_data[indx])
        model_data_s = np.append(model_data_s, data)

    return time, obs_data_s, model_data_s

def import_FRF_Transect(fname):
    """
    This function import a FRF transect csv file
    Comma Separated Value (CSV) ASCII data.  Column Header (not included in the file):
    Locality,Profile,SurveyNumber,Latitude,Longitude,Northing,Easting,FRF_offshore,
    FRF_longshore,Elevation,Ellipsoid,Date,Time,time_sec.

        Data Column headings Explanation:

        1) Locality Code (character) - "39" represents Duck, NC
        2) Profile Number (integer)
        3) Survey Number (integer) - sequential since initial survey
        4) Latitude (decimal degrees)
        5) Longitude (decimal degrees)
        6) Northing (decimal meters) - NAD83 North Carolina State Plane Coordinate
        7) Easting (decimal meters) - NAD83 North Carolina State Plane Coordinates
        8) FRF offshore coordinate: distance offshore from the local baseline
        9) FRF longshore coordinate: distance alongshore from the local origin
        10) Elevation - relative to the North American Vertical Datum in meters (NAVD88)
        11) Ellipsoid - optional, if used this is the geographic ellipsoid of the point
        12) Date - date data point was collected (YYYYMMDD)
        13) Time (hhmmss): Time data point was collected - in 24-hr Eastern Standard Time (EST)
        14) Time_sec: Time (EST) data point was collected in seconds past midnight

        Filename - the filename includes 7 fields with details about the contents of the file:

        1) Location (FRF)
        2) Survey Job (FRF)
        3) Survey Date (YYYYMMDD)
        4) Survey Number
        5) Vertical Datum (NAVD88)
        6) Vessel used (LARC)
        7) Survey Instrument (GPS)

    :param fname: name/location of a FRF measured bathymetry transect file
    :return: dictionary of all fields listed above
    """
    import pytz
    c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13 = [], [], [], [], [], [], [], [], [], [], [], [], [], []
    with open(fname, 'rb') as csvfile:  # opening the file
        reader = csv.reader(csvfile, delimiter=',')  # reading the lines
        for row in reader:
            c0.append(str(row[0]))  # Locality Code
            c1.append(int(row[1]))  # Profile number
            c2.append(int(row[2]))  # survey number
            c3.append(float(row[3]))  # Latitude
            c4.append(float(row[4]))  # Longitude
            c5.append(float(row[5]))  # Northing
            c6.append(float(row[6]))  # Easting
            c7.append(float(row[7]))  # FRF xshore
            c8.append(float(row[8]))  # FRF Long shore
            c9.append(float(row[9]))  # Elevation (NAVD88)
            c10.append(float(row[10]))  # Ellipsoid
            c11.append(row[11])  # YYMMDD
            c12.append(row[12])  # hhmmss
            #    c13.append(row[13])   seconds past midnight
    # convert EST to UTC

    time = []
    for ii in range(0, len(c12)):
        EST = DT.datetime(int(c11[ii][0:4]), int(c11[ii][4:6]), int(c11[ii][6:]),
                          int(c12[ii][:2]), int(c12[ii][2:4]), int(c12[ii][4:]),
                          tzinfo=pytz.timezone('EST'))
        time.append(
            EST.astimezone(pytz.timezone('UTC')).replace(tzinfo=None))  # converting to UTC, and removing UTC metadata

    bathyDict = {'Locality_Code': np.array(c0),
                 'Profile_number': np.array(c1),
                 'Survey_number': np.array(c2),
                 'Latitude': np.array(c3),
                 'Longitude': np.array(c4),
                 'Northing': np.array(c5),
                 'Easting': np.array(c6),
                 'FRF_Xshore': np.array(c7),
                 'FRF_Yshore': np.array(c8),
                 'Elevation': np.array(c9),
                 'Ellipsoid': np.array(c10),
                 'time': np.array(time),  # datetime object
                 'meta': 'date and Time has been converted to a UTC datetimeta object, elevation is in NAVD88',
                 }
    return bathyDict

def vectorRotation(vector, theta=90, axis='z'):
    """
    This function does a vector rotation of the vector input in vector, rotated by theta
    ASSUMES CCW

    :param vector: 2d or 3d vector you want rotated... [x, y, z]
    :param axis: axis you want it rotated about 'x' = [1, 0, 0], 'y' = [0, 1, 0], 'z' = [0, 0, 1]
    :param theta: angle in decimal degrees
    
    :return: vector rotated CCW theta degrees about axis, uses Euler-Rodrigues formula  
    """

    vector = np.asarray(vector)
    assert -360 <= theta <= 360, 'your angle must be a decimal degree value -360  and 360 degrees'
    assert len(vector) >= 2, "You must hand this function a 2D or 3D vector!"
    assert len(vector) <= 3, "You must hand this function a 2D or 3D vector!"

    if len(vector) == 2:
        vector = np.append(vector, 0) # this just converts it to a 2D vector
        ret = '2d'
    else:
        ret = '3d'

    if type(axis) == str:
        assert axis in ['x', 'y', 'z'], 'Acceptable axis inputs are x, y, z, or a 3D vector'

        if axis == 'x':
            axis = [1, 0, 0]
        elif axis == 'y':
            axis = [0, 1, 0]
        elif axis == 'z':
            axis = [0, 0, 1]
        else:
            pass
        axis = np.asarray(axis)
    else:
        axis = np.asarray(axis)
        assert len(axis) == 3, 'Acceptable axis inputs are x, y, z, or a 3D vector'

    theta = 2*math.pi*(theta/360.0)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    r_mat = np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

    r_vector = np.dot(r_mat, vector)

    if ret == '3d':
        return r_vector
    else:
        return r_vector[0:2]


