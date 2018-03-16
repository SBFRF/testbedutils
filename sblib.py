# -*- coding: utf-8 -*-
"""
Author: spicer bak
Contact: spicer.bak@usace.army.mil
Association: USACE CHL Field Research Facility

A base level utility function library that don't fit into specific other places

documented 12/2/17
"""
import numpy as np
import datetime as DT
import warnings
import netCDF4 as nc
import anglesLib

class Bunch(object):
    """allows user to access dictionary data from 'object'
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

def whatIsYesterday(now=DT.date.today(), string=1, days=1):
    """this function finds what yesterday's date string is in the format
    of yyyy-mm-dd

    :param s: now:: the date to start counting backwards from
    
        string:: (1) output is in stiring format (default)
                 (2) output is in datetime format
        days:: how many days to count backwards from
           default = 1
    :return
        the date of yesterday ( "days" previous to input 'now')
    :param now:  (Default value = DT.date.today()

    """

    yesterday = now - DT.timedelta(days)
    if string == 1:
        yesterday = DT.date.strftime(yesterday, '%Y-%m-%d')
    return yesterday

def baseRound(x, base=5):
    """This function will round any value to a multiple of the base,

    :param x: values to be rounded:
    :param base: this is the value by which x is rounded to a multiple of
        ie base = 10  x = [4, 8, 2, 12]  returns [0,10,0,10] (Default value = 5)
    :returns: np.array of floating point numbers rounded to a multiple of base

    """
    x = np.array(x, dtype=float)
    return base * np.round(x/base)

def roundtime(timeIn=None, roundTo=60):
    """"
    Round a datetime object to any time lapse in seconds

    :param dt: datetime.datetime object, default now.
    :param roundTo: Closest number of seconds to round to, default 1 minute.
    
       Author: Thierry Husson 2012 - Use it as you want but don't blame me.
    
       modified by SB to include lists of datetime dataList,
       returned as a list if it came in as a list, if it came in as a datetime object
       it is returned as such
    :return
        list/array of times rounded to 'round to value'
    :param timeIn:  (Default value = None)

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

def createDateList(start, end, delta):
    """creates a generator of dates

    :param start: date to start (date time object ... I think)
    :param end: date to end (datetime object .. I think)
    :param delta: 

    """
    curr = start
    while curr <= end:
        yield curr
        curr += delta

def statsBryant(observations, models):
    """This function does Non-Directional Statsistics
    These statistics are from the Bryant Wave stats CHETN - I - 91

    :param observations: array of observational data
    :param models: array of model data
    :returns: dictiornay
        :key  'bias' average of residuals
        :key  'RMSEdemeaned': RMSEdemeaned,
        :key  'RMSE': R oot mean square error
        :key  'RMSEnorm': normalized root Mean square error also Percent RMSE
        :key  'scatterIndex': ScatterIndex
        :key  'symSlope': symetrical slope
        :key  'corr': R^2 --- coefficient of determination
        :key  'PscoreWilmont': performance score developed by Wilmont
        :key  'PscoreIMEDS':  see Hanson 2007 for description
        :key  'residuals': model - observations
        :

    """
    obsNaNs = np.argwhere(np.isnan(observations)).squeeze()
    modNaNs = np.argwhere(np.isnan(models)).squeeze()
    if type(observations) == np.ma.masked_array or type(models) == np.ma.masked_array:
        raise NotImplementedError('this handles masked arrays poorly, fix or remove before use')
    if len(obsNaNs) > 0:
        warnings.warn('warning found nans in bryant stats')
        observations = np.delete(observations, obsNaNs)
        models = np.delete(models, obsNaNs)  # removing corresponding model data, that cannot be compared
        modNaNs = np.argwhere(np.isnan(models))
        if len(modNaNs) > 0:
            observations = np.delete(observations, modNaNs)
            models = np.delete(models, modNaNs)  # removing c
    elif len(modNaNs) > 0:
        warnings.warn('warning found nans in bryant stats')
        models = np.delete(models, modNaNs)
        observations = np.delete(observations, modNaNs,0)
        obsNaNs = np.argwhere(np.isnan(observations))
        if len(obsNaNs) > 0:
            observations = np.delete(observations, obsNaNs)
            models = np.delete(models, obsNaNs)  # removing cor
    assert len(observations) == len(models), 'these data must be the same length'

    residuals = models - observations
    bias = np.nansum(residuals) / len(residuals)

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
    symr = np.sqrt((models ** 2).sum() / (observations ** 2).sum())
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

def makegif(flist, ofname, size=None, dt=0.5):
    """This function uses imageio to create gifs from a list of images
    
    kwargs for mimwrite http://imageio.readthedocs.org/en/latest/format_gif.html#gif

    :param flist: a sorted list of files to be made into gifs (including path)
    :param ofname: output gif filename (including path)
    :param size: size of pictures (default not resized)
    :param loop: number of loops to do, 0 is default and infinite
    :param dt:  (Default value = 0.5)
    :returns: will write a gif to ofname location

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

def running_mean(data, window):
    """found running mean function on the internet, untested

    :param data: data to run mean
    :param window: window over which to take mean
    :returns: meaned data

    """
    cumsum = np.cumsum(np.insert(data, 0, 0))
    return (cumsum[window:] - cumsum[:-window]) / window

def find_nearest(array, value):
    """Function looks for value in array and returns the closest array value
    (to 'value') and index of that value

    :param array: array to to find things in
    :param value: value to search against
    :return returns a number from the array value, that is closest to value input

    """
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx

def findUniqueFromTuple(a, axis=1):
    """This function finds the unique values of a multi dimensional tuple (quickly)

    :param a: an array of multidimensional size
    :param axis: this is the axis that it looks for unique values using (default is horizontal)
    :returns: array of unique tuples
        warning, values are not sorted

    """
    if axis == 0:
        a = a.T

    b = np.ascontiguousarray(a).view(np.dtype((np.void, a.dtype.itemsize * a.shape[1])))
    _, idx = np.unique(b, return_index=True)

    unique_a = a[idx]
    return unique_a

def weightedAvg(toBeAveraged, weights, avgAxis=0, inverseWeight=False):
    """This function does a weighted average on a multidimensional array

    :param toBeAveraged: values to be averaged (array)
    :param weights: values to be used as weights (does not have to be normalized)
    :param avgAxis: axis over which to average (Default value = 0)
    :param inverseWeight: if true will invert the weights, high values weighted higher with marked FALSE
            to weight high values lower, inverseWeight must be True (Default value = False)
    :returns: an array of weighted average

    """
    if inverseWeight == True:
        weights = 1/weights
    assert toBeAveraged.shape == weights.shape, 'data and weights need to be the same shapes to be averaged'
    averagedData = np.sum(weights * toBeAveraged, axis=avgAxis) / (weights).sum(axis=avgAxis)
    return averagedData

def findbtw(data, lwth, upth, type=0):
    """This function finds both values and indicies of a list values between two values
    :TODO: probably could be improved by using boolean compares

    :param upth: upper level threshold
    :param lwth: lower level threshold
    :param data: list (or numpy array?)
    :param type: 0 = non inclusive  ie. lwth < list <  upth
        1 = low incluisve  ie. lwth <=list <  upth
        2 = high inclusive ie. lwth < list <= upth
        3 = all inclusive  ie  lwth <=list <= upth
    :return
        indicies returns idices that meet established criteria
        values   returns associated values from da (Default value = 0)

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

def timeMatch(obs_time, obs_data, model_time, model_data):
    """This is the time match function from the IMEDs lite version created by ASA
    This has been removed from the IMEDS package to simplify use.
    This method returns the matching model data to the closest obs point.
    
    similar to time match imeds
    
    Time Matching is done by creating a threshold by taking the median of the difference of each time
       then taking the minimum of the difference between the two input times divided by 2.
       a small, arbitrary (as far as I know) factor is then subtracted from that minimum to remove the possiblity
       of matching a time that is exactly half of the sampling interval.

    :param obs_time: observation times, in
    :param obs_data: matching observation data, any shape
    :param model_time: modeling time
    :param model_data: modeling data (any shape)
    :returns: time, (as float)
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

def timeMatch_altimeter(altTime, altData, modTime, modData, window=30 * 60):
    """this function will loop though variable modTim and find the closest value
    then look to see if it's within a window (default, 30 minutes),
     return altimeter data, and matched time
    
     Note: this might be slower than imeds time match or time match,
           which is based on the imeds time match

    :param altTime: altimeter time - tested as epoch (might work in datetime)
    :param altData: altimeter data, some/any floating (int?) value
    :param modTime: base time to match
    :param modData: data to be paired (could be indicies)
    :param window: time in seconds (or time delta, if input as datetimes) (Default value = 30 * 60)
    :returns: timeout: matched time in same format put in
        dataoutt: time matched altimeter data from variable altData
        modout: time matched model data

    """

    # try to convert it from datetime to epochtime
    # this will fail if it already is in epochtime, so it wont do anything.
    dt_check = False
    try:
        timeunits = 'seconds since 1970-01-01 00:00:00'
        altTime_n = nc.date2num(altTime, timeunits)
        del altTime
        altTime = altTime_n
        del altTime_n
    except:
        pass
    try:
        timeunits = 'seconds since 1970-01-01 00:00:00'
        modTime_n = nc.date2num(modTime, timeunits)
        del modTime
        modTime = modTime_n
        del modTime_n
        dt_check = True
    except:
        pass

    assert type(altTime[0]) != DT.datetime, 'time in must be numeric, try epoch!'
    assert type(modTime[0]) != DT.datetime, 'time in must be numeric, try epoch!'

    dataout, timeout, modout = [], [], []
    if type(altData) == np.ma.MaskedArray:
        altTime = altTime[~altData.mask]
        altData = altData[~altData.mask]
    assert len(altData) == len(altTime), 'Altimeter time and data length must be the same'
    for tt, time in enumerate(modTime):
        idx = np.argmin(np.abs(altTime - time))
        if altTime[idx] - time < window:
            # now append data
            dataout.append(altData[idx])
            timeout.append(modTime[tt])
            modout.append(modData[tt])

    if dt_check:
        timeunits = 'seconds since 1970-01-01 00:00:00'
        tOut = np.array(timeout)
        tOutN = nc.num2date(tOut, timeunits)
        del tOut
        tOut = tOutN
        del tOutN

    return tOut, np.array(dataout), np.array(modout)

def reduceDict(dictIn, idxToKeep, exemptList=None):
    """This function will take a dictionary and  reduce it by given indicies, making a COPY of the array
    it assumes that all things of multidimensions have time in the first dimension
    
    WARNING: This could potentially be dangerous as it works by checking each key to see if its the lenth of the
        variable dictIn['time'].  If it is then it will reduce the variable to keep the idxToKeep.  If not it will
        skip the key
    
    This function is useful if trying to reduce a dictionary to specific indicies of interest eg with a time matched index
    Improvement to be safer is welcome, or use with caution.

    :param dictIn: dictionary with no limit to how many keys assumes one key is 'time'
    :param idxToKeep: indices to reduce to, must be shorter than key 'time'
    :param exemptList: this is a list of variables to exclude
            default values are 'time', 'name', 'xFRF', 'yFRF'
    :returns: a dictionary with the same keys that were input
            all keys that came in with the same length as 'time' have been reduced to use the indices of idxToKeep

    """
    assert 'time' in dictIn, 'This function must have a variable "time"'
    if exemptList == None:
        exemptList = ['time', 'name', 'xFRF', 'yFRF']
    idxToKeep = np.array(idxToKeep, dtype=int) # force data to integer type
    dictOut = dictIn.copy()
    for key in dictIn:
        # if things are longer than the indicies of interest and not 'time'
        # print 'key %s size %d' %(key, np.size(dictIn[key], axis=0))
        try:
            if key not in exemptList and dictIn[key].dtype.kind not in ['U','S'] and np.size(dictIn[key], axis=0) == len(dictIn['time']):
                # then reduce
                dictOut[key] = dictIn[key][idxToKeep]  # reduce variable
                # print 'key %s Made it past test and new size %d' %(key, len(dictIn[key]))
        except IndexError:
            pass  # this passes for single number (not arrays)
    dictOut['time'] = dictIn['time'][idxToKeep]  # # once the rest are done finally reduce 'time'
    return dictOut


def waveStat(spec, dirbins, frqbins, lowFreq=0.05, highFreq=0.5):
    """this function will calculate the mean direction from a full spectrum
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

    :param spec: 
    :param dirbins: 
    :param frqbins: 
    :param lowFreq:  (Default value = 0.05)
    :param highFreq:  (Default value = 0.5)
    :returns: Code Translated by Spicer Bak from: fd2BulkStats.m written by Kent Hathaway

    """
    raise NotImplementedError('This function is depricated, The development should be moved to sb.waveLib version!!!!')
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
    Dm = anglesLib.angle_correct(Dm, rad=False)  # fixing directions above or below 360
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
    vavgdir = anglesLib.angle_correct(vavgdir)
    # assert vavgdir == Dm, 'Dm is calculated wrong ... at least once'
    # Mean direction at the peak frequency
    Dmp = np.rad2deg(np.arctan2(np.sum(np.sin(Drad) * Dsp, axis=1),
                                np.sum(np.cos(Drad) * Dsp, axis=1)))  # converting back to degrees
    Dmp = anglesLib.angle_correct(Dmp)
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

def FRFcoord(p1, p2):
    """place older

    :param p1: 
    :param p2: 

    """
    raise ImportError('please use the function in sblib.geoprocess')
