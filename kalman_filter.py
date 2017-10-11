import numpy as np
import datetime as dt

def extract_time(data,index):
    """
    This function takes a dictionary [data] and pulles out all of the keys at specific index [index]

    :param data: dictionary
    :param index: index to be removed
    :return: new dictionary with only the indexs selected returned
    """
    vars = data.keys()
    new = {}
    for vv in vars:
        if vv is 'xm' or vv is 'ym':
            new[vv] = data[vv]
        else:
            new[vv] = data[vv][index]
    return new

def kalman_filter(new, prior, waveHs):
    """
    This function does a kalman filter designed for implmeneting wave height thresholds into the cbathy
    algorithm, this operates on a single time step only!!!


    :param new: a dictionary with keys associated with get data
    :param prior:  a saved dictionary with bathys derived from times when wave heights were below the threshold of choice

    :param waveHs:
    :return:
    """


    n = 2.0
    Cq = 0.067
    sigmax = 100.
    x0 = 150.

    xlen = len(new['xm'])
    ylen = len(new['ym'])
    temp = new['time'] - prior['time']
    delt = temp.days + temp.seconds/(24.*3600.)
    ## maybe we flatten the arrays so we don't have to loop
    # need to make a new x array
    xarray = np.tile(new['xm'],(1,len(new['ym']))).T
   # for ix in range(len(new['xm'])):
    Q = Cq*waveHs**n*np.exp(-((xarray - x0)/sigmax)**2)

   #     for iy in range(len(new['ym'])):
    Pkm = prior['P'].flatten() + Q[:,0]*delt

    K = Pkm/(Pkm + new['depthfCError'].flatten()**2)
    
    hk = prior['depthKF'].flatten() + K*(new['depthfC'].flatten() - prior['depthKF'].flatten())

    Pk = (1-K)*Pkm

    #hnew = np.empty((hk.shape))
    #hnew[:] = np.NAN
    #pnew = np.empty((hk.shape))
    #pnew[:] = np.NAN
    
    #iin = (np.isnan(new['depthfC'].flatten())) #or (np.isnan(new['depthfCError'].flatten()))
    #hk[iin] = prior['depthKF'].flatten()[iin]
    #Pk[iin] = Pkm[iin]

    #iip = np.isnan(prior['depthKF'].flatten())
    #hk[iin] = new['depthfC'].flatten()[iip]
    #Pk[iin] = new['depthfCError'].flatten()[iip]**2.
    
# tests show below is the same
    # new['depthKF'][:] = 0.0
    # new['depthKFError'][:] = 0.0
    # new['P'][:] = 0.0
#    new['Q'][:] = 0.0
#     new['depthKF'] = hk.reshape((ylen,xlen))
# #    new['depthKFError'] = np.sqrt(Pk.reshape((ylen,xlen)))
#     new['P'] = Pk.reshape((ylen,xlen))
# #    new['Q'] = Q.reshape((ylen,xlen))
#     ## fill the new file with the old values when they have no values
#     idd = np.ma.getmask(new['depthKF'])
#     new['depthKF'][idd] = prior['depthKF'][idd]
#     new['P'][idd] = prior['P'][idd]
#     new['depthKFError'] = np.sqrt(new['P'])

    idx = np.argwhere(hk.mask).squeeze()       # find idx of missing points in new array
    hk[idx] = prior['depthKF'].flatten()[idx]  # fill kalman filtered depth estimates with old values when missing
    Pk[idx] = prior['P'].flatten()[idx]        # fill error variance with old values when missing
    # package for departure
    new['depthKF'] = hk.reshape((ylen,xlen))
    new['P'] = Pk.reshape((ylen,xlen))
    new['depthKFError'] = np.sqrt(new['P'])


    return new
