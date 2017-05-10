import numpy as np
from matplotlib import pyplot as plt
import sblib as sb
def ADV(p, u, v, t, waterDepth, gaugeDepth):
    """
    Takes below parameters and makes them into magical wave data

    :param p:  pressure in cm
    :param u: u velocities
    :param v: v velocities
    :param t: time stamp
    :param waterDepth:  positive in meters
    :param gaugeDepth:  negative in meters
    :return:
    """
    from scipy.signal import welch, csd

    Fsamp = 1 / np.median(np.diff(t))  # sample frequency (not period)
    nseg = 512  # 512
    overlap = nseg / 4
    bave = 5  # number of spectral bands to average (taken from kent)

    # Calculate spectra in each componant
    # x is north, y is east in kents code (used for  a's b's)
    pDemeaned = (p - np.mean(p)) / 100  # demean pressure and convert to meters
    [f1, Sp] = welch(x=pDemeaned, window='hanning', fs=Fsamp, nperseg=nseg, noverlap=overlap, nfft=None,
                     return_onesided=True, detrend='linear')
    [f1, Su] = welch(x=u, window='hanning', fs=Fsamp, nperseg=nseg, noverlap=overlap, nfft=None, return_onesided=True,
                     detrend='linear')
    [f1, Sv] = welch(x=v, window='hanning', fs=Fsamp, nperseg=nseg, noverlap=overlap, nfft=None, return_onesided=True,
                     detrend='linear')
    # create cross spectral power across the 3
    [f1, CrossPU] = csd(y=pDemeaned, x=u, fs=Fsamp, nperseg=nseg, noverlap=overlap, return_onesided=True,
                        window='hann')
    [f1, CrossPV] = csd(y=pDemeaned, x=v, fs=Fsamp, nperseg=nseg, noverlap=overlap, return_onesided=True,
                        window='hann')
    [f1, CrossUV] = csd(y=v, x=u, fs=Fsamp, nperseg=nseg, noverlap=overlap, return_onesided=True,
                        window='hann')
    ## now Do some frequency band Averaging
    # first remove first index of array  "remove DC compoonents"
    f1 = f1[1:]
    Sp = np.real(Sp[1:])
    Su = np.real(Su[1:])
    Sv = np.real(Sv[1:])
    CrossPU = CrossPU[1:]
    CrossPV = CrossPV[1:]
    CrossUV = CrossUV[1:]
    # correcting Pressure spectra to surface spectra
    L, _, _ = dispersion(waterDepth, T=1 / f1)
    prF = prFunc(L, waterDepth, gaugeDepth)  #
    Sp = Sp / prF ** 2  # getting surface corrected energy spectrum

    #  Beginning to Setup For Band Averaging
    dk = int(np.floor(bave / 2))

    ki = 0
    freq, SpAve, SuAve, SvAve, CrossPUave, CrossPVave, CrossUVave = [], [], [], [], [], [], []
    for kk in range(dk + 1, len(Sp) - dk, bave):
        avgIdxs = np.linspace(kk - dk, kk + dk, num=int((kk + dk) - (kk - dk) + 1), endpoint=True, dtype=int)
        freq = np.append(freq, f1[int(kk)])  # taking the first frequency for the bin label
        # bin averaging frequency spectra across bave frequency bins
        SpAve.append(np.sum(Sp[avgIdxs]) / bave)
        SuAve.append(np.sum(Su[avgIdxs]) / bave)
        SvAve.append(np.sum(Sv[avgIdxs]) / bave)
        # bin averaging cross correlations
        CrossPUave.append(np.sum(CrossPU[avgIdxs]) / bave)  # Press - FRF X
        CrossPVave.append(np.sum(CrossPV[avgIdxs]) / bave)  # Press - FRF Y
        CrossUVave.append(np.sum(CrossUV[avgIdxs]) / bave)  # FRF X - FRF Y
    # convert back to Numpy array
    freq = np.array(freq)
    SpAve = np.array(SpAve)
    SuAve = np.array(SuAve)
    SvAve = np.array(SvAve)
    CrossPUave = np.array(CrossPUave)
    CrossPVave = np.array(CrossPVave)
    CrossUVave = np.array(CrossUVave)

    # Extracting as and bs
    a1 = np.imag(CrossPUave) / np.sqrt(SpAve * (SuAve + SvAve))
    b1 = -np.imag(CrossPVave) / np.sqrt(SpAve * (SuAve + SvAve))
    a2 = (SuAve - SvAve) / (SvAve + SuAve)
    b2 = -2 * np.real(CrossUVave) / (SvAve + SuAve)

    wfreq = np.arange(0.04, 0.5, 0.0075)  # frequencies to interp to
    wDeg = np.arange(0, 360, 5)
    # interpolating to reasonable frequency banding
    import scipy.interpolate as interp
    f = interp.interp1d(freq, a1, kind='linear')
    a1interp = f(wfreq)
    f = interp.interp1d(freq, a2, kind='linear')
    a2interp = f(wfreq)
    f = interp.interp1d(freq, b1, kind='linear')
    b1interp = f(wfreq)
    f = interp.interp1d(freq, b2, kind='linear')
    b2interp = f(wfreq)
    f = interp.interp1d(freq, SpAve, kind='linear')
    fspec = f(wfreq)

    # mlmSpec = mlm(freqs=freq, dirs=np.arange(0,360,5), c11=SpAve, c22=SuAve, c33=SvAve, c23=np.real(CrossUVave), q12=np.imag(CrossPUave), q13=np.imag(CrossPVave))

    return fspec, a1interp, b1interp, a2interp, b2interp


def qkhfs(w, h):
    """
    Quick iterative calculation of kh in gravity-wave dispersion relationship
    kh = qkhfs(w, h )

    Input
        w - angular wave frequency = 2*pi/T where T = wave period [1/s]
        h - water depth [m]
    Returns
        kh - wavenumber * depth [ ]
    Orbital velocities from kh are accurate to 3e-12 !
    RL Soulsby (2006) \"Simplified calculation of wave orbital velocities\"
    HR Wallingford Report TR 155, February 2006
    Eqns. 12a - 14
    """
    g = 9.81
    x = w ** 2.0 * h / g
    y = np.sqrt(x) * (x < 1.) + x * (x >= 1.)
    # is this appalling code faster than a loop in Python? It is in Matlab.
    t = np.tanh(y)
    y = y - ((y * t - x) / (t + y * (1.0 - t ** 2.0)))
    t = np.tanh(y)
    y = y - ((y * t - x) / (t + y * (1.0 - t ** 2.0)))
    t = np.tanh(y)
    y = y - ((y * t - x) / (t + y * (1.0 - t ** 2.0)))
    kh = y
    return kh


def makeMLMspecFromAsBs(a0, a1, b1, a2, b2, waterDepth, freqs, dirs):
    """

    :param a1:
    :param b1:
    :param a2:
    :param b2:
    :return:
    """
    waveFreqBins = freqs
    waveDirBins = dirs

    L, _, _ = dispersion(waterDepth, 1 / waveFreqBins)  # calculating wave number
    xk = 2 * np.pi / L  # calculating wave number over depth
    ak = 1 / np.tanh(xk * 250)
    c11 = a0
    c22 = 0.5 * (1 + a2) * (ak ** 2) * a0  # %Cnn
    c33 = 0.5 * (1 - a2) * (ak ** 2) * a0  # % Cww
    c23 = 0.5 * b2 * (ak ** 2) * a0  # % Cnw
    q12 = a1 * ak * a0  # % Qvn
    q13 = b1 * ak * a0  # % Qvw
    spec = mlm(waveFreqBins, waveDirBins, c11, c22, c33, c23, q12, q13)

    return spec


def mlm(freqs, dirs, c11, c22, c33, c23, q12, q13):
    """
        From kent hathaway's code translated from fortran, chuck long -> hanson
        mlm.m
    % c  --------------------------------------------------------------------
    % c     Maxumum Likelihood Method (MLM) directional wave analysis
    % c     Reference:
    % c     ---------
    % c     J. Oltman-Shay & R.T. Guza, A Data-Adaptive Ocean Wave Directional
    % c     Spectrum Estimator for Pitch and Roll Type Measurements, Journal
    % c     of Physical Oceanography, Vol. 14, pp. 1800-1810, 1984.
    % c  --------------------------------------------------------------------

    :param freq:  frequencies
    :param dirs:  direction bins
    :param c11:  cross spectrum of 1st channel - pressure (or welch specturm)
    :param c22:  cross spectrum of 2nd channel - u        (or welch specturm)
    :param c33:  cross spectrum of 3rd channel - v        (or welch specturm)
    :param c23:  real part of cross spectrum of 2nd and 3rd channel u v
    :param q12:  imaginary part of 1st and 2nd channel (p u)
    :param q13:  imaginary part of 1st and 3rd channel (p v)
    :return:  2D directional spectrum using MLM method (non-iterative)
            IN  per RADIAN UNITS  -> must be converted back to degrees
    """
    # calculating sin/cos terms from direction angles
    angRad = np.deg2rad(dirs);  # directions in radians (for sin/cos calcs)
    ainc = np.median(np.abs(np.diff(dirs)))  # delta degrees of directions
    sint = np.sin(angRad);
    cost = np.cos(angRad);
    sin2t = np.sin(2. * angRad);
    cos2t = np.cos(2. * angRad);

    dsprdMLM = np.ones((len(freqs), len(dirs)))  # pre allocated space for spectra out
    eftMLM = np.ones((len(freqs), len(dirs)))
    for n1 in range(0, len(freqs)):
        # find spectral values below energy threshold
        if c11[n1] < 1e-5:
            dsprdMLM[n1, :] = 0.0
        else:
            # 1.1 wave number from spectra
            xk = np.sqrt((c22[n1] + c33[n1]) / c11[n1])
            d = c11[n1] * c22[n1] * c33[n1] - c11[n1] * c23[n1] ** 2 - q12[n1] ** 2 * c33[n1] \
                - q13[n1] ** 2 * c22[n1] + 2 * q12[n1] * q13[n1] * c23[n1]
            tmp = 2 * (c22[n1] * c33[n1] - c23[n1] ** 2) / d

            a0 = tmp + (c11[n1] * c22[n1] + c11[n1] * c33[n1] - q12[n1] ** 2 - q13[n1] ** 2) * xk ** 2 / d
            a1 = -2 * (q12[n1] * c33[n1] - q13[n1] * c23[n1]) * xk / d
            b1 = 2 * (q12[n1] * c23[n1] - q13[n1] * c22[n1]) * xk / d
            a2 = 0.5 * (c11[n1] * c33[n1] - c11[n1] * c22[n1] - q13[n1] ** 2 + q12[n1] ** 2) * xk ** 2 / d
            b2 = -(c11[n1] * c23[n1] - q12[n1] * q13[n1]) * xk ** 2 / d
            # 1.2 Directional Spread
            denom = 0.5 * a0 + a1 * cost + b1 * sint + a2 * cos2t + b2 * sin2t;
            # mn = min(denom)
            # if mn <= 6:
            #    denom = denom - mn + 6
            # dd2[n1] = denom
            # dd[n1] = d
            dsprdMLM[n1, :] = np.abs(1 / denom)
            asum = np.sum(dsprdMLM[n1, :])
            #
            # 2.0 Normalize directional spreading function to unity
            #    NOTE: put back the 1/rad units!!!
            asum = asum * ainc
            dsprdMLM[n1, :] = dsprdMLM[n1, :] / asum
    ## 3.0 Generate The 2D spectrum
    for fcount in range(0, len(freqs)):
        eftMLM[fcount, :] = c11[fcount] * dsprdMLM[fcount, :]

    return eftMLM


def prFunc( L, d, z):
    """
    Surface correction for pressure data.  if working with an energy
    spectrum DIVIDE by the returned prf**2

    :param L: wave length  (column)  (m)
    :param d: water depth     ( positve scalar)  (m)
    :param z: gage location below surface (negative) (m)
    :return: Pressure Response (column)
    """
    maxCorrectionFactor = 10
    # checks
    assert len(np.unique(
        z)) == 1, 'pressure response function cannot handle varying depth right now, shouldn''t be hard to modify'
    assert d > 0, 'water depth must be positive (meters)'
    if np.size(z) == 1:
        z = np.array(z)
    elif len(z) > 1:
        z = np.array(z)[0]
    assert z < 0, 'gauge location must be negative (meters)'
    # ## # now calcs
    arg1 = np.cosh(2 * np.pi * (d + z) / L);
    arg2 = np.cosh(2 * np.pi * (d) / L);
    prf = arg1 / arg2
    # limit the amount
    if (prf <= 1 / maxCorrectionFactor).any():
        print 'correcting extreme attenuation'
        prf[np.argwhere(prf <= 1 / maxCorrectionFactor)] = 1

    return prf


def dispersion( h, T):
    """
    Linear Dispersion Relationship
    omega^2 = gk tanh kh
    approximation and iterative method taken from
    ib svendsen "introduction to nearshore hydrodynamics" p 68

    :param h: this is depth in meters
    :param T: this is wave period in seconds
    """
    assert h > 0, 'Water depth must be >0, positive downward convention'
    # set vars
    g = 9.8  # gravity
    residThresh = 0.0001  # threshold for convergence  # hold to 3 decimal places
    omega = 2 * np.pi / T  # wave orbital frequency
    kh_guess = omega ** 2 * h / g  # inital guess at kh
    kh_next = 0  # initalization value
    while (np.abs(kh_next - kh_guess) > residThresh).any():
        kh_next = ((omega ** 2 * h) / g) * 1 / np.tanh(kh_guess)  # used 1/tanh instead of coth
        kh_guess = ((omega ** 2 * h) / g) * 1 / np.tanh(kh_next)  # used 1/tanh instead of coth

    k = kh_guess / h  # calculate wave number
    L = 2 * np.pi / k  # calculate wave length
    c = omega / k  # calculate wave speed
    n = 0.5 * (1 + (2 * k * h) / np.sinh(2 * k * h))  # group/ wave speed ratio

    return L, c, n


def findPeakDirs( spec, dirBin):
    """
    This function finds the peak directions for each frequency at every timestep
    :param spec: spectra in 3 dim (time, freq, dir)
    :param freqBin:  frequency bin associated with 2nd dimension of spec
    :param dirBin:   directions associated with 3rd dimension of spec
    :return:  directions s
    """
    peakDirs = np.zeros((spec.shape[0], spec.shape[1]))
    for t in range(spec.shape[0]):  # looping over time
        for idx in range(spec.shape[1]):  # looping over frequencies
            peakDirs[t, idx] = dirBin[np.argmax(spec[t, idx, :])]  # picking max direction for each frequency
    # for every frequency a direction returned
    return peakDirs


def findPeakFrqs( spec, dirBin):
    """

    :param spec: spectra in 3 dim (time, freq, dir)
    :param freqBin:  frequency bin associated with 2nd dimension of spec
    :param dirBin:   directions associated with 3rd dimension of spec
    :return:  directions s
    """
    peakFreqs = np.zeros((spec.shape[0], spec.shape[2]))
    for t in range(spec.shape[0]):  # looping over time
        for idx in range(spec.shape[1]):  # looping over directions
            peakFreqs[t, idx] = dirBin[np.argmax(spec[t, :, idx])]  # picking max frequency for each direction
    # for every frequency a direction returned
    return peakFreqs


def findTp( spec, wavefreqbin):
    """
    This function finds the Tp of a spectra
    :param spec: 2/3  dimensional spectrum [t, freq, dir]
    :param wavefreqbin:  takes both freq bin or period bin
    :return:
    """
    try:
        assert spec.ndim == 3, 'This function assumes a 3 dimensional spectrum [t, freq, dir]'
    except AssertionError:
        spec = np.expand_dims(spec, axis=0)
    assert wavefreqbin.shape[0] == spec.shape[1], 'This function expects frequency in the 2nd dimension'

    fspec = np.sum(spec, axis=2)  # frequency spec
    freqIdx = np.argmax(fspec, axis=1)  # finding the freq/period of max energy
    fOut = wavefreqbin[freqIdx]  # associating the ind with the freq/period

    return fOut


def findDpAtTp( spec, wavedirbin):
    """
    finds the peak direction at peak frequency
    :param spec:
    :param Tp:
    :param wavedirbin:
    :return:
    """
    try:
        assert spec.ndim == 3, 'This function assumes a 3 dimensional spectrum [t, freq, dir]'
    except AssertionError:
        spec = np.expand_dims(spec, axis=0)  # mkaind 3 dimensions
    assert wavedirbin.shape[0] == spec.shape[2], 'This function expects direction as the 3rd dim'
    TpIdx = np.argmax(spec.sum(axis=2), axis=1)  # finding the Peak period Idx
    DpOut = np.zeros_like(TpIdx)
    for tt in range(TpIdx.shape[0]):
        DpOut[tt] = wavedirbin[np.argmax(spec[tt, TpIdx[tt], :])]  # finding the max idx for each time

    return DpOut


def seaAndSwell2D( specTime, spec, wavefreqbin, wavedirbin, windSpeed, windDirTn, plot=False, depth=26):
    """
    This function will make sea and swell spectra given the below varibles.  The wind Directions must be in the same
    coordinate system as the wave spectra (including wave direction bins)
    :param specTime:  this is a time stamp for the data dimensioned by (time_
    :param spec:  this is the spectral wave energy data dimensioned by (time, wave direction, wave frequency)
    :param wavefreqbin:  # this are the frequencies associated with the spectral wave data
    :param wavedirbin:   # these are the direction bins associated with the spectra, same angle convention and
    :param windSpeed:  # the speed of the wind
    :param windDirTn:   the wind direction of the
    :param depth: depth is assumed to be 26 meter wave rider
    :return:  sea, swell spectra the same size as the input spectra
    """
    assert specTime.shape[0] == spec.shape[0], " data doesn't line up, try using winds from model"
    assert specTime.shape[0] == windSpeed.shape[0], "time must match between wind and wave data"
    angleWindow = 45  # this is the wind window +/-
    # 1  input periods into dispersion relationship  # get back speeds for each frequency
    _, wavespeeds, _ = dispersion(depth, 1 / wavefreqbin)

    # 2  find peak directions for each frequency band
    peakDirsTN = findPeakDirs(spec=spec, dirBin=wavedirbin)

    # 3 now seperate
    swellSea, windSea = np.ones_like(spec) * 1e-6, np.ones_like(spec) * 1e-6
    for tt in range(specTime.shape[0]):
        for ff in range(wavefreqbin.shape[0]):  # loop through frequency bin
            if ff in range(wavefreqbin.shape[0])[-3:]:
                windSea[tt, ff, :] = spec[tt, ff, :]
            if (windDirTn[tt] + angleWindow >= 360):  # if the upper bound of wind direction wraps to upper directions
                # upper bound: wave angle has to be between the 360 and the lower bound angle direction limit
                # lower bound: wave angle has to be above the wind upper bound direction limit
                if (0 <= peakDirsTN[tt, ff] <= sb.angle_correct(windDirTn[tt] + angleWindow)) | (
                                360 > peakDirsTN[tt, ff] >= windDirTn[tt] - angleWindow):
                    if windSpeed[tt] > wavespeeds[ff]:
                        # the wind is positively reinforcing wave generation --> Sea
                        windSea[tt, ff, :] = spec[tt, ff, :]
                    else:
                        swellSea[tt, ff, :] = spec[tt, ff, :]
                else:
                    swellSea[tt, ff, :] = spec[tt, ff, :]
            elif (windDirTn[tt] - 45 < 0):  # if the lower bound of wind direction wraps to negative directions
                # upper bound: the wave angle is smaller than the upper bound wind angle
                # lower bound (with wrap around 360): wave direction is either greater than zero and upper bound wind direction limit
                #   or wave direction is between the lower bound wind direction and 360
                if (0 <= peakDirsTN[tt, ff] <= sb.angle_correct(windDirTn[tt] + angleWindow)) | (
                                sb.angle_correct(windDirTn[tt] - angleWindow) <= peakDirsTN[tt, ff] < 360):
                    if windSpeed[tt] > wavespeeds[ff]:
                        # the wind is positively reinforcing wave generation --> Sea
                        windSea[tt, ff, :] = spec[tt, ff, :]
                    else:
                        swellSea[tt, ff, :] = spec[tt, ff, :]
                else:
                    swellSea[tt, ff, :] = spec[tt, ff, :]
            # if wind direction bounds do not wrap0
            elif (peakDirsTN[tt, ff] <= windDirTn[tt] + angleWindow) & (
                        peakDirsTN[tt, ff] >= windDirTn[tt] - angleWindow):
                # the speeds are worth checking because they are within 45 degrees of the wave propagation
                if windSpeed[tt] > wavespeeds[ff]:
                    # the wind is positively reinforcing wave generation --> Sea
                    windSea[tt, ff, :] = spec[tt, ff, :]
                else:
                    swellSea[tt, ff, :] = spec[tt, ff, :]
            else:
                # the wind direction is not reinforcing the wave direction
                swellSea[tt, ff, :] = spec[tt, ff, :]

        if plot == True:
            from matplotlib import pyplot as plt

            cbarMin = np.min(spec[tt])
            cbarMax = np.max(spec[tt])
            fname = 'figures/SeaSwellcompare_%s.png' % specTime[tt].strftime('%Y%m%dT%H%M%SZ')
            plt.figure(figsize=(12, 12))
            plt.suptitle('Swell And Sea Seperator \n %s' % (specTime[tt].strftime('%Y%m%dT%H%M%SZ')))

            plt.subplot(221)
            plt.title('Total spectra')
            plt.ylabel('wave Direction')
            plt.xlabel('wave Frequency')
            plt.pcolor(wavefreqbin, wavedirbin, spec[tt].T)
            plt.plot(wavefreqbin, peakDirsTN[tt], 'w.-', label='peak wave Direction')
            plt.plot(wavefreqbin, np.tile(sb.angle_correct(windDirTn[tt] + angleWindow), len(wavefreqbin)), 'w--')
            plt.plot(wavefreqbin, np.tile(sb.angle_correct(windDirTn[tt] - angleWindow), len(wavefreqbin)), 'w--')
            plt.plot([0, 0.5], [windDirTn[tt], windDirTn[tt]], 'w-', label='Wind direction')
            # plt.legend(loc='upper right')

            ax1 = plt.subplot(222)
            plt.title('wind and wave speeds')
            plt.xlabel('wave Frequency')
            plt.plot(wavefreqbin, wavespeeds, 'b-', label='wavespeed')
            plt.plot(wavefreqbin, np.tile(windSpeed[tt], len(wavefreqbin)), 'k-', label='Windspeed +')
            ax1.set_ylabel('along wind Componant of WaveSpeed', color='b')
            ax1.text(.2, (windSpeed[tt] + 2) / 2, 'Wind Sea')
            ax1.text(.2, (windSpeed[tt] + 2), 'Swell Sea')
            ax1.legend(loc='best')

            plt.subplot(223)
            plt.title('WindSea')
            plt.pcolor(wavefreqbin, wavedirbin, windSea[tt].T)  # , vmin=cbarMin, vmax=cbarMax)
            plt.ylabel('wave Direction')
            plt.xlabel('wave Frequency')
            plt.colorbar()

            plt.subplot(224)
            plt.title('Swell Sea')
            plt.pcolor(wavefreqbin, wavedirbin, swellSea[tt].T, vmin=cbarMin, vmax=cbarMax)
            plt.ylabel('wave Direction')
            plt.xlabel('wave Frequency')
            plt.colorbar()

            plt.tight_layout(pad=.1, rect=[.03, .03, .95, .9])
            plt.savefig(fname)
            plt.close()

    return windSea, swellSea


def seaAndSwell1D(specTime, spec, wavefreqbin, truncate=0.1):
    """

    :param specTime:
    :param spec:
    :param wavefreqbin:
    :param truncate: value used in save munits as the wave freqbin
    :return:
    """
    assert spec.shape[-1] == len(wavefreqbin), '1D stats need a 1 d spectra'
    fspec = spec

    # 3 now seperate
    swellSea, windSea = np.ones_like(fspec) * 1e-6, np.ones_like(fspec) * 1e-6

    idx = np.argmax(wavefreqbin >= truncate)

    if fspec.ndim == 1:
        windSea[idx:] = fspec[idx:]
        swellSea[:idx] = fspec[:idx]
    else:
        windSea[:, idx:] = fspec[:, idx:]
        swellSea[:, :idx] = fspec[:, :idx]

    return windSea, swellSea


def stats1D(fspec, frqbins, dspec=None, lowFreq=0.05, highFreq=0.5):
    """

    :param fspec: frequency spectra
    :param frqbins: frequency bins associated with the 1d spectra
    :param dspec:  this is the mean direction as a function of frequency ie arctan2(b1,a1)
    :param lowFreq: low frequency cut off for analysis
    :param highFreq: high frequency cutoff for analysis
    :return: a dictionary with statistics
    """
    assert fspec.shape[-1] == len(frqbins), '1D stats need a 1 d spectra'
    fspec = fspec
    frqbins = np.array(frqbins)

    df = np.diff(np.append(frqbins[0], frqbins), n=1)

    # truncating spectra to sea/swell band
    [idx, vals] = sb.findbtw(frqbins, lowFreq, highFreq, type=3)

    m0 = np.sum(fspec * df, axis=1)  # 0th momment
    m1 = np.sum(fspec[:, idx] * df[idx] * frqbins[idx], axis=1)  # 1st moment
    m2 = np.sum(fspec[:, idx] * df[idx] * frqbins[idx] ** 2, axis=1)  # 2nd moment
    m3 = np.sum(fspec[:, idx] * df[idx] * frqbins[idx] ** 3, axis=1)  # 3rd moment
    m4 = np.sum(fspec[:, idx] * df[idx] * frqbins[idx] ** 4, axis=1)  # 4th moment
    m11 = np.sum(fspec[:, idx] * df[idx] * frqbins[idx] ** -1, axis=1)  # negitive one moment
    # wave height
    Hm0 = 4 * np.sqrt(m0)

    ipf = fspec.argmax(axis=1)  # indix of max frequency
    Tp = 1 / frqbins[ipf]  # peak period
    Tm02 = np.sqrt(m0 / m2)  # mean period
    Tm01 = m0 / m1  # average period - cmparible to TS Tm
    Tm10 = m11 / m0
    sprdF = (m0 * m4 - m2 ** 2) / (m0 * m4)

    meta = 'Tp - peak period, Tm - mean period, Tave - average period, comparable to Time series mean period, sprdF - frequency spread, sprdD - directional spread'
    stats = {'Hm0': Hm0,
             'Tp': Tp,
             'Tm': Tm02,
             'Tave': Tm01,
             'sprdF': sprdF,
             'Tm10': Tm10,
             'meta': meta
             }

    return stats


def waveStat(spec, frqbins, dirbins, lowFreq=0.05, highFreq=0.5):
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
    assert type(frqbins) in [np.ndarray, np.ma.MaskedArray], 'the input frqeuency bins must be a numpy array'
    assert type(dirbins) in [np.ndarray, np.ma.MaskedArray], 'the input DIRECTION bins must be a numpy array'
    assert np.array(spec).ndim == 3, 'Spectra must be a 3 dimensional array'
    try:
        assert (spec != 0).all() is not True, 'Spectra must have energy to calculate statistics, all values are 0'
    except AssertionError:
        return 0
    assert spec.shape[1] == frqbins.shape[0], 'The spectra must be dimensioned [t,freq, direciton]'
    assert spec.shape[2] == dirbins.shape[0], 'The spectra must be dimensioned [t,freq, direciton]'
    frqbins = np.array(frqbins)
    dirbins = np.array(dirbins)

    # finding delta freqeucny (may change as in CDIP spectra)
    # frq = np.array(np.zeros(len(frqbins) + 1))  # initializing frqbin bucket
    # # frq[0] = frqbins[0]
    # # frq[1:] = frqbins
    # # df = np.diff(frq, n=1)  # dhange in frequancy banding
    df = np.diff(np.append(frqbins[0], frqbins), n=1)
    dd = np.abs(np.median(np.diff(dirbins)))  # dirbins[2] - dirbins[1]  # assume constant directional bin size
    # finding delta degrees
    # frequency spec
    fspec = np.sum(spec, axis=2) * dd  # fd spectra - sum across the frequcny bins to leave 1 x n-frqbins
    # doing moments over 0.05 to 0.33 Hz (3-20s waves) (mainly for m4 sake)
    [idx, vals] = sb.findbtw(frqbins, lowFreq, highFreq, type=3)

    m0 = np.sum(fspec * df, axis=1)  # 0th momment
    m1 = np.sum(fspec[:, idx] * df[idx] * frqbins[idx], axis=1)  # 1st moment
    m2 = np.sum(fspec[:, idx] * df[idx] * frqbins[idx] ** 2, axis=1)  # 2nd moment
    # m3 = np.sum(fSpecOut[:, idx] * df[idx] * frqbins[idx] ** 3, axis=1)  # 3rd moment
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
    Ds = np.sum(spec * np.tile(df, (len(dirbins), 1)).T, axis=1)  # directional spectra
    Dsp = []
    for ii in range(0, len(ipf)):
        Dsp.append(spec[ii, ipf[ii], :])  # direction spectra at peak-f
    Dsp = np.array(Dsp)
    idp = Dsp.argmax(axis=1)  # index of direction at peak frquency
    Dp = dirbins[idp]  # peak direction

    Drad = np.deg2rad(dirbins)  # making a radian degree bin
    # mean wave direction (e.g. Kuik 1988, USACE WIS)
    Xcomp = np.sum(np.sin(Drad) * Ds, axis=1) / np.sum(Ds * dirbins, axis=1)
    Ycomp = np.sum(np.cos(Drad) * Ds, axis=1) / np.sum(Ds * dirbins, axis=1)
    # Dm = np.rad2deg(np.arctan2(np.sum(np.sin(Drad) * Ds * dirbins, axis=1),
    #                           np.sum(np.cos(Drad) * Ds * dirbins, axis=1)))  # converting back to degrees
    Dm = np.rad2deg(np.arctan2(Xcomp, Ycomp))
    for ii in range(0, np.size(Dm, axis=0)):
        if Dm[ii] >= 360:
            Dm[ii] = 360 - Dm[ii]
        elif Dm[ii] < 0:
            Dm[ii] = 360 + Dm[ii]
            # Vector Dm (Hesser)
    sint = np.sin(Drad)  # sine of dirbins
    cost = np.cos(Drad)  # cosine of dirbins

    sint2 = np.tile(sint, [len(frqbins), 1])  # 2d diretion size of the spectra
    cost2 = np.tile(cost, [len(frqbins), 1])
    delsq = np.tile(df, [len(dirbins), 1]).T

    xsum = np.zeros(np.size(spec, axis=0))
    ysum = np.zeros(np.size(spec, axis=0))
    # summing across all directions
    for ii in range(0, np.size(spec, axis=0)):
        xsum[ii] = sum(np.sum(cost2 * delsq * spec[ii, :, :], axis=1))  # summing along directions, then
        ysum[ii] = sum(np.sum(sint2 * delsq * spec[ii, :, :], axis=1))

    vavgdir = np.arctan2(ysum, xsum)
    vavgdir = np.rad2deg(vavgdir)
    vavgdir = sb.angle_correct(vavgdir)

    # Mean direction at the peak frequency
    Dmp = np.rad2deg(np.arctan2(np.sum(np.sin(Drad) * Dsp * dirbins, axis=1),
                                np.sum(np.cos(Drad) * Dsp * dirbins, axis=1)))  # converting back to degrees
    for ii in range(0, np.size(Dmp, axis=0)):
        if Dmp[ii] >= 360:
            Dmp[ii] = 360 - Dmp[ii]
        elif Dmp[ii] < 0:
            Dmp[ii] = 360 + Dmp[ii]
    # f-spec spread
    sprdF = (m0 * m4 - m2 ** 2) / (m0 * m4)

    # fd-spec spread
    sprdD = np.rad2deg(np.sqrt(2.0 * (1.0 - np.sqrt(Xcomp ** 2 + Ycomp ** 2))))

    ##### Exceprt from Kent's code for spreading - not sure how to handle
    # fd-spec spread, do a linear interp to get closer to half-power
    # % from the delta-deg increments
    # hp = np.max(Dsp)/2;
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
             'Tm10': Tm10,
             'meta': meta
             }
    # print meta
    return stats


def fSpecPeaksValleys(spec1d, wavefreqbin):
    """
    This function takes a 1 dimensional frequency spectra and finds
    :param spec1d: 1 d frequency spectra
    :param plotfname:  file name for figure to be made
    :return: peakindexes, valley indexes
    """
    # import peakutils
    # from scipy import signal
    # # first smooth the spectra to find appropriate peaks
    # # smoothSpec1d = signal.savgol_filter(spec1d, window_length=5, polyorder=2)  # smoothing the spectra
    # smoothSpec1d = spec1d
    # mindist = 3  # this is the minimum number of bins between peaks to be identified
    # threshold = 0.3  # normalized peak threshold
    # peakIdxs = np.squeeze(signal.argrelmax(smoothSpec1d, order=mindist))

    # find indicies of peaks
    # peakIdxs = peakutils.indexes(smoothSpec1d, thres=threshold, min_dist=mindist)

    import sys
    sys.path.append('/home/spike/repos/peakdetect')
    import peakdetect  # https://gist.github.com/sixtenbe/1178136
    lookahead = 3  # len(spec1d) / 15
    delta = 0.08  # minimum difference between a peak and a valley
    peakIdxs, valIdxs = peakdetect.peakdetect(spec1d, wavefreqbin, lookahead=lookahead, delta=delta)

    #
    # # are there more than one peak in the spectra
    # if len(peakIdxs) == 1:  # there's only one frequency in the spectra
    #     valIdxs = None
    # else:  # find the valleys for associated peaks
    #     valIdxs = [] # np.zeros(len(peakIdxs)-1, dtype=int)  # valley indicies
    #     deleteThese = []
    #     for ii in range(0,len(peakIdxs)-1):
    #         neworder = 1
    #         comp = np.squeeze(signal.argrelmin(smoothSpec1d[peakIdxs[ii]:peakIdxs[ii+1]], order=neworder))
    #         while np.size(comp) > 1:  # find 1 valley between 2 peaks
    #             neworder += 1 # increase the order until one valley is found
    #             comp = np.squeeze(signal.argrelmin(smoothSpec1d[peakIdxs[ii]:peakIdxs[ii+1]], order=neworder))
    #             if neworder > 10:
    #                 raise OverflowError, 'bad solution to finding peak'
    #         if np.size(comp) == 0:
    #             deleteThese.append(ii + 1)
    #         else:
    #             valIdxs.append(comp + peakIdxs[ii]) # valley indices are relative to peaks
    #     # remove bad found peaks
    #     peakIdxs = np.delete(peakIdxs, deleteThese)
    #

    ################################################################
    # # begin QA QC for peak Check

    if len(peakIdxs) == 0:  # if no peaks were found (small wave conditions )
        peakIdxs = [[wavefreqbin[np.argmax(spec1d)], spec1d[np.argmax(spec1d)]]]
    elif len(peakIdxs) == len(valIdxs):
        valIdxs.remove(valIdxs[-1])  # remove last valley found, should be odd number
    # # first check for half power
    #     deletedPeaks = []  # deleted indices that don't meet half power criterai
    #     for ff in range(len(peakIdxs)-1):
    #         print 'checking the peaks'
    #         val = spec1d[valIdxs[ff]]    # energy at the valley
    #         peak = spec1d[peakIdxs[ff]]  # energy at the peak
    #         if peak/2 < val:  # if valley is larger than half peak energy delete it
    #             deletedPeaks.append(ff)
    #         # if peak
    # now delete all selected peaks
    # peakIdxs = np.delete(peakIdxs, deletedPeaks)
    # valIdxs = np.delete(valIdxs, deletedPeaks)
    # if len(peakIdxs) > 3:
    #     print 'Too Many Peaks Found!!!'

    ################################################################
    # plot solution
    # if plotfname != None:
    #     plt.figure()
    #     plt.subplot(211)
    #     plt.plot(wavefreqbin, spec1d, label='Spectra')
    #     # plt.plot(wavefreqbin, smoothSpec1d, label='Smoothed')
    #     for idx in peakIdxs:
    #         plt.plot(idx[0], idx[1], 'rd', label='peaks')
    #     for idx in valIdxs:
    #         plt.plot(idx[0], idx[1], 'go', label='valleys' )
    #     plt.ylabel('Energy $m^2/Hz$')
    #     plt.legend()
    #     plt.subplot(212)
    #     plt.plot(wavefreqbin, dspec)
    #     for idx in peakIdxs:
    #         plt.plot(idx[0], dspec[np.squeeze(np.argwhere(idx[0] == wavefreqbin))], 'rd', label='peaks')
    #     plt.ylabel('Direction [deg TN]')
    #     plt.savefig(plotfname)
    #     plt.close()
    return peakIdxs, valIdxs


def isThisWindSea(waveDirectionAtPeak, waveSpeedAtPeak, windDir, windSpeed):
    """

    :rtype: object
    :param waveDirectionAtPeak:  This is a wave direction in degrees (must be same coordinate sys as wind dir)
    :param waveSpeedAtPeak:      This is speed of wave at peak (speed must be same units as wind Dir)
    :param windDir:              this is wind Direction (time matched to Wave direction)
    :param windSpeed:            this is wind speed, Time matched to wave Direction
    :return: True or False value
    """
    windowAngle = 45
    windAngleMax = sb.angle_correct(windDir + windowAngle)
    windAngleMin = sb.angle_correct(windDir - windowAngle)
    # if there is no wrap
    windSea = False  # setting
    if waveSpeedAtPeak <= windSpeed:  # if speeds match up
        if (windDir > 360 - windowAngle) and ((360 > waveDirectionAtPeak >= windAngleMin) or (
                        0 <= waveDirectionAtPeak <= windAngleMax)):  # wind window wraps above 360
            # upper bound: wave angle has to be between the 360 and the lower bound angle direction limit
            # lower bound: wave angle has to be above the wind upper bound direction limit
            windSea = True
        elif (windDir < windowAngle) and ((0 <= waveDirectionAtPeak <= windAngleMax) or (
                        windAngleMin < waveDirectionAtPeak <= 360)):  # wind window wraps below 0
            # upper bound: the wave angle is smaller than the upper bound wind angle
            # lower bound (with wrap around 360): wave direction is either greater than zero and upper bound wind direction limit
            #   or wave direction is between the lower bound wind direction and 360
            windSea = True
        elif (waveDirectionAtPeak > windAngleMin) and (
                    waveDirectionAtPeak < windAngleMax):  # if angles match up with no wrap
            windSea = True

    return windSea


def findidxSeaSwell(spec1d, dspec, wavefreqbin, windSpeed, windDir, depth, plotfname=None):
    """
    This function separates sea and swell starting with the 1d frequency spectra and the 1d directional spectra

    :param spec1d: 1 d frequency spectrum
    :param wavefreqbin: frequency spectrum associated with wave
    :param dspec: 1d directional spectrum (arctan(b1,a1))
    :param waveDirectionAtPeak: the direction of the wave from
    :param plotfname:
    :return: an index corresponding to
    """
    # body of function
    assert (dspec >= 0).all(), 'please correct wave angles between 0 and 360'
    assert len(windDir) == spec1d.shape[0], 'Wind records are not the same length as the frequnecy spectra'
    spec1d = np.expand_dims(spec1d, axis=0) if spec1d.ndim == 1 else spec1d  # expanding dims as necessicary
    dspec = np.expand_dims(dspec, axis=0) if dspec.ndim == 1 else dspec  # expanding the dimesnions is correct
    # loop through how ever many spectra there are
    outidxs = []  # initialization out output
    for i in range(spec1d.shape[0]):
        swellfspec, seafspec = np.ones_like(spec1d) * 1e-6, np.ones_like(spec1d) * 1e-6
        peaks, valleys = fSpecPeaksValleys(spec1d[i], wavefreqbin)  # finding peaks and valleys in fspec
        if valleys == None or len(valleys) == 0:  # there are not 2 wave system componants in the spectra
            waveDirectionAtPeak = dspec[i, np.argwhere(peaks[0][0] == wavefreqbin).squeeze()]
            _, waveSpeedAtPeak, _ = dispersion(depth.squeeze(), 1 / peaks[0][0])
            # find out whether it's sea or swell, (wave direction vs wind direction +/- 45) and log it appropriately
            isit = isThisWindSea(waveDirectionAtPeak, waveSpeedAtPeak, windDir=windDir[i], windSpeed=windSpeed[i])
            if isit == True:
                # this is sea
                idx2separate = 0  # this is wind Sea
            else:
                idx2separate = -1  # ththis is swell
            waveDirectionAtPeak = [waveDirectionAtPeak]
        else:
            # seperate through sea and swell
            waveDirectionAtPeak = []  # initialization
            for ii in range(len(peaks)):  # looping though peaks finding wave direction at each
                waveDirectionAtPeak.append(dspec[i, np.argwhere(peaks[ii][0] == wavefreqbin).squeeze()])
            for ii in range(len(peaks)):  # looping through each valley starting at
                # waveDirectionAtPeak= dspec[np.argwhere(peaks[ii][0] == wavefreqbin).squeeze()]  #
                _, waveSpeedAtPeak, _ = dispersion(depth.squeeze(), 1 / peaks[ii][0])
                isit = isThisWindSea(waveDirectionAtPeak[ii], waveSpeedAtPeak, windDir=windDir[i],
                                          windSpeed=windSpeed[i])
                if isit == True:
                    # its windsea, pick the index and don't go any further
                    idx2separate = np.argwhere(valleys[ii - 1][0] == wavefreqbin).squeeze()  # separation frequency
                    break
                # if isit is False, then we have swell, move to the next peak to find sea
                else:
                    idx2separate = -1  # this is swell
        if plotfname != None:

            if waveSpeedAtPeak > windSpeed[i]:
                title = 'waveSpeed > windSpeed - Swell %s' % plotfname[i].split('/')[-1][5:-1]
            elif waveSpeedAtPeak <= windSpeed[i]:
                title = 'waveSpeed <= windSpeed - MaybeSea %s' % plotfname[i].split('/')[-1][5:-1]

            print plotfname[i].split('/')[-1][5:-1]
            plt.figure()
            plt.suptitle(title)
            plt.subplot(211)
            plt.plot(wavefreqbin, spec1d[i], label='Frequency Spectra')
            # plt.plot(wavefreqbin, smoothSpec1d, label='Smoothed')
            for idx in peaks:
                plt.plot(idx[0], idx[1], 'gx')
            for idx in valleys:
                plt.plot(idx[0], idx[1], 'gx')
            plt.plot(wavefreqbin[idx2separate], spec1d[i, idx2separate], 'rd', label='wave separation frequanc6y')
            plt.ylabel('Energy $m^2/Hz$')
            plt.legend()
            plt.subplot(212)
            plt.plot(wavefreqbin, dspec[i])
            plt.plot([wavefreqbin[0], wavefreqbin[-1]], [windDir[i], windDir[i]], '--', label='Wind Direction')
            plt.plot(wavefreqbin[idx2separate], dspec[i, idx2separate], 'rd', label='Separation Point')
            for ii in range(len(peaks)):
                plt.plot(peaks[ii][0], waveDirectionAtPeak[ii], 'gx')
                plt.text(peaks[ii][0], waveDirectionAtPeak[ii], '  diff = %i' % min(
                    360 - np.abs(waveDirectionAtPeak[ii] - windDir[i]), np.abs(waveDirectionAtPeak[ii] - windDir[i])))
            plt.ylabel('Direction [deg TN]')
            plt.xlabel('Wave Frequency')
            plt.legend()
            plt.savefig(plotfname[i])
            plt.close()

        #
        # swellfspec[:idx2separate] = spec1d[:idx2separate]
        # seafspec[idx2separate:] = spec1d[idx2separate:]
        # Sea.append(seafspec)
        # Swell.append(swellfspec)
        outidxs.append(idx2separate)
    return outidxs


#def decompose2Dspec(spec, wavefreqbin):

    fspec = spec.sum(axis=-1)  # integrating across all directions

    a1, b1 = np.zeros_like(fspec), np.zeros_like(fspec)
    for freq in wavefreqbin:
        a1 = np.arccos(fspec)

        b1 = np.arcsin(fspec)