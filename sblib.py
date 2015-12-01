# -*- coding: utf-8 -*-
"""
Created on Sat Oct 10 10:18:53 2015

@author: u4hncasb
"""
import numpy as np
import datetime as DT


'''
my own library with useful functions and tips from which i found to be helpful
on various general functions that don't fit into specific other codes
'''
__name__
def find_nearest(array,value):
	'''
	Function looks for value in array and returns the closest array value 
	(to 'value') and index of that value 
	
	''' 
	idx = (np.abs(array-value)).argmin()
	return array[idx],idx
def SBcleanangle(directions,deg=360):
      '''
	This function cleans an array of angles (in degrees) to all positive 
	values ranging from 0 to 360
	
	Currently is designed for only degree angles
	'''
      for ii in range(0,len(directions)):
		if directions[ii]>=360:            
			directions[ii]=directions[ii]-360
		elif directions[ii]<0:
			directions[ii]=directions[ii]+360
      return directions
def FRFcoord(p1,p2):
    '''
    #  function [ALat, ALon, spN, spE, Y, X] = frfCoord(p1, p2)
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
    '''
    r2d = 180.0 / np.pi;

    Eom=901951.6805;               #% E Origin State Plane
    Nom=274093.1562;               #% N Origin State Plane
    #ALat0=10.65583950;            % Origin Lat minutes
    #ALon0=44.9811435;             % Origin Lon minutes
    ALat0=36.1775975;              #% Origin Lat minutes
    ALon0=75.7496860;              #% Origin Lon minutes
    DegLat = 110963.35726;         #% m/deg Lat
    DegLon = 89953.36413;          #% m/deg long
    GridAngle=18.1465 /r2d;
    spAngle = (90 - 69.974707831) / r2d

	# Determine Data type
    if np.floor(abs(p1))==75 and np.floor(p2)==36: # lat/lon input
		# to FRF coords 
		ALat=p1
		ALon=p2  # DESIGNATING LAT/LON VARS
		if p1<0:
			p1=-p1
		ALatLeng = (p2 - ALat0) * DegLat
		ALatLeng = -(p1 - ALon0) * DegLon
		R = np.sqrt(ALatLeng**2 + ALonLeng**2)
		Ang1 = np.arctan2(ALongLeng, ALatLeng)
		Ang2 = Ang1 + GridAngle;
		X = R * np.sin(Ang2)
		Y = R * np.cos(Ang2)
		# to StatePlane
		Ang2= Ang2 - spAngle
		AspN = R * np.cos(Ang2)
		AspE = R * np.sin(Ang2)
		spN = AspN + Nom
		spE = AspE + Eom

    elif (p1 > 800000) and p2 >200000: # state plane input
	   spE=p1
	   spN=p2  # designating stateplane vars
	   # to FRF coords 
	   spLengE = p1 - Eom
	   spLengN = p2 - Nom
	   R = np.sqrt(spLengE**2 + spLengN**2)
	   Ang1 = np.arctan2(spLengE,spLengN)
	   Ang2 = Ang1 + spAngle
	   X = R * np.sin(Ang2)
	   Y = R * np.sin(Ang2)
	   # to Lat Lon            
	   Ang2 = Ang1 - (GridAngle-spAngle)   #         % 
	   ALatLeng = R * np.cos(Ang2)
	   ALonLeng = R * np.sin(-Ang2)         #% neg to go west
	   ALat = ALatLeng/DegLat + ALat0    #% was 60 * ALatLeng./DegLat + ALat0;
	   ALon = ALonLeng/DegLon + ALon0

    elif (p1 > -10000 & p1 < 10000) and (p2 > -10000 & p2 < 10000):  # FRF input 
	   X=p1
	   Y=p2;
	   R = np.sqrt(p1**2 + p2**2);
	   Ang1 = np.arctan2(p1, p2);    #      % CW from Y
	   Ang2 = Ang1 - GridAngle;#            % 
	   ALatLeng = R * np.cos(Ang2);
	   ALonLeng = R * np.sin(-Ang2); #     % neg to go west
	   ALat = ALatLeng/DegLat + ALat0;#     % was 60 * ALatLeng./DegLat + ALat0;
	   ALon = ALonLeng/DegLon + ALon0;

		 #  to state plane 
	   Ang2 = Ang1 - spAngle;
	   AspN = R * np.cos(Ang2);
	   AspE = R * np.sin(Ang2);             
	   spN = AspN + Nom;
	   spE = AspE + Eom;


    else:   
          print '<EE> Cound not determine input type, returning NaNs'
          ALat=float('NaN'); ALon=float('NaN'); spN=float('NaN'); spE=float('NaN'); Y=float('NaN'); X=float('NaN');
    coords={'StateplaneE':spE,'StateplaneN':spN,'FRF_Y':Y,'FRF_X':X,'Lat':ALat,'Lon':ALon}
    return coords


def sbfindbtw(data,lwth,upth,type=0):
	'''
	This function finds both values and indicies of a list values between two values
	upth = upper level threshold
	lwth = lower level threshold
	list = list (or numpy array?)
	type: 
		0 = non inclusive  ie. lwth < list <  upth
		1 = low incluisve  ie. lwth <=list <  upth 
		2 = high inclusive ie. lwth < list <= upth
		3 = all inclusive  ie  lwth <=list <= upth
	'''
	indices=[]
	vals=[]
	shp=np.shape(data)
	if len(shp)==2:
		for i,  in enumerate(data):
			for j, elem in enumerate(range):
				if type==0:
					if elem < upth and elem > lwth:
						indices.append((i, j))
						vals.append(elem)
				elif type==1:
					if elem < upth and elem >= lwth:
						indices.append((i, j))
						vals.append(elem)
				elif type==2:
					if elem <= upth and elem > lwth:
						indices.append((i, j))
						vals.append(elem)
				elif type==3:
					if elem <= upth and elem >= lwth:
						indices.append((i, j))
						vals.append(elem)
	if len(shp)==1:
		for j, elem in enumerate(data):
			if type==0:
				if elem < upth and elem > lwth:
				   indices.append((j))
				   vals.append(elem)
			elif type==1:
				if elem < upth and elem >= lwth:
				   indices.append((j))
				   vals.append(elem)
			elif type==2:
			   if elem <= upth and elem > lwth:
				   indices.append((j))
				   vals.append(elem)
			elif type==3:
			   if elem <= upth and elem >= lwth:
				   indices.append((j))
				   vals.append(elem)               
					
	return indices, vals
def roundtime(dt=None, roundTo=60):
	   """Round a datetime object to any time laps in seconds
	   dt : datetime.datetime object, default now.
	   roundTo : Closest number of seconds to round to, default 1 minute.
	   Author: Thierry Husson 2012 - Use it as you want but don't blame me.
	   """
	   if dt == None : dt = DT.datetime.now()
	   seconds = (dt - dt.min).seconds
	   # // is a floor division, not a comment on following line:
	   rounding = (seconds+roundTo/2) // roundTo * roundTo
	   return dt + DT.timedelta(0,rounding-seconds,-dt.microsecond)  
	   
def cart2pol(x,y):
	''' 
	this translates from cartesian coords to polar coordinates (radians)
	'''
	r = np.sqrt(x**2 + y**2)
	theta = np.arctan2(y, x)
	return r, theta
	
def pol2cart( r, theta):
	''' 
	this translates from polar coords (radians) to polar coordinates
	'''
	x = r *np.cos(theta)
	y = r *np.sin(theta) 
	return x, y
def wavestat(spec,dirbins,frqbins):
	"""
	this function will calculate the mean direction from a full spectrum
	only calculates on one 2D spectrum at a time 
	
	Input:
	%     spec  Frequency-direction spectra (2D)       shape(record,frqbin,dirbin)
	%  frqbins  Frequency vector (not assumed constant)
	%  dirbins  Direction vector (assumed constant)
	%
	% Outputs (MKS, Hz, degrees, degrees CW from true north):
	%   Hmo   Significant wave height
	%    Tp   Period of the peak energy in the frequency spectra, (1/Fp).  AKA Tpd, not to be
	%           confused with parabolic fit to spectral period
	%    Tm   Mean spectral period (Tm0,2, from moments 0 & 2), sqrt(m0/m2)
	%  Tave   Average period, frequency sprectra weighted, from first moment (Tm0,1)
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
		
	Code Translated by Spicer Bak from: fd2BulkStats.m written by Kent Hathaway
		
	"""
#        from prep_data import prep_data as pd
#        pd=pd()
	# finding delta freqeucny (may change as in CDIP spectra)
	frq=np.array(np.zeros(len(frqbins)+1)) #initializing frqbin bucket
	frq[0]=frqbins[0]
	frq[1:]=frqbins
	df=np.diff(frq,n=1)
	# finding delta degrees
#        dtheta=np.median(pd.dtfun(dirbins)) #the delta degrees assuing constant
	#frequency spec
	fspec=np.sum(spec,axis=2) #fd spectra - sum across the frequcny bins to leave 1 x n-frqbins
	# doing moments over 0.05 to 0.33 Hz (3-20s waves) (mainly for m4 sake)
	[idx,vals]=sblib.findbtw(frqbins,0.05,0.5,type=3)        
		
	m0=np.sum(fspec*df,axis=1)
	m1=np.sum(fspec[:,idx] * df[idx] * frqbins[idx],axis=1)
	m2=np.sum(fspec[:,idx] * df[idx] * frqbins[idx]**2,axis=1)
	m3=np.sum(fspec[:,idx] * df[idx] * frqbins[idx]**3,axis=1)
	m4=np.sum(fspec[:,idx] * df[idx] * frqbins[idx]**4,axis=1)
	
	
	#sigwave height
	Hm0=4*np.sqrt(m0)
	# period stuff
	ipf=fspec.argmax(axis=1) # indix of max frequency
	Tp=1/frqbins[ipf]        # peak period
	Tm=np.sqrt(m0/m2)        # mean period
	Tave=m0/m1               # average period - cmparible to TS Tm
	
	# directional stuff
	Ds=np.sum(spec*np.tile(df,(len(dirbins),1)).T,axis=1) # directional spectra
	Dsp=[]
	for ii in range(0,len(ipf)):        
		Dsp.append(spec[ii,ipf[ii],:])  #direction spectra at peak-f
	Dsp=np.array(Dsp)
	idp=Dsp.argmax(axis=1)              #index of direction at peak frquency
	Dp=dirbin[idp]                      #peak direction
	
	Drad=np.deg2rad(dirbins)  # making a radian degree bin
	# mean wave direction (e.g. Kuik 1988, USACE WIS)
	Xcomp=np.sum(  np.sin(Drad) * Ds * dirbins,axis=1) / np.sum(Ds * dirbins,axis=1)
	Ycomp=np.sum(  np.sin(Drad) * Ds * dirbins,axis=1) / np.sum(Ds * dirbins,axis=1)
	Dm= np.rad2deg(np.arctan2( np.sum(np.sin(Drad) * Ds * dirbins,axis=1) , np.sum(  np.sin(Drad) * Ds * dirbins,axis=1))) # converting back to degrees
	
	
	# Mean direction at the peak frequency
	Dmp= np.rad2deg(np.arctan2( np.sum(np.sin(Drad) * Dsp * dirbins,axis=1) , np.sum(  np.sin(Drad) * Dsp * dirbins,axis=1))) # converting back to degrees
	
	# f-spec spread
	sprdF= (m0*m4 - m2**2)/(m0*m4)
	
	#fd-spec spread
	sprdD=np.rad2deg( np.sqrt(2.0* (1.0-sqrt(Xcomp**2 + Ycomp**2) )))
	
	#fd-spec spread, do a linear interp to get closer to half-power from the
	# delta-deg increments
	##### Exceprt from Kent's code for spreading - not sure how to handle
	#        % fd-spec spread, do a linear interp to get closer to half-power
	#% from the delta-deg increments
	#hp = max(Dsp)/2;
	#ihp=find(Dsp > hp);
	#
	#  % Left (d1) and right (d2) interps: Y=Dir, X=E   
	#d1=interp1([Dsp(ihp(1)-1) Dsp(ihp(1)+1)], [dwdir(ihp(1)-1) dwdir(ihp(1)+1)], hp);  
	#d2=interp1([Dsp(ihp(end)-1) Dsp(ihp(end)+1)], [dwdir(ihp(end)-1) dwdir(ihp(end)+1)], hp);  
	#sprdDhp = d2 - d1;
	return Hm0, Tp, Tm, Tave,  Dp, Dm, Dmp, sprdF, sprdD