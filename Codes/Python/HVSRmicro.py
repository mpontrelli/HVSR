# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 18:04:24 2020

@author: mpontr01

1/16/2020 added everything beyond window the data and developed code wavav
"""
    
##filename = 'C:\\Users\\mpontr01\\Desktop\\Stations\\Tufts University\\2019_10_26\\Data\\bd_mseed_A00L_1572134400'

## inputs

ftype = '.sac' # filetype, this is to support a bunch of filtypes
sav = 'no' # do you want to save your figures?
outpath = 'C:\\Users\\mpontr01\\Desktop\\Stations\\Tufts University\\2019_10_26\\Figures' # if you want to save your figures use outpath
fs = 100 # I think we'll change this to read in from the file
## filenames, these are for .sac
Vfname = 'C:\\Users\\mpontr01\\Desktop\\boston_site_response\\field_deployments\\Boston area\\908DP\\Data\\2000004171152005_908DP_1_1.sac'
NSfname = 'C:\\Users\\mpontr01\\Desktop\\boston_site_response\\field_deployments\\Boston area\\908DP\\Data\\2000004171152005_908DP_1_2.sac'
EWfname = 'C:\\Users\\mpontr01\\Desktop\\boston_site_response\\field_deployments\\Boston area\\908DP\\Data\\2000004171152005_908DP_1_3.sac'
statname = '908DP'

Allplots = 'yes'
Timeplot = 'no'
HVSRplottog = 'no'
## defaults for windowing the time series
numwin = 10
windowlen = 40
fs = 100
windis = 25
sampnum = windowlen*fs
windisnum = windis *fs

## Filter inputs
lowcut = 0.5
highcut = fs/2 - 1
order = 4

## bounds for plotting and statistcs calculation
upbound = 49
lowbound = 0.1


if ftype in '.sac':
    
    ## import the data
    stv = read(Vfname, debug_headers=True)
    xV = stv[0]
    stNS = read(NSfname, debug_headers=True)
    xNS = stNS[0]
    stEW = read(EWfname, debug_headers=True)
    xEW = stEW[0]
    
    ## now detrend the data

    xV = xV-np.mean(xV)
    xNS = xNS-np.mean(xNS)
    xEW = xEW-np.mean(xEW)

    ## extract values from the metadata, we only need to do this for 1 channel
    station = stv[0].stats.station
    # components
    channel = stv[0].stats.channel
    # fs
    fs = stv[0].stats.sampling_rate
    # length of vector
    npts = stv[0].stats.npts
    # delta t
    delta = stv[0].stats.delta
    # endtime of record for title
    endtime = stv[0].stats.endtime
    # make a time vector
    T = npts/fs
    t = np.linspace(0, T, npts, endpoint=False)/60
    
    ## Now filter the data
    xV = butter_bandpass_filter(xV, lowcut, highcut, fs, order=order)
    xNS = butter_bandpass_filter(xNS, lowcut, highcut, fs, order=order)
    xEW = butter_bandpass_filter(xEW, lowcut, highcut, fs, order=order)
    
    ## now plot the time series if you want
    if Allplots in 'yes' or Timeplot in 'yes':
        timeseriesplot(xNS, xV, xEW, fs, sav, outpath)

    ## Window the data
    xVmatrix=np.zeros([numwin, int(windowlen*fs -1)])
    xNSmatrix=np.zeros([numwin, int(windowlen*fs -1)])
    xEWmatrix=np.zeros([numwin, int(windowlen*fs -1)])
    k = (1,fs)
    for i in range(numwin):
        xVmatrix[i,:] = xV[(k[0]):int((k[1] * windowlen + k[0])-1)]
        xNSmatrix[i,:] = xNS[(k[0]):int((k[1] * windowlen + k[0])-1)]
        xEWmatrix[i,:] = xEW[(k[0]):int((k[1] * windowlen + k[0])-1)]
        k = (k[0] + sampnum + windisnum, fs)
    xHmatrix = xNSmatrix + 1j * xEWmatrix

## apply a Hanning window to the data
    win = np.hanning(sampnum-1)
    for ii in range(numwin):
        xVmatrix[ii,:] = xVmatrix[ii,:]*win
        xHmatrix[ii,:] = xHmatrix[ii,:]*win

    # compute unfiltered magnitude response
    for iii in range(numwin):
        xVmatrix[iii,:] = np.abs(np.fft.fft(xVmatrix[iii,:]))/sampnum
        xHmatrix[iii,:] = np.abs(np.fft.fft(xHmatrix[iii,:]))/sampnum
    
    # compute the frequency axis
    N = sampnum
    faxbinsN = np.linspace(0,fs, N - 1)
    N_2 = math.ceil(N/2)
    fax_HzN = faxbinsN[0:int(N_2)]
    
    ## translate upbound and lowbound in from Hz to sample index
    lowbound = np.argmin(np.abs(fax_HzN - lowbound))
    upbound = np.argmax(np.abs(fax_HzN - lowbound))

    # now cut the magnitude responses in half
    
    xV = []
    xH = []
    xV = xVmatrix[0:N_2]
    xH = xHmatrix[0:N_2]
    xH = np.real(xH)

    ## do the HVSR
    H = xH/xV
    
    ## and compute the Thomspon statistics
    [ahatf, sigma, confinthigh, confintlow] = wavav(H)
    
    ## and cut to the half magnitude
    ahatf = ahatf[0:N_2]
    sigma = sigma[0:N_2]
    confinthigh = confinthigh[0:N_2]
    confintlow = confintlow[0:N_2]
    
    ## now plot the HVSR
    if Allplots in 'yes' or HVSRplottog in 'yes':
        HVSRplot(fax_HzN, ahatf, confinthigh, confintlow, statname, lowbound, upbound, outpath, sav)


