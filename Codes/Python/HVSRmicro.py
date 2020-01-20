# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 18:04:24 2020

@author: mpontr01

1/16/2020 added everything beyond window the data and developed code wavav
"""
    

## inputs

ftype = '.sac' # filetype, this is to support a bunch of filtypes
sav = 'no' # do you want to save your figures?
outpath = 'C:\\Users\\mpontr01\\Desktop\\Stations\\Tufts University\\2019_10_26\\Figures' # if you want to save your figures use outpath

## filenames, these are for .sac
Vfname = 'B:\\Erkan Yilar\\Ambient Noise Data\\Files\\Ambient Noise Data Set\\08.15.2014\\Dan\\Danehy Park Cambridge\\Ground\\Original Files\\2014227134758005_T4260_1_1.sac'
NSfname = 'B:\\Erkan Yilar\\Ambient Noise Data\\Files\\Ambient Noise Data Set\\08.15.2014\\Dan\\Danehy Park Cambridge\\Ground\\Original Files\\2014227134758005_T4260_1_2.sac'
EWfname = 'B:\\Erkan Yilar\\Ambient Noise Data\\Files\\Ambient Noise Data Set\\08.15.2014\\Dan\\Danehy Park Cambridge\\Ground\\Original Files\\2014227134758005_T4260_1_3.sac'

# filename, this is a test for miniseed
filename = 'C:\\Users\\mpontr01\\Box\\2019_2_summer\\Projects\\Boston\\Data\\bd\\Waltham\\bd_data_A23L_1560124800\\bd_mseed_A23L_1560124800'

Allplots = 'no'
Timeplot = 'no'
HVSRplottog = 'yes'
IUMagplot = 'no'
AUMagplot = 'no'
IFMagplot = 'no'
AFMagplot = 'no'

## defaults for windowing the time series
numwin = 10
windowlen = 40
fs = 100
windis = 25
sampnum = windowlen*fs
windisnum = windis *fs
width = 0.5

## Filter inputs
lowcut = 0.5
highcut = fs/2 - 1
order = 4

## bounds for plotting and statistics calculation
upbound = 49
lowbound = 0.5


if ftype in '.sac':
    
    ## import the data
    stv = read(Vfname, debug_headers=True)
    xV = stv[0]
    stNS = read(NSfname, debug_headers=True)
    xNS = stNS[0]
    stEW = read(EWfname, debug_headers=True)
    xEW = stEW[0]
    
if ftype in '.msd':
    stv = read(filename, debug_headers=True)
    xV = stv[0]
    xNS = stv[1]
    xEW = stv[2]
    
## now detrend the data

xV = xV-np.mean(xV)
xNS = xNS-np.mean(xNS)
xEW = xEW-np.mean(xEW)

## extract values from the metadata, we only need to do this for 1 channel
statname = stv[0].stats.station
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
    xVmatrix[iii,:] = 4*np.abs(np.fft.fft(xVmatrix[iii,:]))/sampnum
    xHmatrix[iii,:] = 4*np.abs(np.fft.fft(xHmatrix[iii,:]))/sampnum
    
# compute the frequency axis
N = sampnum
faxbinsN = np.linspace(0,fs, N - 1)
N_2 = math.ceil(N/2)
fax_HzN = faxbinsN[0:int(N_2)]
xV = []
xH = []
xV = xVmatrix[:,0:N_2]
xH = xHmatrix[:,0:N_2]
xH = np.real(xH)

## translate upbound and lowbound in from Hz to sample index
lowbound = np.argmin(np.abs(fax_HzN - lowbound))
upbound = np.argmin(np.abs(fax_HzN - upbound))

# plot individual unfiltered magnitude responses 
if Allplots in 'yes' or IUMagplot in 'yes':
    individmagrespplot(fax_HzN, xH, xV, fs, N_2, lowbound, upbound, outpath, sav)  

# Average the un-smoothed magnitude responses
[ahatfhorz, sigmahorz, confinthighhorz, confintlowhorz] = wavav(xH)
[ahatfvert, sigmavert, confinthighvert, confintlowvert] = wavav(xV)  

# plot averaged unfiltered magnitude responses
if Allplots in 'yes' or AUMagplot in 'yes':
    averagedmagrespplot(fax_HzN, ahatfhorz, ahatfvert, fs,confinthighhorz, confintlowhorz, confinthighvert, confintlowvert, lowbound, outpath, sav)

# compute smoothed magnitude responses
window = int(np.ceil((N/fs)*width))
# check to see if even becasue must be odd
if window % 2 == 0:
    window = window -1
else:
    window = window
XVsmooth = np.zeros([numwin, N_2])
XHsmooth = np.zeros([numwin, N_2])
for i in range(numwin):
    XVsmooth[i, :] = smooth(xV[i,:],window)
    XHsmooth[i,:] = smooth(xH[i,:],window)

# plot individual, smoothed magnitude responses
if Allplots in 'yes' or IFMagplot in 'yes':
    individmagrespplot(fax_HzN, XHsmooth, XVsmooth, fs, N_2, lowbound, upbound, outpath, sav)  

# average the smoothed magnitude responses
[ahatfhorz, sigmahorz, confinthighhorz, confintlowhorz] = wavav(XHsmooth)
[ahatfvert, sigmavert, confinthighvert, confintlowvert] = wavav(XVsmooth)     

# plot averaged, smoothed magnitude responses
if Allplots in 'yes' or AFMagplot in 'yes':
    averagedmagrespplot(fax_HzN, ahatfhorz, ahatfvert, fs,confinthighhorz, confintlowhorz, confinthighvert, confintlowvert, lowbound, outpath, sav)

## do the HVSR
H = XHsmooth/XVsmooth
    
## and compute the Thomspon statistics
[ahatf, sigma, confinthigh, confintlow] = wavav(H)
    

    
## now plot the HVSR
if Allplots in 'yes' or HVSRplottog in 'yes':
    HVSRplot(fax_HzN, ahatf, confinthigh, confintlow, statname, lowbound, upbound, outpath, sav)


