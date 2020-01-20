# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 18:04:24 2020

@author: mpontr01

HISTORY...
MP 1/16/2020: added everything beyond window the data and developed code wavav 
JS 1/18/2020: created working copy; added argparse for command line input 
"""

from scipy.signal import butter, filtfilt
from zipfile import ZipFile
from scipy import signal
import ftplib
from obspy import read
from obspy import UTCDateTime
from numpy import hanning
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import freqz
import math



## inputs

ftype = '.msd' # filetype, this is to support a bunch of filtypes
sav = 'yes' # do you want to save your figures?
outpath = '/Users/jeremy/Desktop'

## filenames, these are for .sac
#Vfname = 'C:\\Users\\mpontr01\\Desktop\\boston_site_response\\field_deployments\\Boston area\\908DP\\Data\\2000004171152005_908DP_1_1.sac'
#NSfname = 'C:\\Users\\mpontr01\\Desktop\\boston_site_response\\field_deployments\\Boston area\\908DP\\Data\\2000004171152005_908DP_1_2.sac'
#EWfname = 'C:\\Users\\mpontr01\\Desktop\\boston_site_response\\field_deployments\\Boston area\\908DP\\Data\\2000004171152005_908DP_1_3.sac'

# filename, this is a test for miniseed
filename = '/Users/jeremy/tufts_2019/research/waves/mseed/2018/258/bd_A23L_1536969600.msd'

Allplots = 'yes'
Timeplot = 'no'
HVSRplottog = 'no'
## defaults for windowing the time series
numwin = 100
windowlen = 40
fs = 200
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

'''
if ftype in '.sac':
    
    ## import the data
    stv = read(Vfname, debug_headers=True)
    xV = stv[0]
    stNS = read(NSfname, debug_headers=True)
    xNS = stNS[0]
    stEW = read(EWfname, debug_headers=True)
    xEW = stEW[0]
'''   
    
if ftype in '.msd':
    stv = read(filename, debug_headers=True)
    xV = stv[0]
    xNS = stv[1]
    xEW = stv[2]
    
## now detrend the data

xV = xV-np.mean(xV) #NOTE: These obspy traces are automatically turned into numpy arrays after this operation; CONSIDER DETRENDING AS OBSPY TRACE FOR PERFORMANCE, THEN `st[0].data` FOR NUMPY ARRAY
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
t = np.linspace(0, T, npts, endpoint=False)/60 #JS to MP: why divide by 60??

"""BEGIN FILTERING"""
def butter_bandpass(lowcut, highcut, fs, order=5): #JS to MP: why are there two butter bandpass functions? Can we consolidate?
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = filtfilt(b, a, data) # I edited this from the recipe (it was lfilt)
    return y


## Now filter the data
xV = butter_bandpass_filter(xV, lowcut, highcut, fs, order=order)
xNS = butter_bandpass_filter(xNS, lowcut, highcut, fs, order=order)
xEW = butter_bandpass_filter(xEW, lowcut, highcut, fs, order=order)

 
#plt.plot(t, xV, label='Filtered signal')

#plt.show() 

"""END FILTERING"""



"""BEGIN WINDOWING"""
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

print(len(H))
print(H)

def wavav(H):
    size1 = H.shape
    len1 = size1[0]
    q = np.log(H)
    ahatf = np.exp(np.nansum(q, axis=0)/len1)
    for i in range(len1):
        q[i,:] = (np.log(H[i,:]) - np.log(ahatf))**2
    sigma = np.sqrt(np.nansum(q, axis = 0)/len1)
    confinthigh = np.exp(np.log(ahatf)+1.96*sigma)
    confintlow = np.exp(np.log(ahatf)-1.96*sigma)
    return ahatf, sigma, confinthigh, confintlow

def HVSRplot(fax_HzN, ahatf, confinthigh, confintlow, statname, lowbound, upbound, outpath, sav):
    fig2 = plt.figure()
    plt.yscale("log")
    plt.xscale("log")
    plt.rcParams.update({'font.size': 12})
    plt.rcParams['font.family'] = "Times New Roman"
    plt.fill_between(fax_HzN, confintlow, confinthigh, color= [0.9, 0.9, 0.9])
    plt.plot(fax_HzN, ahatf, color = [0, 0.30196, 0.6588], linewidth = 2)
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Amplification')
    plt.title(statname)
    plt.grid(True)
    plt.ylim((0.1, 100))
    plt.xlim((fax_HzN[lowbound], fax_HzN[upbound]))
   
    ## save if sav is toggled on
    if sav in 'yes':
        fig.savefig(outpath  + '\HVSR.jpg', dpi=100)

wavav(H)
HVSRplot()

"""END WINDOWING"""



"""    
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
"""
