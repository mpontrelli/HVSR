# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 09:31:40 2019

@author: mpontr01
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

# INPUTS
    # datatype = either bd or ws, ground motion or weather respectively
    # station - the station you want to pull data from 
    # outpath - path where the mseed datafile will go
    # file that you want to extract, you have to go into the ftp directory to pick the one you want
    # pw - must be a string

# OUTPUTS
    # This writes a datafile from the D.O data ftp server to the outpath
def ftp_connect(datatype, station, outpath, file, pw):
    #outpath = 'C:/Users/mpontr01/Desktop/Stations/Tufts University/2019_10_26/Data'
    ftp = ftplib.FTP("ftp.aftac.gov")
    ftp.login()
    ftp.cwd('outgoing/TuftsUniversity/BostonSeismicStudy/' + datatype + '/' + station + '/')
    ftp.retrlines('LIST') 
    with open(file, "wb") as gFile:
        ftp.retrbinary('RETR'+ ' ' +  file, gFile.write)
    ZipFile(file).extractall(path = outpath, members = None, pwd=bytes(pw, 'utf-8'))

def mag_response(x, fs):
    N = len(x)
    X = np.fft.fft(x)
    X_mag = abs(X)/(2*N)
    fax_HzN1 = np.linspace(0, fs, N)
    N_2 = math.ceil(N//2) - fs
    fax_HzN = fax_HzN1[1:int(N_2)]
    X_mag = X_mag[1:int(N_2)]
    return X_mag, fax_HzN

def complex_time(xNS, xEW):
    xH = xNS + 1j*xEW
    return xH
    
    

# INPUTS
    # filename - name of ground motion file to be processed
    # outpath - folder into which the program writes the figures
    # sav - if toggled on, then save the plots
    
# OUTPUTS
    # plots of time series, magnitude response, HVSR, spectrograms of each component
def mseed_timeseries(filename, outpath, sav):
    st = read(filename, debug_headers=True)
    print(st)
    print(st[0].stats)
    xV = st[0]
    xNS = st[1]
    xEW = st[2]

    # detrend the data
    xV = xV-np.mean(xV)
    xNS = xNS-np.mean(xNS)
    xEW = xEW-np.mean(xEW)

    # extract values from the metadata
    # station
    station = st[0].stats.station
    # components
    vertical = st[0].stats.channel
    northsouth = st[1].stats.channel
    eastwest = st[2].stats.channel
    # fs
    fs = st[0].stats.sampling_rate
    # length of vector
    npts = st[0].stats.npts
    # delta t
    delta = st[0].stats.delta
    # endtime of record for title
    endtime = st[0].stats.endtime
    # make a time vector
    T = npts/fs
    t = np.linspace(0, T, npts, endpoint=False)/3600

    # Now design the filter
    lowcut = 0.5
    highcut = fs/2 - 1
    order = 4

    # you can plot the filter if you want
    # Plot the frequency response.
#    plt.figure(1)
#    plt.clf()
#    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
#    w, h = freqz(b, a, worN=2000)
#    plt.plot((fs * 0.5 / np.pi) * w, abs(h), label="order = %d" % order)
#    plt.xlabel('Frequency (Hz)')
#    plt.ylabel('Gain')
#    plt.grid(True)
#    plt.legend(loc='best')

    # now filter the data
    xV = butter_bandpass_filter(xV, lowcut, highcut, fs, order=order)
    xNS = butter_bandpass_filter(xNS, lowcut, highcut, fs, order=order)
    xEW = butter_bandpass_filter(xEW, lowcut, highcut, fs, order=order)

    # compute the maximum value, d is used to create y limits
    a = max(abs(xV))
    b = max(abs(xNS))
    c = max(abs(xEW))
    d = max(a,b,c)

    # now plot the data
    # Vertical

    plt.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0.2, hspace=0.2)
    plt.rcParams.update({'font.size': 12})
    plt.rcParams['font.family'] = "Times New Roman"
    fig = plt.figure(2)
    fig.set_size_inches(8, 8)
    ax1 = fig.add_subplot(3,1,1)
    ax1.plot(t, xV, label='Filtered signal')
    plt.ylabel('counts')
    plt.title(str(station) + ' ' + str(vertical) + ' ' + str(endtime.month) + '/' + str(endtime.day)+ '/'+ str(endtime.year))
    plt.grid(True)
    plt.ylim((-d, d))
    plt.xlim((0, 24))

    # North-south
    ax2 = fig.add_subplot(3,1,2)
    ax2.plot(t, xNS, label='Filtered signal')
    plt.ylabel('counts')
    plt.title(str(northsouth))
    plt.grid(True)
    plt.ylim((-d, d))
    plt.xlim((0, 24))

    # East-west
    ax3 = fig.add_subplot(3,1,3)
    ax3.plot(t, xEW, label='Filtered signal')
    plt.xlabel('time (hours)')
    plt.ylabel('counts')
    plt.title(str(eastwest))
    plt.grid(True)
    plt.ylim((-d, d))
    plt.xlim((0, 24))

    # and save
    if sav == 'yes':
        fig.savefig(outpath + str(station)+str(endtime.month)+str(endtime.day)+str(endtime.year)+ '.jpg', dpi=100)
    
    return xV, xNS, xEW
    # Now perform the HVSR
    
def sac_timeseries(filename, outpath, sav):
    st = read(filename, debug_headers=True)
    print(st)
    print(st[0].stats)
    x = st[0]


    # detrend the data
    x = x-np.mean(x)


    # extract values from the metadata
    # station
    station = st[0].stats.station
    # components
    channel = st[0].stats.channel

    # fs
    fs = st[0].stats.sampling_rate
    # length of vector
    npts = st[0].stats.npts
    # delta t
    delta = st[0].stats.delta
    # endtime of record for title
    endtime = st[0].stats.endtime
    # make a time vector
    T = npts/fs
    t = np.linspace(0, T, npts, endpoint=False)/60

    # Now design the filter
    lowcut = 0.5
    highcut = fs/2 - 1
    order = 4

    # you can plot the filter if you want
    # Plot the frequency response.
#    plt.figure(1)
#    plt.clf()
#    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
#    w, h = freqz(b, a, worN=2000)
#    plt.plot((fs * 0.5 / np.pi) * w, abs(h), label="order = %d" % order)
#    plt.xlabel('Frequency (Hz)')
#    plt.ylabel('Gain')
#    plt.grid(True)
#    plt.legend(loc='best')

    # now filter the data
    x = butter_bandpass_filter(x, lowcut, highcut, fs, order=order)


    # compute the maximum value, d is used to create y limits
    a = max(abs(x))


    # now plot the data
    # Vertical

    plt.rcParams.update({'font.size': 12})
    plt.rcParams['font.family'] = "Times New Roman"
    fig = plt.figure(2)
    fig.set_size_inches(8, 8)
  
    plt.plot(t, x, label='Filtered signal')
    plt.ylabel('counts')
    plt.title(str(station) + ' ' + str(channel) + ' ' + str(endtime.month) + '/' + str(endtime.day)+ '/'+ str(endtime.year))
    plt.grid(True)
    plt.ylim((-a, a))
    #plt.xlim((0, t(len[t]-1)))



    # and save
    if sav == 'yes':
        fig.savefig(outpath + str(station)+str(endtime.month)+str(endtime.day)+str(endtime.year)+ '.jpg', dpi=100)
    
    return x


def HVSR(xV, xH):
        HV_1 = xH/xV
        return HV_1
    
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

## FILTERS

## Butterworth filter
    # butterworth filter designed for processing microtremor data
    # documentation comes from https://scipy-cookbook.readthedocs.io/items/ButterworthBandpass.html
    
    # INPUTS
    
    # OUTPUTS
    
# Author: Marshall Pontrelli
# Date: 11/4/2019
def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a


def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = filtfilt(b, a, data) # I edited this from the recipe (it was lfilt)
    return y




def run():
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.signal import freqz

    # Sample rate and desired cutoff frequencies (in Hz).
    fs = 5000.0
    lowcut = 500.0
    highcut = 1250.0

    # Plot the frequency response for a few different orders.
    plt.figure(1)
    plt.clf()
    for order in [3, 6, 9]:
        b, a = butter_bandpass(lowcut, highcut, fs, order=order)
        w, h = freqz(b, a, worN=2000)
        plt.plot((fs * 0.5 / np.pi) * w, abs(h), label="order = %d" % order)

    plt.plot([0, 0.5 * fs], [np.sqrt(0.5), np.sqrt(0.5)],
             '--', label='sqrt(0.5)')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Gain')
    plt.grid(True)
    plt.legend(loc='best')

    # Filter a noisy signal.
    T = 0.05
    nsamples = T * fs
    t = np.linspace(0, T, nsamples, endpoint=False)
    a = 0.02
    f0 = 600.0
    x = 0.1 * np.sin(2 * np.pi * 1.2 * np.sqrt(t))
    x += 0.01 * np.cos(2 * np.pi * 312 * t + 0.1)
    x += a * np.cos(2 * np.pi * f0 * t + .11)
    x += 0.03 * np.cos(2 * np.pi * 2000 * t)
    plt.figure(2)
    plt.clf()
    plt.plot(t, x, label='Noisy signal')

    y = butter_bandpass_filter(x, lowcut, highcut, fs, order=6)
    plt.plot(t, y, label='Filtered signal (%g Hz)' % f0)
    plt.xlabel('time (seconds)')
    plt.hlines([-a, a], 0, T, linestyles='--')
    plt.grid(True)
    plt.axis('tight')
    plt.legend(loc='upper left')

    plt.show()


run()

# this is a basic moving average filter. Python may have a better one but I just
# wanted to copy the MATLAB smooth function https://www.mathworks.com/help/curvefit/smooth.html 
# window is in samples
# I pulled this from stackoverflow: https://stackoverflow.com/questions/40443020/matlabs-smooth-implementation-n-point-moving-average-in-numpy-python
def smooth(a,window):
    # a: NumPy 1-D array containing the data to be smoothed
    # WSZ: smoothing window size needs, which must be odd number,
    # as in the original MATLAB implementation
    out0 = np.convolve(a,np.ones(window,dtype=int),'valid')/window    
    r = np.arange(1,window-1,2)
    start = np.cumsum(a[:window-1])[::2]/r
    stop = (np.cumsum(a[:-window:-1])[::2]/r)[::-1]
    q = np.concatenate((  start , out0, stop  ))
    return q
## PLOTS

def timeseriesplot(xNS, xV, xEW, fs, sav, outpath):
    # find length of xV
    npts = len(xV)
    # find the maximum value for bounds for plotting
    a = np.max(np.abs(xV))
    b = np.max(np.abs(xNS))
    c = np.max(np.abs(xEW))
    d = np.array([a,b,c])
    d = np.max(d)
    
    # create a time vector and plot time series
    time = np.linspace(0, npts/fs, npts, endpoint=False)/60
    
    plt.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0.2, hspace=0.2)
    plt.rcParams.update({'font.size': 12})
    plt.rcParams['font.family'] = "Times New Roman"
    fig = plt.figure()
    fig.set_size_inches(8, 8)
    ax1 = fig.add_subplot(3,1,1)
    ax1.plot(time, xV, label='Time (mins)')
    plt.ylabel('counts')
    plt.title('V')
    plt.grid(True)
    plt.ylim((-d, d))
    plt.xlim((0, npts/(fs*60)))

    # North-south
    ax2 = fig.add_subplot(3,1,2)
    ax2.plot(time, xNS, label='Time (mins)')
    plt.ylabel('counts')
    plt.title('EW')
    plt.grid(True)
    plt.ylim((-d, d))
    plt.xlim((0, npts/(fs*60)))

    # East-west
    ax3 = fig.add_subplot(3,1,3)
    ax3.plot(time, xEW, label='Filtered signal')
    plt.xlabel('time (hours)')
    plt.ylabel('counts')
    plt.title('NS')
    plt.grid(True)
    plt.ylim((-d, d))
    plt.xlim((0, npts/(fs*60)))

    # and save
    if sav in 'yes':
        fig.savefig(outpath  + '\timeseries.jpg', dpi=100)
        
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

def individmagrespplot(fax_HzN, xH, xV, fs,N_2, lowbound, upbound, outpath, sav):
    fig3 = plt.figure()
    ax1 = fig3.add_subplot(2,1,1)
    plt.yscale("log")
    plt.xscale("log")
    plt.rcParams.update({'font.size': 12})
    plt.rcParams['font.family'] = "Times New Roman"
    ax1.plot(fax_HzN, xH[0:N_2, :].T, color = [0.8, 0.8, 0.8])
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Mag resp')
    plt.title('Horizontal')
    plt.grid(True)
    plt.ylim((0.1, 1000))
    plt.xlim((fax_HzN[lowbound], fax_HzN[upbound]))
    
    
<<<<<<< HEAD
    
=======
    fig4 = plt.figure()
    ax1 = fig4.add_subplot(2,1,2)
    plt.yscale("log")
    plt.xscale("log")
    plt.rcParams.update({'font.size': 12})
    plt.rcParams['font.family'] = "Times New Roman"
    ax1.plot(fax_HzN, xV[0:N_2, :].T, color = [0.8, 0.8, 0.8])
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Mag resp')
    plt.title('Vertical')
    plt.grid(True)
    plt.ylim((0.1, 1000))
    plt.xlim((fax_HzN[lowbound], fax_HzN[upbound]))
    
    ## save if sav is toggled on
    if sav in 'yes':
        fig.savefig(outpath  + '\HVSR.jpg', dpi=100)
        
def averagedmagrespplot(fax_HzN, ahatfhorz, ahatfvert, fs,confinthighhorz, confintlowhorz, confinthighvert, confintlowvert, lowbound, outpath, sav):    
    fig5 = plt.figure()
    ax1 = fig5.add_subplot(2,1,1)
    plt.yscale("log")
    plt.xscale("log")
    plt.rcParams.update({'font.size': 12})
    plt.rcParams['font.family'] = "Times New Roman"
    plt.fill_between(fax_HzN, confintlowhorz, confinthighhorz, color= [0.9, 0.9, 0.9])
    ax1.plot(fax_HzN, ahatfhorz, color = [0, 0.30196, 0.6588])
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Mag resp')
    plt.title('Horizontal')
    plt.grid(True)
    plt.ylim((0.1, 1000))
    plt.xlim((fax_HzN[lowbound], fax_HzN[upbound]))
    
    
    fig6 = plt.figure()
    ax1 = fig6.add_subplot(2,1,2)
    plt.yscale("log")
    plt.xscale("log")
    plt.rcParams.update({'font.size': 12})
    plt.rcParams['font.family'] = "Times New Roman"
    plt.fill_between(fax_HzN, confintlowvert, confinthighvert, color= [0.9, 0.9, 0.9])
    ax1.plot(fax_HzN, ahatfvert, color = [0, 0.30196, 0.6588])
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Mag resp')
    plt.title('Vertical')
    plt.grid(True)
    plt.ylim((0.1, 1000))
    plt.xlim((fax_HzN[lowbound], fax_HzN[upbound]))
    
    ## save if sav is toggled on
    if sav in 'yes':
        fig.savefig(outpath  + '\HVSR.jpg', dpi=100)
        
        
>>>>>>> b2fba09c062ebdb8dd3bfe5ff5b1662417c678bc
