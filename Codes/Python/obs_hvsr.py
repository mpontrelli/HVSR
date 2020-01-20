#!/usr/bin/env python3

#for testing: to load script in py3 command prompt: `exec(open('obs_hvsr.py').read())`

"""BEGIN IMPORTS"""
from obspy import read
from obspy import UTCDateTime
import numpy as np
import math
import scipy as sp
import matplotlib.pyplot as plt
"""END IMPORTS"""


filename = '/Users/jeremy/tufts_2019/research/waves/mseed/2018/258/tufts_trimmed.msd' #3600 seconds

#windows vars
window_len = 100
fs = 200
win_dis = 100
sampnum = window_len*fs
windisnum = win_dis *fs

#filter vars
lowpass=0.5
highpass=fs/2 - 1
order=4

#begin 
st = read(filename)

start = st[0].stats.starttime
end = st[0].stats.endtime
st_len = np.abs(end - start)


#detrend
st.detrend(type='demean')

#bandpass filter
st.filter('bandpass', freqmin=lowpass, freqmax=highpass, corners=order, zerophase=False)

#window
windows = []
tapered_windows = []
num_win = 10

for window in st.slide(window_length=window_len, step=win_dis):
    windows.append(window)

#print(windows)

for window in windows:
    tapered_windows.append(window.taper(type='hann', max_percentage=0.05)) #NOTE: must taper AFTER detrending/filtering


vert = []
for each in tapered_windows:
    vert.append(each[0].data)

horz = []
for each in tapered_windows:
    horz.append(each[1].data + 1j * each[2].data)


#freq axe [from Pontrelli]
N = len(horz[0])
faxbinsN = np.linspace(0, fs, N - 1)
N_2 = math.ceil(N/2)
fax_HzN = faxbinsN[0:int(N_2)]

vert_array = np.array(vert) #turn into np array, could do as `vert = np.array(vert)` but want to show steps w/ new vars
horz_array = np.array(horz)

vert_fft = np.abs(np.fft.fft(vert_array)) / sampnum #need to confirm why normalize by `sampnum`, simple fft tutorial may explain
horz_fft = np.abs(np.fft.fft(horz_array)) / sampnum 

h_v = horz_fft / vert_fft #divides horz np array by vert np array into new HVSR array, simply called `h_v`

#holder holds half of the values in each window, MP does this in `HVSR_micro.m`, review tutorial he gave last year, due to inherent FFT properties
holder = []
for i in h_v:
    half_h_v = np.vstack(i[0:N_2])
    holder.append(half_h_v) #w/out will just keep writing over each iteration, half of values, i.e. 1-100 Hz ? maybe 101 or 99 actually? in a LIST

fig,ax = plt.subplots()

plt.loglog(fax_HzN, holder[0])

"""
for i in holder:
    plt.loglog(fax_HzN, holder[i])
"""
plt.show()





"""
#mag resp
mV, mH = [], []
#script slow here
for i in vert:
    mV.append(np.abs(sp.fftpack.fft(i)) / sampnum)

for i in horz:
    mH.append(np.abs(sp.fftpack.fft(i)) / sampnum)

half_mV, half_mH = [], []

for each in mV:
    half_mV.append(mV[each][0:N_2])

for each in mH:
    half_mH.append(mH[each][0:N_2])


"""









