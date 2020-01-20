#!/usr/bin/env python3

#for testing: to load script in py3 command prompt: `exec(open('obs_hvsr.py').read())`
#NOTE: PROVIDE CUT 3-COMPONENT WAVEFORMS ONLY (WORKING ON SINGLE CHANNEL VERSION)


"""BEGIN IMPORTS"""
from obspy import read
from obspy import UTCDateTime
import numpy as np
import math
import scipy as sp
import matplotlib.pyplot as plt
import argparse
import sys
"""END IMPORTS"""

def main():
    
    """BEGIN ARGUMENT PARSER"""
    parser = argparse.ArgumentParser(
    	formatter_class=argparse.RawDescriptionHelpFormatter,
    	description="Calculate HVSR using 3-Component Seismic Data\nPontrelli and Salerno, 2020\n\nRun `python obs_hvsr.py --help`")
    
    parser.add_argument("-f", "--file", type=str,
        required=True, help="Input file path, format: `/path/to/file.msd`")
    
    parser.add_argument("-sr","--sample_rate", type=int, 
        required=False, default=200, help="Sample Rate in Hz, ex: `100`, default: 200")
    
    parser.add_argument("-wl", "--window_length", type=int,
        required=False, default=40, help="Window Length in Sec, ex: `50`, default: 40")

    parser.add_argument("-wd", "--window_distance", type=int,
        required=False, default=25, help="Window Distance in Sec, ex: `100`, default: 25")
    
    parser.add_argument("-v", "--verbose", action="count",
        default=0, help="increase spewage")
    
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)   # <--If script is executed w/out args, it will auto show help
        sys.exit(1) # ...and exit out without running anything

    args = parser.parse_args()
    fyle = args.file
    window_len = args.window_length
    fs = args.sample_rate
    win_dis = args.window_distance
    """END ARGUMENT PARSER"""
     
    #windows vars
    sampnum = window_len*fs
    windisnum = win_dis *fs

    #filter vars
    lowpass=0.5
    highpass=fs/2 - 1
    order=4

    """BEGIN READ, DETREND, FILTER""" 
    st = read(fyle)

    start = st[0].stats.starttime
    end = st[0].stats.endtime
    st_len = np.abs(end - start)


    #detrend
    st.detrend(type='demean')

    #bandpass filter
    st.filter('bandpass', freqmin=lowpass, freqmax=highpass, corners=order, zerophase=False)
    """END READ, DETREND, FILTER"""

    """BEGIN WINDOWING"""
    windows = []
    tapered_windows = []
    num_win = 10

    
    for window in st.slide(window_length=window_len, step=win_dis): #window_length = seconds of window, win_dis = offset of 2nd win from start of 1st win!!
        windows.append(window) #so, these are overlapping if window_len = 40 and win_dis = 25, UNLIKE SESAME, compare w/ MP
    "END WINDOWING"""
  
    """BEGIN TAPERING OF WINDOWS"""
    for window in windows:
        tapered_windows.append(window.taper(type='hann', max_percentage=0.05)) #NOTE: must taper AFTER detrending/filtering
    """END TAPERING OF WINDOWS"""

    """BEGIN EXTRACT VERTICAL AND HORIZONTAL COMPONENTS FROM THE WINDOWS"""
    vert = []
    for each in tapered_windows:
        vert.append(each[0].data)

    horz = []
    for each in tapered_windows:
        horz.append(each[1].data + 1j * each[2].data)
    """END EXTRACT VERTICAL AND HORIZONTAL COMPONENTS FROM THE WINDOWS"""

    """BEGIN COMPUTE FREQ AXIS"""
    #freq axe [from MP HVSR_micro.py]
    N = len(horz[0])
    faxbinsN = np.linspace(0, fs, N - 1)
    N_2 = math.ceil(N/2)
    fax_HzN = faxbinsN[0:int(N_2)]
    """END COMPUTE FREQ AXIS"""
    
    """BEGIN COMPUTE FFT AND H/V"""
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
    """END COMPUTE FFT AND H/V"""

    """BEGIN HALF-ASSED PLOTTING"""
    fig,ax = plt.subplots()

    plt.loglog(fax_HzN, holder[0])
    plt.show()

    plt.show()
    """END HALF-ASSED PLOTTING"""

if __name__ == '__main__':
    main()
