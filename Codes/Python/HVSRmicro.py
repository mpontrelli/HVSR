# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 18:04:24 2020

@author: mpontr01

1/16/2020 added everything beyond window the data and developed code wavav
"""

sav = 'no'
#filename = 'C:\\Users\\mpontr01\\Desktop\\Stations\\Tufts University\\2019_10_26\\Data\\bd_mseed_A00L_1572134400'
Vfname = 'C:\\Users\\mpontr01\\Desktop\\boston_site_response\\field_deployments\\Boston area\\908DP\\Data\\2000004171152005_908DP_1_1.sac'
NSfname = 'C:\\Users\\mpontr01\\Desktop\\boston_site_response\\field_deployments\\Boston area\\908DP\\Data\\2000004171152005_908DP_1_2.sac'
EWfname = 'C:\\Users\\mpontr01\\Desktop\\boston_site_response\\field_deployments\\Boston area\\908DP\\Data\\2000004171152005_908DP_1_3.sac'
outpath = 'C:\\Users\\mpontr01\\Desktop\\Stations\\Tufts University\\2019_10_26\\Figures'
#[xV, xNS, xEW] = mseed_timeseries(filename, outpath, sav)
xV = sac_timeseries(Vfname, outpath, sav)
xNS = sac_timeseries(NSfname, outpath, sav)
xEW = sac_timeseries(EWfname, outpath, sav)
numwin = 15
windowlen = 40
fs = 100
windis = 25
sampnum = windowlen*fs
windisnum = windis *fs
# Window the data
xVmatrix=np.zeros([numwin, windowlen*fs -1])
xNSmatrix=np.zeros([numwin, windowlen*fs -1])
xEWmatrix=np.zeros([numwin, windowlen*fs -1])
k = (1,fs)
for i in range(numwin):
    xVmatrix[i,:] = xV[(k[0]):(k[1] * windowlen + k[0])-1]
    xNSmatrix[i,:] = xNS[(k[0]):(k[1] * windowlen + k[0])-1]
    xEWmatrix[i,:] = xEW[(k[0]):(k[1] * windowlen + k[0])-1]
    k = (k[0] + sampnum + windisnum, fs)
xHmatrix = xNSmatrix + 1j * xEWmatrix

# window the data
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

xV = []
xH = []
# now cut the magnitude responses in half
xV = xVmatrix[0:N_2]
xH = xHmatrix[0:N_2]
xH = np.real(xH)


H = xH/xV
size1 = H.shape
len1 = size1[0]
q = np.log(H)
ahatf = np.exp(np.nansum(q, axis=0)/len1)
for i in range(len1):
    q[i,:] = (np.log(H[i,:]) - np.log(ahatf))**2
sigma = np.sqrt(np.nansum(q, axis = 0)/len1)
confinthigh = np.exp(np.log(ahatf)+1.96*sigma)
confintlow = np.exp(np.log(ahatf)-1.96*sigma)

ahatf = ahatf[0:N_2]
sigma = sigma[0:N_2]
confinthigh = confinthigh[0:N_2]
confintlow = confintlow[0:N_2]


plt.yscale("log")
plt.xscale("log")


plt.fill_between(fax_HzN, confintlow, confinthigh, color= [0.9, 0.9, 0.9])
plt.plot(fax_HzN, ahatf)
#plt.plot(fax_HzN, confintlow)

HV_1 = HV_1[0:N_2]
plt.figure
plt.yscale("log")
plt.xscale("log")
plt.plot(fax_HzN, HV_1)

