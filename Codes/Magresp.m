function [N_2, fax_HzN, XH_magfilt, XV_magfilt]=  Magresp(xNS, xV, xEW, fs)

%Triangular filter
N = length(xNS); %length North_South_Component
width = .01; %width for triangle moving average filter in hz
q = ceil((N/fs)*width); %width for triangle moving average filter in samples

w = hamming(q);
a=5; %second filter coefficient

xH = xNS + 1i*xEW; %complex time

XH = fft(xH); %fft North_South_Component x
XV = fft(xV); %fft Vertical_Component x


XH_mag = abs(XH)/N; %normalized magnitude spectra North_South_Component x
XV_mag = abs(XV)/N; %normalized magnitude spectra Vertical_Component x
%YNS_mag = abs(YNS); %magnitude spectra North_South_Component y
%YV_mag = abs(YV); %magnitude spectra Vertical_Component y
%YEW_mag = abs(YEW); %magnitude spectra East_West_Component y

XV_magfilt1=filtfilt(w,a,XV_mag); %filtered magnitude spectra North_South_Component x 
XH_magfilt1=filtfilt(w,a,XH_mag); %filtered magnitude spectra North_South_Component x 

fax_binsN = [0 : N-1]; %samples in NS component
fax_HzN1 = fax_binsN*fs/N; %frequency axis NS (Hz)

if fs == 100
    ff = 5;
end
if fs == 200
    ff = 10;
end
N_2 = ceil(N/ff); %half magnitude spectrum
fax_HzN = fax_HzN1(1 : N_2);
% XH_magfilt1 = kohmachi(XH_mag,fax_HzN1,30);
% XV_magfilt1 = kohmachi(XV_mag,fax_HzN1,30);
XH_magfilt = XH_magfilt1(1: N_2);
XV_magfilt = XV_magfilt1(1: N_2);



% XH_magfilt = XH_mag(1: N_2);
% XV_magfilt = XV_mag(1: N_2);
% XH_magfilt = kohmachi(XH_magfilt,fax_HzN,30);
% XV_magfilt = kohmachi(XV_magfilt,fax_HzN,30);
end