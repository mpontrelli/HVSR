%% Analyze Mexico City RACM station HVSR, SSR and TTF using correlation

% Author Marshall Pontrelli
% Date: 12/6/2019
close all
clear all

%% Inputs
statname = 'CE32';
%% load data
% HVSR
HVload = load(strcat('C:\Users\mpontr01\Box\2019_3_fall\AGU\Data\HVSRs\',statname, '.mat'));
HV = HVload.ahatfHV;
HVfreq = HVload.freq;

% SSR
SSRload = load(strcat('C:\Users\mpontr01\Box\2019_3_fall\AGU\Data\SSRs\',statname, '.mat'));
SSR = SSRload.ahatfSSR;
SSRfreq = SSRload.newfaxhz;

% TTF
TTFload = load(strcat('C:\Users\mpontr01\Box\2019_3_fall\AGU\Data\TTFs\',statname, '.mat'));
TTF = TTFload.amps;
TTFfreq = TTFload.freq;
%% now interpolate onto SSR 
HV = interp1(HVfreq, HV, SSRfreq);
TTF = interp1(TTFfreq, TTF, SSRfreq);

%% now load the waveform

[xNS,xV,xEW, fs] = readfile1('C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Data\TP13\TP1320170919181440');

[xNS, xV, xEW] = Butter(xNS, xV, xEW, fs); %filter the data
N = length(xNS) + length(TTF) - 1;
TTF = [TTF zeros(1,N-length(TTF))];
factor = length(TTF)/length(xEW);
timevec = 0:1/fs: (length(xNS)-1)/fs;
timevec2 = 0:1/(factor*fs): (length(xNS)-1)/fs;
xNS = interp1(timevec, xNS, timevec2);
figure
plot(timevec2, xNS)
TTF = TTF(1:length(xNS));
XNS = fft(xNS);


q = complex(TTF);
waveform = ifft(TTF.*XNS);

%% Now load waveform on the ground
[xNS,xV,xEW, fs] = readfile1('C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Data\CE32\CE3220170919181440');
[xNS, xV, xEW] = Butter(xNS, xV, xEW, fs); %filter the data
figure
factor = length(waveform)/length(xEW);
time = 1:1/(factor*fs):(length(xNS)-1)/fs;
waveform = waveform(1:length(time));
plot(time, (waveform))
hold on
time2 = 0:1/fs:(length(xNS) -1)/fs ;
xNS = interp1(time2, xNS, time);
plot(time, xNS)
