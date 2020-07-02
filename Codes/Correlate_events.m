%% Correlate events

% Author: Marshall Pontrelli
% Date: 5/15/2020

close all
clear all

%% CE32
filename = 'C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Processed_data2\NZ31\NZ3120170919181440';
CE32 = load(filename);
CE32ID = CE32.data.meta.event.ID;
CE32az = CE32.data.meta.event.azimuth;
CE32dist = CE32.data.meta.event.epi_dist;
CE32magrespNS = CE32.data.processing.filtereddata.acceleration.NS.mag_resps.smooth;
CE32wave = CE32.data.processing.filtereddata.acceleration.NS.waveform_orig;
CE32time = CE32.data.processing.filtereddata.time_orig;
% The next 3 lines normalize the waveform by the sqrt of its total energy
waveform1sq = CE32wave.^2;
waveform1norm = sqrt(sum(waveform1sq));
CE32wave = CE32wave/waveform1norm;

%% AU11
filename = 'C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Processed_data2\NZ20\NZ2020170919181440';

AU11 = load(filename);
AU11ID = AU11.data.meta.event.ID;
AU11az = AU11.data.meta.event.azimuth;
AU11dist = AU11.data.meta.event.epi_dist;
AU11magrespNS = AU11.data.processing.filtereddata.acceleration.NS.mag_resps.smooth;
AU11wave = AU11.data.processing.filtereddata.acceleration.NS.waveform_orig;
AU11time = AU11.data.processing.filtereddata.time_orig;
% The next 3 lines normalize the waveform by the sqrt of its total energy
waveform1sq = AU11wave.^2;
waveform1norm = sqrt(sum(waveform1sq));
AU11wave = AU11wave/waveform1norm;


%% Now correlate
[c12,lags12]=xcorr(CE32wave,AU11wave);
[maxc12,indc12]=max(c12);
%lagstime1=lags12(indc12)*delta;
figure
subplot(3,1,1)
plot(CE32time,CE32wave)
xlim([0 350])
ylim([-0.05, 0.05])
subplot(3,1,2)
plot(AU11time,AU11wave)
xlim([0 350])
ylim([-0.05, 0.05])
subplot(3,1,3)
plot(lags12/100,c12)