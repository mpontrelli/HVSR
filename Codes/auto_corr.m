%% auto_corr
% Auto_correlate ground motion data


%% Author: Marshall Pontrelli
% Date: 5/1/2020

%% Input filenames
Vfname = 'C:\Users\mpontr01\Box\People\Marshall and Jeremy\NEHRP\DATA\pratt_stations_zn\zn_net.492117\DC09.ZN.mseed';
EWfname = 'C:\Users\mpontr01\Box\People\Marshall and Jeremy\NEHRP\DATA\pratt_stations_zn\zn_net.492117\DC09.ZN.mseed';
NSfname = 'C:\Users\mpontr01\Box\People\Marshall and Jeremy\NEHRP\DATA\pratt_stations_zn\zn_net.492117\DC09.ZN.mseed';

Vfname2 = 'C:\Users\mpontr01\Box\People\Marshall and Jeremy\NEHRP\DATA\pratt_stations_zn\zn_net.492117\DC07.ZN.mseed';
EWfname2 = 'C:\Users\mpontr01\Box\People\Marshall and Jeremy\NEHRP\DATA\pratt_stations_zn\zn_net.492117\DC07.ZN.mseed';
NSfname2 = 'C:\Users\mpontr01\Box\People\Marshall and Jeremy\NEHRP\DATA\pratt_stations_zn\zn_net.492117\DC07.ZN.mseed';

%% read DMC data
Station = 'DC09';
Network = 'ZN';
Component = 'EHE';
time1 = '2015-06-07 06:07:00';
time2 = '2015-06-08 06:07:00';

[sampletimes,trace1,waveformfilt] = getDMCData(time1,time2,Station,Network,Component);

Station = 'DC30';
Network = 'ZN';
Component = 'EHE';
time1 = '2015-06-07 06:07:00';
time2 = '2015-06-08 06:07:00';

[sampletimes1,trace2,waveformfilt1] = getDMCData(time1,time2,Station,Network,Component);
%% necessary inputs
LowCorner = 3;
HighCorner = 5;
Npoles = 4;
fs = 100;
Filterplot = 'no';

%%
[waveformfilt] = Butter2(waveformfilt, fs, 'LowCorner', LowCorner, 'HighCorner', HighCorner, 'Npoles', Npoles , 'Filterplot', Filterplot);
[waveformfilt1] = Butter2(waveformfilt1, fs, 'LowCorner', LowCorner, 'HighCorner', HighCorner, 'Npoles', Npoles , 'Filterplot', Filterplot);

%% read files and convert into vectors
    
[~,~,ext] = fileparts(Vfname);
%.sacBinary
if strcmp(ext, '.sac') == 1 % use "rdmseed" if file is in miniseed format
[V] = ReadSacBinaryFile(Vfname); %vertical
[NS] = ReadSacBinaryFile(NSfname); %North-south
[EW] = ReadSacBinaryFile(EWfname); %East-West
    
%miniseed
elseif strcmp(ext, '.msd') == 1 % use "rdmseed" if file is in miniseed format
    x = rdmseed(NSfname);
    [NS, EW, V] = openmseed(x);
    
elseif strcmp(ext, '.mseed') == 1 % use "rdmseed" if file is in miniseed format
    x = rdmseed(NSfname);
    [NS, EW, V] = openmseed(x);
end

[~,~,ext] = fileparts(Vfname2);
%.sacBinary
if strcmp(ext, '.sac') == 1 % use "rdmseed" if file is in miniseed format
[V2] = ReadSacBinaryFile(Vfname2); %vertical
[NS2] = ReadSacBinaryFile(NSfname2); %North-south
[EW2] = ReadSacBinaryFile(EWfname2); %East-West
    
%miniseed
elseif strcmp(ext, '.msd') == 1 % use "rdmseed" if file is in miniseed format
    x = rdmseed(NSfname2);
    [NS2, EW2, V2] = openmseed(x);
    
elseif strcmp(ext, '.mseed') == 1 % use "rdmseed" if file is in miniseed format
    x = rdmseed(NSfname2);
    [NS2, EW2, V2] = openmseed(x);
end
%% Filter
[V] = Butter2(V, fs, 'LowCorner', LowCorner, 'HighCorner', HighCorner, 'Npoles', Npoles , 'Filterplot', Filterplot);
Filterplot = 'no'; % toggle off filter plot so it doesn't plot response three times
[NS] = Butter2(NS, fs, 'LowCorner', LowCorner, 'HighCorner', HighCorner, 'Npoles', Npoles , 'Filterplot', Filterplot);
[EW] = Butter2(EW, fs, 'LowCorner', LowCorner, 'HighCorner', HighCorner, 'Npoles', Npoles , 'Filterplot', Filterplot);


[V2] = Butter2(V2, fs, 'LowCorner', LowCorner, 'HighCorner', HighCorner, 'Npoles', Npoles , 'Filterplot', Filterplot);
Filterplot = 'no'; % toggle off filter plot so it doesn't plot response three times
[NS2] = Butter2(NS2, fs, 'LowCorner', LowCorner, 'HighCorner', HighCorner, 'Npoles', Npoles , 'Filterplot', Filterplot);
[EW2] = Butter2(EW2, fs, 'LowCorner', LowCorner, 'HighCorner', HighCorner, 'Npoles', Npoles , 'Filterplot', Filterplot);

%% Now auto correlate
% NS = NS(1:end-1);
% EW =EW(1:end-1);
% V = V(1:end-1);
%waveformfilt = waveformfilt(1:end-1);
[c,lags] = xcorr(waveformfilt, waveformfilt1,'normalized');

%%
figure
plot(lags/(60*fs),c)