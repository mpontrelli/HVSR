%function HVSR(path, datapath, individ, varargin)
close all
clear all
path = 'C:\Users\mpontr01\Box\2019_2_summer\Projects\HVSR';
datapath = 'C:\Users\mpontr01\Box\2019_2_summer\Projects\Mexico City\Data';
warning('off','all') %The warnings are from the triangular filter which is 
cd(datapath)
stationlist = dir;
stationlist = stationlist(3:length(stationlist));
for eee = 1:length(stationlist)
    station = stationlist(eee);
    statname = station.name;
    station = strcat(station.folder, '\', statname);
    %go into data directory and build structure of all files in it
    cd(station)
    %cd 'C:\Users\Marshall\Box Sync\tFall_2018\Research\Mexico_City\Data\AE02';
    files = dir;
    files = files(3:length(files));
    %create empty matrix that gets filled with all H/V values, counter is used
    %to index this matrix
    XH_final_matrix = [];
    XV_final_matrix = [];
    peakfreq = [];
    peakamp = [];
    counter = 0;
    %change directory back to codes to access functions needed 
    %cd 'C:\Users\Marshall\Box Sync\tFall_2018\Research\Mexico_City\Codes';
    d = strcat(path, '\HVSR\Codes');
    cd(d)
    for file = files'
        filename = strcat(station,'\',file.name);
        [xNS,xV,xEW, fs] = readfile1(filename);
        [PGANS,PGAV,PGAEW] = PGA(xNS,xV,xEW);
        if PGANS < 0.1 && PGAV < 0.1 && PGAEW < 0.1
            counter = counter + 1;
            %fs = station.fs; %sampling frequency in hz
            [xNS, xV, xEW] = Butter(xNS, xV, xEW, fs); %filter the data
            [N_2, fax_HzN, XH_magfilt,XV_magfilt] =  Magresp(xNS, xV, xEW, fs); %Compute mag responses and run through triangular filter
            
            %perform H/V
            %[H_V1] = HV(XH_magfilt,XV_magfilt);
            %make Hz vector and linear interpolate all H/V ETFs to this vector
            newfaxhz = 0:0.001:20;
            newXH = interp1(fax_HzN, XH_magfilt, newfaxhz);
            newXV = interp1(fax_HzN, XV_magfilt, newfaxhz);
            XH_final_matrix(counter, :) = newXH; 
            XV_final_matrix(counter, :) = newXV; 
            else
            [xNS, xV, xEW] = Butter(xNS, xV, xEW, fs); %filter the data
            [N_2, fax_HzN, XH_magfilt,XV_magfilt] =  Magresp(xNS, xV, xEW, fs); %Compute mag responses and run through triangular filter
            %perform H/V
            [H_V1] = HV(XH_magfilt,XV_magfilt);
            %make Hz vector and linear interpolate all H/V ETFs to this vector
            newfaxhz = 0:0.001:20;
            newH_V1 = interp1(fax_HzN, H_V1, newfaxhz);
            clear newfaxhz
            clear newH_V1
            clear H_V1 
        end
    end
newfaxhz = 0:0.001:20;  
%statistics per Thompson et al 2012 page 34
%compute maximum likelihood estimator of median
[magH, sigmaH, confinthighH, confintlowH] = HVSRavg(XH_final_matrix);
[magV, sigmaV, confinthighV, confintlowV] = HVSRavg(XV_final_matrix);
ahatf = HV(magH, magV);
confinthigh = HV(confinthighH,confinthighV);
confintlow = HV(confintlowH, confintlowV);

%plot
HVSRplot(ahatf, newfaxhz, confinthigh, confintlow, statname);

N = length(ahatf); %length North_South_Component
width = .1; %width for triangle moving average filter in hz
q = ceil((N/20)*width); %width for triangle moving average filter in samples
e = smooth(ahatf, q, 'moving');
%% Determine if peak is a peak
count = 0;
[amps,amplocs] = findpeaks(e);
for hh = 1:length(amps)
    freqpeak(hh) = newfaxhz(amplocs(hh));
end
[trough,troughloc] = findpeaks(-1*e);
for hh = 1:length(trough)
    freqtrough(hh) = newfaxhz(troughloc(hh));
end
if freqtrough(1) > freqpeak(1)
    Y = 1;
    X = 1;
    trough = horzcat(Y,trough);
    troughloc = horzcat(X,troughloc);
end
if freqpeak(length(freqpeak)) > freqtrough(length(freqtrough))
    Y = -1;
    X = 19991;
    trough = horzcat(trough', Y);
    troughloc = horzcat(troughloc',X);
end
for hh = 1:length(trough)
    freqtrough(hh) = newfaxhz(troughloc(hh));
end
for k = 1:length(amps)
    height = amps(k)/sqrt(2);
    lefttrough = -1 * trough(k);
    righttrough = -1 * trough(k + 1);
    if lefttrough < height && righttrough < height
        count  = count + 1;
        peakfreq(count) = newfaxhz(amplocs(k));
        peakamp(count) = amps(k); 
        amplocs2(count) = amplocs(k);
    end
end
%points = plot(peakfreq,peakamp,'o','MarkerSize',15,'MarkerEdgeColor','g','MarkerFaceColor','g');
%% Compute desired statistics
Taxstat = [];
for f = 1:length(peakamp)
    taxstat(f,1) = f; 
    taxstat(f,2) = peakfreq(f); 
    A = peakamp(f);
    taxstat(f,3) = A;
    amploc2 = amplocs2(f);
    [I1, I2, f1, f2, hpb] =  HalfPowerBand2(A, amploc2, newfaxhz, ahatf); 
    taxstat(f,4) = hpb;
    a = sigma(I1:I2);
    sigmai = median(a);
    taxstat(f,5) = sigmai;
end
newfaxhz = 0:0.001:20;
if nargin == 2
    continue
elseif strcmp(individ, 'yes') == 1
    individplot(HV_final_matrix, newfaxhz, statname)
end
fclose('all')
end
%end