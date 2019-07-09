function HVSR(path, datapath, individ, varargin)
disp(nargin)
%function HVSR(path, datapath, individ)
%close all
%path = 'C:\Users\mpontr01\Box\2019_2_summer\Projects\HVSR';
%datapath = 'C:\Users\mpontr01\Box\2019_2_summer\Projects\Mexico City\Data';
%individ = 'yes';
warning('off','all') %The warnings are from the triangular filter which is 
%still a piece of the code, though it can be removed. 

%% AU11
 stationlist = {'AU11'};%,'AL01', 'AO24', 'AP68', 'AR14', 'AU11', 'AU46', 'BA49',...
%     'BL45', 'BO39', 'CA20', 'CA59', 'CB43', 'CC55',...
%     'CE18', 'CE23', 'CE32', 'CH84', 'CI05', 'CJ03', 'CJ04', 'CO47', 'CO56',...
%     'CP28', 'CS78', 'CT64', 'CU80', 'DM12', 'DR16', 'DX37',...
%     'EO30', 'ES57', 'EX08', 'EX09', 'EX12', 'FJ74', 'GA62', 'GC38', 'GR27',...
%     'HA41', 'HJ72', 'IB22', 'IM40', 'JA43', 'JC54', 'LI33', 'LI58', 'LV17',...
%     'ME52', 'MI15', 'MY19', 'NZ20', 'NZ31', 'PA34', 'PD42', 'PE10', 'RI76',...
%     'RM48', 'SI53', 'SP51', 'TE07', 'TH35', 'TL08',...
%     'TL55', 'UC44', 'VG09', 'VM29', 'XP06'};
for eee = 1:length(stationlist)
    station = stationlist{eee};
    disp(station)
    d = strcat(datapath,'\',station);
    %go into data directory and build structure of all files in it
    cd(d)
    %cd 'C:\Users\Marshall\Box Sync\tFall_2018\Research\Mexico_City\Data\AE02';
    files = dir;
    files = files(3:length(files));
    
    %create empty matrix that gets filled with all H/V values, counter is used
    %to index this matrix
    HV_final_matrix = [];
    peakfreq = [];
    peakamp = [];
    counter = 0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %function

    %change directory back to codes to access functions needed 
    %cd 'C:\Users\Marshall\Box Sync\tFall_2018\Research\Mexico_City\Codes';
    d = strcat(path, '\HVSR\Codes');
    cd(d)
    for file = files'
        filename = strcat(datapath,'\', station,'\',file.name);
        %filename = strcat('C:\Users\Marshall\Box Sync\tFall_2018\Research\Mexico_City\Data\AE02\',file.name);
        [xNS,xV,xEW, fs] = readfile1(filename);
        [PGANS,PGAV,PGAEW] = PGA(xNS,xV,xEW);
        if PGANS < 0.1 && PGAV < 0.1 && PGAEW < 0.1
            counter = counter + 1;
            %fs = station.fs; %sampling frequency in hz
            [xNS, xV, xEW] = Butter(xNS, xV, xEW, fs); %filter the data
            [N_2, fax_HzN, XH_magfilt,XV_magfilt] =  Magresp(xNS, xV, xEW, fs); %Compute mag responses and run through triangular filter
    
            %perform H/V
            [H_V1] = HV(XH_magfilt,XV_magfilt);

            %make Hz vector and linear interpolate all H/V ETFs to this vector
            newfaxhz = 0:0.001:20;
            newH_V1 = interp1(fax_HzN, H_V1, newfaxhz);
            HV_final_matrix(counter, :) = newH_V1; 
            clear H_V1
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

[ahatf, sigma, confinthigh, confintlow] = HVSRavg(HV_final_matrix);
HVSRplot(ahatf, newfaxhz, confinthigh, confintlow, station);

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
    individplot(HV_final_matrix, newfaxhz, station)
end
fclose('all')
end
%end