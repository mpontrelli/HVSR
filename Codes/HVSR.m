function [newfaxhz, ahatf, sigma, lowbound, taxstat, statsfinal, fsmin, recmax, varargout] = HVSR(path, datapath, varargin)
%create Input Parser object
p = inputParser;
%add inputs to the scheme
defaultxbound = ([0, 20]);
defaultybound = ([0.0001, 100]);
addRequired(p,'path',@ischar);
addRequired(p,'datapath',@ischar);
addParameter(p, 'wavecut', 'No', @ischar);
addParameter(p, 'HVSR', 'No', @ischar);
addParameter(p, 'individHVSR', 'No', @ischar);
addParameter(p, 'magresps', 'No', @ischar);
addParameter(p, 'magxbounds', defaultxbound, @isnumeric);
addParameter(p, 'magybounds', defaultybound, @isnumeric);
%parse the inputs
parse(p, path, datapath, varargin{:})
%set varibales from the parse
wavecut = p.Results.wavecut;
HVSR = p.Results.HVSR;
individHVSR = p.Results.individHVSR;
magresps = p.Results.magresps;
magxbounds = p.Results.magxbounds;
magybounds = p.Results.magybounds;
%start function
cd(datapath)
stationlist = dir;
stationlist = stationlist(3:length(stationlist));
statsfinal = [];
for eee = 1:length(stationlist)
    station = stationlist(eee);
    statname = station.name;
    station = strcat(datapath, '\', statname);
    disp(station)
    %go into data directory and build structure of all files in it
    cd(station)
    files = dir;
    files = files(3:length(files));
    %create empty matrix that gets filled with all H/V values, counter is used
    %to index this matrix
    HV_final_matrix = [];
    XH_final_matrix = [];
    XV_final_matrix = [];
    lowbound_matrix = [];
    counter = 0;
    %change directory back to codes to access functions needed 
    d = strcat(path, '\HVSR\Codes');
    cd(d)
    %Loop through the records once to get the number of samples in each
    %record for the x-axis linear interpolation and the sampling freqeuncy 
    %of each record. The number of samples value will be the
    %sample number of the record with the maximum number of samples. The
    %sampling freqeuncy is required for the freqeuncy axis high bound
    [fsmin, recmax, recmin] = statrecinfo(files, station);
    for file = files'
        filename = strcat(station,'\',file.name);
        [xNS,xV,xEW, fs] = readfile1(filename);
        [PGANS,PGAV,PGAEW] = PGA(xNS,xV,xEW);
        if PGANS < 0.1 && PGAV < 0.1 && PGAEW < 0.1
            if strcmp(wavecut, 'yes') == 1   
                [xNS,xV,xEW] = waveformcut(xNS,xV,xEW, recmin);
                counter = counter + 1;
                %fs = station.fs; %sampling frequency in hz
                [xNS, xV, xEW] = Butter(xNS, xV, xEW, fs); %filter the data
                [N_2, fax_HzN, XH_magfilt, XV_magfilt, XH_mag, XV_mag, lowbound] =  Magresp(xNS, xV, xEW, fs, fsmin); %Compute mag responses and run through triangular filter
                [H_V1] = HV(XH_magfilt, XV_magfilt);
                XH_final_matrix(counter, :) = XH_mag;
                XV_final_matrix(counter, :) = XV_mag; 
                HV_final_matrix(counter, :) = H_V1; 
                clear H_V1
            else
            counter = counter + 1;
            %fs = station.fs; %sampling frequency in hz
            [xNS, xV, xEW] = Butter(xNS, xV, xEW, fs); %filter the data
            [N_2, fax_HzN, XH_magfilt, XV_magfilt, XH_mag, XV_mag, lowbound] =  Magresp(xNS, xV, xEW, fs, fsmin); %Compute mag responses and run through triangular filter
            %perform H/V
            [H_V1] = HV(XH_magfilt, XV_magfilt);
            %make Hz vector and linear interpolate all H/V ETFs to this vector
            newfaxhz = 0: (1/ (ceil(recmax/2) - fsmin))*(fsmin/2 - 1): (fsmin/2 - 1);
            [~, lowindex] = min(abs(newfaxhz - lowbound));
            lowbound_matrix(counter, :) = lowindex;
            %mag resp matrix build
            newXH_mag = interp1(fax_HzN, XH_mag, newfaxhz);
            XH_final_matrix(counter, :) = newXH_mag; 
            newXV_mag = interp1(fax_HzN, XV_mag, newfaxhz);
            XV_final_matrix(counter, :) = newXV_mag; 
            newH_V1 = interp1(fax_HzN, H_V1, newfaxhz);
            HV_final_matrix(counter, :) = newH_V1; 
            clear H_V1
            end
            
%             [xNS, xV, xEW] = Butter(xNS, xV, xEW, fs); %filter the data
%             [N_2, fax_HzN, XH_magfilt, XV_magfilt, XH_mag, XV_mag, lowbound] =  Magresp(xNS, xV, xEW, fs, fsmin); %Compute mag responses and run through triangular filter
%             %perform H/V
%             [H_V1] = HV(XH_magfilt,XV_magfilt);
%             %make Hz vector and linear interpolate all H/V ETFs to this vector
%             newfaxhz = 0: (1/ (ceil(recmax/2) - fsmin))*(fsmin/2 - 1): (fsmin/2 - 1);
%             [~, lowindex] = min(abs(newfaxhz - lowbound));
%             newH_V1 = interp1(fax_HzN, H_V1, newfaxhz); 
%             clear newfaxhz
%             clear newH_V1
%             clear H_V1 
        end
    end
if strcmp(wavecut, 'yes') == 1  
    newfaxhz = fax_HzN;
    [~, lowbound] = min(abs(newfaxhz - lowbound));
else
    newfaxhz = 0: (1/ (ceil(recmax/2) - fsmin))*(fsmin/2 - 1) : (fsmin/2 - 1);
    lowbound = max(lowbound_matrix);
end
% lowbound comes from lowbound _matrix and corresponds to the lowest
% frequency that can be resolved at the shortest time series record in the
% station database
%HVSR
if strcmp(HVSR, 'yes') == 1   
    [ahatf, sigma, confinthigh, confintlow] = wavav(HV_final_matrix);
    HVSRplot(ahatf, newfaxhz, confinthigh, confintlow, lowbound, statname);  
    [matrix, matrix1, peakind,ahatf1,newfaxhz1] = peakiden(ahatf, newfaxhz, lowbound, fsmin);
    [taxstat] = specratstat(peakind, matrix, matrix1, ahatf1, newfaxhz1, sigma, statname,lowbound);
    statsfinal = vertcat(statsfinal, taxstat);
end

if strcmp(individHVSR, 'yes') == 1   
    individplot(HV_final_matrix, newfaxhz, statname);
end
%Mag responses
if strcmp(magresps, 'yes') == 1     
    [ahatfXH, sigma, confinthighXH, confintlowXH] = wavav(XH_final_matrix);
    [ahatfXV, sigma, confinthighXV, confintlowXV] = wavav(XV_final_matrix);
    magrespplot(ahatfXH, newfaxhz, confinthighXH, confintlowXH, statname, magxbounds, magybounds, ' horizontal', taxstat);
    magrespplot(ahatfXV, newfaxhz, confinthighXV, confintlowXV, statname, magxbounds, magybounds, ' vertical', taxstat);
end
fclose('all')
end
end