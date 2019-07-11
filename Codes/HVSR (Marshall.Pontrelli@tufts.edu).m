function [statsfinal, varargout] = HVSR(path, datapath, varargin)
%testing github Yay!
%create Input Parser object
p = inputParser;
%add inputs to the scheme
defaultxbound = ([0, 20]);
defaultybound = ([0.0001, 100]);
addRequired(p,'path',@ischar);
addRequired(p,'datapath',@ischar);
addParameter(p, 'HVSR', 'No', @ischar);
addParameter(p, 'individHVSR', 'No', @ischar);
addParameter(p, 'magresps', 'No', @ischar);
addParameter(p, 'magxbounds', defaultxbound, @isnumeric);
addParameter(p, 'magybounds', defaultybound, @isnumeric);
%parse the inputs
parse(p, path, datapath, varargin{:})
%set varibales from the parse
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
    station = strcat(station.folder, '\', statname);
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
    for file = files'
        filename = strcat(station,'\',file.name);
        [xNS,xV,xEW, fs] = readfile1(filename);
        [PGANS,PGAV,PGAEW] = PGA(xNS,xV,xEW);
        if PGANS < 0.1 && PGAV < 0.1 && PGAEW < 0.1
            counter = counter + 1;
            %fs = station.fs; %sampling frequency in hz
            [xNS, xV, xEW] = Butter(xNS, xV, xEW, fs); %filter the data
            [N_2, fax_HzN, XH_magfilt, XV_magfilt, XH_mag, XV_mag, lowbound] =  Magresp(xNS, xV, xEW, fs); %Compute mag responses and run through triangular filter
            %perform H/V
            [H_V1] = HV(XH_magfilt, XV_magfilt);
            %make Hz vector and linear interpolate all H/V ETFs to this vector
            newfaxhz = 0:0.001:20;
            %mag resp matrix build
            newXH_mag = interp1(fax_HzN, XH_mag, newfaxhz);
            XH_final_matrix(counter, :) = newXH_mag; 
            newXV_mag = interp1(fax_HzN, XV_mag, newfaxhz);
            XV_final_matrix(counter, :) = newXV_mag; 
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

%HVSR
if strcmp(HVSR, 'yes') == 1   
    [ahatf, sigma, confinthigh, confintlow] = wavav(HV_final_matrix);
    HVSRplot(ahatf, newfaxhz, confinthigh, confintlow, statname);  
    [peakamp, peakfreq, amplocs2] = peakiden(ahatf, newfaxhz);
    [taxstat] = specratstat(peakamp, peakfreq, amplocs2, ahatf, newfaxhz, sigma, statname);
    statsfinal = vertcat(statsfinal, taxstat);
end

if strcmp(individHVSR, 'yes') == 1   
    individplot(HV_final_matrix, newfaxhz, statname);
end
%Mag responses
if strcmp(magresps, 'yes') == 1     
    [ahatfXH, sigma, confinthighXH, confintlowXH] = wavav(XH_final_matrix);
    [ahatfXV, sigma, confinthighXV, confintlowXV] = wavav(XV_final_matrix);
    magrespplot(ahatfXH, newfaxhz, confinthighXH, confintlowXH, statname, magxbounds, magybounds, ' horizontal');
    magrespplot(ahatfXV, newfaxhz, confinthighXV, confintlowXV, statname, magxbounds, magybounds, ' vertical');
end
fclose('all')
end
end