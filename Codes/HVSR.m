%% HVSR
% HVSR is designed to process many earthquake ground motion records at one
% site and process them using Nakamura's HVSR technique. There are several
% potential output plots including the site averaged HVSR, the HVSRs at a
% site all on one plot before averaging) and the magnitude responses. 

    %% data structure
    % data should be stored as follows: a master folder contains
    % subfolders with individual stations, each containing waveforms of the
    % same format. Naturally, the station folders should be the name of the
    % station

    %% INPUTS
    
    % path - File directory to accompanying HVSR software folder
    
    % datapath - path where the data are stored. This is a master folder
    % containing folders of stations containing ground motions at that
    % station
    
    % varargin 
        % wavecut - When toggled on ('yes'), this conditional cuts all 
        % waveforms to the length of the shortest waveform in the dataset
        % instead of using linear interpolation
        
        % HVSR - When toggled on ('yes'), this takes the HVSR matrix and
        % averages them together using equations 2-4 from Thompson et al.
        % 2012, then it plots them. It then identifies significant peaks
        % in the HVSR and computes the peak statistics on them.
        
        % individHVSR - When toggled on ('yes'), this plots the individual
        % HVSRs all on 1 plot
        
        % magresps - When toggle on ('yes'), this plots the horizontal and vertical
        % averaged magnitude responses 
        
        % magxbound - x bounds for the magnitude responses. The default xbound is ([0, 20]);
        
        % magybounds - y bounds for the magnitude responses. The default ybound is ([0.0001, 100]);
    
    
    %% OUTPUTS
    
    % newfaxhz - Interpolated frequency vector
    
    % ahatf - Average HVSR from Thompson et al. 2012 equations 2 and 3
    
    % sigma - Standard deviation from Thompson et al. 2012 equation 4
    
    % lowbound - Lowest frequency that can be resolved with the length of
    % the data
    
    % taxstat - Output matrix of HVSR statistics
    
    % statsfinal - Output matrix of HVSR statistics (this gets updated if
    % you want to run multiple stations to create a matrix of statistics on
    % many stations
    
    % fsmin - minimum sampling frequency in your dataset. This was added to
    % account for multiple sampling frequencies in the Mexico City RACM 
    % dataset (100 and 200 Hz)
    
    % recmax - maximum record length in the station dataset
   
%% Author: Marshall Pontrelli
% Date: developed between September, 2017 and August, 2019
%% Start
function [newfaxhz, ahatf, sigma, lowbound, taxstat, statsfinal, fsmin, recmax, varargout] = HVSR(path, datapath, varargin)
%% parse inputs
% create Input Parser object
p = inputParser;
% add inputs to the scheme
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
% parse the inputs
parse(p, path, datapath, varargin{:})
% set varibales from the parse
wavecut = p.Results.wavecut;
HVSR = p.Results.HVSR;
individHVSR = p.Results.individHVSR;
magresps = p.Results.magresps;
magxbounds = p.Results.magxbounds;
magybounds = p.Results.magybounds;
%% start function

% ACCESSING THE DATA
% go into the data folder and get a list of stations
cd(datapath)
stationlist = dir;
stationlist = stationlist(3:length(stationlist));
statsfinal = [];
% start the for loop that goes through all the station folders
for eee = 1:length(stationlist)
    station = stationlist(eee); % read the folder info
    statname = station.name; % extract the station name (this is the name of folder)
    station = strcat(datapath, '\', statname); % create the path name of the station folder
    % disp(station % you can uncheck this if you want to double check that the paths are working
    % go into data directory and build structure of all files in it
    cd(station)
    files = dir;
    files = files(3:length(files));
    % now you're just looking at the waveform files of the station
    % create empty matrix that gets filled with all H/V values, counter is used
    % to index this matrix
    HV_final_matrix = [];
    XH_final_matrix = [];
    XV_final_matrix = [];
    lowbound_matrix = [];
    counter = 0;
    % change directory back to codes to access functions needed 
    d = strcat(path, '\HVSR\Codes');
    cd(d)
    % Loop through the records once to get the number of samples in each
    % record for the x-axis linear interpolation and the sampling freqeuncy 
    % of each record. The number of samples value will be the
    % sample number of the record with the maximum number of samples. The
    % lowest sampling freqeuncy (fsmin) is required for the frequency axis
    % high bond and is used as the nyquist frequency for the whole analysis
    [fsmin, recmax, recmin] = statrecinfo(files, station);
    
    % Loop through all the ground motion records
    for file = files'
        filename = strcat(station,'\',file.name);
        [xNS,xV,xEW, fs] = readfile1(filename);
        N = length(xNS); % length of the signal
        %Low bound calculation
        lowbound = 1/(N/fs); % this is the lowest frequency that we can image
        [PGANS,PGAV,PGAEW] = PGA(xNS,xV,xEW);
        filename = strcat(station,'\',file.name); % pull out entire file name
        [xNS,xV,xEW, fs] = readfile1(filename); % read NS, V, EW and fs
        [PGANS,PGAV,PGAEW] = PGA(xNS,xV,xEW); %extract PGA from the record
        if PGANS < 0.1 && PGAV < 0.1 && PGAEW < 0.1
            if strcmp(wavecut, 'yes') == 1   
                [xNS,xV,xEW] = waveformcut(xNS,xV,xEW, recmin);
                counter = counter + 1;
                
                % filter the data
                [xNS, xV, xEW] = Butter(xNS, xV, xEW, fs);
                
                % compute the complex time series
                [xH] =  complex_time(xNS, xEW);
                
                %compute magnitude responses
                ff = 2; % half magnitude spectra
                % if you have a sampling frequency greater than your lowest
                % sampling frequency, dived by more (for example, in Mexico
                % City RACM dataset, there are records with both 100 and
                % 200 Hz. If a record has 200 Hz fs, split response spectra
                % up by 4, not 2, then all the averaged outputs will be 0
                % to 50 Hz)
                if fs > fsmin 
                    ff = ff * (fs / fsmin);
                end
                [XH_mag, fax_HzN] =  Magresp(xH, fs, ff);
                [XV_mag] =  Magresp(xV, fs, ff);
                
                % smooth the magnitude response
                width = .1; % width for smoothing filter in hz
                q = ceil((N/fs)*width); % width for smoothing filter in samples
                XV_magfilt=smooth(XV_mag,q); 
                XH_magfilt=smooth(XH_mag,q);  
                [H_V1] = HV(XH_magfilt, XV_magfilt);
                XH_final_matrix(counter, :) = XH_magfilt;
                XV_final_matrix(counter, :) = XV_magfilt; 
                HV_final_matrix(counter, :) = H_V1; 
                clear H_V1
            else
            counter = counter + 1;
            
            % filter the data
            [xNS, xV, xEW] = Butter(xNS, xV, xEW, fs); %filter the data
            
            % compute the complex time series
            [xH] =  complex_time(xNS, xEW);
            
            %compute magnitude responses
            ff = 2; % half magnitude spectra
            % if you have a sampling frequency greater than your lowest
            % sampling frequency, dived by more (for example, in Mexico
            % City RACM dataset, there are records with both 100 and
            % 200 Hz. If a record has 200 Hz fs, split response spectra
            % up by 4, not 2, then all the averaged outputs will be 0
            % to 50 Hz)
            if fs > fsmin 
                ff = ff * (fs / fsmin);
            end
            [XH_mag, fax_HzN] =  Magresp(xH, fs, ff);
            [XV_mag] =  Magresp(xV, fs, ff);
                
            % smooth the magnitude response
            width = .1; % width for smoothing filter in hz
            q = ceil((N/fs)*width); % width for smoothing filter in samples
            XV_magfilt=smooth(XV_mag,q); 
            XH_magfilt=smooth(XH_mag,q);  
 
            %perform H/V
            [H_V1] = HV(XH_magfilt, XV_magfilt);
            
            %make Hz vector and linear interpolate all H/V ETFs to this vector
            newfaxhz = 0: (1/ (ceil(recmax/2) - fsmin))*(fsmin/2 - 1): (fsmin/2 - 1);
            [~, lowindex] = min(abs(newfaxhz - lowbound));
            lowbound_matrix(counter, :) = lowindex;
            %mag resp matrix build
            newXH_mag = interp1(fax_HzN, XH_magfilt, newfaxhz);
            XH_final_matrix(counter, :) = newXH_mag; 
            newXV_mag = interp1(fax_HzN, XV_magfilt, newfaxhz);
            XV_final_matrix(counter, :) = newXV_mag; 
            newH_V1 = interp1(fax_HzN, H_V1, newfaxhz);
            HV_final_matrix(counter, :) = newH_V1; 
            clear H_V1
            end
        end
    end
if strcmp(wavecut, 'yes') == 1  
    newfaxhz = fax_HzN;
    [~, lowbound] = min(abs(newfaxhz - lowbound));
else
    newfaxhz = 0: (1/ (ceil(recmax/2) - fsmin))*(fsmin/2 - 1) : (fsmin/2 - 1);
    lowbound = max(lowbound_matrix);
end
%  comes from lowbound _matrix and corresponds to the lowest
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