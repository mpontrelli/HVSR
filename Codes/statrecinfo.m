%% statrecinfo

% statrec info is designed specifically for the Mexico City RACM dataset.
% It loops through all the records at a station and pulls out the longest
% record, shortest record and lowest sampling frequency used to record any
% of the records. It is only needed if wavecut is toggled on. If you want
% to use HVSR on another dataset and use wavecut, then you must reconfigure
% statrecinfo to get the same information from the other dataset. Likely,
% this just involves changing readfile1. It was developed in response to
% questions over upsampling of data onto the longest record using linear 
% interpolation.

    % INPUTS
    
    % files - list of files in the folder
    
    % station - folder path that contains the waveform files
    
    % OUTPUTS
    
    % fsmin - lowest sampling frequency used to record any of the data
    
    % recmax = maximum number of samples in the recording
    
    % minimum number of samples in the recording
    
%% Author: Marshall Pontrelli
% Date: developed summer, 2019 
%% 
function [fsmin, recmax, recmin] = statrecinfo(files, station)
counter1 = 0;
lengthvec = [];
fsvec = [];
for file = files'
    counter1 = counter1 + 1;
    filename = strcat(station,'\',file.name);
    [xNS,xV,xEW, fs] = readfile1(filename);
    lengthvec(counter1, :) = length(xNS);
    fsvec(counter1, :) = fs;
end
fsmin = min(fsvec);
recmax = max(lengthvec);
recmin = min(lengthvec);

end