%% Rosetta

% Rosetta reads and extracts stations and event info metadata from the Mexico City 
% RACM event textfile. It goes into the textfile and reads the sampling 
% frequency and finds and translates the soil type to be plotted on the 
% final plot that HVSR outputs. If the user wishes to use HVSR for another
% dataset, Rosetta must be written such that readfile1 outputs the
% necessary information.

    % INPUTS
    
    % filename - path of the file containing the ground motion info and
    % metadata
    
    %OUTPUTS 
    
    % fs - sampling frequency
    
    % soil - soil type on which the station sits

%% Author: Marshall Pontrelli
% Date: developed between September, 2017 and August, 2019
%% 
function [fs, soil] = Rosetta(filename)
[fid3,~] = fopen(filename,'r');

% Now read the file
%
% Skip the first 26 lines because they only contain header information.
% This section retrieves the soil type by reading enough letters to
% be able to categorize the site into one of 4 categories: Lake Bed,
% Transition, Compact or Structure. There is only one structure station. 
for jjj = 1:26
    line = fgetl(fid3);
end

b = length(line);
% Soil reads the first three letters of the line 
soil = line(42:44);
% soil1 reads the first letter of the second word so in the instance of
%"terrano", which occurs for both Lake Zone and Compact, Rosetta can use
%the first letter of the second word to determine what soil type the
%station is in. 
soil1 = line(50);

%LAKE
if strcmp(soil,'ARC') == 1
    soil = 'Lake Zone';

elseif  strcmp(soil,'Alt') == 1
        soil = 'Lake Zone';
        
elseif  strcmp(soil,'ALT') == 1
soil = 'Lake Zone';

elseif strcmp(soil,'Ter') == 1 && strcmp(soil1, 'b') == 1
        soil = 'Lake Zone';
elseif strcmp(soil, 'TER') == 1 && strcmp(soil1, 'B') == 1
    soil = 'Lake Zone';
elseif strcmp(soil,'Arc') == 1
soil = 'Lake Zone';
%LAKE END

%TRANSITION
elseif strcmp(soil,'Tra') == 1
    soil = 'Transition';

elseif strcmp(soil,'TRA') == 1
soil = 'Transition';
%TRANSITION END

%COMPACT
elseif strcmp(soil,'Are') == 1
    soil = 'Compact';
elseif strcmp(soil,'Ter') == 1
    soil = 'Compact';
elseif strcmp(soil,'ARE') == 1
    soil = 'Compact';
elseif strcmp(soil, 'TER') == 1 && strcmp(soil1, 'E') == 1
    soil = 'Compact';

%COMPACT END

%STRUCTURE
elseif strcmp(soil,'EST') == 1
    soil = 'Structure';
%STRUCTURE END
end

%Now we read the sampling frequency which varies from station to station

for jjj = 1:13
    line = fgetl(fid3);
end
a = length(line);
fs = str2num(line(a-2:a));
end