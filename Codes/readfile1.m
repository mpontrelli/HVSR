%% readfile1

% readfile1 is designed for the event files in the RACM network. Each file
% has 109 lines of heading and then three columns which contain the NS, V
% and EW components of the associated seismometer. It reads these and
% outputs them as a vector. It also reads the sampling frequency of the
% station and the soil type of the station using the subrouting "Rosetta".
% If the user wishes to run HVSR on another ground motion dataset, they
% must rewrite readfile1 such that it pulls the same information from their
% dataset. This includes rewriting "Rosetta" which is also taylored to the
% RACM dataset.

    % INPUTS
    
    % filename - Filename of the ground motion datafile
    
    % OUTPUTS
    
    % xNS - vector of north-south ground motion
    
    % xV - vector of vertical ground motion
    
    % xEW - vector of east-west ground motion
    
    % fs - sampling frequency of record
    
    % soil - soil type as indicated in the RACM network database

%% Author: Marshall Pontrelli
% Date: developed summer, 2019 
%% 
function [xNS,xV,xEW, fs, soil] = readfile1(filename)
%Set space deliminatorIt
deliminator='';
%Data starts in row 109, column 1
R=109; %Row data start text file
C=0; %column data start text file
data=dlmread(filename,deliminator,R,C); %retrieve data
xNS=data(:,1); %North_South_Component
xV=data(:,2); %Vertical_Component
xEW=data(:,3); %East_West_Component
end 


