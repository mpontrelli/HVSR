%% Butter
% Butter generates poles and zeros of a 4th order butterworth filter with
% highcorner equal to the nyquist frequency of the record. It then zero phase filters
% the data.

    % INPUTS
    
    % xNS - north south component of the ground motion record

    % xV - vertical component of the ground motion record
        
    % xEW - east west component of the ground motion record
    
    % fs - the sampling frequency of the record
    
    % OUTPUTS
    
    % xNS - filtered north south component of the ground motion record

    % xV - filtered vertical component of the ground motion record
        
    % xEW - filtered east west component of the ground motion record

%% Author: Marshall Pontrelli
% Date: developed between September, 2017 and August, 2019
%% Start   
function [xNS, xV, xEW] = Butter(xNS, xV, xEW, fs)
% Here are the filter parameters
HighCorner = (fs / 2) - 1.01;
Npoles = 4;  % Corner for 1 pass of the two-pass filter
% Filter Design
fN = (fs / 2) - 1;
Highcut = HighCorner / fN;
[bb1, aa1] = butter(Npoles, Highcut);
% Now filter 
xNS = filtfilt(bb1,aa1,xNS - mean(xNS));
xV = filtfilt(bb1,aa1,xV - mean(xV));
xEW=filtfilt(bb1,aa1,xEW - mean(xEW));
end