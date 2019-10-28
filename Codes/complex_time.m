%% Complex time
% Compute the complex time series which combines the horizontal and
% vertical components of the record. Based of Steidl et al. 1994.

    % INPUTS
    
    % xNS - north south component
    
    % xEW - east west component
    
    % OUTPUTS 
    
    % xH - complex time series of xEW and xNS
%% Author: Marshall Pontrelli
% Date: developed between September, 2017 and August, 2019
%% Start
function [xH] =  complex_time(xNS, xEW)
    xH = xNS + 1i*xEW;
end
