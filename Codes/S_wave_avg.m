%% S_wave_avg
% compute the average shear wave velocity profile

    % INPUTS
    % d - a vector of depths
    
    % v - a  vector of velocities (same length as d)
    
    % OUTPUTS
    
    % v_avg - the average shear wave velocity of the profile
    
    
%% Author: Marshall Pontrelli
% Date: 11/25/2020

function [v_avg] =  S_wave_avg(d,v)
    den = zeros(1,length(d));
    for i = 1:length(d)
        den(i) = d(i)/v(i);
    end
    den = sum(den);
    v_avg = sum(d)/den;
end