%% PGA
% PGA picks out the peak ground acceleration of each of the three
% components of the input ground motion. The factor "/981" is to convert to
% g's from gals. This must be modified if the input ground motion is in
% units other than gals. 

    % INPUTS
    
    % xNS - north south component of the ground motion record

    % xV - vertical component of the ground motion record
        
    % xEW - east west component of the ground motion record
    
    % OUTPUTS
    
    % PGANS - peak ground acceleration of the NS component
    
    % PGAV - peak ground acceleration of the V component
    
    % PGAEW - peak ground acceleration of the EW component
    
%% Author: Marshall Pontrelli
% Date: developed between September, 2017 and August, 2019   

%% Start
function [PGANS,PGAV,PGAEW] = PGA(xNS,xV,xEW)
PGANS = max(abs(xNS))/981;
PGAV = max(abs(xV))/981;
PGAEW = max(abs(xEW)/981);
end 