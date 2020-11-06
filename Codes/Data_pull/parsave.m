%% parsave
% save functions from within a parallel for loop
    % INPUTS
    
    % OUTPUTS
    
    
%% Author: Marshall Pontrelli
% Date: 4/17/2020

function parsave(fname, data1,samplerate, sensitivity, sensunits)
  save(fname, 'data1','samplerate', 'sensitivity', 'sensunits')
end