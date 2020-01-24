%% Butter2
  % a Bandpass butterworth filter to detrend groundmotion microtremor data
    %INPUTS
    %x - waveform
       
    %OUTPUTS
    %y - butterworth filtered waveform
%% Author: Marshall Pontrelli
% Date: developed between September, 2017 and August, 2019

% Updates: 
% January 2020 - added inputs LowCorner, Highcorner, Npoles, fs and
% filterplot to allow for filter parameter input in HVSRmicro. It is also
% now used in Mexico_City_processing for processing the time domain
% waveforms.

%% Start
function [y] = Butter2(x, LowCorner, HighCorner, Npoles, fs, Filterplot)
% Here are the filter parameters


fN=fs/2;
Lowcut=LowCorner/fN;
Highcut=HighCorner/fN;
[bb1, aa1]=butter(Npoles,[Lowcut Highcut]);
%freqz(bb1,aa1)
y=filtfilt(bb1,aa1,(x-mean(x)));

if strcmp(Filterplot, 'yes') == 1
    freqz(bb1,aa1)
end

end