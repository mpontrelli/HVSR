%% Butter2
  % a Bandpass butterworth filter to detrend groundmotion microtremor data
    %INPUTS
    %x - waveform
       
    %OUTPUTS
    %y - butterworth filtered waveform

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