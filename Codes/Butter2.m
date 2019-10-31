%% Butter2
  % a Bandpass butterworth filter to detrend groundmotion microtremor data
    %INPUTS
    %x - waveform
       
    %OUTPUTS
    %y - butterworth filtered waveform

function [y] = Butter2(x)
% Here are the filter parameters
LowCorner=0.5;
%HighCorner=9;
HighCorner=49;
Npoles=4;  % Corner for 1 pass of the two-pass filter

fN=100/2;
Lowcut=LowCorner/fN;
Highcut=HighCorner/fN;
[bb1 aa1]=butter(Npoles,[Lowcut Highcut]);
%freqz(bb1,aa1)
y=filtfilt(bb1,aa1,(x-mean(x)));
end