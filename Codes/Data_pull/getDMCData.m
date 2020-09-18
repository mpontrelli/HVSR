function [sampletimes,trace1,waveformfilt, samplerate, sensitivity, sensunits] = getDMCData(time1,time2,Station,Network,Component, LowCorner, HighCorner)
%Function used to pull data from IRIS 

% Before running this file the first time, execute the following statement
javaaddpath IRIS-WS-2.0.18.jar
% in the Matlat command window

% Set the filter parameters
LowCorner = LowCorner;
HighCorner = HighCorner;
Npoles=4;  % Corner for 1 pass of the two-pass filter

% Set the data time window and station
%time1='2018-04-22 20:35:51';
%time2='2018-04-23 24:00:00';
%Station='PKME';
%Network='US';
%Component='BHZ';

% Set the time axis tick interval
xticks=10;

% Get the data
mytrace=irisFetch.Traces(Network,Station,'*',Component,time1,time2);

%Conditional Statement: if the the Station did not record data at that time
%then output nothing
%Check to see if the Struc is empty
Check = isempty(mytrace);
if Check == 1
    sampletimes = 0; trace1 = 0; waveformfilt = 0;
    return
end

trace1 = mytrace(1);
sensitivity = trace1.sensitivity;
sensunits = trace1.sensitivityUnits;
samplerate = trace1.sampleRate;

% Filter the data
fN=trace1.sampleRate/2;
Lowcut=LowCorner/fN;
Highcut=HighCorner/fN;
[bb1, aa1]=butter(Npoles,[Lowcut Highcut]);

waveformfilt=filtfilt(bb1,aa1,double(trace1.data-mean(trace1.data)));

% Next prepare and plot the data
sampletimes=linspace(trace1.startTime,trace1.endTime,trace1.sampleCount);
end 