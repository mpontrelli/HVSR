

% Before running this file the first time, execute the following statement
javaaddpath IRIS-WS-2.0.18.jar
% in the Matlat command window

% Set the filter parameters
LowCorner=.25;
%HighCorner=9;
HighCorner=9;
Npoles=4;  % Corner for 1 pass of the two-pass filter

% Set the data time window and station
time1='2020-11-07 14:10:00';
time2='2020-11-07 14:15:00';
Station='WES';
%Station='ACCN';
Network='*';
Component='BHZ';
% Component='HHZ';

% Set the time axis tick interval
xticks=10;


% Get the data
mytrace=irisFetch.Traces(Network,Station,'*',Component,time1,time2);
trace1 = mytrace(1);

% Filter the data
fN=trace1.sampleRate/2;
Lowcut=LowCorner/fN;
Highcut=HighCorner/fN;
[bb1 aa1]=butter(Npoles,[Lowcut Highcut]);

waveformfilt=filtfilt(bb1,aa1,double(trace1.data-mean(trace1.data)));

% Next prepare and plot the data
sampletimes=linspace(trace1.startTime,trace1.endTime,trace1.sampleCount);
figure
plot(sampletimes,trace1.data);
%xtickvec=[trace1.startTime:xticks/(60*60*24):trace1.endTime];

datetick('x',13,'keeplimits','keepticks');
title(trace1.station)
xEW = trace1.data;
%set(gca,'XTick',xtickvec);
%set(gca,'XTick',sampletimes);