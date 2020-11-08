clear all

% Before running this file the first time, execute the following statement
javaaddpath IRIS-WS-2.0.18.jar
% in the Matlat command window

time1='2020-11-08 14:10:00';
time2='2020-11-08 14:15:00';

Channel='BHZ';
Network='*';

% Here are the filter parameters
LowCorner=3;
%LowCorner=.01;
HighCorner=15;
Npoles=2;  % Corner for 1 pass of the two-pass filter

N4Stations={
'D62A'
'E63A'
'F64A'
'G65A'
'E62A'
'F63A'
'F62A'
'G62A'
'H62A'
'I63A'
'I62A'
'J61A'
'K62A'
'L64A'
'L61B'
'M63A'
'N62A'
'J59A'
'J58A'
'J57A'
'J56A'
'K58A'
'L59A'
'K57A'
'P60A'
};

% Loop over all of the stations to get the waveforms
for jjkk=1:length(N4Stations)
    Station=N4Stations{jjkk};

% Get the data
mytrace=irisFetch.Traces(Network,Station,'*',Channel,time1,time2);

% Test to see if the data exist
if isempty(mytrace)==0
    % The data exist for this station, so process it
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
plot(sampletimes,waveformfilt);
datetick('x',13,'keeplimits')
title(Station)
end  % if isempty(mytrace)==0

end  % for jjkk=1:length(N4Stations)
