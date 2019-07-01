% readfile1 is designed for the event files in the RACM network. Each file
%has 109 lines of heading and then three columns which contain the NS, V
%and EW components of the associated seismometer. It reads these and
%outputs them as a vector. It also reads the sampling frequency of the
%station and the soil type of the station using the subrouting "Rosetta"

function [xNS,xV,xEW, fs, soil] = readfile1(filename)
%Set space deliminatorIt
deliminator='';
%Data starts in row 109, column 1
R=109; %Row data start text file
C=0; %column data start text file
data=dlmread(filename,deliminator,R,C); %retrieve data
[fs, soil] = Rosetta(filename);

xNS=data(:,1); %North_South_Component
xV=data(:,2); %Vertical_Component
xEW=data(:,3); %East_West_Component

end 


