%% Rosetta

% Rosetta reads and extracts stations and event info metadata from the Mexico City 
% RACM event textfile. It goes into the textfile and reads the sampling 
% frequency and finds and translates the soil type to be plotted on the 
% final plot that HVSR outputs. If the user wishes to use HVSR for another
% dataset, Rosetta must be written such that readfile1 outputs the
% necessary information. The style is baed off of Dr. John Ebel's codes
% which extract information from similar datafiles. It's where I learned
% it.

    % INPUTS
    
    % filename - path of the file containing the ground motion info and
    % metadata
    
    %OUTPUTS 
    
    % fs - sampling frequency
    
    % soil - soil type on which the station sits

%% Author: Marshall Pontrelli
% Date: developed between September, 2017 and August, 2019

% Update - expanded 2/12/2020 to read all metadata, create a structure out
% of it, incorporate readfile 1 and process the data all in 1 go.
%% start 
function [data] = Rosetta(filename, varargin)
    [fid3,~] = fopen(filename,'r');
    % Now read the file
    %
    % Skip the first 16 lines because they only contain header information.
    %% STATION NAME then read the 17th line which contains the station name
    for jjj = 1:17
        line = fgetl(fid3);
    end
    b = length(line);
    data.meta.station.name = line(42:b);
    
    %% LATITUDE Now skip 5 lines and read the 6th containing latitude of the station
    for jjj = 1:6
        line = fgetl(fid3);
    end
    
    % This code came form https://www.mathworks.com/matlabcentral/answers/113957-how-to-separate-numbers-and-text-from-a-string
    s = string(line);
    data.meta.station.lat = str2double(regexp(s,'\d+(\.\d+)?|\.\d+','match'));
    
    %% LONGITUDE now skip a line and read longitude
    line = fgetl(fid3);
    s = string(line);
    data.meta.station.lon = str2double(regexp(s,'\d+(\.\d+)?|\.\d+','match'));

    
    %% ALTITUDE now skip a line and read altitude
    line = fgetl(fid3);
    s = string(line);
    data.meta.station.alt = str2double(regexp(s,'\d+(\.\d+)?|\.\d+','match'));    
    
    %% SOIL Now skip a line and read in the soil type by reading enough letters to
    % be able to categorize the site into one of 4 categories: Lake Bed,
    % Transition, Compact or Structure. There is only one structure station. 
    line = fgetl(fid3);

    % Soil reads the first three letters of the line 
    soil = line(42:44);
    % soil1 reads the first letter of the second word so in the instance of
    %"terrano", which occurs for both Lake Zone and Compact, Rosetta can use
    %the first letter of the second word to determine what soil type the
    %station is in. 
    soil1 = line(50);

    % LAKE
    if strcmp(soil,'ARC') == 1
        data.meta.station.soil = 'Lake Zone';

    elseif  strcmp(soil,'Alt') == 1
        data.meta.station.soil = 'Lake Zone';
        
    elseif  strcmp(soil,'ALT') == 1
        data.meta.station.soil = 'Lake Zone';

    elseif strcmp(soil,'Ter') == 1 && strcmp(soil1, 'b') == 1
        data.meta.station.soil = 'Lake Zone';
    elseif strcmp(soil, 'TER') == 1 && strcmp(soil1, 'B') == 1
        data.meta.station.soil = 'Lake Zone';
    elseif strcmp(soil,'Arc') == 1
        data.meta.station.soil = 'Lake Zone';
    % LAKE END

    % TRANSITION
    elseif strcmp(soil,'Tra') == 1
        data.meta.station.soil = 'Transition';

    elseif strcmp(soil,'TRA') == 1
        data.meta.station.soil = 'Transition';
    %TRANSITION END

    % COMPACT
    elseif strcmp(soil,'Are') == 1
        data.meta.station.soil = 'Compact';
    elseif strcmp(soil,'Ter') == 1
        data.meta.station.soil = 'Compact';
    elseif strcmp(soil,'ARE') == 1
        data.meta.station.soil = 'Compact';
    elseif strcmp(soil, 'TER') == 1 && strcmp(soil1, 'E') == 1
        data.meta.station.soil = 'Compact';

    % COMPACT END

    % STRUCTURE
    elseif strcmp(soil,'EST') == 1
        data.meta.station.soil = 'Structure';
    % STRUCTURE END
    end
    
    %% INSTRUMENT
    % MODEL skip 7 lines and read the 8th containing instrument information
    for jjj = 1:8
        line = fgetl(fid3);
    end
    a = string(line);
    b = strsplit(a);
    data.meta.instrument.model = char(b(end));
    
    % SERIAL NUMBER now skip a line and get the instrument serial number
    line = fgetl(fid3);
    s = string(line);
    data.meta.instrument.serial_number = str2double(regexp(s,'\d+(\.\d+)?|\.\d+','match'));    

    % SAMPLING FREQUENCY now skip 4 lines and get the sampling frequency
    for jjj = 1:4
        line = fgetl(fid3);
    end    
    s = string(line);
    a = str2double(regexp(s,'\d+(\.\d+)?|\.\d+','match'));  
    data.meta.instrument.fs = a(end);
    fs = a(end);
    
    %% SENSOR NATURAL FREQUENCY now skip 4 lines and get the sensor natural frequency 
    for jjj = 1:4
        line = fgetl(fid3);
    end    
    s = string(line);
    a = str2double(regexp(s,'\d+(\.\d+)?|\.\d+','match'));  
    data.meta.instrument.fn = a(end);
    
    %% SENSOR DAMPING now skip 2 lines and get the sensor damping
    for jjj = 1:2
        line = fgetl(fid3);
    end    
    s = string(line);
    a = str2double(regexp(s,'\d+(\.\d+)?|\.\d+','match'));  
    data.meta.instrument.damping.EW = a(3);
    data.meta.instrument.damping.NS = a(4);
    data.meta.instrument.damping.V = a(end);
    
    %% SENSOR Trigger now skip 4 lines and get the sensor trigger threshold
    for jjj = 1:4
        line = fgetl(fid3);
    end    
    s = string(line);
    a = str2double(regexp(s,'\d+(\.\d+)?|\.\d+','match'));  
    data.meta.instrument.trigger.EW = a(3);
    data.meta.instrument.trigger.NS = a(4);
    data.meta.instrument.trigger.V = a(end);

    %% RECORD now skip 2 lines and read the preevent and post event memory 
    for jjj = 1:2
        line = fgetl(fid3);
    end    
    s = string(line);
    a = str2double(regexp(s,'\d+(\.\d+)?|\.\d+','match'));  
    data.meta.record.pre = a;
    line = fgetl(fid3);
    s = string(line);
    a = str2double(regexp(s,'\d+(\.\d+)?|\.\d+','match'));  
    data.meta.record.post = a;
    
    %% Now skip 5 lines and read the event data
    %% DATE
    for jjj = 1:5
        line = fgetl(fid3);
    end    
    s = string(line);
    a = strsplit(s);  
    data.meta.event.date = a(end);
    
    %% TIME
    line = fgetl(fid3);
    s = string(line);
    a = strsplit(s);  
    data.meta.event.time = a(end);    
    
    %% MAGNITUDE
    line = fgetl(fid3);
    s = string(line);
    a = str2double(regexp(s,'\d+(\.\d+)?|\.\d+','match'));  
    data.meta.event.mag = mean(a);    
    
    %% LATITUDE
    line = fgetl(fid3);
    s = string(line);
    a = str2double(regexp(s,'\d+(\.\d+)?|\.\d+','match'));  
    data.meta.event.lat = a;    
    
    %% LONGITUDE
    line = fgetl(fid3);
    s = string(line);
    a = str2double(regexp(s,'\d+(\.\d+)?|\.\d+','match'));  
    data.meta.event.lon = a; 
    
    %% DEPTH
    line = fgetl(fid3);
    s = string(line);
    a = str2double(regexp(s,'\d+(\.\d+)?|\.\d+','match'));  
    data.meta.event.depth = a; 
    
    %% Now back to the record
    % Skip 6 lines and read the time of first sample
    for jjj = 1:6
        line = fgetl(fid3);
    end    
    s = string(line);
    a = strsplit(s);  
    data.meta.record.time_first_sample = a(end); 
    
    %% duration of record
    for jjj = 1:2
        line = fgetl(fid3);
    end    
    s = string(line);
    a = str2double(regexp(s,'\d+(\.\d+)?|\.\d+','match'));  
    data.meta.record.duration = a(end); 
    
    %% record number of samples
    for jjj = 1:2
        line = fgetl(fid3);
    end    
    s = string(line);
    a = str2double(regexp(s,'\d+(\.\d+)?|\.\d+','match'));  
    data.meta.record.samples = a(end); 
    
    %% record PGAs
    for jjj = 1:2
        line = fgetl(fid3);
    end    
    s = string(line);
    a = str2double(regexp(s,'\d+(\.\d+)?|\.\d+','match'));  
    data.meta.record.PGA.NS = a(3); 
    data.meta.record.PGA.EW = a(end); 
    data.meta.record.PGA.V = a(4); 
    
    %% record PGA indices
    line = fgetl(fid3);
    s = string(line);
    a = str2double(regexp(s,'\d+(\.\d+)?|\.\d+','match'));  
    data.meta.record.PGA.NS_index = a(3); 
    data.meta.record.PGA.EW_index = a(end); 
    data.meta.record.PGA.V_index = a(4); 
    
    %% record units
    for jjj = 1:3
        line = fgetl(fid3);
    end    
    s = string(line);
    a = strsplit(s);  
    b = strcat(a(end-1), a(end));
    data.meta.record.units = b;
    
    fclose(fid3);
      
    
    %% Ok so we read in all the metadata, now let's read in the data
    % Set space deliminator
    deliminator = '';
    % Data starts in row 109, column 1
    R = 109; % Row data start text file
    C = 0; % column data start text file
    data1 = dlmread(filename,deliminator,R,C); % retrieve data
    NS = data1(:,1)/980; % North_South_Component
    EW = data1(:,3)/980; % East_West_Component
    V = data1(:,2)/980; % Vertical_Component
    time = (1:length(V))/fs;
    data.processing.rawdata.NS = NS; %North_South_Component
    data.processing.rawdata.EW = EW; %East_West_Component
    data.processing.rawdata.V = V; %Vertical_Component
    data.processing.rawdata.time = time;

    %% Now time to start processing the data. Start by running a Butterworth filter over it.
    
    LowCorner = 0.1;
    HighCorner = fs/2 - 1;
    Npoles = 4;

    [V] = (Butter2(V, fs, 'LowCorner', LowCorner, 'HighCorner', HighCorner, 'Npoles', Npoles));
    [NS] = (Butter2(NS, fs, 'LowCorner', LowCorner, 'HighCorner', HighCorner, 'Npoles', Npoles));
    [EW] = (Butter2(EW, fs, 'LowCorner', LowCorner, 'HighCorner', HighCorner, 'Npoles', Npoles));

    data.processing.filtereddata.NS = NS; %North_South_Component
    data.processing.filtereddata.EW = EW; %East_West_Component
    data.processing.filtereddata.V = V; %Vertical_Component
    data.processing.filtereddata.time = time;
    data.processing.filtereddata.filter.lowcorner = LowCorner;
    data.processing.filtereddata.filter.highcorner = HighCorner;
    data.processing.filtereddata.filter.npoles = Npoles;
    
    %% Now do Arias intensity
    [Iaval, D5D95, D5D75, rate_arias, Ianorm] =  Arias(time, NS, fs);
    data.processing.filtereddata.arias.NS.intensity = Iaval;
    data.processing.filtereddata.arias.NS.D595 = D5D95;
    data.processing.filtereddata.arias.NS.D575 = D5D75;
    data.processing.filtereddata.arias.NS.rate = rate_arias;
    data.processing.filtereddata.arias.NS.normalized = Ianorm;
    
    [Iaval, D5D95, D5D75, rate_arias, Ianorm] =  Arias(time, EW, fs);
    data.processing.filtereddata.arias.EW.intensity = Iaval;
    data.processing.filtereddata.arias.EW.D595 = D5D95;
    data.processing.filtereddata.arias.EW.D575 = D5D75;
    data.processing.filtereddata.arias.EW.rate = rate_arias;
    data.processing.filtereddata.arias.EW.normalized = Ianorm;
    
    [Iaval, D5D95, D5D75, rate_arias, Ianorm] =  Arias(time, V, fs);
    data.processing.filtereddata.arias.V.intensity = Iaval;
    data.processing.filtereddata.arias.V.D595 = D5D95;
    data.processing.filtereddata.arias.V.D575 = D5D75;
    data.processing.filtereddata.arias.V.rate = rate_arias;
    data.processing.filtereddata.arias.V.normalized = Ianorm;
    %% Response spectra
    [y, displacement, velocity, period, max1, per_max, resp_02_secs, resp_05_secs, resp_1_secs, resp_2_secs, resp_5_secs] = Response_Spectra(NS, fs, 5);    
    data.processing.filtereddata.spectra.NS.accel = y;
    data.processing.filtereddata.spectra.NS.velocity = velocity;
    data.processing.filtereddata.spectra.NS.velocity = displacement;
    data.processing.filtereddata.spectra.period_vec = period;
    
    
end