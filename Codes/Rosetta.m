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
    % Skip the first 9 lines because they only contain header information.

    %% EVENT ID then read the 9th line which contains the event ID
    for jjj = 1:9
        line = fgetl(fid3);
    end
    b = length(line);
    data.meta.event.ID = line(46:b); 
    %% STATION NAME then skip 8 lines and read the 17th line which contains the station name
    for jjj = 1:8
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
    data.meta.station.lon = -str2double(regexp(s,'\d+(\.\d+)?|\.\d+','match'));

    
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
    elseif strcmp(soil, 'TER') == 1 && strcmp(soil1, 'F') == 1
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
    if length(a) < 3
        data.meta.instrument.trigger.EW = [];
        data.meta.instrument.trigger.NS = [];
        data.meta.instrument.trigger.V = [];
    else
        data.meta.instrument.trigger.EW = a(3);
        data.meta.instrument.trigger.NS = a(4);
        data.meta.instrument.trigger.V = a(end);
    end

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
    data.meta.event.lon = -a; 
    
    %% DEPTH
    line = fgetl(fid3);
    s = string(line);
    a = str2double(regexp(s,'\d+(\.\d+)?|\.\d+','match'));  
    data.meta.event.depth = a; 
    
    %% azimuth
    az = azimuth(data.meta.station.lat,data.meta.station.lon,data.meta.event.lat, data.meta.event.lon);
    data.meta.event.azimuth = az;
    
    %% epcentral distance
    dist = deg2km(distance(data.meta.station.lat,data.meta.station.lon,data.meta.event.lat, data.meta.event.lon));
    data.meta.event.epi_dist = dist;
    
    %% hypocentral distance
    
    
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
    if length(a) < 3
        data.meta.record.PGA.NS_index = []; 
        data.meta.record.PGA.EW_index = []; 
        data.meta.record.PGA.V_index = [];  
    else
    data.meta.record.PGA.NS_index = a(3); 
    data.meta.record.PGA.EW_index = a(end); 
    data.meta.record.PGA.V_index = a(4); 
    end
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
    data.processing.rawdata.time = time';

    %% Now time to start processing the data. Start by running a Butterworth filter over it.
    
    LowCorner = 0.01;
    HighCorner = fs/2 - 1;
    Npoles = 4;

    [NS, V, EW] = Butter(NS, V, EW, fs);
%     [V] = (Butter2(V, fs, 'LowCorner', LowCorner, 'HighCorner', HighCorner, 'Npoles', Npoles));
%     [NS] = (Butter2(NS, fs, 'LowCorner', LowCorner, 'HighCorner', HighCorner, 'Npoles', Npoles));
%     [EW] = (Butter2(EW, fs, 'LowCorner', LowCorner, 'HighCorner', HighCorner, 'Npoles', Npoles));
    %% rotate the horizontals
    for i = 1:length(NS)
        rot(i) = NS(i)*cos(pi/4) + EW(i)*sin(pi/4);
    end
    rot = rot';
    
    %% save vectors that aren't zero padded for later and make them in m/s
    NS_orig = NS * 9.8;
    EW_orig = EW * 9.8;
    V_orig = V * 9.8;
    rot_orig = rot * 9.8;
   
    time_orig = time;
    data.processing.filtereddata.acceleration.NS.waveform_orig = NS_orig; %North_South_Component
    data.processing.filtereddata.acceleration.EW.waveform_orig = EW_orig; %East_West_Component
    data.processing.filtereddata.acceleration.V.waveform_orig = V_orig; %Vertical_Component
    data.processing.filtereddata.acceleration.rotated.waveform_orig = rot_orig; %rotated
    data.processing.filtereddata.time_orig = time_orig; %rotated
    %% Zero pad on either side to equal 100000 or 200000 depending on fs
    if fs == 100
        len = length(V);
        f = (100000 - len)/2;
        if mod(f,1) == 0
            c1 = zeros(1,f);
            c2 = zeros(1,f);
            loc1 = f;
            loc2 = 100000 - f;
        else
            zz = randi([1, 2], 1);
            z = [0.5, -0.5];
            zzz = z(zz);
            q = (f + zzz);
            q1 = (f - zzz);
            c1 = zeros(1, q);
            c2 = zeros(1, q1);
            loc1 = q;
            loc2 = 100000 - q1;
        end 
    end
    if fs == 200
        len = length(V);
        f = (200000 - len)/2;
        if mod(f,1) == 0
            c1 = zeros(1,f);
            c2 = zeros(1,f);
        else
            zz = randi([1, 2], 1);
            z = [0.5, -0.5];
            zzz = z(zz);
            q = (f + zzz);
            q1 = (f - zzz);
            c1 = zeros(1, q);
            c2 = zeros(1, q1); 
        end   
    end
    time = vertcat(c1',time',c2');
    NS = vertcat(c1',NS,c2');
    EW = vertcat(c1',EW,c2');
    V = vertcat(c1',V,c2');
    rot = vertcat(c1',rot,c2');
    
    %

    data.processing.filtereddata.acceleration.NS.waveform = NS; %North_South_Component
    data.processing.filtereddata.acceleration.EW.waveform = EW; %East_West_Component
    data.processing.filtereddata.acceleration.V.waveform = V; %Vertical_Component
    data.processing.filtereddata.acceleration.rotated.waveform = rot; %rotated
    data.processing.filtereddata.time = time';
    data.processing.filtereddata.filter.lowcorner = LowCorner;
    data.processing.filtereddata.filter.highcorner = HighCorner;
    data.processing.filtereddata.filter.npoles = Npoles;
    

    
     
    %% North-south
    %% Now integrate the waveform for velocity and displacement
    [PGA1, PGV1, PGD1, v, d] =  waveform_integrate(NS, fs);
    data.processing.filtereddata.acceleration.NS.PGA = PGA1;
    data.processing.filtereddata.velocity.NS.waveform = v;
    data.processing.filtereddata.velocity.NS.PGV = PGV1;
    data.processing.filtereddata.displacement.NS.waveform = d;
    data.processing.filtereddata.displacement.NS.PGD = PGD1;
    
        %% Now do Arias intensity
    [Iaval, D5D95, D5D75, rate_arias, Ianorm] =  Arias(time_orig, NS_orig, fs);
    data.processing.filtereddata.acceleration.NS.arias.intensity = Iaval;
    data.processing.filtereddata.acceleration.NS.arias.D595 = D5D95;
    data.processing.filtereddata.acceleration.NS.arias.D575 = D5D75;
    data.processing.filtereddata.acceleration.NS.arias.rate = rate_arias;
    data.processing.filtereddata.acceleration.NS.arias.normalized = Ianorm;
    
   %% Now do response spectra
    [y, displacement, velocity, period, max1, max1V, max1D, per_maxA, per_maxV,...
        per_maxD, resp_02_secs, resp_02_secsV, resp_02_secsD ,resp_05_secs ,...
        resp_05_secsV,resp_05_secsD, resp_1_secs, resp_1_secsV,resp_1_secsD,...
        resp_2_secs,resp_2_secsV, resp_2_secsD, resp_5_secs,resp_5_secsV,...
        resp_5_secsD] = Response_Spectra(NS_orig, fs, 5);
    data.processing.filtereddata.acceleration.NS.spectra.waveform = y;
    data.processing.filtereddata.acceleration.NS.spectra.max = max1;
    data.processing.filtereddata.acceleration.NS.spectra.max_period = per_maxA;
    data.processing.filtereddata.acceleration.NS.spectra.seconds_02 = resp_02_secs;
    data.processing.filtereddata.acceleration.NS.spectra.seconds_05 = resp_05_secs;
    data.processing.filtereddata.acceleration.NS.spectra.seconds_1 = resp_1_secs;
    data.processing.filtereddata.acceleration.NS.spectra.seconds_2 = resp_2_secs;
    data.processing.filtereddata.acceleration.NS.spectra.seconds_5 = resp_5_secs;
    
    data.processing.filtereddata.velocity.NS.spectra.waveform = velocity;
    data.processing.filtereddata.velocity.NS.spectra.max = max1V;
    data.processing.filtereddata.velocity.NS.spectra.max_period = per_maxV;
    data.processing.filtereddata.velocity.NS.spectra.seconds_02 = resp_02_secsV;
    data.processing.filtereddata.velocity.NS.spectra.seconds_05 = resp_05_secsV;
    data.processing.filtereddata.velocity.NS.spectra.seconds_1 = resp_1_secsV;
    data.processing.filtereddata.velocity.NS.spectra.seconds_2 = resp_2_secsV;
    data.processing.filtereddata.velocity.NS.spectra.seconds_5 = resp_5_secsV;
    
    data.processing.filtereddata.displacement.NS.spectra.waveform = displacement;
    data.processing.filtereddata.displacement.NS.spectra.max = max1D;
    data.processing.filtereddata.displacement.NS.spectra.max_period = per_maxD;
    data.processing.filtereddata.displacement.NS.spectra.seconds_02 = resp_02_secsD;
    data.processing.filtereddata.displacement.NS.spectra.seconds_05 = resp_05_secsD;
    data.processing.filtereddata.displacement.NS.spectra.seconds_1 = resp_1_secsD;
    data.processing.filtereddata.displacement.NS.spectra.seconds_2 = resp_2_secsD;
    data.processing.filtereddata.displacement.NS.spectra.seconds_5 = resp_5_secsD;
    
    %% Now do magnitude response
    [NS_mag_a, NS_mag_smooth_a] =  magresp_ros(NS, fs, len, c1, c2);
    data.processing.filtereddata.acceleration.NS.mag_resps.unfiltered = NS_mag_a;
    data.processing.filtereddata.acceleration.NS.mag_resps.smooth = NS_mag_smooth_a;
    [NS_mag, NS_mag_smooth] =  magresp_ros(v', fs, len, c1, c2);
    data.processing.filtereddata.velocity.NS.mag_resps.unfiltered = NS_mag;
    data.processing.filtereddata.velocity.NS.mag_resps.smooth = NS_mag_smooth;
    [NS_mag, NS_mag_smooth] =  magresp_ros(d', fs, len, c1, c2);
    data.processing.filtereddata.displacement.NS.mag_resps.unfiltered = NS_mag;
    data.processing.filtereddata.displacement.NS.mag_resps.smooth = NS_mag_smooth;
    %% Compute the frequency - axis
    N = length(NS);
    fax_binsN = (0 : N-1); %samples in NS component
    fax_HzN1 = fax_binsN*fs/N; %frequency axis NS (Hz)
    N_2 = floor(N/2); %half magnitude spectrum
    fax_HzN = fax_HzN1(1 : N_2);
    data.processing.filtereddata.freq_vec = fax_HzN';
    
    %% East - West
    %% Now integrate the waveform for velocity and displacement
    [PGA1, PGV1, PGD1, v, d] =  waveform_integrate(EW, fs);
    data.processing.filtereddata.acceleration.EW.PGA = PGA1;
    data.processing.filtereddata.velocity.EW.waveform = v;
    data.processing.filtereddata.velocity.EW.PGV = PGV1;
    data.processing.filtereddata.displacement.EW.waveform = d;
    data.processing.filtereddata.displacement.EW.PGD = PGD1;
    
    %% Now do Arias intensity
    [Iaval, D5D95, D5D75, rate_arias, Ianorm] =  Arias(time_orig, EW_orig, fs);
    data.processing.filtereddata.acceleration.EW.arias.intensity = Iaval;
    data.processing.filtereddata.acceleration.EW.arias.D595 = D5D95;
    data.processing.filtereddata.acceleration.EW.arias.D575 = D5D75;
    data.processing.filtereddata.acceleration.EW.arias.rate = rate_arias;
    data.processing.filtereddata.acceleration.EW.arias.normalized = Ianorm;
    
    %% Now do response spectra
   [y, displacement, velocity, period, max1, max1V, max1D, per_maxA, per_maxV,...
        per_maxD, resp_02_secs, resp_02_secsV, resp_02_secsD ,resp_05_secs ,...
        resp_05_secsV,resp_05_secsD, resp_1_secs, resp_1_secsV,resp_1_secsD,...
        resp_2_secs,resp_2_secsV, resp_2_secsD, resp_5_secs,resp_5_secsV,...
        resp_5_secsD] = Response_Spectra(EW_orig, fs, 5);
    data.processing.filtereddata.acceleration.EW.spectra.waveform = y;
    data.processing.filtereddata.acceleration.EW.spectra.max = max1;
    data.processing.filtereddata.acceleration.EW.spectra.max_period = per_maxA;
    data.processing.filtereddata.acceleration.EW.spectra.seconds_02 = resp_02_secs;
    data.processing.filtereddata.acceleration.EW.spectra.seconds_05 = resp_05_secs;
    data.processing.filtereddata.acceleration.EW.spectra.seconds_1 = resp_1_secs;
    data.processing.filtereddata.acceleration.EW.spectra.seconds_2 = resp_2_secs;
    data.processing.filtereddata.acceleration.EW.spectra.seconds_5 = resp_5_secs;
    
    data.processing.filtereddata.velocity.EW.spectra.waveform = velocity;
    data.processing.filtereddata.velocity.EW.spectra.max = max1V;
    data.processing.filtereddata.velocity.EW.spectra.max_period = per_maxV;
    data.processing.filtereddata.velocity.EW.spectra.seconds_02 = resp_02_secsV;
    data.processing.filtereddata.velocity.EW.spectra.seconds_05 = resp_05_secsV;
    data.processing.filtereddata.velocity.EW.spectra.seconds_1 = resp_1_secsV;
    data.processing.filtereddata.velocity.EW.spectra.seconds_2 = resp_2_secsV;
    data.processing.filtereddata.velocity.EW.spectra.seconds_5 = resp_5_secsV;
    
    data.processing.filtereddata.displacement.EW.spectra.waveform = displacement;
    data.processing.filtereddata.displacement.EW.spectra.max = max1D;
    data.processing.filtereddata.displacement.EW.spectra.max_period = per_maxD;
    data.processing.filtereddata.displacement.EW.spectra.seconds_02 = resp_02_secsD;
    data.processing.filtereddata.displacement.EW.spectra.seconds_05 = resp_05_secsD;
    data.processing.filtereddata.displacement.EW.spectra.seconds_1 = resp_1_secsD;
    data.processing.filtereddata.displacement.EW.spectra.seconds_2 = resp_2_secsD;
    data.processing.filtereddata.displacement.EW.spectra.seconds_5 = resp_5_secsD;

    %% Now do magnitude response
    [EW_mag_a, EW_mag_smooth_a] =  magresp_ros(EW, fs, len, c1, c2);
    data.processing.filtereddata.acceleration.EW.mag_resps.unfiltered = EW_mag_a;
    data.processing.filtereddata.acceleration.EW.mag_resps.smooth = EW_mag_smooth_a;
    [EW_mag, EW_mag_smooth] =  magresp_ros(v', fs, len, c1, c2);
    data.processing.filtereddata.velocity.EW.mag_resps.unfiltered = EW_mag;
    data.processing.filtereddata.velocity.EW.mag_resps.smooth = EW_mag_smooth;
    [EW_mag, EW_mag_smooth] =  magresp_ros(d', fs, len, c1, c2);
    data.processing.filtereddata.displacement.EW.mag_resps.unfiltered = EW_mag;
    data.processing.filtereddata.displacement.EW.mag_resps.smooth = EW_mag_smooth;
    %% Vertical
    %% Integrate the waveform
    [PGA1, PGV1, PGD1, v, d] =  waveform_integrate(V, fs);
    data.processing.filtereddata.acceleration.V.PGA = PGA1;
    data.processing.filtereddata.velocity.V.waveform = v;
    data.processing.filtereddata.velocity.V.PGV = PGV1;
    data.processing.filtereddata.displacement.V.waveform = d;
    data.processing.filtereddata.displacement.V.PGD = PGD1;
    
    %% Now do Arias intensity    
    [Iaval, D5D95, D5D75, rate_arias, Ianorm] =  Arias(time_orig, V_orig, fs);
    data.processing.filtereddata.acceleration.V.arias.intensity = Iaval;
    data.processing.filtereddata.acceleration.V.arias.D595 = D5D95;
    data.processing.filtereddata.acceleration.V.arias.D575 = D5D75;
    data.processing.filtereddata.acceleration.V.arias.rate = rate_arias;
    data.processing.filtereddata.acceleration.V.arias.normalized = Ianorm;
    
    %% Now response spectra
    [y, displacement, velocity, period, max1, max1V, max1D, per_maxA, per_maxV,...
        per_maxD, resp_02_secs, resp_02_secsV, resp_02_secsD ,resp_05_secs ,...
        resp_05_secsV,resp_05_secsD, resp_1_secs, resp_1_secsV,resp_1_secsD,...
        resp_2_secs,resp_2_secsV, resp_2_secsD, resp_5_secs,resp_5_secsV,...
        resp_5_secsD] = Response_Spectra(V_orig, fs, 5);
    data.processing.filtereddata.acceleration.V.spectra.waveform = y;
    data.processing.filtereddata.acceleration.V.spectra.max = max1;
    data.processing.filtereddata.acceleration.V.spectra.max_period = per_maxA;
    data.processing.filtereddata.acceleration.V.spectra.seconds_02 = resp_02_secs;
    data.processing.filtereddata.acceleration.V.spectra.seconds_05 = resp_05_secs;
    data.processing.filtereddata.acceleration.V.spectra.seconds_1 = resp_1_secs;
    data.processing.filtereddata.acceleration.V.spectra.seconds_2 = resp_2_secs;
    data.processing.filtereddata.acceleration.V.spectra.seconds_5 = resp_5_secs;
    
    data.processing.filtereddata.velocity.V.spectra.waveform = velocity;
    data.processing.filtereddata.velocity.V.spectra.max = max1V;
    data.processing.filtereddata.velocity.V.spectra.max_period = per_maxV;
    data.processing.filtereddata.velocity.V.spectra.seconds_02 = resp_02_secsV;
    data.processing.filtereddata.velocity.V.spectra.seconds_05 = resp_05_secsV;
    data.processing.filtereddata.velocity.V.spectra.seconds_1 = resp_1_secsV;
    data.processing.filtereddata.velocity.V.spectra.seconds_2 = resp_2_secsV;
    data.processing.filtereddata.velocity.V.spectra.seconds_5 = resp_5_secsV;
    
    data.processing.filtereddata.displacement.V.spectra.waveform = displacement;
    data.processing.filtereddata.displacement.V.spectra.max = max1D;
    data.processing.filtereddata.displacement.V.spectra.max_period = per_maxD;
    data.processing.filtereddata.displacement.V.spectra.seconds_02 = resp_02_secsD;
    data.processing.filtereddata.displacement.V.spectra.seconds_05 = resp_05_secsD;
    data.processing.filtereddata.displacement.V.spectra.seconds_1 = resp_1_secsD;
    data.processing.filtereddata.displacement.V.spectra.seconds_2 = resp_2_secsD;
    data.processing.filtereddata.displacement.V.spectra.seconds_5 = resp_5_secsD;    
    
    %% Now do magnitude response
    [V_mag_a, V_mag_smooth_a] =  magresp_ros(V, fs, len, c1, c2);
    data.processing.filtereddata.acceleration.V.mag_resps.unfiltered = V_mag_a;
    data.processing.filtereddata.acceleration.V.mag_resps.smooth = V_mag_smooth_a;
    [V_mag, V_mag_smooth] =  magresp_ros(v', fs, len, c1, c2);
    data.processing.filtereddata.velocity.V.mag_resps.unfiltered = V_mag;
    data.processing.filtereddata.velocity.V.mag_resps.smooth = V_mag_smooth;
    [V_mag, V_mag_smooth] =  magresp_ros(d', fs, len, c1, c2);
    data.processing.filtereddata.displacement.V.mag_resps.unfiltered = V_mag;
    data.processing.filtereddata.displacement.V.mag_resps.smooth = V_mag_smooth;    
    
    %% Rotated
    %5 Integrate the waveform
    [PGA, PGV, PGD, v, d] =  waveform_integrate(rot, fs);
    data.processing.filtereddata.acceleration.rotated.PGA = PGA;
    data.processing.filtereddata.velocity.rotated.waveform = v;
    data.processing.filtereddata.velocity.rotated.PGV = PGV;
    data.processing.filtereddata.displacement.rotated.waveform = d;
    data.processing.filtereddata.displacement.rotated.PGD = PGD;
    
    %% Now do Arias intensity
    [Iaval, D5D95, D5D75, rate_arias, Ianorm] =  Arias(time_orig, rot_orig, fs);
    data.processing.filtereddata.acceleration.rotated.arias.intensity = Iaval;
    data.processing.filtereddata.acceleration.rotated.arias.D595 = D5D95;
    data.processing.filtereddata.acceleration.rotated.arias.D575 = D5D75;
    data.processing.filtereddata.acceleration.rotated.arias.rate = rate_arias;
    data.processing.filtereddata.acceleration.rotated.arias.normalized = Ianorm;

    %% Now do response spectra
    [y, displacement, velocity, period, max1, max1V, max1D, per_maxA, per_maxV,...
        per_maxD, resp_02_secs, resp_02_secsV, resp_02_secsD ,resp_05_secs ,...
        resp_05_secsV,resp_05_secsD, resp_1_secs, resp_1_secsV,resp_1_secsD,...
        resp_2_secs,resp_2_secsV, resp_2_secsD, resp_5_secs,resp_5_secsV,...
        resp_5_secsD] = Response_Spectra(rot_orig, fs, 5);
    data.processing.filtereddata.acceleration.rotated.spectra.waveform = y;
    data.processing.filtereddata.acceleration.rotated.spectra.max = max1;
    data.processing.filtereddata.acceleration.rotated.spectra.max_period = per_maxA;
    data.processing.filtereddata.acceleration.rotated.spectra.seconds_02 = resp_02_secs;
    data.processing.filtereddata.acceleration.rotated.spectra.seconds_05 = resp_05_secs;
    data.processing.filtereddata.acceleration.rotated.spectra.seconds_1 = resp_1_secs;
    data.processing.filtereddata.acceleration.rotated.spectra.seconds_2 = resp_2_secs;
    data.processing.filtereddata.acceleration.rotated.spectra.seconds_5 = resp_5_secs;
    
    data.processing.filtereddata.velocity.rotated.spectra.waveform = velocity;
    data.processing.filtereddata.velocity.rotated.spectra.max = max1V;
    data.processing.filtereddata.velocity.rotated.spectra.max_period = per_maxV;
    data.processing.filtereddata.velocity.rotated.spectra.seconds_02 = resp_02_secsV;
    data.processing.filtereddata.velocity.rotated.spectra.seconds_05 = resp_05_secsV;
    data.processing.filtereddata.velocity.rotated.spectra.seconds_1 = resp_1_secsV;
    data.processing.filtereddata.velocity.rotated.spectra.seconds_2 = resp_2_secsV;
    data.processing.filtereddata.velocity.rotated.spectra.seconds_5 = resp_5_secsV;
    
    data.processing.filtereddata.displacement.rotated.spectra.waveform = displacement;
    data.processing.filtereddata.displacement.rotated.spectra.max = max1D;
    data.processing.filtereddata.displacement.rotated.spectra.max_period = per_maxD;
    data.processing.filtereddata.displacement.rotated.spectra.seconds_02 = resp_02_secsD;
    data.processing.filtereddata.displacement.rotated.spectra.seconds_05 = resp_05_secsD;
    data.processing.filtereddata.displacement.rotated.spectra.seconds_1 = resp_1_secsD;
    data.processing.filtereddata.displacement.rotated.spectra.seconds_2 = resp_2_secsD;
    data.processing.filtereddata.displacement.rotated.spectra.seconds_5 = resp_5_secsD;    
    
    data.processing.filtereddata.period_vec = period;
    
     %% Now do magnitude response
    [rot_mag_a, rot_mag_smooth_a] =  magresp_ros(rot, fs, len, c1, c2);
    data.processing.filtereddata.acceleration.rotated.mag_resps.unfiltered = rot_mag_a;
    data.processing.filtereddata.acceleration.rotated.mag_resps.smooth = rot_mag_smooth_a;
    [rot_mag, rot_mag_smooth] =  magresp_ros(v', fs, len, c1, c2);
    data.processing.filtereddata.velocity.rotated.mag_resps.unfiltered = rot_mag;
    data.processing.filtereddata.velocity.rotated.mag_resps.smooth = rot_mag_smooth;
    [rot_mag, rot_mag_smooth] =  magresp_ros(d', fs, len, c1, c2);
    data.processing.filtereddata.displacement.rotated.mag_resps.unfiltered = rot_mag;
    data.processing.filtereddata.displacement.rotated.mag_resps.smooth = rot_mag_smooth;  
    
    %% Now onto the HVSRs
    %% create upbound and lowbound values
    lowbound = 0.1;
    upbound = 5;
    [~, lowbound] = min(abs(fax_HzN - lowbound));
    [~, upbound] = min(abs(fax_HzN - upbound));
    %% first do the complex combination
    H = NS + 1i*EW;   
    [comp_mag, comp_mag_smooth] =  magresp_ros(H, fs, len, c1, c2);
    data.processing.filtereddata.acceleration.complex.mag_resps.unfiltered = comp_mag;
    data.processing.filtereddata.acceleration.complex.mag_resps.smooth = comp_mag_smooth;
    [H_V] = HV(comp_mag,V_mag_a);
    data.processing.filtereddata.acceleration.complex.HVSR.unfilt = H_V;
    [H_V] = HV(comp_mag_smooth,V_mag_smooth_a);
    data.processing.filtereddata.acceleration.complex.HVSR.smooth.HV = H_V;
    [M, I] = max(H_V(lowbound:upbound));
    freq_max = fax_HzN(I + lowbound);
    data.processing.filtereddata.acceleration.complex.HVSR.smooth.Amp = M;
    data.processing.filtereddata.acceleration.complex.HVSR.smooth.fn = freq_max;
    data.processing.filtereddata.acceleration.complex.HVSR.smooth.fn_index = I + lowbound;
    %% NS
    [H_V] = HV(NS_mag_a,V_mag_a);
    data.processing.filtereddata.acceleration.NS.HVSR.unfilt = H_V;
    [H_V] = HV(NS_mag_smooth_a,V_mag_smooth_a);
    data.processing.filtereddata.acceleration.NS.HVSR.smooth.HV = H_V;
    [M, I] = max(H_V(lowbound:upbound));
    freq_max = fax_HzN(I + lowbound);
    data.processing.filtereddata.acceleration.NS.HVSR.smooth.Amp = M;
    data.processing.filtereddata.acceleration.NS.HVSR.smooth.fn = freq_max;
    data.processing.filtereddata.acceleration.NS.HVSR.smooth.fn_index = I + lowbound;
    %% EW
    [H_V] = HV(EW_mag_a,V_mag_a);
    data.processing.filtereddata.acceleration.EW.HVSR.unfilt = H_V;
    [H_V] = HV(EW_mag_smooth_a,V_mag_smooth_a);
    data.processing.filtereddata.acceleration.EW.HVSR.smooth.HV = H_V;
    [M, I] = max(H_V(lowbound:upbound));
    freq_max = fax_HzN(I + lowbound);
    data.processing.filtereddata.acceleration.EW.HVSR.smooth.Amp = M;
    data.processing.filtereddata.acceleration.EW.HVSR.smooth.fn = freq_max;
    data.processing.filtereddata.acceleration.EW.HVSR.smooth.fn_index = I + lowbound;
    
    %% rot
    [H_V] = HV(rot_mag_a,V_mag_a);
    data.processing.filtereddata.acceleration.rotated.HVSR.unfilt = H_V;
    [H_V] = HV(rot_mag_smooth_a,V_mag_smooth_a);
    data.processing.filtereddata.acceleration.rotated.HVSR.smooth.HV = H_V;
    [M, I] = max(H_V(lowbound:upbound));
    freq_max = fax_HzN(I + lowbound);
    data.processing.filtereddata.acceleration.rotated.HVSR.smooth.Amp = M;
    data.processing.filtereddata.acceleration.rotated.HVSR.smooth.fn = freq_max;
    data.processing.filtereddata.acceleration.rotated.HVSR.smooth.fn_index = I + lowbound;
end