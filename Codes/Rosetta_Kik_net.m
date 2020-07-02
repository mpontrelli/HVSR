%% Read the Kik-net data files and create structure 

% Author: Marshall Pontrelli
% Date: Winter 2020, alot of work done on 5/19/2020


function [data] = Rosetta_Kik_net(folder,event)
codepath = 'C:\Users\mpontr01\Desktop\HVSR\Codes';
filename = strcat(folder,event,'.EW1');
%% EW1
[fid3,~] = fopen(filename,'r');

% Now read the file
%
%% Start by getting the event origin time
line = fgetl(fid3);
b = length(line);
data.meta.event.origin_time = line(19:b);

%% Now get event latitude
line = fgetl(fid3);
s = string(line);
data.meta.event.lat = str2double(regexp(s,'\d+(\.\d+)?|\.\d+','match'));

%% Now event longitude
line = fgetl(fid3);
s = string(line);
data.meta.event.lon = str2double(regexp(s,'\d+(\.\d+)?|\.\d+','match'));

%% Now depth
line = fgetl(fid3);
s = string(line);
data.meta.event.depth = str2double(regexp(s,'\d+(\.\d+)?|\.\d+','match'));

%% Now magnitude
line = fgetl(fid3);
s = string(line);
data.meta.event.mag = str2double(regexp(s,'\d+(\.\d+)?|\.\d+','match'));

%% Now station name
line = fgetl(fid3);
b = length(line);
data.meta.station.name = line(19:b);

%% Now station lat
line = fgetl(fid3);
s = string(line);
data.meta.station.lat = str2double(regexp(s,'\d+(\.\d+)?|\.\d+','match'));

%% Now station lon
line = fgetl(fid3);
s = string(line);
data.meta.station.lon = str2double(regexp(s,'\d+(\.\d+)?|\.\d+','match'));

%% Now station elevation
line = fgetl(fid3);
s = string(line);
data.meta.station.elev = str2double(regexp(s,'\d+(\.\d+)?|\.\d+','match'));

%% Now record time
line = fgetl(fid3);
b = length(line);
data.meta.record.start_time = line(19:b);

%% Now sampling frequency
line = fgetl(fid3);
s = string(line);
fs = str2double(regexp(s,'\d+(\.\d+)?|\.\d+','match'));
data.meta.instrument.fs = fs;

%% Now record duration
line = fgetl(fid3);
s = string(line);
data.meta.record.duration = str2double(regexp(s,'\d+(\.\d+)?|\.\d+','match'));

%% Now dir, I don't know what that is, but let's save it
line = fgetl(fid3);
s = string(line);
data.meta.dir = str2double(regexp(s,'\d+(\.\d+)?|\.\d+','match'));

%% Now scale factor
line = fgetl(fid3);
s = string(line);
q = str2double(regexp(s,'\d+(\.\d+)?|\.\d+','match'));
scale_factor = q(1)/q(2);
data.meta.instrument.EW1.scale_factor = scale_factor;

%% Now max acceleration, in gals
line = fgetl(fid3);
s = string(line);
data.meta.record.max_accel = str2double(regexp(s,'\d+(\.\d+)?|\.\d+','match'));

%% Now onto the EW1 data

deliminator='';
%Data starts in row 109, column 1
R = 17; %Row data start text file
C = 0; %column data start text file
data1 = dlmread(filename,deliminator,R,C); %retrieve data
data1 = data1';
EW1 = data1(:);
EW1 = scale_factor*(EW1 - mean(EW1));
data.processing.rawdata.EW1 = EW1;

%% Now onto the EW2 data
filename = strcat(folder,event,'.EW2');
[fid3,~] = fopen(filename,'r');

% All info is the same, but we need to get the scale factor, skip 13 lines
for jjj = 1:14
    line = fgetl(fid3);
end
s = string(line);
q = str2double(regexp(s,'\d+(\.\d+)?|\.\d+','match'));
scale_factor = q(1)/q(2);
data.meta.instrument.EW2.scale_factor = scale_factor;

%% Now get the EW2 data

deliminator='';
%Data starts in row 109, column 1
R = 17; %Row data start text file
C = 0; %column data start text file
data1 = dlmread(filename,deliminator,R,C); %retrieve data
data1 = data1';
EW2 = data1(:);
EW2 = scale_factor*(EW2 - mean(EW2));
data.processing.rawdata.EW2 = EW2;

%% Now onto the NS1 data
filename = strcat(folder,event,'.NS1');
[fid3,~] = fopen(filename,'r');

% All info is the same, but we need to get the scale factor, skip 13 lines
for jjj = 1:14
    line = fgetl(fid3);
end
s = string(line);
q = str2double(regexp(s,'\d+(\.\d+)?|\.\d+','match'));
scale_factor = q(1)/q(2);
data.meta.instrument.NS1.scale_factor = scale_factor;

%% Now get the NS1 data

deliminator='';
%Data starts in row 109, column 1
R = 17; %Row data start text file
C = 0; %column data start text file
data1 = dlmread(filename,deliminator,R,C); %retrieve data
data1 = data1';
NS1 = data1(:);
NS1 = scale_factor*(NS1 - mean(NS1));
data.processing.rawdata.NS1 = NS1;

%% Now onto the NS2 data
filename = strcat(folder,event,'.NS2');
[fid3,~] = fopen(filename,'r');

% All info is the same, but we need to get the scale factor, skip 13 lines
for jjj = 1:14
    line = fgetl(fid3);
end
s = string(line);
q = str2double(regexp(s,'\d+(\.\d+)?|\.\d+','match'));
scale_factor = q(1)/q(2);
data.meta.instrument.NS2.scale_factor = scale_factor;

%% Now get the NS2 data

deliminator='';
%Data starts in row 109, column 1
R = 17; %Row data start text file
C = 0; %column data start text file
data1 = dlmread(filename,deliminator,R,C); %retrieve data
data1 = data1';
NS2 = data1(:);
NS2 = scale_factor*(NS2 - mean(NS2));
data.processing.rawdata.NS2 = NS2;

%% Now onto the UD1 data
filename = strcat(folder,event,'.UD1');
[fid3,~] = fopen(filename,'r');

% All info is the same, but we need to get the scale factor, skip 13 lines
for jjj = 1:14
    line = fgetl(fid3);
end
s = string(line);
q = str2double(regexp(s,'\d+(\.\d+)?|\.\d+','match'));
scale_factor = q(1)/q(2);
data.meta.instrument.UD1.scale_factor = scale_factor;

%% Now get the UD1 data

deliminator='';
%Data starts in row 109, column 1
R = 17; %Row data start text file
C = 0; %column data start text file
data1 = dlmread(filename,deliminator,R,C); %retrieve data
data1 = data1';
UD1 = data1(:);
UD1 = scale_factor*(UD1 - mean(UD1));
data.processing.rawdata.UD1 = UD1;

%% Now onto the UD2 data
filename = strcat(folder,event,'.UD2');
[fid3,~] = fopen(filename,'r');

% All info is the same, but we need to get the scale factor, skip 13 lines
for jjj = 1:14
    line = fgetl(fid3);
end
s = string(line);
q = str2double(regexp(s,'\d+(\.\d+)?|\.\d+','match'));
scale_factor = q(1)/q(2);
data.meta.instrument.UD2.scale_factor = scale_factor;

%% Now get the UD1 data

deliminator='';
%Data starts in row 109, column 1
R = 17; %Row data start text file
C = 0; %column data start text file
data1 = dlmread(filename,deliminator,R,C); %retrieve data
data1 = data1';
UD2 = data1(:);
UD2 = scale_factor*(UD2 - mean(UD2));
data.processing.rawdata.UD2 = UD2;

%% OK, now all the data is loaded, let's process it

% Declare some filter parameters  
LowCorner = 0.01;
HighCorner = fs/2 - 1;
Npoles = 4;
cd(codepath)
[EW1] = (Butter2(EW1, fs, 'LowCorner', LowCorner, 'HighCorner', HighCorner, 'Npoles', Npoles));
[EW2] = (Butter2(EW2, fs, 'LowCorner', LowCorner, 'HighCorner', HighCorner, 'Npoles', Npoles));
[NS1] = (Butter2(NS1, fs, 'LowCorner', LowCorner, 'HighCorner', HighCorner, 'Npoles', Npoles));
[NS2] = (Butter2(NS2, fs, 'LowCorner', LowCorner, 'HighCorner', HighCorner, 'Npoles', Npoles));
[UD1] = (Butter2(UD1, fs, 'LowCorner', LowCorner, 'HighCorner', HighCorner, 'Npoles', Npoles));
[UD2] = (Butter2(UD2, fs, 'LowCorner', LowCorner, 'HighCorner', HighCorner, 'Npoles', Npoles));

%% Rotate the horizontals
for i = 1:length(EW1)
    rot1(i) = NS1(i)*cos(pi/4) + EW1(i)*sin(pi/4);
end
rot1 = rot1';

for i = 1:length(EW2)
    rot2(i) = NS2(i)*cos(pi/4) + EW2(i)*sin(pi/4);
end
rot2 = rot2';

%% Now make a time vector
time = (1:length(EW1))/fs;

%% Now save
NS1_orig = NS1;
NS2_orig = NS2;
EW1_orig = EW1;
EW2_orig = EW2;
UD1_orig = UD1;
UD2_orig = UD2;
rot1_orig = rot1';
rot2_orig = rot2';
time_orig = time';

data.processing.filtereddata.acceleration.NS1.waveform_orig = NS1_orig; %North_South_Component
data.processing.filtereddata.acceleration.NS2.waveform_orig = NS2_orig; %North_South_Component
data.processing.filtereddata.acceleration.EW1.waveform_orig = EW1_orig; %East_West_Component
data.processing.filtereddata.acceleration.EW2.waveform_orig = EW2_orig; %East_West_Component
data.processing.filtereddata.acceleration.UD1.waveform_orig = UD1_orig; %Vertical_Component
data.processing.filtereddata.acceleration.UD2.waveform_orig = UD2_orig; %Vertical_Component
data.processing.filtereddata.acceleration.rotated1.waveform_orig = rot1_orig; %rotated
data.processing.filtereddata.acceleration.rotated2.waveform_orig = rot2_orig; %rotated
data.processing.filtereddata.time_orig = time_orig; %rotated

%% Now zero pad
if fs == 100
    len = length(EW1);
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
    len = length(EW1);
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
NS1 = vertcat(c1',NS1,c2');
NS2 = vertcat(c1',NS2,c2');
EW1 = vertcat(c1',EW1,c2');
EW2 = vertcat(c1',EW2,c2');
UD1 = vertcat(c1',UD1,c2');
UD2 = vertcat(c1',UD2,c2');
rot1 = vertcat(c1',rot1,c2');
rot2 = vertcat(c1',rot2,c2');

%% Now save zero padded
data.processing.filtereddata.acceleration.NS1.waveform = NS1; %North_South_Component
data.processing.filtereddata.acceleration.NS2.waveform = NS2; %North_South_Component
data.processing.filtereddata.acceleration.EW1.waveform = EW1; %East_West_Component
data.processing.filtereddata.acceleration.EW2.waveform = EW2; %East_West_Component
data.processing.filtereddata.acceleration.UD1.waveform = UD1; %Vertical_Component
data.processing.filtereddata.acceleration.UD2.waveform = UD2; %Vertical_Component
data.processing.filtereddata.acceleration.rotated1.waveform = rot1; %rotated
data.processing.filtereddata.acceleration.rotated1.waveform = rot2; %rotated
data.processing.filtereddata.time = time;
data.processing.filtereddata.filter.lowcorner = LowCorner;
data.processing.filtereddata.filter.highcorner = HighCorner;
data.processing.filtereddata.filter.npoles = Npoles;
    
%% Start with NS1
%% Now integrate the waveforms
[PGA1, PGV1, PGD1, v, d] =  waveform_integrate(NS1, fs);
data.processing.filtereddata.acceleration.NS1.PGA = PGA1;
data.processing.filtereddata.velocity.NS1.waveform = v;
data.processing.filtereddata.velocity.NS1.PGV = PGV1;
data.processing.filtereddata.displacement.NS1.waveform = d;
data.processing.filtereddata.displacement.NS1.PGD = PGD1;

%% Now do Arias intensity
[IaX, D5, D75, D95, Iaval, D5D95, D5D75, rate_arias, Ianorm] =  Arias(time_orig, NS1_orig, fs);
data.processing.filtereddata.acceleration.NS1.arias.no_norm = IaX;
data.processing.filtereddata.acceleration.NS1.arias.intensity = Iaval;
data.processing.filtereddata.acceleration.NS1.arias.D5 = D5;
data.processing.filtereddata.acceleration.NS1.arias.D75 = D75;
data.processing.filtereddata.acceleration.NS1.arias.D95 = D95;
data.processing.filtereddata.acceleration.NS1.arias.D595 = D5D95;
data.processing.filtereddata.acceleration.NS1.arias.D575 = D5D75;
data.processing.filtereddata.acceleration.NS1.arias.rate = rate_arias;
data.processing.filtereddata.acceleration.NS1.arias.normalized = Ianorm;

%% Now do response spectra
[y, displacement, velocity, period, max1, max1V, max1D, per_maxA, per_maxV,...
    per_maxD, resp_02_secs, resp_02_secsV, resp_02_secsD ,resp_05_secs ,...
    resp_05_secsV,resp_05_secsD, resp_1_secs, resp_1_secsV,resp_1_secsD,...
    resp_2_secs,resp_2_secsV, resp_2_secsD, resp_5_secs,resp_5_secsV,...
    resp_5_secsD] = Response_Spectra(NS1_orig, fs, 5);
data.processing.filtereddata.acceleration.NS1.spectra.waveform = y;
data.processing.filtereddata.acceleration.NS1.spectra.max = max1;
data.processing.filtereddata.acceleration.NS1.spectra.max_period = per_maxA;
data.processing.filtereddata.acceleration.NS1.spectra.seconds_02 = resp_02_secs;
data.processing.filtereddata.acceleration.NS1.spectra.seconds_05 = resp_05_secs;
data.processing.filtereddata.acceleration.NS1.spectra.seconds_1 = resp_1_secs;
data.processing.filtereddata.acceleration.NS1.spectra.seconds_2 = resp_2_secs;
data.processing.filtereddata.acceleration.NS1.spectra.seconds_5 = resp_5_secs;

data.processing.filtereddata.velocity.NS1.spectra.waveform = velocity;
data.processing.filtereddata.velocity.NS1.spectra.max = max1V;
data.processing.filtereddata.velocity.NS1.spectra.max_period = per_maxV;
data.processing.filtereddata.velocity.NS1.spectra.seconds_02 = resp_02_secsV;
data.processing.filtereddata.velocity.NS1.spectra.seconds_05 = resp_05_secsV;
data.processing.filtereddata.velocity.NS1.spectra.seconds_1 = resp_1_secsV;
data.processing.filtereddata.velocity.NS1.spectra.seconds_2 = resp_2_secsV;
data.processing.filtereddata.velocity.NS1.spectra.seconds_5 = resp_5_secsV;
    
data.processing.filtereddata.displacement.NS1.spectra.waveform = displacement;
data.processing.filtereddata.displacement.NS1.spectra.max = max1D;
data.processing.filtereddata.displacement.NS1.spectra.max_period = per_maxD;
data.processing.filtereddata.displacement.NS1.spectra.seconds_02 = resp_02_secsD;
data.processing.filtereddata.displacement.NS1.spectra.seconds_05 = resp_05_secsD;
data.processing.filtereddata.displacement.NS1.spectra.seconds_1 = resp_1_secsD;
data.processing.filtereddata.displacement.NS1.spectra.seconds_2 = resp_2_secsD;
data.processing.filtereddata.displacement.NS1.spectra.seconds_5 = resp_5_secsD;

%% Now do magnitude response
[NS1_mag_a, NS1_mag_smooth_a] =  magresp_ros(NS1, fs, len, c1, c2);
data.processing.filtereddata.acceleration.NS1.mag_resps.unfiltered = NS1_mag_a;
data.processing.filtereddata.acceleration.NS1.mag_resps.smooth = NS1_mag_smooth_a;
[NS1_mag, NS1_mag_smooth] =  magresp_ros(v', fs, len, c1, c2);
data.processing.filtereddata.velocity.NS1.mag_resps.unfiltered = NS1_mag;
data.processing.filtereddata.velocity.NS1.mag_resps.smooth = NS1_mag_smooth;
[NS1_mag, NS1_mag_smooth] =  magresp_ros(d', fs, len, c1, c2);
data.processing.filtereddata.displacement.NS1.mag_resps.unfiltered = NS1_mag;
data.processing.filtereddata.displacement.NS1.mag_resps.smooth = NS1_mag_smooth;

%% Compute the frequency - axis
N = length(NS1);
fax_binsN = (0 : N-1); %samples in NS component
fax_HzN1 = fax_binsN*fs/N; %frequency axis NS (Hz)
N_2 = floor(N/2); %half magnitude spectrum
fax_HzN = fax_HzN1(1 : N_2);
data.processing.filtereddata.freq_vec = fax_HzN';

%% Now NS2
%% Now integrate the waveforms
[PGA1, PGV1, PGD1, v, d] =  waveform_integrate(NS2, fs);
data.processing.filtereddata.acceleration.NS2.PGA = PGA1;
data.processing.filtereddata.velocity.NS2.waveform = v;
data.processing.filtereddata.velocity.NS2.PGV = PGV1;
data.processing.filtereddata.displacement.NS2.waveform = d;
data.processing.filtereddata.displacement.NS2.PGD = PGD1;

%% Now do Arias intensity
[IaX, D5, D75, D95, Iaval, D5D95, D5D75, rate_arias, Ianorm] =  Arias(time_orig, NS2_orig, fs);
data.processing.filtereddata.acceleration.NS2.arias.no_norm = IaX;
data.processing.filtereddata.acceleration.NS2.arias.intensity = Iaval;
data.processing.filtereddata.acceleration.NS2.arias.D5 = D5;
data.processing.filtereddata.acceleration.NS2.arias.D75 = D75;
data.processing.filtereddata.acceleration.NS2.arias.D95 = D95;
data.processing.filtereddata.acceleration.NS2.arias.D595 = D5D95;
data.processing.filtereddata.acceleration.NS2.arias.D575 = D5D75;
data.processing.filtereddata.acceleration.NS2.arias.rate = rate_arias;
data.processing.filtereddata.acceleration.NS2.arias.normalized = Ianorm;

%% Now do response spectra
[y, displacement, velocity, period, max1, max1V, max1D, per_maxA, per_maxV,...
    per_maxD, resp_02_secs, resp_02_secsV, resp_02_secsD ,resp_05_secs ,...
    resp_05_secsV,resp_05_secsD, resp_1_secs, resp_1_secsV,resp_1_secsD,...
    resp_2_secs,resp_2_secsV, resp_2_secsD, resp_5_secs,resp_5_secsV,...
    resp_5_secsD] = Response_Spectra(NS2_orig, fs, 5);
data.processing.filtereddata.acceleration.NS2.spectra.waveform = y;
data.processing.filtereddata.acceleration.NS2.spectra.max = max1;
data.processing.filtereddata.acceleration.NS2.spectra.max_period = per_maxA;
data.processing.filtereddata.acceleration.NS2.spectra.seconds_02 = resp_02_secs;
data.processing.filtereddata.acceleration.NS2.spectra.seconds_05 = resp_05_secs;
data.processing.filtereddata.acceleration.NS2.spectra.seconds_1 = resp_1_secs;
data.processing.filtereddata.acceleration.NS2.spectra.seconds_2 = resp_2_secs;
data.processing.filtereddata.acceleration.NS2.spectra.seconds_5 = resp_5_secs;

data.processing.filtereddata.velocity.NS2.spectra.waveform = velocity;
data.processing.filtereddata.velocity.NS2.spectra.max = max1V;
data.processing.filtereddata.velocity.NS2.spectra.max_period = per_maxV;
data.processing.filtereddata.velocity.NS2.spectra.seconds_02 = resp_02_secsV;
data.processing.filtereddata.velocity.NS2.spectra.seconds_05 = resp_05_secsV;
data.processing.filtereddata.velocity.NS2.spectra.seconds_1 = resp_1_secsV;
data.processing.filtereddata.velocity.NS2.spectra.seconds_2 = resp_2_secsV;
data.processing.filtereddata.velocity.NS2.spectra.seconds_5 = resp_5_secsV;
    
data.processing.filtereddata.displacement.NS2.spectra.waveform = displacement;
data.processing.filtereddata.displacement.NS2.spectra.max = max1D;
data.processing.filtereddata.displacement.NS2.spectra.max_period = per_maxD;
data.processing.filtereddata.displacement.NS2.spectra.seconds_02 = resp_02_secsD;
data.processing.filtereddata.displacement.NS2.spectra.seconds_05 = resp_05_secsD;
data.processing.filtereddata.displacement.NS2.spectra.seconds_1 = resp_1_secsD;
data.processing.filtereddata.displacement.NS2.spectra.seconds_2 = resp_2_secsD;
data.processing.filtereddata.displacement.NS2.spectra.seconds_5 = resp_5_secsD;

%% Now do magnitude response
[NS2_mag_a, NS2_mag_smooth_a] =  magresp_ros(NS2, fs, len, c1, c2);
data.processing.filtereddata.acceleration.NS2.mag_resps.unfiltered = NS2_mag_a;
data.processing.filtereddata.acceleration.NS2.mag_resps.smooth = NS2_mag_smooth_a;
[NS2_mag, NS2_mag_smooth] =  magresp_ros(v', fs, len, c1, c2);
data.processing.filtereddata.velocity.NS2.mag_resps.unfiltered = NS2_mag;
data.processing.filtereddata.velocity.NS2.mag_resps.smooth = NS2_mag_smooth;
[NS2_mag, NS2_mag_smooth] =  magresp_ros(d', fs, len, c1, c2);
data.processing.filtereddata.displacement.NS2.mag_resps.unfiltered = NS2_mag;
data.processing.filtereddata.displacement.NS2.mag_resps.smooth = NS2_mag_smooth;

    
%% Now EW1
%% Now integrate the waveforms
[PGA1, PGV1, PGD1, v, d] =  waveform_integrate(EW1, fs);
data.processing.filtereddata.acceleration.EW1.PGA = PGA1;
data.processing.filtereddata.velocity.EW1.waveform = v;
data.processing.filtereddata.velocity.EW1.PGV = PGV1;
data.processing.filtereddata.displacement.EW1.waveform = d;
data.processing.filtereddata.displacement.EW1.PGD = PGD1;

%% Now do Arias intensity
[IaX, D5, D75, D95, Iaval, D5D95, D5D75, rate_arias, Ianorm] =  Arias(time_orig, EW1_orig, fs);
data.processing.filtereddata.acceleration.EW1.arias.no_norm = IaX;
data.processing.filtereddata.acceleration.EW1.arias.intensity = Iaval;
data.processing.filtereddata.acceleration.EW1.arias.D5 = D5;
data.processing.filtereddata.acceleration.EW1.arias.D75 = D75;
data.processing.filtereddata.acceleration.EW1.arias.D95 = D95;
data.processing.filtereddata.acceleration.EW1.arias.D595 = D5D95;
data.processing.filtereddata.acceleration.EW1.arias.D575 = D5D75;
data.processing.filtereddata.acceleration.EW1.arias.rate = rate_arias;
data.processing.filtereddata.acceleration.EW1.arias.normalized = Ianorm;

%% Now do response spectra
[y, displacement, velocity, period, max1, max1V, max1D, per_maxA, per_maxV,...
    per_maxD, resp_02_secs, resp_02_secsV, resp_02_secsD ,resp_05_secs ,...
    resp_05_secsV,resp_05_secsD, resp_1_secs, resp_1_secsV,resp_1_secsD,...
    resp_2_secs,resp_2_secsV, resp_2_secsD, resp_5_secs,resp_5_secsV,...
    resp_5_secsD] = Response_Spectra(EW1_orig, fs, 5);
data.processing.filtereddata.acceleration.EW1.spectra.waveform = y;
data.processing.filtereddata.acceleration.EW1.spectra.max = max1;
data.processing.filtereddata.acceleration.EW1.spectra.max_period = per_maxA;
data.processing.filtereddata.acceleration.EW1.spectra.seconds_02 = resp_02_secs;
data.processing.filtereddata.acceleration.EW1.spectra.seconds_05 = resp_05_secs;
data.processing.filtereddata.acceleration.EW1.spectra.seconds_1 = resp_1_secs;
data.processing.filtereddata.acceleration.EW1.spectra.seconds_2 = resp_2_secs;
data.processing.filtereddata.acceleration.EW1.spectra.seconds_5 = resp_5_secs;

data.processing.filtereddata.velocity.EW1.spectra.waveform = velocity;
data.processing.filtereddata.velocity.EW1.spectra.max = max1V;
data.processing.filtereddata.velocity.EW1.spectra.max_period = per_maxV;
data.processing.filtereddata.velocity.EW1.spectra.seconds_02 = resp_02_secsV;
data.processing.filtereddata.velocity.EW1.spectra.seconds_05 = resp_05_secsV;
data.processing.filtereddata.velocity.EW1.spectra.seconds_1 = resp_1_secsV;
data.processing.filtereddata.velocity.EW1.spectra.seconds_2 = resp_2_secsV;
data.processing.filtereddata.velocity.EW1.spectra.seconds_5 = resp_5_secsV;
    
data.processing.filtereddata.displacement.EW1.spectra.waveform = displacement;
data.processing.filtereddata.displacement.EW1.spectra.max = max1D;
data.processing.filtereddata.displacement.EW1.spectra.max_period = per_maxD;
data.processing.filtereddata.displacement.EW1.spectra.seconds_02 = resp_02_secsD;
data.processing.filtereddata.displacement.EW1.spectra.seconds_05 = resp_05_secsD;
data.processing.filtereddata.displacement.EW1.spectra.seconds_1 = resp_1_secsD;
data.processing.filtereddata.displacement.EW1.spectra.seconds_2 = resp_2_secsD;
data.processing.filtereddata.displacement.EW1.spectra.seconds_5 = resp_5_secsD;

%% Now do magnitude response
[EW1_mag_a, EW1_mag_smooth_a] =  magresp_ros(EW1, fs, len, c1, c2);
data.processing.filtereddata.acceleration.EW1.mag_resps.unfiltered = EW1_mag_a;
data.processing.filtereddata.acceleration.EW1.mag_resps.smooth = EW1_mag_smooth_a;
[EW1_mag, EW1_mag_smooth] =  magresp_ros(v', fs, len, c1, c2);
data.processing.filtereddata.velocity.EW1.mag_resps.unfiltered = EW1_mag;
data.processing.filtereddata.velocity.EW1.mag_resps.smooth = EW1_mag_smooth;
[EW1_mag, EW1_mag_smooth] =  magresp_ros(d', fs, len, c1, c2);
data.processing.filtereddata.displacement.EW1.mag_resps.unfiltered = EW1_mag;
data.processing.filtereddata.displacement.EW1.mag_resps.smooth = EW1_mag_smooth;

%% Now EW2
%% Now integrate the waveforms
[PGA1, PGV1, PGD1, v, d] =  waveform_integrate(EW2, fs);
data.processing.filtereddata.acceleration.EW2.PGA = PGA1;
data.processing.filtereddata.velocity.EW2.waveform = v;
data.processing.filtereddata.velocity.EW2.PGV = PGV1;
data.processing.filtereddata.displacement.EW2.waveform = d;
data.processing.filtereddata.displacement.EW2.PGD = PGD1;

%% Now do Arias intensity
[IaX, D5, D75, D95, Iaval, D5D95, D5D75, rate_arias, Ianorm] =  Arias(time_orig, EW2_orig, fs);
data.processing.filtereddata.acceleration.EW2.arias.no_norm = IaX;
data.processing.filtereddata.acceleration.EW2.arias.intensity = Iaval;
data.processing.filtereddata.acceleration.EW2.arias.D5 = D5;
data.processing.filtereddata.acceleration.EW2.arias.D75 = D75;
data.processing.filtereddata.acceleration.EW2.arias.D95 = D95;
data.processing.filtereddata.acceleration.EW2.arias.D595 = D5D95;
data.processing.filtereddata.acceleration.EW2.arias.D575 = D5D75;
data.processing.filtereddata.acceleration.EW2.arias.rate = rate_arias;
data.processing.filtereddata.acceleration.EW2.arias.normalized = Ianorm;

%% Now do response spectra
[y, displacement, velocity, period, max1, max1V, max1D, per_maxA, per_maxV,...
    per_maxD, resp_02_secs, resp_02_secsV, resp_02_secsD ,resp_05_secs ,...
    resp_05_secsV,resp_05_secsD, resp_1_secs, resp_1_secsV,resp_1_secsD,...
    resp_2_secs,resp_2_secsV, resp_2_secsD, resp_5_secs,resp_5_secsV,...
    resp_5_secsD] = Response_Spectra(EW2_orig, fs, 5);
data.processing.filtereddata.acceleration.EW2.spectra.waveform = y;
data.processing.filtereddata.acceleration.EW2.spectra.max = max1;
data.processing.filtereddata.acceleration.EW2.spectra.max_period = per_maxA;
data.processing.filtereddata.acceleration.EW2.spectra.seconds_02 = resp_02_secs;
data.processing.filtereddata.acceleration.EW2.spectra.seconds_05 = resp_05_secs;
data.processing.filtereddata.acceleration.EW2.spectra.seconds_1 = resp_1_secs;
data.processing.filtereddata.acceleration.EW2.spectra.seconds_2 = resp_2_secs;
data.processing.filtereddata.acceleration.EW2.spectra.seconds_5 = resp_5_secs;

data.processing.filtereddata.velocity.EW2.spectra.waveform = velocity;
data.processing.filtereddata.velocity.EW2.spectra.max = max1V;
data.processing.filtereddata.velocity.EW2.spectra.max_period = per_maxV;
data.processing.filtereddata.velocity.EW2.spectra.seconds_02 = resp_02_secsV;
data.processing.filtereddata.velocity.EW2.spectra.seconds_05 = resp_05_secsV;
data.processing.filtereddata.velocity.EW2.spectra.seconds_1 = resp_1_secsV;
data.processing.filtereddata.velocity.EW2.spectra.seconds_2 = resp_2_secsV;
data.processing.filtereddata.velocity.EW2.spectra.seconds_5 = resp_5_secsV;
    
data.processing.filtereddata.displacement.EW2.spectra.waveform = displacement;
data.processing.filtereddata.displacement.EW2.spectra.max = max1D;
data.processing.filtereddata.displacement.EW2.spectra.max_period = per_maxD;
data.processing.filtereddata.displacement.EW2.spectra.seconds_02 = resp_02_secsD;
data.processing.filtereddata.displacement.EW2.spectra.seconds_05 = resp_05_secsD;
data.processing.filtereddata.displacement.EW2.spectra.seconds_1 = resp_1_secsD;
data.processing.filtereddata.displacement.EW2.spectra.seconds_2 = resp_2_secsD;
data.processing.filtereddata.displacement.EW2.spectra.seconds_5 = resp_5_secsD;

%% Now do magnitude response
[EW2_mag_a, EW2_mag_smooth_a] =  magresp_ros(EW2, fs, len, c1, c2);
data.processing.filtereddata.acceleration.EW2.mag_resps.unfiltered = EW2_mag_a;
data.processing.filtereddata.acceleration.EW2.mag_resps.smooth = EW2_mag_smooth_a;
[EW2_mag, EW2_mag_smooth] =  magresp_ros(v', fs, len, c1, c2);
data.processing.filtereddata.velocity.EW2.mag_resps.unfiltered = EW2_mag;
data.processing.filtereddata.velocity.EW2.mag_resps.smooth = EW2_mag_smooth;
[EW2_mag, EW2_mag_smooth] =  magresp_ros(d', fs, len, c1, c2);
data.processing.filtereddata.displacement.EW2.mag_resps.unfiltered = EW2_mag;
data.processing.filtereddata.displacement.EW2.mag_resps.smooth = EW2_mag_smooth;


%% Now UD1
%% Now integrate the waveforms
[PGA1, PGV1, PGD1, v, d] =  waveform_integrate(UD1, fs);
data.processing.filtereddata.acceleration.UD1.PGA = PGA1;
data.processing.filtereddata.velocity.UD1.waveform = v;
data.processing.filtereddata.velocity.UD1.PGV = PGV1;
data.processing.filtereddata.displacement.UD1.waveform = d;
data.processing.filtereddata.displacement.UD1.PGD = PGD1;

%% Now do Arias intensity
[IaX, D5, D75, D95, Iaval, D5D95, D5D75, rate_arias, Ianorm] =  Arias(time_orig, UD1_orig, fs);
data.processing.filtereddata.acceleration.UD1.arias.no_norm = IaX;
data.processing.filtereddata.acceleration.UD1.arias.intensity = Iaval;
data.processing.filtereddata.acceleration.UD1.arias.D5 = D5;
data.processing.filtereddata.acceleration.UD1.arias.D75 = D75;
data.processing.filtereddata.acceleration.UD1.arias.D95 = D95;
data.processing.filtereddata.acceleration.UD1.arias.D595 = D5D95;
data.processing.filtereddata.acceleration.UD1.arias.D575 = D5D75;
data.processing.filtereddata.acceleration.UD1.arias.rate = rate_arias;
data.processing.filtereddata.acceleration.UD1.arias.normalized = Ianorm;

%% Now do response spectra
[y, displacement, velocity, period, max1, max1V, max1D, per_maxA, per_maxV,...
    per_maxD, resp_02_secs, resp_02_secsV, resp_02_secsD ,resp_05_secs ,...
    resp_05_secsV,resp_05_secsD, resp_1_secs, resp_1_secsV,resp_1_secsD,...
    resp_2_secs,resp_2_secsV, resp_2_secsD, resp_5_secs,resp_5_secsV,...
    resp_5_secsD] = Response_Spectra(UD1_orig, fs, 5);
data.processing.filtereddata.acceleration.UD1.spectra.waveform = y;
data.processing.filtereddata.acceleration.UD1.spectra.max = max1;
data.processing.filtereddata.acceleration.UD1.spectra.max_period = per_maxA;
data.processing.filtereddata.acceleration.UD1.spectra.seconds_02 = resp_02_secs;
data.processing.filtereddata.acceleration.UD1.spectra.seconds_05 = resp_05_secs;
data.processing.filtereddata.acceleration.UD1.spectra.seconds_1 = resp_1_secs;
data.processing.filtereddata.acceleration.UD1.spectra.seconds_2 = resp_2_secs;
data.processing.filtereddata.acceleration.UD1.spectra.seconds_5 = resp_5_secs;

data.processing.filtereddata.velocity.UD1.spectra.waveform = velocity;
data.processing.filtereddata.velocity.UD1.spectra.max = max1V;
data.processing.filtereddata.velocity.UD1.spectra.max_period = per_maxV;
data.processing.filtereddata.velocity.UD1.spectra.seconds_02 = resp_02_secsV;
data.processing.filtereddata.velocity.UD1.spectra.seconds_05 = resp_05_secsV;
data.processing.filtereddata.velocity.UD1.spectra.seconds_1 = resp_1_secsV;
data.processing.filtereddata.velocity.UD1.spectra.seconds_2 = resp_2_secsV;
data.processing.filtereddata.velocity.UD1.spectra.seconds_5 = resp_5_secsV;
    
data.processing.filtereddata.displacement.UD1.spectra.waveform = displacement;
data.processing.filtereddata.displacement.UD1.spectra.max = max1D;
data.processing.filtereddata.displacement.UD1.spectra.max_period = per_maxD;
data.processing.filtereddata.displacement.UD1.spectra.seconds_02 = resp_02_secsD;
data.processing.filtereddata.displacement.UD1.spectra.seconds_05 = resp_05_secsD;
data.processing.filtereddata.displacement.UD1.spectra.seconds_1 = resp_1_secsD;
data.processing.filtereddata.displacement.UD1.spectra.seconds_2 = resp_2_secsD;
data.processing.filtereddata.displacement.UD1.spectra.seconds_5 = resp_5_secsD;

%% Now do magnitude response
[UD1_mag_a, UD1_mag_smooth_a] =  magresp_ros(UD1, fs, len, c1, c2);
data.processing.filtereddata.acceleration.UD1.mag_resps.unfiltered = UD1_mag_a;
data.processing.filtereddata.acceleration.UD1.mag_resps.smooth = UD1_mag_smooth_a;
[UD1_mag, UD1_mag_smooth] =  magresp_ros(v', fs, len, c1, c2);
data.processing.filtereddata.velocity.UD1.mag_resps.unfiltered = UD1_mag;
data.processing.filtereddata.velocity.UD1.mag_resps.smooth = UD1_mag_smooth;
[UD1_mag, UD1_mag_smooth] =  magresp_ros(d', fs, len, c1, c2);
data.processing.filtereddata.displacement.UD1.mag_resps.unfiltered = UD1_mag;
data.processing.filtereddata.displacement.UD1.mag_resps.smooth = UD1_mag_smooth;

%% Now UD2
%% Now integrate the waveforms
[PGA1, PGV1, PGD1, v, d] =  waveform_integrate(UD2, fs);
data.processing.filtereddata.acceleration.UD2.PGA = PGA1;
data.processing.filtereddata.velocity.UD2.waveform = v;
data.processing.filtereddata.velocity.UD2.PGV = PGV1;
data.processing.filtereddata.displacement.UD2.waveform = d;
data.processing.filtereddata.displacement.UD2.PGD = PGD1;

%% Now do Arias intensity
[IaX, D5, D75, D95, Iaval, D5D95, D5D75, rate_arias, Ianorm] =  Arias(time_orig, UD2_orig, fs);
data.processing.filtereddata.acceleration.UD2.arias.no_norm = IaX;
data.processing.filtereddata.acceleration.UD2.arias.intensity = Iaval;
data.processing.filtereddata.acceleration.UD2.arias.D5 = D5;
data.processing.filtereddata.acceleration.UD2.arias.D75 = D75;
data.processing.filtereddata.acceleration.UD2.arias.D95 = D95;
data.processing.filtereddata.acceleration.UD2.arias.D595 = D5D95;
data.processing.filtereddata.acceleration.UD2.arias.D575 = D5D75;
data.processing.filtereddata.acceleration.UD2.arias.rate = rate_arias;
data.processing.filtereddata.acceleration.UD2.arias.normalized = Ianorm;

%% Now do response spectra
[y, displacement, velocity, period, max1, max1V, max1D, per_maxA, per_maxV,...
    per_maxD, resp_02_secs, resp_02_secsV, resp_02_secsD ,resp_05_secs ,...
    resp_05_secsV,resp_05_secsD, resp_1_secs, resp_1_secsV,resp_1_secsD,...
    resp_2_secs,resp_2_secsV, resp_2_secsD, resp_5_secs,resp_5_secsV,...
    resp_5_secsD] = Response_Spectra(UD2_orig, fs, 5);
data.processing.filtereddata.acceleration.UD2.spectra.waveform = y;
data.processing.filtereddata.acceleration.UD2.spectra.max = max1;
data.processing.filtereddata.acceleration.UD2.spectra.max_period = per_maxA;
data.processing.filtereddata.acceleration.UD2.spectra.seconds_02 = resp_02_secs;
data.processing.filtereddata.acceleration.UD2.spectra.seconds_05 = resp_05_secs;
data.processing.filtereddata.acceleration.UD2.spectra.seconds_1 = resp_1_secs;
data.processing.filtereddata.acceleration.UD2.spectra.seconds_2 = resp_2_secs;
data.processing.filtereddata.acceleration.UD2.spectra.seconds_5 = resp_5_secs;

data.processing.filtereddata.velocity.UD2.spectra.waveform = velocity;
data.processing.filtereddata.velocity.UD2.spectra.max = max1V;
data.processing.filtereddata.velocity.UD2.spectra.max_period = per_maxV;
data.processing.filtereddata.velocity.UD2.spectra.seconds_02 = resp_02_secsV;
data.processing.filtereddata.velocity.UD2.spectra.seconds_05 = resp_05_secsV;
data.processing.filtereddata.velocity.UD2.spectra.seconds_1 = resp_1_secsV;
data.processing.filtereddata.velocity.UD2.spectra.seconds_2 = resp_2_secsV;
data.processing.filtereddata.velocity.UD2.spectra.seconds_5 = resp_5_secsV;
    
data.processing.filtereddata.displacement.UD2.spectra.waveform = displacement;
data.processing.filtereddata.displacement.UD2.spectra.max = max1D;
data.processing.filtereddata.displacement.UD2.spectra.max_period = per_maxD;
data.processing.filtereddata.displacement.UD2.spectra.seconds_02 = resp_02_secsD;
data.processing.filtereddata.displacement.UD2.spectra.seconds_05 = resp_05_secsD;
data.processing.filtereddata.displacement.UD2.spectra.seconds_1 = resp_1_secsD;
data.processing.filtereddata.displacement.UD2.spectra.seconds_2 = resp_2_secsD;
data.processing.filtereddata.displacement.UD2.spectra.seconds_5 = resp_5_secsD;

%% Now do magnitude response
[UD2_mag_a, UD2_mag_smooth_a] =  magresp_ros(UD2, fs, len, c1, c2);
data.processing.filtereddata.acceleration.UD2.mag_resps.unfiltered = UD2_mag_a;
data.processing.filtereddata.acceleration.UD2.mag_resps.smooth = UD2_mag_smooth_a;
[UD2_mag, UD2_mag_smooth] =  magresp_ros(v', fs, len, c1, c2);
data.processing.filtereddata.velocity.UD2.mag_resps.unfiltered = UD2_mag;
data.processing.filtereddata.velocity.UD2.mag_resps.smooth = UD2_mag_smooth;
[UD2_mag, UD2_mag_smooth] =  magresp_ros(d', fs, len, c1, c2);
data.processing.filtereddata.displacement.UD2.mag_resps.unfiltered = UD2_mag;
data.processing.filtereddata.displacement.UD2.mag_resps.smooth = UD2_mag_smooth;

%% Now rot1
%% Now integrate the waveforms
[PGA1, PGV1, PGD1, v, d] =  waveform_integrate(rot1, fs);
data.processing.filtereddata.acceleration.rotated1.PGA = PGA1;
data.processing.filtereddata.velocity.rotated1.waveform = v;
data.processing.filtereddata.velocity.rotated1.PGV = PGV1;
data.processing.filtereddata.displacement.rotated1.waveform = d;
data.processing.filtereddata.displacement.rotated1.PGD = PGD1;

%% Now do Arias intensity
[IaX, D5, D75, D95, Iaval, D5D95, D5D75, rate_arias, Ianorm] =  Arias(time_orig, rot1_orig, fs);
data.processing.filtereddata.acceleration.rotated1.arias.no_norm = IaX;
data.processing.filtereddata.acceleration.rotated1.arias.intensity = Iaval;
data.processing.filtereddata.acceleration.rotated1.arias.D5 = D5;
data.processing.filtereddata.acceleration.rotated1.arias.D75 = D75;
data.processing.filtereddata.acceleration.rotated1.arias.D95 = D95;
data.processing.filtereddata.acceleration.rotated1.arias.D595 = D5D95;
data.processing.filtereddata.acceleration.rotated1.arias.D575 = D5D75;
data.processing.filtereddata.acceleration.rotated1.arias.rate = rate_arias;
data.processing.filtereddata.acceleration.rotated1.arias.normalized = Ianorm;

%% Now do response spectra
[y, displacement, velocity, period, max1, max1V, max1D, per_maxA, per_maxV,...
    per_maxD, resp_02_secs, resp_02_secsV, resp_02_secsD ,resp_05_secs ,...
    resp_05_secsV,resp_05_secsD, resp_1_secs, resp_1_secsV,resp_1_secsD,...
    resp_2_secs,resp_2_secsV, resp_2_secsD, resp_5_secs,resp_5_secsV,...
    resp_5_secsD] = Response_Spectra(rot1_orig, fs, 5);
data.processing.filtereddata.acceleration.rotated1.spectra.waveform = y;
data.processing.filtereddata.acceleration.rotated1.spectra.max = max1;
data.processing.filtereddata.acceleration.rotated1.spectra.max_period = per_maxA;
data.processing.filtereddata.acceleration.rotated1.spectra.seconds_02 = resp_02_secs;
data.processing.filtereddata.acceleration.rotated1.spectra.seconds_05 = resp_05_secs;
data.processing.filtereddata.acceleration.rotated1.spectra.seconds_1 = resp_1_secs;
data.processing.filtereddata.acceleration.rotated1.spectra.seconds_2 = resp_2_secs;
data.processing.filtereddata.acceleration.rotated1.spectra.seconds_5 = resp_5_secs;

data.processing.filtereddata.velocity.rotated1.spectra.waveform = velocity;
data.processing.filtereddata.velocity.rotated1.spectra.max = max1V;
data.processing.filtereddata.velocity.rotated1.spectra.max_period = per_maxV;
data.processing.filtereddata.velocity.rotated1.spectra.seconds_02 = resp_02_secsV;
data.processing.filtereddata.velocity.rotated1.spectra.seconds_05 = resp_05_secsV;
data.processing.filtereddata.velocity.rotated1.spectra.seconds_1 = resp_1_secsV;
data.processing.filtereddata.velocity.rotated1.spectra.seconds_2 = resp_2_secsV;
data.processing.filtereddata.velocity.rotated1.spectra.seconds_5 = resp_5_secsV;
    
data.processing.filtereddata.displacement.rotated1.spectra.waveform = displacement;
data.processing.filtereddata.displacement.rotated1.spectra.max = max1D;
data.processing.filtereddata.displacement.rotated1.spectra.max_period = per_maxD;
data.processing.filtereddata.displacement.rotated1.spectra.seconds_02 = resp_02_secsD;
data.processing.filtereddata.displacement.rotated1.spectra.seconds_05 = resp_05_secsD;
data.processing.filtereddata.displacement.rotated1.spectra.seconds_1 = resp_1_secsD;
data.processing.filtereddata.displacement.rotated1.spectra.seconds_2 = resp_2_secsD;
data.processing.filtereddata.displacement.rotated1.spectra.seconds_5 = resp_5_secsD;

%% Now do magnitude response
[rot1_mag_a, rot1_mag_smooth_a] =  magresp_ros(rot1, fs, len, c1, c2);
data.processing.filtereddata.acceleration.rotated1.mag_resps.unfiltered = rot1_mag_a;
data.processing.filtereddata.acceleration.rotated1.mag_resps.smooth = rot1_mag_smooth_a;
[rot1_mag, rot1_mag_smooth] =  magresp_ros(v', fs, len, c1, c2);
data.processing.filtereddata.velocity.rotated1.mag_resps.unfiltered = rot1_mag;
data.processing.filtereddata.velocity.rotated1.mag_resps.smooth = rot1_mag_smooth;
[rot1_mag, rot1_mag_smooth] =  magresp_ros(d', fs, len, c1, c2);
data.processing.filtereddata.displacement.rotated1.mag_resps.unfiltered = rot1_mag;
data.processing.filtereddata.displacement.rotated1.mag_resps.smooth = rot1_mag_smooth;

%% Now rot2
%% Now integrate the waveforms
[PGA1, PGV1, PGD1, v, d] =  waveform_integrate(rot2, fs);
data.processing.filtereddata.acceleration.rotated2.PGA = PGA1;
data.processing.filtereddata.velocity.rotated2.waveform = v;
data.processing.filtereddata.velocity.rotated2.PGV = PGV1;
data.processing.filtereddata.displacement.rotated2.waveform = d;
data.processing.filtereddata.displacement.rotated2.PGD = PGD1;

%% Now do Arias intensity
[IaX, D5, D75, D95, Iaval, D5D95, D5D75, rate_arias, Ianorm] =  Arias(time_orig, rot2_orig, fs);
data.processing.filtereddata.acceleration.rotated2.arias.no_norm = IaX;
data.processing.filtereddata.acceleration.rotated2.arias.intensity = Iaval;
data.processing.filtereddata.acceleration.rotated2.arias.D5 = D5;
data.processing.filtereddata.acceleration.rotated2.arias.D75 = D75;
data.processing.filtereddata.acceleration.rotated2.arias.D95 = D95;
data.processing.filtereddata.acceleration.rotated2.arias.D595 = D5D95;
data.processing.filtereddata.acceleration.rotated2.arias.D575 = D5D75;
data.processing.filtereddata.acceleration.rotated2.arias.rate = rate_arias;
data.processing.filtereddata.acceleration.rotated2.arias.normalized = Ianorm;

%% Now do response spectra
[y, displacement, velocity, period, max1, max1V, max1D, per_maxA, per_maxV,...
    per_maxD, resp_02_secs, resp_02_secsV, resp_02_secsD ,resp_05_secs ,...
    resp_05_secsV,resp_05_secsD, resp_1_secs, resp_1_secsV,resp_1_secsD,...
    resp_2_secs,resp_2_secsV, resp_2_secsD, resp_5_secs,resp_5_secsV,...
    resp_5_secsD] = Response_Spectra(rot2_orig, fs, 5);
data.processing.filtereddata.acceleration.rotated2.spectra.waveform = y;
data.processing.filtereddata.acceleration.rotated2.spectra.max = max1;
data.processing.filtereddata.acceleration.rotated2.spectra.max_period = per_maxA;
data.processing.filtereddata.acceleration.rotated2.spectra.seconds_02 = resp_02_secs;
data.processing.filtereddata.acceleration.rotated2.spectra.seconds_05 = resp_05_secs;
data.processing.filtereddata.acceleration.rotated2.spectra.seconds_1 = resp_1_secs;
data.processing.filtereddata.acceleration.rotated2.spectra.seconds_2 = resp_2_secs;
data.processing.filtereddata.acceleration.rotated2.spectra.seconds_5 = resp_5_secs;

data.processing.filtereddata.velocity.rotated2.spectra.waveform = velocity;
data.processing.filtereddata.velocity.rotated2.spectra.max = max1V;
data.processing.filtereddata.velocity.rotated2.spectra.max_period = per_maxV;
data.processing.filtereddata.velocity.rotated2.spectra.seconds_02 = resp_02_secsV;
data.processing.filtereddata.velocity.rotated2.spectra.seconds_05 = resp_05_secsV;
data.processing.filtereddata.velocity.rotated2.spectra.seconds_1 = resp_1_secsV;
data.processing.filtereddata.velocity.rotated2.spectra.seconds_2 = resp_2_secsV;
data.processing.filtereddata.velocity.rotated2.spectra.seconds_5 = resp_5_secsV;
    
data.processing.filtereddata.displacement.rotated2.spectra.waveform = displacement;
data.processing.filtereddata.displacement.rotated2.spectra.max = max1D;
data.processing.filtereddata.displacement.rotated2.spectra.max_period = per_maxD;
data.processing.filtereddata.displacement.rotated2.spectra.seconds_02 = resp_02_secsD;
data.processing.filtereddata.displacement.rotated2.spectra.seconds_05 = resp_05_secsD;
data.processing.filtereddata.displacement.rotated2.spectra.seconds_1 = resp_1_secsD;
data.processing.filtereddata.displacement.rotated2.spectra.seconds_2 = resp_2_secsD;
data.processing.filtereddata.displacement.rotated2.spectra.seconds_5 = resp_5_secsD;

%% Now do magnitude response
[rot2_mag_a, rot2_mag_smooth_a] =  magresp_ros(rot2, fs, len, c1, c2);
data.processing.filtereddata.acceleration.rotated2.mag_resps.unfiltered = rot2_mag_a;
data.processing.filtereddata.acceleration.rotated2.mag_resps.smooth = rot2_mag_smooth_a;
[rot2_mag, rot2_mag_smooth] =  magresp_ros(v', fs, len, c1, c2);
data.processing.filtereddata.velocity.rotated2.mag_resps.unfiltered = rot2_mag;
data.processing.filtereddata.velocity.rotated2.mag_resps.smooth = rot2_mag_smooth;
[rot2_mag, rot2_mag_smooth] =  magresp_ros(d', fs, len, c1, c2);
data.processing.filtereddata.displacement.rotated2.mag_resps.unfiltered = rot2_mag;
data.processing.filtereddata.displacement.rotated2.mag_resps.smooth = rot2_mag_smooth;


%% Now onto the HVSRs and BSRs
%% create upbound and lowbound values
lowbound = 0.1;
upbound = 10;
[~, lowbound] = min(abs(fax_HzN - lowbound));
[~, upbound] = min(abs(fax_HzN - upbound));
%% first do the complex combination

%% comp1
H = NS1 + 1i*EW1;   
[comp_mag1, comp_mag_smooth1] =  magresp_ros(H, fs, len, c1, c2);
data.processing.filtereddata.acceleration.complex1.mag_resps.unfiltered = comp_mag1;
data.processing.filtereddata.acceleration.complex1.mag_resps.smooth = comp_mag_smooth1;
[H_V] = HV(comp_mag1,UD1_mag_a);
data.processing.filtereddata.acceleration.complex1.HVSR.unfilt = H_V;
[H_V] = HV(comp_mag_smooth1,UD1_mag_smooth_a);
data.processing.filtereddata.acceleration.complex1.HVSR.smooth.HV = H_V;


%% comp2
H = NS2 + 1i*EW2;   
[comp_mag2, comp_mag_smooth2] =  magresp_ros(H, fs, len, c1, c2);
data.processing.filtereddata.acceleration.complex2.mag_resps.unfiltered = comp_mag2;
data.processing.filtereddata.acceleration.complex2.mag_resps.smooth = comp_mag_smooth2;
[H_V] = HV(comp_mag2,UD2_mag_a);
data.processing.filtereddata.acceleration.complex2.HVSR.unfilt = H_V;
[H_V] = HV(comp_mag_smooth2,UD2_mag_smooth_a);
data.processing.filtereddata.acceleration.complex2.HVSR.smooth.HV = H_V;

%% Complex BSR
BSR = comp_mag_smooth2./comp_mag_smooth1;
data.processing.filtereddata.acceleration.BSR.complex = BSR;

%% Now geometric mean
%% geo_mean1
H = sqrt(NS1_mag_a.*EW1_mag_a);
H_smooth = sqrt(NS1_mag_smooth_a.*EW1_mag_smooth_a);
data.processing.filtereddata.acceleration.geo_mean1.mag_resps.unfiltered = H;
data.processing.filtereddata.acceleration.geo_mean1.mag_resps.smooth = H_smooth;
[H_V] = HV(H,UD1_mag_a);
data.processing.filtereddata.acceleration.geo_mean1.HVSR.unfilt = H_V;
[H_V] = HV(H_smooth,UD1_mag_smooth_a);
data.processing.filtereddata.acceleration.geo_mean1.HVSR.smooth.HV = H_V;

%% geo_mean2
H2 = sqrt(NS2_mag_a.*EW2_mag_a);
H2_smooth = sqrt(NS2_mag_smooth_a.*EW2_mag_smooth_a);
data.processing.filtereddata.acceleration.geo_mean2.mag_resps.unfiltered = H2;
data.processing.filtereddata.acceleration.geo_mean2.mag_resps.smooth = H2_smooth;
[H_V] = HV(H2,UD2_mag_a);
data.processing.filtereddata.acceleration.geo_mean2.HVSR.unfilt = H_V;
[H_V] = HV(H2_smooth,UD2_mag_smooth_a);
data.processing.filtereddata.acceleration.geo_mean2.HVSR.smooth.HV = H_V;

%% Geomean BSR
BSR = H2_smooth./H_smooth;
data.processing.filtereddata.acceleration.BSR.geo_mean = BSR;

%% Now NS1
[H_V] = HV(NS1_mag_a,UD1_mag_a);
data.processing.filtereddata.acceleration.NS1.HVSR.unfilt = H_V;
[H_V] = HV(NS1_mag_smooth_a,UD1_mag_smooth_a);
data.processing.filtereddata.acceleration.NS1.HVSR.smooth.HV = H_V;

%% Now NS2
[H_V] = HV(NS2_mag_a,UD2_mag_a);
data.processing.filtereddata.acceleration.NS2.HVSR.unfilt = H_V;
[H_V] = HV(NS2_mag_smooth_a,UD2_mag_smooth_a);
data.processing.filtereddata.acceleration.NS2.HVSR.smooth.HV = H_V;

%% Now NS BSR
BSR = NS2_mag_smooth_a./NS1_mag_smooth_a;
data.processing.filtereddata.acceleration.BSR.NS = BSR;
%% Now EW1
[H_V] = HV(EW1_mag_a,UD1_mag_a);
data.processing.filtereddata.acceleration.EW1.HVSR.unfilt = H_V;
[H_V] = HV(EW1_mag_smooth_a,UD1_mag_smooth_a);
data.processing.filtereddata.acceleration.EW1.HVSR.smooth.HV = H_V;

%% Now EW2
[H_V] = HV(EW2_mag_a,EW2_mag_a);
data.processing.filtereddata.acceleration.EW2.HVSR.unfilt = H_V;
[H_V] = HV(EW2_mag_smooth_a,UD2_mag_smooth_a);
data.processing.filtereddata.acceleration.EW2.HVSR.smooth.HV = H_V;

%% Now EW BSR
BSR = EW2_mag_smooth_a./EW1_mag_smooth_a;
data.processing.filtereddata.acceleration.BSR.EW = BSR;

%% Now rot1
[H_V] = HV(rot1_mag_a,UD1_mag_a);
data.processing.filtereddata.acceleration.rotated1.HVSR.unfilt = H_V;
[H_V] = HV(rot1_mag_smooth_a,UD1_mag_smooth_a);
data.processing.filtereddata.acceleration.rotated1.HVSR.smooth.HV = H_V;

%% Now EW2
[H_V] = HV(rot2_mag_a,rot2_mag_a);
data.processing.filtereddata.acceleration.rotated2.HVSR.unfilt = H_V;
[H_V] = HV(rot2_mag_smooth_a,UD2_mag_smooth_a);
data.processing.filtereddata.acceleration.rotated2.HVSR.smooth.HV = H_V;

%% Now EW BSR
BSR = rot2_mag_smooth_a./rot1_mag_smooth_a;
data.processing.filtereddata.acceleration.BSR.rot= BSR;


fclose('all')
clear rot1
clear rot2
clear NS1
clear NS2
clear EW1
clear EW2
clear UD1
clear UD2
end