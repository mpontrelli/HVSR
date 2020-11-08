%% Iris_pull
% Pull data off the IRIS database and filter. Gets read by GUI


%% Author: Marshall Pontrelli
% Date: 9/18/2020

%% Notes
% filter code: Broadband = 1, Short-period = 2, Long-period = 3
close all
clear all

%% filter inputs (Hz)
broad_low = 0.1;
broad_high = 19;
short_low = 5;
short_high = 19;
long_low = 0.1;
long_high = 5;
%% start infinite loop
i = 3;
while i > 2
    t_now = datetime;
    t_now = datestr(t_now);
    hour_min = t_now(13:16);
    if strcmp(hour_min, '00:0')
        disp(t_now)
    else
        data = cell(1,2);
    for j = 1:5
    %% Pull Weston data
    disp(j)
    if j ==1
    Station = 'WES';
    Network = 'NE';
   
    d = date;
    time1 = strcat(d,{' '}, '00:00:00');
    time1 = time1{1};
    time2 = datetime('now') + hours(3) + minutes(50);

    % Broadband
    LowCorner = broad_low;
    HighCorner = broad_high;
    [sampletimes,trace1,data{j}{3,1}, samplerate, sensitivity, sensunits] ...
        = getDMCData(time1,time2,Station,Network,'HHZ', LowCorner, HighCorner);
    [~,~,data{j}{1,1}, ~, ~, ~] ...
        = getDMCData(time1,time2,Station,Network,'HH1', LowCorner, HighCorner);
    [~,~,data{j}{2,1}, ~, ~, ~] ...
        = getDMCData(time1,time2,Station,Network,'HH2', LowCorner, HighCorner);
   
    
    % Short-period
    LowCorner = short_low;
    HighCorner = short_high;
    [~,~,data{j}{3,2}, ~, ~, ~] ...
        = getDMCData(time1,time2,Station,Network,'HHZ', LowCorner, HighCorner);
    [~,~,data{j}{1,2}, ~, ~, ~] ...
        = getDMCData(time1,time2,Station,Network,'HH1', LowCorner, HighCorner);
    [~,~,data{j}{2,2}, ~, ~, ~] ...
        = getDMCData(time1,time2,Station,Network,'HH2', LowCorner, HighCorner);
       
    % Long-period
    LowCorner = long_low;
    HighCorner = long_high;
    [~,~,data{j}{3, 3}, ~, ~, ~] ...
        = getDMCData(time1,time2,Station,Network,'HHZ', LowCorner, HighCorner);
    [~,~,data{j}{1,3}, ~, ~, ~] ...
        = getDMCData(time1,time2,Station,Network,'HH1', LowCorner, HighCorner);
    [sampletimes,trace1,data{j}{2,3}, samplerate, sensitivity, sensunits] ...
        = getDMCData(time1,time2,Station,Network,'HH2', LowCorner, HighCorner);
   data1 = data{j};
    parsave(strcat('C:\Users\',getenv('username'),'\Desktop\Data_pull\',Network,'_',Station,'.mat'),data1,...
        samplerate, sensitivity, sensunits)
    end
    %% Pull F64A
    if j==2
    Station = 'F64A';
    Network = 'N4';
   
    d = date;
    time1 = strcat(d,{' '}, '00:00:00');
    time1 = time1{1};
    time2 = datetime('now') + hours(3) + minutes(50);

    % Broadband
    LowCorner = broad_low;
    HighCorner = broad_high;
    [sampletimes,trace1,data{j}{3,1}, samplerate, sensitivity, sensunits] ...
        = getDMCData(time1,time2,Station,Network,'HHZ', LowCorner, HighCorner);
    [~,~,data{j}{1,1}, ~, ~, ~] ...
        = getDMCData(time1,time2,Station,Network,'HH1', LowCorner, HighCorner);
    [~,~,data{j}{2,1}, ~, ~, ~] ...
        = getDMCData(time1,time2,Station,Network,'HH2', LowCorner, HighCorner);
   
    
    % Short-period
    LowCorner = short_low;
    HighCorner = short_high;
    [~,~,data{j}{3,2}, ~, ~, ~] ...
        = getDMCData(time1,time2,Station,Network,'HHZ', LowCorner, HighCorner);
    [~,~,data{j}{1,2}, ~, ~, ~] ...
        = getDMCData(time1,time2,Station,Network,'HH1', LowCorner, HighCorner);
    [~,~,data{j}{2,2}, ~, ~, ~] ...
        = getDMCData(time1,time2,Station,Network,'HH2', LowCorner, HighCorner);
       
    % Long-period
    LowCorner = long_low;
    HighCorner = long_high;
    [~,~,data{j}{3, 3}, ~, ~, ~] ...
        = getDMCData(time1,time2,Station,Network,'HHZ', LowCorner, HighCorner);
    [~,~,data{j}{1,3}, ~, ~, ~] ...
        = getDMCData(time1,time2,Station,Network,'HH1', LowCorner, HighCorner);
    [sampletimes,trace1,data{j}{2,3}, samplerate, sensitivity, sensunits] ...
        = getDMCData(time1,time2,Station,Network,'HH2', LowCorner, HighCorner);
   data1 = data{j};
    parsave(strcat('C:\Users\',getenv('username'),'\Desktop\Data_pull\',Network,'_',Station,'.mat'),data1,...
        samplerate, sensitivity, sensunits)
    end
    %% Pull M65A
    if j==3
    Station = 'M65A';
    Network = 'TA';
   
    d = date;
    time1 = strcat(d,{' '}, '00:00:00');
    time1 = time1{1};
    time2 = datetime('now') + hours(3) + minutes(50);

    % Broadband
    LowCorner = broad_low;
    HighCorner = broad_high;
    [sampletimes,trace1,data{j}{3,1}, samplerate, sensitivity, sensunits] ...
        = getDMCData(time1,time2,Station,Network,'HHZ', LowCorner, HighCorner);
    [~,~,data{j}{1,1}, ~, ~, ~] ...
        = getDMCData(time1,time2,Station,Network,'HHE', LowCorner, HighCorner);
    [~,~,data{j}{2,1}, ~, ~, ~] ...
        = getDMCData(time1,time2,Station,Network,'HHN', LowCorner, HighCorner);
   
    
    % Short-period
    LowCorner = short_low;
    HighCorner = short_high;
    [~,~,data{j}{3,2}, ~, ~, ~] ...
        = getDMCData(time1,time2,Station,Network,'HHZ', LowCorner, HighCorner);
    [~,~,data{j}{1,2}, ~, ~, ~] ...
        = getDMCData(time1,time2,Station,Network,'HHE', LowCorner, HighCorner);
    [~,~,data{j}{2,2}, ~, ~, ~] ...
        = getDMCData(time1,time2,Station,Network,'HHN', LowCorner, HighCorner);
       
    % Long-period
    LowCorner = long_low;
    HighCorner = long_high;
    [~,~,data{j}{3, 3}, ~, ~, ~] ...
        = getDMCData(time1,time2,Station,Network,'HHZ', LowCorner, HighCorner);
    [~,~,data{j}{1,3}, ~, ~, ~] ...
        = getDMCData(time1,time2,Station,Network,'HHE', LowCorner, HighCorner);
    [sampletimes,trace1,data{j}{2,3}, samplerate, sensitivity, sensunits] ...
        = getDMCData(time1,time2,Station,Network,'HHN', LowCorner, HighCorner);
   data1 = data{j};
    parsave(strcat('C:\Users\',getenv('username'),'\Desktop\Data_pull\',Network,'_',Station,'.mat'),data1,...
        samplerate, sensitivity, sensunits)
    end
    %% Pull PKME
    if j==4
    Station = 'PKME';
    Network = 'US';
   
    d = date;
    time1 = strcat(d,{' '}, '00:00:00');
    time1 = time1{1};
    time2 = datetime('now') + hours(3) + minutes(50);

    % Broadband
    LowCorner = broad_low;
    HighCorner = broad_high;
    [sampletimes,trace1,data{j}{3,1}, samplerate, sensitivity, sensunits] ...
        = getDMCData(time1,time2,Station,Network,'BHZ', LowCorner, HighCorner);
    [~,~,data{j}{1,1}, ~, ~, ~] ...
        = getDMCData(time1,time2,Station,Network,'BH1', LowCorner, HighCorner);
    [~,~,data{j}{2,1}, ~, ~, ~] ...
        = getDMCData(time1,time2,Station,Network,'BH2', LowCorner, HighCorner);
   
    
    % Short-period
    LowCorner = short_low;
    HighCorner = short_high;
    [~,~,data{j}{3,2}, ~, ~, ~] ...
        = getDMCData(time1,time2,Station,Network,'BHZ', LowCorner, HighCorner);
    [~,~,data{j}{1,2}, ~, ~, ~] ...
        = getDMCData(time1,time2,Station,Network,'BH1', LowCorner, HighCorner);
    [~,~,data{j}{2,2}, ~, ~, ~] ...
        = getDMCData(time1,time2,Station,Network,'BH2', LowCorner, HighCorner);
       
    % Long-period
    LowCorner = long_low;
    HighCorner = long_high;
    [~,~,data{j}{3, 3}, ~, ~, ~] ...
        = getDMCData(time1,time2,Station,Network,'BHZ', LowCorner, HighCorner);
    [~,~,data{j}{1,3}, ~, ~, ~] ...
        = getDMCData(time1,time2,Station,Network,'BH1', LowCorner, HighCorner);
    [sampletimes,trace1,data{j}{2,3}, samplerate, sensitivity, sensunits] ...
        = getDMCData(time1,time2,Station,Network,'BH2', LowCorner, HighCorner);
   data1 = data{j};
    parsave(strcat('C:\Users\',getenv('username'),'\Desktop\Data_pull\',Network,'_',Station,'.mat'),data1,...
        samplerate, sensitivity, sensunits)
    end
    %% Pull G62A
    if j==5
    Station = 'G62A';
    Network = 'N4';
   
    d = date;
    time1 = strcat(d,{' '}, '00:00:00');
    time1 = time1{1};
    time2 = datetime('now') + hours(3) + minutes(50);

    % Broadband
    LowCorner = broad_low;
    HighCorner = broad_high;
    [sampletimes,trace1,data{j}{3,1}, samplerate, sensitivity, sensunits] ...
        = getDMCData(time1,time2,Station,Network,'HHZ', LowCorner, HighCorner);
    [~,~,data{j}{1,1}, ~, ~, ~] ...
        = getDMCData(time1,time2,Station,Network,'HH1', LowCorner, HighCorner);
    [~,~,data{j}{2,1}, ~, ~, ~] ...
        = getDMCData(time1,time2,Station,Network,'HH2', LowCorner, HighCorner);
   
    
    % Short-period
    LowCorner = short_low;
    HighCorner = short_high;
    [~,~,data{j}{3,2}, ~, ~, ~] ...
        = getDMCData(time1,time2,Station,Network,'HHZ', LowCorner, HighCorner);
    [~,~,data{j}{1,2}, ~, ~, ~] ...
        = getDMCData(time1,time2,Station,Network,'HH1', LowCorner, HighCorner);
    [~,~,data{j}{2,2}, ~, ~, ~] ...
        = getDMCData(time1,time2,Station,Network,'HH2', LowCorner, HighCorner);
       
    % Long-period
    LowCorner = long_low;
    HighCorner = long_high;
    [~,~,data{j}{3, 3}, ~, ~, ~] ...
        = getDMCData(time1,time2,Station,Network,'HHZ', LowCorner, HighCorner);
    [~,~,data{j}{1,3}, ~, ~, ~] ...
        = getDMCData(time1,time2,Station,Network,'HH1', LowCorner, HighCorner);
    [sampletimes,trace1,data{j}{2,3}, samplerate, sensitivity, sensunits] ...
        = getDMCData(time1,time2,Station,Network,'HH2', LowCorner, HighCorner);
   data1 = data{j};
    parsave(strcat('C:\Users\',getenv('username'),'\Desktop\Data_pull\',Network,'_',Station,'.mat'),data1,...
        samplerate, sensitivity, sensunits)
    end
    end
    pause(300)
    end
    
end