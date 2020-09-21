%% Iris_pull
% Pull data off the IRIS database and filter. Gets read by GUI


%% Author: Marshall Pontrelli
% Date: 9/18/2020

%% Notes
% filter code: Broadband = 1, Short-period = 2, Long-period = 3
close all
clear all
%% start infinite loop
i = 3;
while i > 2
    t_now = datetime;
    t_now = datestr(t_now);
    hour_min = t_now(13:16);
    if strcmp(hour_min, '00:0')
        disp(t_now)
    else
        
    % Pull Weston data
    Station = 'WES';
    Network = 'NE';
   
    d = date;
    time1 = strcat(d,{' '}, '00:00:00');
    time1 = time1{1};
    time2 = datetime('now') + hours(3) + minutes(50);

    % Broadband
    LowCorner = 0.1;
    HighCorner = 49;
    [sampletimes,trace1,data{3,1}, samplerate, sensitivity, sensunits] ...
        = getDMCData(time1,time2,Station,Network,'HHZ', LowCorner, HighCorner);
    [~,~,data{1,1}, ~, ~, ~] ...
        = getDMCData(time1,time2,Station,Network,'HH1', LowCorner, HighCorner);
    [~,~,data{2,1}, ~, ~, ~] ...
        = getDMCData(time1,time2,Station,Network,'HH2', LowCorner, HighCorner);
   
    
    % Short-period
    LowCorner = 5;
    HighCorner = 49;
    [~,~,data{3,2}, ~, ~, ~] ...
        = getDMCData(time1,time2,Station,Network,'HHZ', LowCorner, HighCorner);
    [~,~,data{1,2}, ~, ~, ~] ...
        = getDMCData(time1,time2,Station,Network,'HH1', LowCorner, HighCorner);
    [~,~,data{2,2}, ~, ~, ~] ...
        = getDMCData(time1,time2,Station,Network,'HH2', LowCorner, HighCorner);
       
    % Long-period
    LowCorner = 0.1;
    HighCorner = 5;
    [~,~,data{3, 3}, ~, ~, ~] ...
        = getDMCData(time1,time2,Station,Network,'HHZ', LowCorner, HighCorner);
    [~,~,data{1,3}, ~, ~, ~] ...
        = getDMCData(time1,time2,Station,Network,'HH1', LowCorner, HighCorner);
    [sampletimes,trace1,data{2,3}, samplerate, sensitivity, sensunits] ...
        = getDMCData(time1,time2,Station,Network,'HH2', LowCorner, HighCorner);
   
    save(strcat('C:\Users\',getenv('username'),'\Desktop\Data_pull\','WES.mat'),'data',...
        'samplerate', 'sensitivity', 'sensunits')
    
    end
    pause(300)
end