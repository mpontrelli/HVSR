%% Read_reftek
% Reads reftek file and outputs .sac file

    % INPUTS
    
    % Filename - Reftek filename
    
%% Author: Marshall Pontrelli
% Date: 6/23/2021

path = cd;
filename = strcat(path, '\', '152241006_000B8D08');
command = strcat('pas2sac',{' '},filename,{' '},path);
command = command{1};

cd(strcat('C:\Users\',getenv('username'),'\Desktop\HVSR\Codes\reftek'));
[~,~] = system(command);