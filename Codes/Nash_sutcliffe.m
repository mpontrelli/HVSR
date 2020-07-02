
%% Compute Nash sutcliffe on all vs. azimuth
% Author: Marshall Pontrelli
% Date: 5/20/2020

% For HVSR taxonomy USGS proposal 2020
close all
clear all
all = xlsread('C:\Users\mpontr01\Box\2020_1_spring\Research\Proposals\HV_classification\Azimuths.xlsx', 2, 'B1:B59');
az = xlsread('C:\Users\mpontr01\Box\2020_1_spring\Research\Proposals\HV_classification\Azimuths.xlsx', 2, 'D1:D59');

%% now process
Num = sum((az - all).^2);
Den = sum((az - mean(az)).^2);
E = 100*(1 - Num/Den);