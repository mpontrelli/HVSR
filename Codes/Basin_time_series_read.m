% Basin_time_series_read looks at the time series data from the 2017 Puebla
% Mexico event. We are looking to discern between 1D and 3D basin effects
% and started with this event. 


% Author - Marshall Pontrelli
% Date - 1/13/2020
%% start
close all
clear all

UC44 = open('C:\Users\mpontr01\Box\2020_1_spring\SSA\Time domain\Data_mat_files\UC44.mat');
CE32 = open('C:\Users\mpontr01\Box\2020_1_spring\SSA\Time domain\Data_mat_files\CE32.mat');
UC44NS = UC44.xNS;
CE32NS = CE32.xNS;
CE32NS = CE32NS(1:length(UC44NS));

figure
subplot(2,1,1)
plot(CE32NS)
subplot(2,1,2)
plot(UC44NS)