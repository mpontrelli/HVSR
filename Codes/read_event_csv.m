%% read_event_csv
  % Reads the data structure from create_mat_file and turns it into a .csv
  % file
%% Author: Marshall Pontrelli
% Date: 4/23/2020
%%
close all
clear all

%%
T = readtable('C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Final_tables\Events.csv');%, 'HeaderLines',3);

%% azimuth
azimuth = table2array((T(:,6)));
histogram(azimuth)

%% epicentral distance
epi_dist = table2array((T(:,7)));
histogram(epi_dist)

%% distance vs. azimuth
plot(epi_dist, azimuth, '*')


%% CE32 analysis
%% load data
CE32_az = T{1187:1230,{'event_azimuth'}}; % azimuth
CE32_depth = T{1187:1230,{'event_depth'}}; % depth
CE32_EW_freq = T{1187:1230,{'A_EW_HV_freq'}}; % EW_freq
CE32_EW_amp = T{1187:1230,{'A_EW_HV_amp'}}; % EW_amp
CE32_NS_freq = T{1187:1230,{'A_NS_HV_freq'}}; % NS_freq
CE32_NS_amp = T{1187:1230,{'A_NS_HV_amp'}}; % NS_amp
CE32_comp_amp = T{1187:1230,{'A_complex_HV_amp'}}; % Complex_amp
CE32_comp_freq = T{1187:1230,{'A_complex_HV_freq'}}; % Complex_freq
CE32_NS_PGA = T{1187:1230,{'NS_PGA'}}; % NS PGA



Amp_diff = CE32_NS_amp - CE32_EW_amp;

%% plot data
figure
histogram(CE32_comp_freq)
figure
plot(CE32_az, CE32_EW_amp, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
xlabel('Azimuth')
ylabel('Amplitude')
title('EW')
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
grid on
box on

figure
plot(CE32_az, CE32_NS_amp, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
xlabel('Azimuth')
ylabel('Amplitude')
title('NS')
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
grid on
box on

figure
plot(CE32_az, Amp_diff, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
xlabel('Azimuth')
ylabel('Amplitude difference')
title('NS_amp - EW_amp')
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
grid on
box on

figure
plot(CE32_az, CE32_comp_amp, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
xlabel('Azimuth')
ylabel('Amplitude')
title('Complex')
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
grid on
box on

figure
plot(CE32_depth, CE32_comp_amp, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
xlabel('Depth')
ylabel('Amplitude')
title('Complex')
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
grid on
box on

figure
plot(CE32_NS_freq, CE32_NS_PGA, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
xlabel('Fundamental frequency (Hz)')
ylabel('PGA')
title('Complex')
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
grid on
box on

figure
plot(CE32_NS_freq, CE32_EW_freq, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
xlabel('NS_freq (Hz)')
ylabel('EW_freq (Hz)')
title('Frequency difference')
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
grid on
box on

%% TH35 analysis
%% load data
TH35_az = T{4511:4581,{'event_azimuth'}}; % azimuth
TH35_depth = T{4511:4581,{'event_depth'}}; % depth
TH35_EW_freq = T{4511:4581,{'A_EW_HV_freq'}}; % EW_freq
TH35_EW_amp = T{4511:4581,{'A_EW_HV_amp'}}; % EW_amp
TH35_NS_freq = T{4511:4581,{'A_NS_HV_freq'}}; % NS_freq
TH35_NS_amp = T{4511:4581,{'A_NS_HV_amp'}}; % NS_amp
TH35_comp_amp = T{4511:4581,{'A_complex_HV_amp'}}; % Complex_amp
TH35_comp_freq = T{4511:4581,{'A_complex_HV_freq'}}; % Complex_freq
TH35_NS_PGA = T{4511:4581,{'NS_PGA'}}; % NS PGA



TH35_diff = TH35_az_NS_amp - TH35_az_EW_amp;

%% plot data

figure
histogram(TH35_comp_freq)
figure
plot(TH35_az, TH35_EW_amp, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
xlabel('Azimuth')
ylabel('Amplitude')
title('EW')
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
grid on
box on

figure
plot(TH35_az, TH35_NS_amp, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
xlabel('Azimuth')
ylabel('Amplitude')
title('NS')
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
grid on
box on

figure
plot(TH35_az, TH35_diff, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
xlabel('Azimuth')
ylabel('Amplitude difference')
title('NS_amp - EW_amp')
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
grid on
box on

figure
plot(TH35_az, TH35_comp_amp, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
xlabel('Azimuth')
ylabel('Amplitude')
title('Complex')
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
grid on
box on

figure
plot(TH35_depth, TH35_comp_amp, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
xlabel('Depth')
ylabel('Amplitude')
title('Complex')
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
grid on
box on

figure
plot(TH35_NS_freq, TH35_NS_PGA, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
xlabel('Fundamental frequency (Hz)')
ylabel('PGA')
title('Complex')
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
grid on
box on

figure
plot(TH35_NS_freq, TH35_EW_freq, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
xlabel('NS_freq (Hz)')
ylabel('EW_freq (Hz)')
title('Frequency difference')
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
grid on
box on