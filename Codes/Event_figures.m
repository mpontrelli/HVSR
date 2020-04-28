%% Event_figures
%% Author: Marshall Pontrelli
% Date: 4/27/2020

%%
close all
clear all

%%
CE32 = load('C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Processed_data2\CE32\CE3220170919181440.mat');
TP13 = load('C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Processed_data2\TP13\TP1320170919181440.mat');
%% inputs
fs = CE32.data.meta.instrument.fs;
time = CE32.data.processing.filtereddata.time_orig;
timelong = CE32.data.processing.filtereddata.time;
%% time series
% Acceleration CE32

NS = CE32.data.processing.filtereddata.acceleration.NS.waveform_orig;
EW = CE32.data.processing.filtereddata.acceleration.EW.waveform_orig;
V = CE32.data.processing.filtereddata.acceleration.EW.waveform_orig;
timeseriesplot(NS,EW, V, fs,'yaxislabel', 'm/s^2')

% Acceleration TP13
NS = TP13.data.processing.filtereddata.acceleration.NS.waveform_orig;
EW = TP13.data.processing.filtereddata.acceleration.EW.waveform_orig;
V = TP13.data.processing.filtereddata.acceleration.EW.waveform_orig;
timeseriesplot(NS,EW, V, fs,'yaxislabel', 'm/s^2')

%% time series with velocity and displacement CE32

NS = CE32.data.processing.filtereddata.acceleration.NS.waveform_orig;
[PGA, PGV, PGD, v, d] =  waveform_integrate(NS, fs);

% find the maximum value to make bounds for plotting
dd = max(abs(NS));

% Acceleration
subplot(3,1,1)
plot(time,NS)
ylabel('m/s^2')
xlim([0 length(NS)/fs])
ylim([-dd dd])
grid on 
box on
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);

% Velocity
subplot(3,1,2)
plot(time,v)
ylabel('m/s')
xlim([0 length(NS)/fs])
ylim([-dd dd])
grid on 
box on
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);

% displacement
subplot(3,1,3)
plot(time,d);
xlabel('Time (secs)')
ylabel('m')
xlim([0 length(NS)/fs])
ylim([-dd dd])
grid on 
box on
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);
    
%makes figure full screen and font in times
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
%% time series with velocity and displacement TP13
time = TP13.data.processing.filtereddata.time_orig;
NS = TP13.data.processing.filtereddata.acceleration.NS.waveform_orig;
[PGA, PGV, PGD, v, d] =  waveform_integrate(NS, fs);

% find the maximum value to make bounds for plotting
dd = max(abs(NS));

% Acceleration
subplot(3,1,1)
plot(time,NS)
ylabel('m/s^2')
xlim([0 length(NS)/fs])
ylim([-dd dd])
grid on 
box on
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);

% Velocity
subplot(3,1,2)
plot(time,v)
ylabel('m/s')
xlim([0 length(NS)/fs])
ylim([-dd dd])
grid on 
box on
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);

% displacement
subplot(3,1,3)
plot(time,d);
xlabel('Time (secs)')
ylabel('m')
xlim([0 length(NS)/fs])
ylim([-dd dd])
grid on 
box on
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);
    
%makes figure full screen and font in times

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
%% Arias intensity
% CE32
NS_arias_int = CE32.data.processing.filtereddata.acceleration.NS.arias.intensity;
NS_arias_D595 = CE32.data.processing.filtereddata.acceleration.NS.arias.D595;
NS_arias_D575 = CE32.data.processing.filtereddata.acceleration.NS.arias.D575;
Arias_plot(time, Ianorm)
%% Response spectra