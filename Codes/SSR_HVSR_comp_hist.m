%% Hist SSR-HVSR comparison correaltion r values
% Author: Marshall Pontrelli
% Date: 6/6/2020

% For Pontrelli et al. 2020
close all 
clear all
sheet = 1;

%% Load HV data
corrs = xlsread('C:\Users\mpontr01\Box\Pontrelli_et_al_2020\HVSR_SSR_comparison.xlsx', 1, 'B2:B51');

%%
figure
histogram(corrs)
xlabel('Correlation value (r)')
ylabel('Number of stations')
grid on
box on
xlim([0.5,1])
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);