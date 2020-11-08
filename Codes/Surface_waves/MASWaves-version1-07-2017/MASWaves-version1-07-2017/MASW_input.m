%% MASW_input

% Inputs soil profile into raylee_lysmer

%% Author: Marshall Pontrelli
% Date: 11/6/2020

close all
clear all

%% 
%% inputs
min_offset = 2;
max_offset = 24;
N = 24;
x = linspace(min_offset,max_offset, N);
cT_min = 40;
cT_max = 400;
delta_cT = 1;
fs = 500;
datapath = 'C:\Users\mpontr01\Box\Projects\New_England_field_work\Springfield\L62A\';
%% import the data
% ACCESSING THE DATA
% go into the data folder and get a list of stations
path = pwd;
cd(datapath)
tracelist = dir;
tracelist = tracelist(3:length(tracelist));
cd('C:\Users\mpontr01\Desktop\HVSR\Codes\Surface_waves\MASWaves-version1-07-2017\MASWaves-version1-07-2017');
% start the for loop that goes through all the station folders
B = zeros(501, 24);
for eee = 1:length(tracelist)
    trace = tracelist(eee); % read the folder info
    tracename = trace.name;
    filename = strcat(datapath,tracename);
    A = importdata(filename);
    B = A + B;
    q{eee} = importdata(filename);
end

%% Now make dispersion curve
[f,c,A] = MASWaves_dispersion_imaging(B,N,x,fs,cT_min,cT_max,delta_cT);

%% Now plot dispersion curve
resolution = 100;
fmin = 0; % Hz
fmax = 50; % Hz
FigWidth = 7; % cm
FigHeight = 7; % cm
FigFontSize = 8; % pt
figure
[fplot,cplot,Aplot] = MASWaves_plot_dispersion_image_2D(f,c,A,fmin,fmax,...
resolution,FigWidth,FigHeight,FigFontSize);

%% Select dispersion curve
f_receivers = 4.5; % Hz
select = 'numbers';
up_low_boundary = 'yes';
p = 95; % Percentage
[f_curve0,c_curve0,lambda_curve0,...
f_curve0_up,c_curve0_up,lambda_curve0_up,...
f_curve0_low,c_curve0_low,lambda_curve0_low] = ...
MASWaves_extract_dispersion_curve(f,c,A,fmin,fmax,f_receivers,... 
select,up_low_boundary,p);

%% Plot dispersion curve
FigWidth = 9; % cm
FigHeight = 6; % cm
FigFontSize = 8; % pt
type = 'f_c';
up_low_boundary = 'yes';
figure
MASWaves_plot_dispersion_curve(f_curve0,c_curve0,lambda_curve0,...
f_curve0_up,c_curve0_up,lambda_curve0_up,f_curve0_low,c_curve0_low,...
lambda_curve0_low,type,up_low_boundary,FigWidth,FigHeight,FigFontSize);

%% Now forward model
% Repeated use of MASWaves_theoretical_dispersion_curve.m, MASWaves_misfit.m
% and MASWaves_plot_theor_exp_dispersion_curves.m
% (For iteration, the layer parameters should be updated and this code section run % again).

c_test_min = 0; % m/s
c_test_max = 500; % m/s
delta_c_test = 0.5; % m/s
c_test = c_test_min:delta_c_test:c_test_max; % m/s
% Layer parameters
n = 6;

h = [2 2 2 2 2 2 Inf]; % m
beta = [150 165 170 175 190 200 210]; % m/s
alpha = beta*1.87; % m/s
rho = [1850 1850 1850 1850 1850 1850 1850]; % kg/m^3
up_low_boundary = 'yes';
[c_t,lambda_t] = MASWaves_theoretical_dispersion_curve...
(c_test,lambda_curve0,h,alpha,beta,rho,n);
up_low_boundary = 'yes';
FigWidth = 8; % cm
FigHeight = 10; % cm
FigFontSize = 8; % pt
figure
MASWaves_plot_theor_exp_dispersion_curves(c_t,lambda_t,...
c_curve0,lambda_curve0,c_curve0_up,lambda_curve0_up,...
c_curve0_low,lambda_curve0_low,up_low_boundary,...
FigWidth,FigHeight,FigFontSize)
e = MASWaves_misfit(c_t,c_curve0);

%% plot velocity profile
figure
velocity_profile(h, beta)