%% MASW_input

% Inputs soil profile into raylee_lysmer

%% Author: Marshall Pontrelli
% Date: 11/6/2020

close all
clear all

%% 
%% inputs
path = pwd;
min_offset = 2;
max_offset = 24;
N = 24;
x = linspace(min_offset,max_offset, N);
cT_min = 40;
cT_max = 1000;
delta_cT = 1;
fs = 500;
datapath = strcat(path, '\Data\');
%% import the data
% ACCESSING THE DATA
% go into the data folder and get a list of stations
cd(datapath)
tracelist = dir;
tracelist = tracelist(3:length(tracelist));
cd(strcat('C:\Users\',getenv('username'),'\Desktop\HVSR\Codes\Surface_waves\MASWaves-version1-07-2017\MASWaves-version1-07-2017'));
% start the for loop that goes through all the station folders
B = zeros(501, 24);
for eee = 1:length(tracelist)
    trace = tracelist(eee); % read the folder info
    tracename = trace.name;
    filename = strcat(datapath,tracename);
    A = importdata(filename);
    A = A(:,3:26);
    B = A + B;
    q{eee} = importdata(filename);
end

%% figure out how big B is
BB = B;
B = B(:,1:23);
c = size(B);
a = c(2);
b = c(1);
time = 1000*(0:b-1)/(fs);
data.trace = B;

%% Now normalize B
for i = 1:a
    sta = B(:,i);
    max_sta = max(sta);
    col = sta/max_sta;
    B_norm(:,i) = col;
end
data.trace_norm = B_norm;

%% now plot
traces = figure;
for i = 1:a
    sta = B_norm(:,i);
    % Extract positive and negative part
    yp = (sta + abs(sta))/2;
    yn = (sta - abs(sta))/2;
    hold on
    plot(time,yn+i,'b')
    hold on
    fill(time,yp+i,i,'Facecolor','k')
    hold on
end
view([90, 90])
grid on
box on
xlabel('Time (ms)','FontName', 'Times New Roman', 'FontSize', 18,'rotation', 270, 'VerticalAlignment','middle');

yticks([1 a])
yticklabels({num2str(min_offset),num2str(max_offset)})
ax = gca;
ax.YAxisLocation = 'right';


set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);    
xlim([0,500])
ylim([-0.5, i+1])
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.5, 0.96]);
title('Offset (m)','FontName', 'Times New Roman', 'FontSize', 18);

out_name = 'trace';
saveas(traces, strcat(path, '\','Figures','\', strcat(out_name,'.jpg')), 'jpg');
saveas(traces, strcat(path, '\','Figures','\', strcat(out_name,'.fig')), 'fig');

%% Now make dispersion curve
[f,c,A] = MASWaves_dispersion_imaging(BB,N,x,fs,cT_min,cT_max,delta_cT);

%% Now plot dispersion curve
resolution = 100;
fmin = 0; % Hz
fmax = 100; % Hz
FigWidth = 15; % cm
FigHeight = 15; % cm
FigFontSize = 10; % pt
disp_curve = figure;
[fplot,cplot,Aplot] = MASWaves_plot_dispersion_image_2D(f,c,A,fmin,fmax,...
resolution,FigWidth,FigHeight,FigFontSize);

out_name = 'Disp_curve';
saveas(disp_curve, strcat(path, '\','Figures','\', strcat(out_name,'.jpg')), 'jpg');
saveas(disp_curve, strcat(path, '\','Figures','\', strcat(out_name,'.fig')), 'fig');
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
disp_pick = figure;
MASWaves_plot_dispersion_curve(f_curve0,c_curve0,lambda_curve0,...
f_curve0_up,c_curve0_up,lambda_curve0_up,f_curve0_low,c_curve0_low,...
lambda_curve0_low,type,up_low_boundary,FigWidth,FigHeight,FigFontSize);
out_name = 'Disp_picks';
saveas(disp_pick, strcat(path, '\','Figures','\', strcat(out_name,'.jpg')), 'jpg');
saveas(disp_pick, strcat(path, '\','Figures','\', strcat(out_name,'.fig')), 'fig');
%% Now forward model
% Repeated use of MASWaves_theoretical_dispersion_curve.m, MASWaves_misfit.m
% and MASWaves_plot_theor_exp_dispersion_curves.m
% (For iteration, the layer parameters should be updated and this code section run % again).
cd(strcat('C:\Users\',getenv('username'),'\Desktop\HVSR\Codes\Surface_waves\MASWaves-version1-07-2017\MASWaves-version1-07-2017'));

c_test_min = 0; % m/s
c_test_max = 500; % m/s
delta_c_test = 0.5; % m/s
c_test = c_test_min:delta_c_test:c_test_max; % m/s
% Layer parameters
n = 7;
h = [5,5,5,5,5,5, 40]; % m
data.thicknesses = h;
beta = [200,220,230,240,250,250, 250]; % m/s
data.velocities = beta;
alpha = beta*1.87; % m/s
rho = [1250 1250 1250 1850 1850 1850 1850]; % kg/m^3
data.rhos = rho;

% Now model
[c_t,lambda_t] = MASWaves_theoretical_dispersion_curve...
(c_test,lambda_curve0,h,alpha,beta,rho,n);
up_low_boundary = 'yes';

% Now plot
FigWidth = 8; % cm
FigHeight = 10; % cm
FigFontSize = 8; % pt
Theor_disp = figure;
MASWaves_plot_theor_exp_dispersion_curves(c_t,lambda_t,...
c_curve0,lambda_curve0,c_curve0_up,lambda_curve0_up,...
c_curve0_low,lambda_curve0_low,up_low_boundary,...
FigWidth,FigHeight,FigFontSize)
e = MASWaves_misfit(c_t,c_curve0);

out_name = 'Theor_disp';
saveas(Theor_disp, strcat(path, '\','Figures','\', strcat(out_name,'.jpg')), 'jpg');
saveas(Theor_disp, strcat(path, '\','Figures','\', strcat(out_name,'.fig')), 'fig');

% plot velocity profile
Velocity_profile = figure;
velocity_profile(h, beta)
out_name = 'Theor_disp';
saveas(Theor_disp, strcat(path, '\','Figures','\', strcat(out_name,'.jpg')), 'jpg');
saveas(Theor_disp, strcat(path, '\','Figures','\', strcat(out_name,'.fig')), 'fig');
cd(strcat('C:\Users\',getenv('username'),'\Desktop\HVSR\Codes'));

out_name = 'Velocity_profile';
saveas(Velocity_profile, strcat(path, '\','Figures','\', strcat(out_name,'.jpg')), 'jpg');
saveas(Velocity_profile, strcat(path, '\','Figures','\', strcat(out_name,'.fig')), 'fig');

%% average velocity calculation
d = h(1:n);
v = beta(1:n);
[v_avg] =  S_wave_avg(d,v);
data.v_avg = v_avg;
[vs_30] =  Vs_30(h, beta);
data.vs_30 = vs_30;
Vs_disp = strcat('Vs average = ', {' '},num2str(v_avg), {' '},'m/s');
Vs_disp = Vs_disp{1};
disp(Vs_disp)
depth = sum(d);
data.depth = depth;
depth_disp = strcat('Depth = ', {' '},num2str(depth), {' '},'m');
depth_disp = depth_disp{1};
disp(depth_disp)
fn = v_avg/(4*depth);
fn_disp = strcat('f0 = ', {' '},num2str(fn), {' '},'Hz');
fn_disp = fn_disp{1};
disp(fn_disp)
data.f0 = fn;

%% save data structure
%% Now save data
save(strcat(path, '\','Figures','\','data.mat'),'data')