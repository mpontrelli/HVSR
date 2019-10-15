%% HVSRmicro2
  %Read file in either .sac, sacbinary or .mseed and compute the micro
  %tremo HVSR. Options for number of windows, window length and distance
  %are available. This arose out of doing microtremor surveys in Boston, MA
  %and developed over time starting in the summer of 2019
  
    %INPUTS
    %statname - Name of the station, used in plotting titles (string)

    %lowbound - lowcutoff value for the plots (used to get rid of large
    %error bars at low frequencies (samples)
    
    %fs - sampling frequency (Hz)
    
    %windowlen - length of the time series windows for computing the HVSR
    %(seconds)
    
    %numwin - number of windows to be averaged in HVSR computation
    %(number)
    
    %windis - distance between windows (seconds)
    
    %Vfname - filename of vertical component (string)
    
    %NSfname - filename of NS component (string)
    
    %EWfname - filename of EW component (string)
    
    %TTF - if this is toggled on, go into NRATTLE folder and plot the
    %outputs of NRATTLE. This should ONLY be on if you have done a TTF in
    %NRATTLE. ('yes' or 'no')
    
    %outpath - the filepath for the figure outputs (string)
    
    %sav - if toggled on, saves the figures to specified output
    
    %OUTPUTS
    %6 figures:
    
    %1 - time series of each component
    
    %2 - Individual, non-smoothed vertical and complex horizontal magnitude
    %responses
    
    %3 - Averaged, non-smoothed vertical and complex horizontal magnitude
    %responses
    
    %4 - Individual, smoothed vertical and complex horizontal magnitude
    %responses
    
    %5 - Averaged, smoothed vertical and complex horizontal magnitude
    %responses
    
    %6 - HVSR computed from smoothed magnitude responses. If TTF is toggled
    %on, this also plots the Theoretical transfer function computed from
    %NRATTLE
    
    
  %Author: Marshall Pontrelli
  %Co-author Justin Reyes
  %Summer 2019

close all
clear all
%% INPUTS

statname = '';
lowbound = 10;
fs = 100;
fsmin = fs;
windowlen = 50;
numwin = 15;
windis = 10;
Vfname = 'C:\Users\Marshall\Desktop\boston-site-response\field_deployments\radian_demo_deployment\potential_boreholes\Kraft2-borehole-11\Site_response\data\2000049134427006_kraft21_1_1.sac';
NSfname = 'C:\Users\Marshall\Desktop\boston-site-response\field_deployments\radian_demo_deployment\potential_boreholes\Kraft2-borehole-11\Site_response\data\2000049134427006_kraft22_1_2.sac';
EWfname = 'C:\Users\Marshall\Desktop\boston-site-response\field_deployments\radian_demo_deployment\potential_boreholes\Kraft2-borehole-11\Site_response\data\2000049134427006_kraft23_1_3.sac';
TTF = 'no';
outpath = 'C:\Users\Marshall\Desktop\boston-site-response\field_deployments\radian_demo_deployment\potential_boreholes\Kraft2-borehole-11\Site_response\figures';
sav = 'no';

%turn windows into samples for windowing calculations
sampnum = windowlen*fs; 
windisnum = windis*fs;

%% read files and convert into vectors
%.sacBinary
[xV] = ReadSacBinaryFile(Vfname); %vertical
[xNS] = ReadSacBinaryFile(NSfname); %North-south
[xEW] = ReadSacBinaryFile(EWfname); %East-West
%rdsac
% [xV] = rdsac('Vfname');
% [xV] = xV.d;
% [xNS] = rdsac('NSfname');
% [xNS] = xNS.d;
% [xEW] = rdsac('EWfname');
% [xEW] = xEW.d;

%miniseed
%A = rdmseed('fname');
%% Filter
[xV] = Butter2(xV);
[xNS] = Butter2(xNS);
[xEW] = Butter2(xEW);

%% find the maximum value in the microtremor to make bounds for plotting
a = max(abs(xV));
a(2) = max(abs(xNS));
a(3) = max(abs(xEW));
d = max(a);

%% create a time vector and plot time series (OUTPUT 1)
time = (1:length(xV))/fs;
timeseries = figure;
subplot(3,1,1)
plot(time,xV)
title('V')
xlabel('Time (secs)')
ylabel('counts')
xlim([0 length(xV)/fs])
ylim([-d d])
set(gca,'FontSize',20)
grid on 
box on

subplot(3,1,2)
plot(time,xEW)
title('EW')
xlabel('Time (secs)')
ylabel('counts')
xlim([0 length(xV)/fs])
ylim([-d d])
set(gca,'FontSize',20)
grid on 
box on


subplot(3,1,3)
plot(time,xNS);
title('NS')
xlabel('Time (secs)')
ylabel('counts')
xlim([0 length(xV)/fs])
ylim([-d d])
set(gca,'FontSize',20)
grid on 
box on

%makes figure full screen
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

%save
if strcmp(sav, 'yes') == 1
    saveas(timeseries, strcat(outpath, '\', 'timeseries.jpg'), 'jpg');
    saveas(timeseries, strcat(outpath, '\', 'timeseries.fig'), 'fig');
end

%% Window data
%Window the data with 'numwin' non-overlapping windows of 'windowlen' secs and 
% 'windis' secs apart
k = [1,fs];
for iii = 1:numwin
    xVmatrix(iii,:) = xV((k(1)):(k(2)*windowlen+k(1))-1);
    xNSmatrix(iii,:) = xNS((k(1)):(k(2)*windowlen+k(1))-1);
    xEWmatrix(iii,:) = xEW((k(1)):(k(2)*windowlen+k(1))-1);
    k(1) = k(1)+ sampnum + windisnum;
end

%% Compute the complex time series
xHmatrix = xNSmatrix + 1i.*xEWmatrix; 

%% Compute unfiltered magnitude responses
% Compute the fft for each data window
for iii = 1:numwin
    XVmatrix(iii,:) = abs(fft(xVmatrix(iii,:)))/sampnum;
    XHmatrix(iii,:) = abs(fft(xHmatrix(iii,:)))/sampnum;
end

%Computing the frequency -axis
N = sampnum;
fax_binsN = (0 : N-1); %samples in NS component
fax_HzN1 = fax_binsN*fs/N; %frequency axis NS (Hz)
N_2 = ceil(N/2); %half magnitude spectrum
fax_HzN = fax_HzN1(1 : N_2);
for iii = 1:numwin
    XVmat = XVmatrix(iii,:);
    XHmat = XHmatrix(iii, :);
    XVmatrix2(iii,:) = XVmat(1 : N_2);
    XHmatrix2(iii,:) = XHmat(1 : N_2);
end
%% plot individual unfiltered magnitude responses (OUTPUT 2)
individualunfiltered = figure;
subplot(1,2,1)
plot(fax_HzN, XHmatrix2, 'Color',  'k' , 'Linewidth', .5);
title('Horizontal', 'FontSize', 20)
xlabel('Frequency (Hz)','FontSize', 18)
ylabel('Amplification','FontSize', 18)
set(gca,'FontSize',20,'YScale', 'log')
xlim([fax_HzN(10) 40])
ylim([0.1 1000])
xticks([0.1 1 10])
xticklabels({'0.1','1','10'})
yticks([0.1 1 10 100 1000])
yticklabels({'0.1','1','10', '100', '1000'})
grid on 
box on

subplot(1,2,2)
plot(fax_HzN, XVmatrix2, 'Color',  'k' , 'Linewidth', .5);
title('Vertical', 'FontSize', 20)
xlabel('Frequency (Hz)','FontSize', 18)
ylabel('Amplification','FontSize', 18)
set(gca,'FontSize',20,'YScale', 'log')
xlim([fax_HzN(10) 41])
ylim([0.1 1000])
xticks([ 1 10 40])
xticklabels({'1','10', '40'})
yticks([0.1 1 10 100 1000])
yticklabels({'0.1','1','10', '100', '1000'})
grid on 
box on

%makes figure full screen
set(gcf, 'Units', 'Normalized', 'OuterPosition', [1, 0.5, 1, 0.5]);

%save
if strcmp(sav, 'yes') == 1
    saveas(individualunfiltered, strcat(outpath, '\', 'individualunfiltered.jpg'));
    saveas(individualunfiltered, strcat(outpath, '\', 'individualunfiltered.fig'));
end

%% Average the un-smoothed magnitude responses
[ahatfhorz, sigmahorz, confinthighhorz, confintlowhorz] =  wavav(XHmatrix2);
[ahatfvert, sigmavert, confinthighvert, confintlowvert] =  wavav(XVmatrix2);

%% Plot averaged unfiltered magnitude responses (OUTPUT 3)
averageunfiltered = figure;
subplot(1,2,1)
hold on
confidenceinterval=shadedplot(fax_HzN, confinthighhorz, confintlowhorz,[.9,.9,.9],'k');
hold on
plot(fax_HzN, ahatfhorz, 'Color', [0 0.30196 0.6588] , 'Linewidth', 1.5);
xlabel('Frequency (Hz)','FontSize', 18)
ylabel('Amplification','FontSize', 18)
title('Horizontal','FontSize', 18)
set(gca,'FontSize',20,'YScale', 'log')
xlim([fax_HzN(10) 40])
ylim([0.1 1000])
xticks([0.1 1 10])
xticklabels({'0.1','1','10'})
yticks([0.1 1 10 100 1000])
yticklabels({'0.1','1','10', '100', '1000'})
grid on 
box on


subplot(1,2,2)
hold on
confidenceinterval=shadedplot(fax_HzN, confinthighvert, confintlowvert,[.9,.9,.9],'k');
hold on
plot(fax_HzN, ahatfvert, 'Color', [0 0.30196 0.6588] , 'Linewidth', 1.5);
xlabel('Frequency (Hz)','FontSize', 18)
ylabel('Amplification','FontSize', 18)
title('Vertical','FontSize', 18)
set(gca,'FontSize',20,'YScale', 'log')
xlim([fax_HzN(10) 40])
ylim([0.1 1000])
xticks([0.1 1 10])
xticklabels({'0.1','1','10'})
yticks([0.1 1 10 100 1000])
yticklabels({'0.1','1','10', '100', '1000'})
grid on 
box on

%makes figure full screen
set(gcf, 'Units', 'Normalized', 'OuterPosition', [1, 0.5, 1, 0.5]);

%save
if strcmp(sav, 'yes') == 1
    saveas(averageunfiltered, strcat(outpath, '\', 'averageunfiltered.jpg'));
    saveas(averageunfiltered, strcat(outpath, '\', 'averageunfiltered.fig'));
end

%% compute smoothed magnitude responses
width = .1; %width for triangle moving average filter in hz
window = ceil((N/20)*width); %width for triangle moving average filter in samples where 20 is the number of Hz on your x-axis
for iii = 1:numwin
    XVmatrix3(iii,:) = smooth(XVmatrix2(iii,:),window);
    XHmatrix3(iii,:) = smooth(XHmatrix2(iii,:),window);
end

%% Plot individual, smoothed magnitude responses (OUTPUT 4)
individualfiltered = figure;
subplot(1,2,1)
plot(fax_HzN, XHmatrix3, 'Color',  'k' , 'Linewidth', .5);
title('Horizontal', 'FontSize', 20)
xlabel('Frequency (Hz)','FontSize', 18)
ylabel('Amplification','FontSize', 18)
set(gca,'FontSize',20,'YScale', 'log')
xlim([fax_HzN(10) 40])
ylim([0.1 1000])
xticks([0.1 1 10])
xticklabels({'0.1','1','10'})
yticks([0.1 1 10 100 1000])
yticklabels({'0.1','1','10', '100', '1000'})
grid on 
box on

subplot(1,2,2)
plot(fax_HzN, XVmatrix3, 'Color',  'k' , 'Linewidth', .5);
title('Vertical', 'FontSize', 20)
xlabel('Frequency (Hz)','FontSize', 18)
ylabel('Amplification','FontSize', 18)
set(gca,'FontSize',20,'YScale', 'log')
xlim([fax_HzN(10) 40])
ylim([0.1 1000])
xticks([0.1 1 10])
xticklabels({'0.1','1','10'})
yticks([0.1 1 10 100 1000])
yticklabels({'0.1','1','10', '100', '1000'})
grid on 
box on

%makes figure full screen
set(gcf, 'Units', 'Normalized', 'OuterPosition', [1, 0.5, 1, 0.5]);

%save
if strcmp(sav, 'yes') == 1
    saveas(individualfiltered, strcat(outpath, '\', 'individualfiltered.jpg'));
    saveas(individualfiltered, strcat(outpath, '\', 'individualfiltered.fig'));
end

%% Average the smoothed magnitude responses
[ahatfhorz, sigmahorz, confinthighhorz, confintlowhorz] =  wavav(XHmatrix3);
[ahatfvert, sigmavert, confinthighvert, confintlowvert] =  wavav(XVmatrix3);

%% Plot averaged, smoothed magnitude responses (OUTPUT 5)
averagefiltered = figure;
subplot(1,2,1)
hold on
confidenceinterval=shadedplot(fax_HzN, confinthighhorz, confintlowhorz,[.9,.9,.9],'k');
hold on
plot(fax_HzN, ahatfhorz, 'Color', [0 0.30196 0.6588] , 'Linewidth', 1.5);
xlabel('Frequency (Hz)','FontSize', 18)
ylabel('Amplification','FontSize', 18)
title('Horizontal','FontSize', 18)
set(gca,'FontSize',20,'YScale', 'log','XScale', 'log')
xlim([fax_HzN(10) 40])
ylim([0.1 1000])
xticks([0.1 1 10])
xticklabels({'0.1','1','10'})
yticks([0.1 1 10 100 1000])
yticklabels({'0.1','1','10', '100', '1000'})
xticks([0.1 1 10])
xticklabels({'0.1','1','10'})
grid on 
box on

subplot(1,2,2)
hold on
confidenceinterval=shadedplot(fax_HzN, confinthighvert, confintlowvert,[.9,.9,.9],'k');
hold on
plot(fax_HzN, ahatfvert, 'Color', [0 0.30196 0.6588] , 'Linewidth', 1.5);
xlabel('Frequency (Hz)','FontSize', 18)
ylabel('Amplification','FontSize', 18)
title('Vertical','FontSize', 18)
set(gca,'FontSize',20,'YScale', 'log', 'XScale', 'log')
xlim([fax_HzN(10) 40])
ylim([0.1 1000])
xticks([0.1 1 10])
xticklabels({'0.1','1','10'})
yticks([0.1 1 10 100 1000])
yticklabels({'0.1','1','10', '100', '1000'})
xticks([0.1 1 10])
xticklabels({'0.1','1','10'})
grid on 
box on

%makes figure full screen
set(gcf, 'Units', 'Normalized', 'OuterPosition', [1, 0.5, 1, 0.5]);

%save
if strcmp(sav, 'yes') == 1
    saveas(averagefiltered, strcat(outpath, '\', 'averagefiltered.jpg'));
    saveas(averagefiltered, strcat(outpath, '\', 'averagefiltered.fig'));
end

%% Compute the HVSR
for iii = 1:numwin
    [H_V(iii,:)] = HV(XHmatrix3(iii,:),XVmatrix3(iii,:));
end

%% average the HVSR
[ahatf, sigma, confinthigh, confintlow] =  wavav(H_V);

%% Plot the HVSR (OUTPUT 6)
HVSR = figure;
hold on
confidenceinterval=shadedplot(fax_HzN(10:length(fax_HzN)), confinthigh(10:length(fax_HzN)), confintlow(10:length(fax_HzN)),[.9,.9,.9],'k');
hold on
ETF = plot(fax_HzN(10 :length(fax_HzN)), ahatf(10:length(fax_HzN)), 'Color', [0 0.30196 0.6588] , 'Linewidth', 1.5);
xlabel('Frequency (Hz)','FontSize', 18)
ylabel('Amplification','FontSize', 18)
title(strcat(statname), 'FontSize', 18)
set(gca,'FontSize',20,'YScale', 'log')
xlim([fax_HzN(10) 40])
ylim([0.1 100])
xticks([1 10 40])
xticklabels({'1','10', '40'})
yticks([0.1 1 10 100])
yticklabels({'0.1','1','10', '100'})
grid on 
box on
hold on



%% if you want to plot TTF from NRATTLE
if strcmp(TTF, 'yes') == 1
    Read_amps_4_plot
    TTF = plot(fax_HzN(10 :length(fax_HzN)), amps(10:length(fax_HzN)), 'Color', 'k', 'linewidth', 2);
    legend([ETF, TTF], {'HVSR', 'TTF'})
end

%makes figure full screen
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

%save
if strcmp(sav, 'yes') == 1
    saveas(HVSR, strcat(outpath, '\', 'HVSR.jpg'));
    saveas(HVSR, strcat(outpath, '\', 'HVSR.fig'));
end

%% compute statistics on HVSR
[matrix, matrix1, peakind,ahatf1,newfaxhz1] = peakiden(ahatf, fax_HzN, lowbound, fsmin);
[taxstat] = specratstat(peakind, matrix, matrix1, ahatf1, newfaxhz1, sigma, statname, lowbound);



%% End