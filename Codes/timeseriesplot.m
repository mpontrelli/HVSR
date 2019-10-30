%% timeseries
% plot NS, EW and V from three input ground motion vectors. The vectors
% should already be filtered.

    % INPUTS
    
    % xNS - north south component of the ground motion record

    % xV - vertical component of the ground motion record
        
    % xEW - east west component of the ground motion record
    
    % fs - record sampling frequency
    
    % sav - is toggled on ('yes') then save the output figure
    
    % outpath - the output path
    
    % OUTPUTS
    
    % Full screen plot of the ground motion
    
%% Author: Marshall Pontrelli
% Date: 10/30/2019 

%% Start
function timeseriesplot(xNS,xV,xEW, fs, sav, outpath)
%% find the maximum value to make bounds for plotting
a = max(abs(xV));
a(2) = max(abs(xNS));
a(3) = max(abs(xEW));
d = max(a);

%% create a time vector and plot time series
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