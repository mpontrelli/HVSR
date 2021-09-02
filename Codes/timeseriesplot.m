%% timeseriesplot 

% plot NS, EW and V from three input ground motion vectors. The vectors
% should already be filtered.

    % INPUTS
    
    % xNS - north south component of the ground motion record

    % xEW - east west component of the ground motion record
    
    % xV - vertical component of the ground motion record

    % fs - record sampling frequency
    
    % sav - is toggled on ('yes') then save the output figure
    
    % outpath - the output path
    
    % OUTPUTS
    
    % Full screen plot of the ground motion
    
%% Author: Marshall Pontrelli
% Date: 10/30/2019 
% edited - 2/11/2020 - formatting and axis labeling

%% Start
function timeseriesplot(xNS,xEW, xV, fs, varargin)
    %% parse inputs
    % create Input Parser object
    p = inputParser;

    % Required inputs
    addRequired(p, 'xNS',@isnumeric);
    addRequired(p, 'xEW',@isnumeric);
    addRequired(p, 'xV',@isnumeric);
    addRequired(p, 'fs',@isnumeric);
    
    % Optional inputs
    addParameter(p, 'yaxislabel', '', @ischar);
    addParameter(p, 'sav', 'no', @ischar);
    addParameter(p, 'outpath', 'no', @ischar);
    addParameter(p, 'out_name', 'timeseriesplot')
    
    % parse the inputs
    parse(p, xNS, xEW, xV,fs, varargin{:})
    % set varibales from the parse
    yaxis = p.Results.yaxislabel;
    sav = p.Results.sav;
    outpath = p.Results.outpath;
    out_name = p.Results.out_name;
    
    %% find the maximum value to make bounds for plotting
    a = max(abs(xV));
    a(2) = max(abs(xNS));
    a(3) = max(abs(xEW));
    d = max(a);

    %% create a time vector and plot time series
    time = (1:length(xV))/fs;
    timeseries = figure;
    
    % North - South
    subplot(3,1,1)
    plot(time,xNS)
    title('NS')
    ylabel(yaxis)
    xlim([0 length(xNS)/fs])
    ylim([-d d])
    grid on 
    box on
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);

    % East - West
    subplot(3,1,2)
    plot(time,xEW)
    title('EW')
    ylabel(yaxis)
    xlim([0 length(xV)/fs])
    ylim([-d d])
    grid on 
    box on
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);

    % Vertical
    subplot(3,1,3)
    plot(time,xV);
    title('V')
    xlabel('Time (secs)')
    ylabel(yaxis)
    xlim([0 length(xV)/fs])
    ylim([-d d])
    grid on 
    box on
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);
    
    %makes figure full screen and font in times
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

    
    %save
    if strcmp(sav, 'yes') == 1
        saveas(timeseries, strcat(outpath, '\', strcat(out_name,'.jpg')), 'jpg');
        saveas(timeseries, strcat(outpath, '\', strcat(out_name,'.fig')), 'fig');
    end
end