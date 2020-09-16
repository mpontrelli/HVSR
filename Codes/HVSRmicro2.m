%% HVSRmicro2
  % Read file in either .sac, sacbinary or .mseed and compute the micro
  % tremor HVSR. Options for number of windows, window length and distance
  % are available. This arose out of doing microtremor surveys in Boston, MA
  % and developed over time starting in the summer of 2019
  
    % INPUTS 
    % statname - Name of the station, used in plotting titles (string)

    % lowbound - lowcutoff value for the plots, can be any frequency. If a 
    % weird frequency is put in, 3.17 or something, it will find the closest 
    % frequency and use that. The bound between bupbound and lowbound is also
    % used for the statistics calculation, ie. the program only fi nds the peaks 
    % in between the bounds.(used to get rid of large error bars at low frequencies 
    % (must be a frequency between 0 and fs/2, default, 1 which is the first frequency)
    
    % fs - sampling frequency (Hz)
    
    % windowlen - length of the time series windows for computing the HVSR
    % (seconds)
    
    % numwin - number of windows to be averaged in HVSR computation
    % (number)
    
    % windis - distance between windows (seconds)
    
    % Vfname - filename of vertical component (string)
    
    % NSfname - filename of NS component (string)
    
    % EWfname - filename of EW component (string)
    
    % TTF - if this is toggled on, go into NRATTLE folder and plot the
    % outputs of NRATTLE. This should ONLY be on if you have done a TTF in
    % NRATTLE. ('yes' or 'no')
    
    % Allplots - if toggled on, plots everything
    
    % Timeplot - if toggled on, plots timseries
    
    % IUMagplot - if toggled on, plots individual, unfiltered mag responses
    
    % AUMagplot - if toggled on, plots averaged, unfiltered mag responses
    
    % IFMagplot - if toggled on, plots individual, filtered mag responses
    
    % AFMagplot - if toggled on, plots averaged, filtered mag responses
    
    % HVSRplot - if toggled on, plots HVSR
    
    % Filterplot - if toggled on, plots the Butterworth filter used on the
    % time series. 
    
    % outpath - the filepath for the figure outputs (string)
    
    % sav - if toggled on, saves the figures to specified output
    
    % width - the width of the smoothing filter for the magnitude
    % responses in Hz, default is 0.5 Hz.
    
    % OUTPUTS
    % 6 figures:
    
    % 1 - time series of each component
    
    % 2 - Individual, non-smoothed vertical and complex horizontal magnitude
    % responses
    
    % 3 - Averaged, non-smoothed vertical and complex horizontal magnitude
    %responses
    
    % 4 - Individual, smoothed vertical and complex horizontal magnitude
    %responses
    
    % 5 - Averaged, smoothed vertical and complex horizontal magnitude
    %responses
    
    % 6 - HVSR computed from smoothed magnitude responses. If TTF is toggled
    %on, this also plots the Theoretical transfer function computed from
    %NRATTLE
    
    
  %% Author: Marshall Pontrelli
  % Co-authors: Justin Reyes and Jeremy Salerno
  % Summer 2019
  
  % Update, ongoing edits in January 2020. See Github repository HVSR for
  % update descriptions
  
  % 8/26/2020 - Added detrending windows per advice of Jeremy Salerno
  
%% Start
function [ahatf, fax_HzN, datamat, datamatmax, confinthigh, confintlow] = HVSRmicro2(Vfname, NSfname, EWfname, fs, statname, varargin)
    %% parse inputs
    
    % create Input Parser object
    p = inputParser;
    
    % add inputs to the scheme
    defaultwindowlen = 40;
    defaultnumwin = 10;
    defaultwindis = 25;
    defaultlowbound = 0;
    defaultupbound =  fs/2 -1;
    defaultLowCorner = 0.1;
    defaultHighCorner = fs/2 - 1;
    defaultNpoles = 4;
    defaultwidth = 0.5;
    
    % Required inputs
    addRequired(p,'Vfname',@ischar);
    addRequired(p,'NSfname',@ischar);
    addRequired(p,'EWfname',@ischar);
    addRequired(p,'fs', @isnumeric);
    addRequired(p,'statname', @ischar);

    % Optional inputs
    addParameter(p, 'lowbound', defaultlowbound, @isnumeric);
    addParameter(p, 'upbound', defaultupbound, @isnumeric);
    addParameter(p, 'TTF', 'no', @ischar);
    addParameter(p, 'outpath', 'no', @ischar);
    addParameter(p, 'sav', 'no', @ischar);
    addParameter(p, 'windowlen', defaultwindowlen, @isnumeric);
    addParameter(p, 'numwin', defaultnumwin, @isnumeric);
    addParameter(p, 'windis', defaultwindis, @isnumeric);
    addParameter(p, 'Allplots', 'no', @ischar);
    addParameter(p, 'Timeplot', 'no', @ischar);
    addParameter(p, 'IUMagplot', 'no', @ischar);
    addParameter(p, 'AUMagplot', 'no', @ischar);
    addParameter(p, 'IFMagplot', 'no', @ischar);
    addParameter(p, 'AFMagplot', 'no', @ischar);
    addParameter(p, 'HVSRplot', 'no', @ischar);
    addParameter(p, 'Filterplot', 'no', @ischar);
    addParameter(p, 'LowCorner', defaultLowCorner, @isnumeric);
    addParameter(p, 'HighCorner', defaultHighCorner, @isnumeric);
    addParameter(p, 'Npoles', defaultNpoles, @isnumeric);
    addParameter(p, 'width', defaultwidth, @isnumeric);
    % parse the inputs
    parse(p, Vfname, NSfname, EWfname,fs, statname, varargin{:})
    % set varibales from the parse
    windowlen = p.Results.windowlen;
    numwin = p.Results.numwin;
    windis = p.Results.windis;
    TTF = p.Results.TTF;
    outpath = p.Results.outpath;
    sav = p.Results.sav;
    lowbound = p.Results.lowbound;
    upbound = p.Results.upbound;
    Allplots = p.Results.Allplots;
    Timeplot = p.Results.Timeplot;
    IUMagplot = p.Results.IUMagplot;
    AUMagplot = p.Results.AUMagplot;
    IFMagplot = p.Results.IFMagplot;
    AFMagplot = p.Results.AFMagplot;
    HVSRplot = p.Results.HVSRplot;
    Filterplot = p.Results.Filterplot;
    LowCorner = p.Results.LowCorner;
    HighCorner = p.Results.HighCorner;
    Npoles = p.Results.Npoles;
    width = p.Results.width;
    %turn windows into samples for windowing calculations
    sampnum = windowlen*fs; 
    windisnum = windis*fs;

    %% read files and convert into vectors
    
    [~,~,ext] = fileparts(Vfname);
   %.sacBinary
    if strcmp(ext, '.sac') == 1 % use "rdmseed" if file is in miniseed format
        [V] = ReadSacBinaryFile(Vfname); %vertical
        [NS] = ReadSacBinaryFile(NSfname); %North-south
        [EW] = ReadSacBinaryFile(EWfname); %East-West
    
    %rdsac
    %[xV] = rdsac('Vfname');
    %[xV] = xV.d;
    %[xNS] = rdsac('NSfname');
    %[xNS] = xNS.d;
    %[xEW] = rdsac('EWfname');
    %[xEW] = xEW.d;


    %miniseed
    elseif strcmp(ext, '.msd') == 1 % use "rdmseed" if file is in miniseed format
        x = rdmseed(NSfname);
        [NS, EW, V] = openmseed(x);
    
    elseif strcmp(ext, '.mseed') == 1 % use "rdmseed" if file is in miniseed format
        x = rdmseed(NSfname);
        [NS, EW, V] = openmseed(x);
    end
    %% detrend and Filter
    opol = 6;
    t = (1:length(V))/fs;
    t = t';
    [p,~,mu] = polyfit(t,V,opol); % code from https://www.mathworks.com/help/signal/ug/remove-trends-from-data.html
    f_y = polyval(p,t,[],mu);
    V = V - f_y;
    [V] = Butter2(V, fs, 'LowCorner', LowCorner, 'HighCorner', HighCorner, 'Npoles', Npoles , 'Filterplot', Filterplot);
    Filterplot = 'no'; % toggle off filter plot so it doesn't plot response three times
    [p,~,mu] = polyfit(t,NS,opol); % code from https://www.mathworks.com/help/signal/ug/remove-trends-from-data.html
    f_y = polyval(p,t,[],mu);
    NS = NS - f_y;
    [NS] = Butter2(NS, fs, 'LowCorner', LowCorner, 'HighCorner', HighCorner, 'Npoles', Npoles , 'Filterplot', Filterplot);
    [p,~,mu] = polyfit(t,EW,opol); % code from https://www.mathworks.com/help/signal/ug/remove-trends-from-data.html
    f_y = polyval(p,t,[],mu);
    EW = EW - f_y;
    [EW] = Butter2(EW, fs, 'LowCorner', LowCorner, 'HighCorner', HighCorner, 'Npoles', Npoles , 'Filterplot', Filterplot);

    %% Create a time series plot (Output 1)]
    if strcmp(Allplots, 'yes') == 1 || strcmp(Timeplot, 'yes') == 1
        timeseriesplot(NS,EW,V, fs)
    end

    %% Window data
    %Window the data with 'numwin' windows of 'windowlen' secs and 
    % 'windis' secs apart. This does support overlapping windows
    k = [1,fs];
    NSmatrix = zeros(numwin, sampnum);
    EWmatrix = zeros(numwin, sampnum);
    Vmatrix = zeros(numwin, sampnum);
    for iii = 1:numwin
        Vmatrix(iii,:) = V((k(1)):(k(2)*windowlen+k(1))-1);
        NSmatrix(iii,:) = NS((k(1)):(k(2)*windowlen+k(1))-1);
        EWmatrix(iii,:) = EW((k(1)):(k(2)*windowlen+k(1))-1);
        k(1) = k(1)+ sampnum + windisnum;
    end

    %% Compute the complex time series
    %Steidl et al. 1996
    %Hmatrix = NSmatrix + 1i.*EWmatrix; 
    %% detrend the windows
    for i = 1:numwin
        t = (1:length(Vmatrix(i,:)))/fs;
        [p,~,mu] = polyfit(t,Vmatrix(i,:),opol); 
        f_y = polyval(p,t,[],mu);
        Vmatrix(i,:) = Vmatrix(i,:) - f_y;
        [Vmatrix(i,:)] = Butter2(Vmatrix(i,:), fs, 'LowCorner', LowCorner, 'HighCorner', HighCorner, 'Npoles', Npoles , 'Filterplot', Filterplot);
        [p,~,mu] = polyfit(t,NSmatrix(i,:),opol); 
        f_y = polyval(p,t,[],mu);
        NSmatrix(i,:) = NSmatrix(i,:) - f_y;
        NSmatrix(i,:) = Butter2(NSmatrix(i,:), fs, 'LowCorner', LowCorner, 'HighCorner', HighCorner, 'Npoles', Npoles , 'Filterplot', Filterplot);
        [p,~,mu] = polyfit(t,EWmatrix(i,:),opol); 
        f_y = polyval(p,t,[],mu);
        EWmatrix(i,:) = EWmatrix(i,:) - f_y;
        EWmatrix(i,:) = Butter2(EWmatrix(i,:), fs, 'LowCorner', LowCorner, 'HighCorner', HighCorner, 'Npoles', Npoles , 'Filterplot', Filterplot);
    end
    %% window the data
    win = hann(sampnum)';
    for i = 1:numwin
        Vmatrix(i,:) = Vmatrix(i,:).*win;
        NSmatrix(i,:) = NSmatrix(i,:).*win;
        EWmatrix(i,:) = EWmatrix(i,:).*win;
    end
    

    
    %% Compute unfiltered magnitude responses
    % Compute the fft for each data window
    for iii = 1:numwin
        Vmatrix(iii,:) = 4*abs(fft(Vmatrix(iii,:)))/sampnum;
        EWmatrix(iii,:) = 4*abs(fft(EWmatrix(iii,:)))/sampnum;
        NSmatrix(iii,:) = 4*abs(fft(NSmatrix(iii,:)))/sampnum;
    end

    %Computing the frequency -axis
    N = sampnum;
    fax_binsN = (0 : N-1); %samples in NS component
    fax_HzN1 = fax_binsN*fs/N; %frequency axis NS (Hz)
    N_2 = ceil(N/2); %half magnitude spectrum
    fax_HzN = fax_HzN1(1 : N_2);
    NSmatrix2 = zeros(numwin, N_2);
    EWmatrix2 = zeros(numwin, N_2);
    Vmatrix2 = zeros(numwin, N_2);
    for iii = 1:numwin
        Vmat = Vmatrix(iii,:);
        NSmat = NSmatrix(iii, :);
        EWmat = EWmatrix(iii, :);
        Vmatrix2(iii,:) = Vmat(1 : N_2);
        NSmatrix2(iii,:) = NSmat(1 : N_2);
        EWmatrix2(iii,:) = EWmat(1 : N_2);
    end
    
    %% Combine the horizontals
    Hmatrix2 = sqrt(NSmatrix2.*EWmatrix2);

    %% create upbound and lowbound in terms of sample number
    [~, lowbound] = min(abs(fax_HzN - lowbound));
    [~, upbound] = min(abs(fax_HzN - upbound));
    
    %% plot individual unfiltered magnitude responses (OUTPUT 2)
    if strcmp(Allplots, 'yes') == 1 || strcmp(IUMagplot, 'yes') == 1
        individmagrespplot(fax_HzN, Hmatrix2, Vmatrix2, fs, lowbound, outpath, sav)
    end

    %% Average the un-smoothed magnitude responses
    [ahatfhorz, sigmahorz, confinthighhorz, confintlowhorz] =  wavav(Hmatrix2);
    [ahatfvert, sigmavert, confinthighvert, confintlowvert] =  wavav(Vmatrix2);

    %% Plot averaged unfiltered magnitude responses (OUTPUT 3)
    if strcmp(Allplots, 'yes') == 1 || strcmp(AUMagplot, 'yes') == 1
        averagedmagrespplot(fax_HzN, ahatfhorz, ahatfvert, fs,confinthighhorz, confintlowhorz, confinthighvert, confintlowvert, lowbound, outpath, sav)
    end


    %% compute smoothed magnitude responses
    window = ceil((N/fs)*width); %width for smoothing filter in samples where 20 is the number of Hz on your x-axis
    Hmatrix3 = zeros(numwin, N_2);
    Vmatrix3 = zeros(numwin, N_2);
    for iii = 1:numwin
        Vmatrix3(iii,:) = smooth(Vmatrix2(iii,:),window);
        Hmatrix3(iii,:) = smooth(Hmatrix2(iii,:),window);
    end

    %% Plot individual, smoothed magnitude responses (OUTPUT 4)
    if strcmp(Allplots, 'yes') == 1 || strcmp(IFMagplot, 'yes') == 1
        individmagrespplot(fax_HzN, Hmatrix3, Vmatrix3, fs, lowbound, outpath, sav)
    end

    %% Average the smoothed magnitude responses
    [ahatfhorz, sigmahorz, confinthighhorz, confintlowhorz] =  wavav(Hmatrix3);
    [ahatfvert, sigmavert, confinthighvert, confintlowvert] =  wavav(Vmatrix3);

    %% Plot averaged, smoothed magnitude responses (OUTPUT 5)
    if strcmp(Allplots, 'yes') == 1 || strcmp(AFMagplot, 'yes') == 1
        averagedmagrespplot(fax_HzN, ahatfhorz, ahatfvert, fs,confinthighhorz, confintlowhorz, confinthighvert, confintlowvert, lowbound, outpath, sav)
    end

    %% Compute the HVSR
    H_V = zeros(numwin, N_2);
    for iii = 1:numwin
        [H_V(iii,:)] = HV(Hmatrix3(iii,:),Vmatrix3(iii,:));
    end
    %% average the HVSR
    [ahatf, sigma, confinthigh, confintlow] =  wavav(H_V);

    %% Plot the HVSR (OUTPUT 6)
    if strcmp(Allplots, 'yes') == 1 || strcmp(HVSRplot, 'yes') == 1
        HVSRmicroplot(fax_HzN, ahatf, confinthigh, confintlow, statname, lowbound, upbound, outpath, sav, TTF)
    end
    % compute statistics on HVSR
    [peakfreqs, peakamps, hpb, f1s, f2s, Areamat, proms, amps, peakind2, freqs, sigs, I1s, I2s] = peakiden(ahatf, fax_HzN, sigma, lowbound, upbound);

    [~, I] = max(amps);

    datamat = vertcat(freqs,amps,proms,hpb,sigs);
    datamatmax = datamat(:,1);
    datamatmax = datamatmax';
    datamat = datamat';

    %% now plot, make the figure and set the base
    figure
    hold on
    confidenceinterval=shadedplot(fax_HzN(1:length(fax_HzN)), confinthigh(1:length(fax_HzN)), confintlow(1:length(fax_HzN)),[.9,.9,.9],[1,1,1]);
    hold on
    ETF = plot(fax_HzN(1 :length(fax_HzN)), ahatf(1:length(fax_HzN)), 'Color', [0 0.30196 0.6588] , 'Linewidth', 1.5);
    xlabel('Frequency (Hz)','FontSize', 18)
    ylabel('Amplification','FontSize', 18)
    title(strcat(statname), 'FontSize', 18)
    set(gca,'YScale', 'log', 'XScale', 'log','FontName', 'Times New Roman', 'FontSize', 18)
    %makes figure full screen
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    grid on 
    box on
    hold on

    %% now plot fundamental resonance
    line([freqs(I),freqs(I)],[0.1, 100],'LineStyle', '--', 'color','k')
    hold on


    %% Now the hpb-prominence cross
    hold on
    plot([freqs(I),freqs(I)],[amps(I), amps(I) - proms(I)], 'c', 'Linewidth', 2)
    hold on

    freq_ar = peakfreqs{I};
    amp_ar = peakamps{I};
    ampd = amp_ar(I1s(I));
    ampdd = amp_ar(I2s(I));
    plot([f1s(I),f2s(I)],[ampd, ampdd], 'c', 'LineWidth', 2)
    hold on

    a = "Max peak freq =" + " "  + num2str(freqs(I),3);
    b = strcat("Amp =" + " "  + num2str(amps(I),3));
    c = strcat("HPB =" + " "  +num2str(hpb(I),3));
    d = strcat("Prom =" + " "  +num2str(proms(I),3));
    e = strcat("Area =" + " "  +num2str(Areamat(I),3));
    f = strcat("\sigma =" + " "  +num2str(sigs(I),3));
    str = {a, b, c, d, e, f};
    xlim([fax_HzN(lowbound) fax_HzN(upbound)])
    ylim([0.1 40])
    q = find(ahatf == amps(I));
    q = fax_HzN(q);
    text(q - 10,12,str, 'FontName', 'Times New Roman', 'FontSize', 18,'Color', 'black', 'HorizontalAlignment', 'right', 'EdgeColor','k','BackgroundColor', 'w')



end