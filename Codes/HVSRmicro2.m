%% HVSRmicro2
  % Read file in either .sac, sacbinary or .mseed and compute the micro
  % tremor HVSR. Options for number of windows, window length and distance
  % are available. This arose out of doing microtremor surveys in Boston, MA
  % and developed over time starting in the summer of 2019
  
    % INPUTS
    % statname - Name of the station, used in plotting titles (string)

    % lowbound - lowcutoff value for the plots, can be any frequency. If a 
    % weird frequency is put in, 3.17 or something, it will find the closest 
    % frequency and use that. The bound between upbound and lowbound is also
    % used for the statistics calculation, ie. the program only finds the peaks 
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
  
%% Start
function [ahatf, fax_HzN, taxstat] = HVSRmicro2(Vfname, NSfname, EWfname, fs, statname, varargin)
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
    if ext == '.sac' % use "rdmseed" if file is in miniseed format
        [V] = ReadSacBinaryFile(Vfname); %vertical
        [NS] = ReadSacBinaryFile(NSfname); %North-south
        [EW] = ReadSacBinaryFile(EWfname); %East-West
    end
    %rdsac
    %[xV] = rdsac('Vfname');
    %[xV] = xV.d;
    %[xNS] = rdsac('NSfname');
    %[xNS] = xNS.d;
    %[xEW] = rdsac('EWfname');
    %[xEW] = xEW.d;


    %miniseed
    if ext == '.msd' % use "rdmseed" if file is in miniseed format
        x = rdmseed(NSfname);
        [NS, EW, V] = openmseed(x);
    end
    %% Filter
    [V] = Butter2(V, fs, 'LowCorner', LowCorner, 'HighCorner', HighCorner, 'Npoles', Npoles , 'Filterplot', Filterplot);
    Filterplot = 'no'; % toggle off filter plot so it doesn't plot response three times
    [NS] = Butter2(NS, fs, 'LowCorner', LowCorner, 'HighCorner', HighCorner, 'Npoles', Npoles , 'Filterplot', Filterplot);
    [EW] = Butter2(EW, fs, 'LowCorner', LowCorner, 'HighCorner', HighCorner, 'Npoles', Npoles , 'Filterplot', Filterplot);

    %% Create a time series plot (Output 1)]
    if strcmp(Allplots, 'yes') == 1 || strcmp(Timeplot, 'yes') == 1
        timeseriesplot(NS,V,EW, fs)
    end

    %% Window data
    %Window the data with 'numwin' windows of 'windowlen' secs and 
    % 'windis' secs apart. This does support overlapping windows
    k = [1,fs];
    for iii = 1:numwin
        Vmatrix(iii,:) = V((k(1)):(k(2)*windowlen+k(1))-1);
        NSmatrix(iii,:) = NS((k(1)):(k(2)*windowlen+k(1))-1);
        EWmatrix(iii,:) = EW((k(1)):(k(2)*windowlen+k(1))-1);
        k(1) = k(1)+ sampnum + windisnum;
    end

    %% Compute the complex time series
    %Steidl et al. 1996
    Hmatrix = NSmatrix + 1i.*EWmatrix; 

    %% window the data
    win = hann(sampnum)';
    for i = 1:numwin
        Vmatrix(i,:) = Vmatrix(i,:).*win;
        Hmatrix(i,:) = Hmatrix(i,:).*win;
    end
    %% Compute unfiltered magnitude responses
    % Compute the fft for each data window
    for iii = 1:numwin
        Vmatrix(iii,:) = 4*abs(fft(Vmatrix(iii,:)))/sampnum;
        Hmatrix(iii,:) = 4*abs(fft(Hmatrix(iii,:)))/sampnum;
    end

    %Computing the frequency -axis
    N = sampnum;
    fax_binsN = (0 : N-1); %samples in NS component
    fax_HzN1 = fax_binsN*fs/N; %frequency axis NS (Hz)
    N_2 = ceil(N/2); %half magnitude spectrum
    fax_HzN = fax_HzN1(1 : N_2);
    for iii = 1:numwin
        Vmat = Vmatrix(iii,:);
        Hmat = Hmatrix(iii, :);
        Vmatrix2(iii,:) = Vmat(1 : N_2);
        Hmatrix2(iii,:) = Hmat(1 : N_2);
    end

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
    [matrix, matrix1, peakind,ahatf1,newfaxhz1] = peakiden(ahatf, fax_HzN, lowbound, upbound);
    [taxstat] = specratstat(peakind, matrix, matrix1, ahatf1, newfaxhz1, sigma, statname, lowbound, upbound);
    % save
    if strcmp(sav, 'yes') == 1
        taxstat2 = cell2table(taxstat);
        writetable(taxstat2, strcat(outpath, '\', 'statistics.txt'));
    end

end