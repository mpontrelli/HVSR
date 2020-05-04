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
    Hmatrix2 = zeros(numwin, N_2);
    Vmatrix2 = zeros(numwin, N_2);
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
    [matrix, matrix1, peakind,ahatf1,newfaxhz1, peakfreqs, peakamps, Areamat] = peakiden(ahatf, fax_HzN, lowbound, upbound);
    [taxstat,sigma1, sigs] = specratstat(peakind, matrix, matrix1', ahatf1, newfaxhz1, sigma, statname, lowbound, upbound);
    % save
    if strcmp(sav, 'yes') == 1
        taxstat2 = cell2table(taxstat);
        writetable(taxstat2, strcat(outpath, '\', 'statistics.txt'));
    end
    
    %% now plot with filled peaks
    % Compute half power bandwidths
    for f = 1:length(peakind)
        loc = peakind(f);
        A = matrix(f,2);
        [I1, I, f1, f2, hpb] =  HalfPowerBand2(A, loc, newfaxhz1, ahatf1);
        hpb1(f,1) = hpb;
        hpb1(f,2) = f1;
        hpb1(f,3) = f2;
        hpb1(f,4) = I1;
        hpb1(f,5) = I;
    end
    
    % Pull out vectors: amplitude, frequency, left hpb frequency, right
    % hpb frequency, prominence, hpb
    amps = matrix(:,2)'; freqs = matrix(:,1)'; f1s = hpb1(:,2)'; f2s = hpb1(:,3)';
    proms = matrix1; hpbss = hpb1(:,1)';
    
    % Now set some conditions to classify the peak
    [~, amp_class] = sort(amps, 'descend');
    [~, prom_class] = sort(proms, 'descend');
    [~, hpb_class] = sort(hpbss, 'descend');
    [~, area_class] = sort(Areamat, 'descend');
    [~, sig_class] = sort(sigs, 'descend');
    classmatrix = vertcat(amp_class, prom_class, hpb_class, area_class, sig_class);
    % now plot
    figure

    % start with the confidence interval
    confidenceinterval=shadedplot(fax_HzN(lowbound:length(fax_HzN)), confinthigh(lowbound:length(confinthigh)), confintlow(lowbound:length(confintlow)),[.9,.9,.9],[1 1 1]);
    hold on
    
    % Plot filled in peaks color coordinated based on their area
    len = length(peakfreqs);
    red = [1, 0, 0];
    pink = [255, 192, 203]/255;
    colors_p = [linspace(red(1),pink(1),len)', linspace(red(2),pink(2),len)', linspace(red(3),pink(3),len)'];
    col = colormap(colors_p);
    for i = 1:length(peakfreqs)
        freq_ar = peakfreqs{i};
        amp_ar = peakamps{i};
        fill(freq_ar, amp_ar, col(i,:), 'LineStyle','none')
        alpha(.5)
        hold on
    end
%     hold on
%     freq_ar = peakfreqs{1};
%     amp_ar = peakamps{1};
%     fill(freq_ar, amp_ar, 'r', 'LineStyle','none')
%     alpha(.5)
%     hold on
%     freq_ar = peakfreqs{2};
%     amp_ar = peakamps{2};
%     fill(freq_ar, amp_ar, 'b', 'LineStyle','none')
%     alpha(.5)
%     hold on

    plot(fax_HzN,ahatf, 'LineWidth',2, 'Color',[0 0.5 0])
    for i = 1:2%length(f1s)
        dddd = f1s(i);
        ampd = find(fax_HzN == dddd);
        ampd = ahatf(ampd);
        ddddd = f2s(i);
        ampdd = find(fax_HzN == ddddd);
        ampdd = ahatf(ampdd);
        plot([dddd,ddddd],[ampd, ampdd], 'k')
        hold on
    end
    hold on

    % now plot fundamental resonance
    line([freqs(1),freqs(1)],[0.1, 40],'LineStyle', '--', 'color','r')

    % and the second peak
    line([freqs(2),freqs(2)],[0.1, 40], 'LineStyle','--', 'color', 'b')

    hold on
    for i = 1:2%length(amps)
        plot([freqs(i),freqs(i)],[amps(i), amps(i) - proms(i)], 'k')
        hold on
    end

    hold on
    xlim([0.1 20])
    ylim([0.5 40])
    xticks([.1 1 20])
    xticklabels({'0.1', '1', '10'})
    yticks([1 10 100])
    yticklabels({ '1','10', '100'})
    xlabel('Frequency (Hz)')
    ylabel('Amplification')
    title(statname)
    set(gca,'YScale', 'log','FontName', 'Times New Roman', 'FontSize', 14)
    grid on
    box on
    hold on

    str = {strcat('Peak 1 prom = ',{' '},num2str(proms(1))), strcat('Peak 2 prom = ',{' '},num2str(proms(2)))};
    stra = str{1}{1};
    strb = str{2}{1};
    strc = strcat('Peak 1 area = ',num2str(Areamat(1)));
    strd = strcat('Peak 1 hpb = ',num2str(hpbss(1)));
    stre = strcat('Peak 2 area = ',num2str(Areamat(2)));
    strf = strcat('Peak 2 hpb = ',num2str(hpbss(2)));
    strreal = {stra, strd, strc};
    text(.15,20,strreal, 'FontName', 'Times New Roman', 'FontSize', 12,'Color', 'red')
    strreal2 ={strb, strf,stre};
    text(.15,5,strreal2, 'FontName', 'Times New Roman', 'FontSize', 12, 'Color', 'blue')

    % Now the arrow for freq and amp
    freq = strcat('fn = ',num2str(matrix(2,1)));
    amp = strcat('Amp = ',num2str(matrix(2,2)));
    strreal2 = {freq, amp};
    text(freqs(2)-5 ,amps(2) +1,strreal2,'FontName', 'Times New Roman', 'FontSize', 12, 'Color', 'blue' )


end