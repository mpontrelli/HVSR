%% Sac_binary_2_mat
% reads three sac binary files and saves them as .mat files

function [data] = Sac_binary_2_mat(NSfname, EWfname, Vfname, varargin) %name, outpath)

   %% parse inputs
    % create Input Parser object
    p = inputParser;
    
    % Required inputs
    addRequired(p, 'NSfname',@ischar);
    addRequired(p, 'EWfname',@ischar);
    addRequired(p, 'Vfname',@ischar);

    % Optional inputs
    addParameter(p, 'name', 'none', @ischar);
    addParameter(p, 'sav', 'no', @ischar);
    addParameter(p, 'outpath', 'no', @ischar);
    addParameter(p, 'windowlen', 40, @isnumeric);
    addParameter(p, 'numwin', 10, @isnumeric);
    addParameter(p, 'windis', 25, @isnumeric);
    addParameter(p, 'width', 0.5, @isnumeric);
    
    % parse the inputs
    parse(p, NSfname, EWfname, Vfname, varargin{:})
    % set varibales from the parse
    name = p.Results.name;
    sav = p.Results.sav;
    outpath = p.Results.outpath;
    windowlen = p.Results.windowlen;
    numwin = p.Results.numwin;
    windis = p.Results.windis;
    width = p.Results.width;
    data.proc.window_length = windowlen;
    data.proc.number_windows = numwin;
    data.proc.window_distance = windis;
    data.proc.filt_width = width;
    
    %% read the files
    [NS] = ReadSacBinaryFile(NSfname); %North-south
    [EW] = ReadSacBinaryFile(EWfname); %East-West
    [V, delta, stla, stlo] = ReadSacBinaryFile(Vfname); %vertical
    fs = 1/delta;
    sampnum = windowlen*fs; 
    windisnum = windis*fs;
    
    %% filter them
    [NS] = Butter2(NS, fs);
    [EW] = Butter2(EW, fs);
    [V] = Butter2(V,fs);
   
    %% save into a structure
    data.NS = NS;
    data.EW = EW;
    data.V = V;
    data.fs = fs;
    data.lat = stla;
    data.lon = stlo;
    data.npts = length(NS);
    data.name = name;
    
    %% Window data
    %Window the data with 'numwin' non-overlapping windows of 'windowlen' secs and 
    % 'windis' secs apart
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
        NSmatrix(i,:) = NSmatrix(i,:).*win;
        EWmatrix(i,:) = EWmatrix(i,:).*win;
    end
    
    %% Compute unfiltered magnitude responses
    % Compute the fft for each data window
    for iii = 1:numwin
        Vmatrix(iii,:) = 4*abs(fft(Vmatrix(iii,:)))/sampnum;
        Hmatrix(iii,:) = 4*abs(fft(Hmatrix(iii,:)))/sampnum;
        NSmatrix(iii,:) = 4*abs(fft(NSmatrix(iii,:)))/sampnum;
        EWmatrix(iii,:) = 4*abs(fft(EWmatrix(iii,:)))/sampnum;
    end
    
    %% Compute the frequency - axis
    N = sampnum;
    fax_binsN = (0 : N-1); %samples in NS component
    fax_HzN1 = fax_binsN*fs/N; %frequency axis NS (Hz)
    N_2 = ceil(N/2); %half magnitude spectrum
    fax_HzN = fax_HzN1(1 : N_2);
    data.HVSR.freq = fax_HzN;
    
    for iii = 1:numwin
        Vmat = Vmatrix(iii,:);
        Hmat = Hmatrix(iii, :);
        NSmat = NSmatrix(iii, :);
        EWmat = EWmatrix(iii, :);
        Vmatrix2(iii,:) = Vmat(1 : N_2);
        Hmatrix2(iii,:) = Hmat(1 : N_2);
        NSmatrix2(iii,:) = NSmat(1 : N_2);
        EWmatrix2(iii,:) = EWmat(1 : N_2);
    end

    %% compute smoothed magnitude responses
    window = ceil((N/fs)*width); %width for smoothing filter in samples where 20 is the number of Hz on your x-axis
    for iii = 1:numwin
        Vmatrix3(iii,:) = smooth(Vmatrix2(iii,:),window);
        Hmatrix3(iii,:) = smooth(Hmatrix2(iii,:),window);
        NSmatrix3(iii,:) = smooth(NSmatrix2(iii,:),window);
        EWmatrix3(iii,:) = smooth(EWmatrix2(iii,:),window);
    end
    
    %% Compute the HVSR
    for iii = 1:numwin
        [HV_comp(iii,:)] = HV(Hmatrix3(iii,:),Vmatrix3(iii,:));
        [HV_NS(iii,:)] = HV(NSmatrix3(iii,:),Vmatrix3(iii,:));
        [HV_EW(iii,:)] = HV(EWmatrix3(iii,:),Vmatrix3(iii,:));
    end
    
    %% average the HVSR
    [ahatf_comp, sigma_comp, confinthigh_comp, confintlow_comp] =  wavav(HV_comp);
    data.HVSR.complex.ahatf = ahatf_comp;
    data.HVSR.complex.confinthigh = confinthigh_comp;
    data.HVSR.complex.confintlow = confintlow_comp;
    
    [ahatf_NS, sigma_NS, confinthigh_NS, confintlow_NS] =  wavav(HV_NS);
    data.HVSR.NS.ahatf = ahatf_NS;
    data.HVSR.NS.confinthigh = confinthigh_NS;
    data.HVSR.NS.confintlow = confintlow_NS;
    [ahatf_EW, sigma_EW, confinthigh_EW, confintlow_EW] =  wavav(HV_EW);
    data.HVSR.EW.ahatf = ahatf_EW;
    data.HVSR.EW.confinthigh = confinthigh_EW;
    data.HVSR.EW.confintlow = confintlow_EW;
    

    % compute statistics on HVSR
    % complex
    [M, I] = max(ahatf_comp);
    data.HVSR.complex.amp = M;
    data.HVSR.complex.fn = fax_HzN(I);
    
    % North south
    [M, I] = max(ahatf_NS);
    data.HVSR.NS.amp = M;
    data.HVSR.NS.fn = fax_HzN(I);
    
    % East West
    [M, I] = max(ahatf_EW);
    data.HVSR.EW.amp = M;
    data.HVSR.EW.fn = fax_HzN(I);
    
    %% save the structure
    if strcmp(sav, 'yes') == 1
        save(strcat(outpath, '\',name, '.mat'), 'data')
    end
    
end