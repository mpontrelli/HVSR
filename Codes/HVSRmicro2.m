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
    
    
  %% Author: Marshall Pontrelli
  %Co-author Justin Reyes
  %Summer 2019
  
%% Start
function [ahatf, fax_HzN, taxstat] = HVSRmicro2(Vfname, NSfname, EWfname, fs, statname, varargin)
%% parse inputs
% create Input Parser object
p = inputParser;
% add inputs to the scheme
defaultwindowlen = 40;
defaultnumwin = 10;
defaultwindis = 25;
defaultlowbound = 1/defaultwindowlen;
defaultupbound = ceil(fs/2)-1;
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
addParameter(p, 'plots', 'no', @ischar);
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
plots = p.Results.plots;
fsmin = fs;
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

%% Create a time series plot (Output 1)]
if strcmp(plots, 'yes') == 1
    timeseriesplot(xNS,xV,xEW, fs, sav, outpath)
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
%Steidl et al. 1994
xHmatrix = xNSmatrix + 1i.*xEWmatrix; 

%% window the data
win = hann(sampnum)';
for i = 1:numwin
    xVmatrix(i,:) = xVmatrix(i,:).*win;
    xHmatrix(i,:) = xHmatrix(i,:).*win;
end
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
if strcmp(plots, 'yes') == 1
    individmagrespplot(fax_HzN, XHmatrix2, XVmatrix2, fs, lowbound, outpath, sav)
end


%% Average the un-smoothed magnitude responses
[ahatfhorz, sigmahorz, confinthighhorz, confintlowhorz] =  wavav(XHmatrix2);
[ahatfvert, sigmavert, confinthighvert, confintlowvert] =  wavav(XVmatrix2);

%% Plot averaged unfiltered magnitude responses (OUTPUT 3)
if strcmp(plots, 'yes') == 1
    averagedmagrespplot(fax_HzN, ahatfhorz, ahatfvert, fs,confinthighhorz, confintlowhorz, confinthighvert, confintlowvert, lowbound, outpath, sav)
end


%% compute smoothed magnitude responses
width = .5; %width for triangle moving average filter in hz
window = ceil((N/fs)*width); %width for triangle moving average filter in samples where 20 is the number of Hz on your x-axis
for iii = 1:numwin
    XVmatrix3(iii,:) = smooth(XVmatrix2(iii,:),window);
    XHmatrix3(iii,:) = smooth(XHmatrix2(iii,:),window);
end

%% Plot individual, smoothed magnitude responses (OUTPUT 4)
if strcmp(plots, 'yes') == 1
    individmagrespplot(fax_HzN, XHmatrix3, XVmatrix3, fs, lowbound, outpath, sav)
end


%% Average the smoothed magnitude responses
[ahatfhorz, sigmahorz, confinthighhorz, confintlowhorz] =  wavav(XHmatrix3);
[ahatfvert, sigmavert, confinthighvert, confintlowvert] =  wavav(XVmatrix3);

%% Plot averaged, smoothed magnitude responses (OUTPUT 5)
if strcmp(plots, 'yes') == 1
    averagedmagrespplot(fax_HzN, ahatfhorz, ahatfvert, fs,confinthighhorz, confintlowhorz, confinthighvert, confintlowvert, lowbound, outpath, sav)
end


%% Compute the HVSR
for iii = 1:numwin
    [H_V(iii,:)] = HV(XHmatrix3(iii,:),XVmatrix3(iii,:));
end

%% average the HVSR
[ahatf, sigma, confinthigh, confintlow] =  wavav(H_V);

%% Plot the HVSR (OUTPUT 6)
if strcmp(plots, 'yes') == 1
    HVSRmicroplot(fax_HzN, ahatf, fs, confinthigh, confintlow, statname, lowbound, upbound, outpath, sav, TTF)
end
% compute statistics on HVSR
[matrix, matrix1, peakind,ahatf1,newfaxhz1] = peakiden(ahatf, fax_HzN, lowbound, upbound, fsmin);
[taxstat] = specratstat(peakind, matrix, matrix1, ahatf1, newfaxhz1, sigma, statname, lowbound);
% save
if strcmp(sav, 'yes') == 1
    taxstat2 = cell2table(taxstat);
    writetable(taxstat2, strcat(outpath, '\', 'statistics.txt'));
 end

end


%% End