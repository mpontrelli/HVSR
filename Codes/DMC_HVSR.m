close all
clear all


% Before running this file the first time, execute the following statement
%javaaddpath IRIS-WS-2.0.18.jar

network = 'TA';
statname = 'I64A';
starttime = '2015-02-11 15:18:00';
endtime = '2015-02-12 16:19:05';

[sampletimes,trace1,waveformfilt] = getDMCData(starttime,endtime ,statname,network,'HHZ');
fs = trace1.sampleRate;
V = trace1.data;
%V = V(1:length(V)-1);
%%
% INPUTS
windowlen = 40;
numwin = 20;
windis = 25;
TTF = 'no';
outpath = 'no';
sav = 'no';
lowbound = 0.2;
upbound = fs/2 -1;
Allplots = 'yes';
Timeplot = 'no';
IUMagplot = 'no';
AUMagplot = 'no';
IFMagplot = 'no';
AFMagplot = 'no';
HVSRplot = 'no';
Filterplot = 'no';
LowCorner = 0.1;
HighCorner = fs/2 - 1;
Npoles = 4;
width = 0.5;
%turn windows into samples for windowing calculations
sampnum = windowlen*fs; 
windisnum = windis*fs;


[~,trace1,waveformfilt] = getDMCData(starttime,endtime ,statname,network,'HHE');
EW = trace1.data;
%EW = EW(1:length(EW)-1);

[~,trace1,waveformfilt] = getDMCData(starttime,endtime ,statname,network,'HHN');
NS = trace1.data;
%NS = NS(1:length(NS)-1);


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
[~, upbound] = (min(abs(fax_HzN - upbound)));
    
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

% %% Plot individual, smoothed magnitude responses (OUTPUT 4)
% if strcmp(Allplots, 'yes') == 1 || strcmp(IFMagplot, 'yes') == 1
%     individmagrespplot(fax_HzN, Hmatrix3, Vmatrix3, fs, lowbound, outpath, sav)
% end

%% Average the smoothed magnitude responses
[ahatfhorz, sigmahorz, confinthighhorz, confintlowhorz] =  wavav(Hmatrix3);
[ahatfvert, sigmavert, confinthighvert, confintlowvert] =  wavav(Vmatrix3);

%% Plot averaged, smoothed magnitude responses (OUTPUT 5)
% if strcmp(Allplots, 'yes') == 1 || strcmp(AFMagplot, 'yes') == 1
%     averagedmagrespplot(fax_HzN, ahatfhorz, ahatfvert, fs,confinthighhorz, confintlowhorz, confinthighvert, confintlowvert, lowbound, outpath, sav)
% end

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