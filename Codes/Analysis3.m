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
  
close all
clear all
%% Start
path = cd;
nm1 = '2021231085651006';
extrafile = 'yes';
nm = '2021231090000006';
cd(strcat('C:\Users\',getenv('username'),'\Desktop\HVSR\Codes'));
Vfname = strcat(path,'\',nm,'_STATI_1_1.sac');
NSfname = strcat(path,'\',nm,'_STATI_1_2.sac');
EWfname = strcat(path,'\',nm,'_STATI_1_3.sac');
fs = 100;
statname = ''; 


wincut = 'no'; % do you need to cut a window?
which_win = [1]; % which window?
windowlen = 40;
numwin = 20;
windis = 1;
lowbound = 0.1;
upbound =  49;
LowCorner = 0.1;
HighCorner = fs/2 - 1;
Npoles = 4;
width = 0.5;
TTF = 'no';
outpath = strcat(path,'\');
sav = 'yes'; 
Allplots =  'yes';
Timeplot = 'no';
IUMagplot = 'no';
AUMagplot = 'no';
IFMagplot = 'no';
AFMagplot = 'no';
HVSRplot = 'no';
Filterplot = 'no';
sampnum = windowlen*fs; 
windisnum = windis*fs;

%% read files and convert into vectors
    
[~,~,ext] = fileparts(Vfname);
%.sacBinary
if strcmp(ext, '.sac') == 1 % use "rdmseed" if file is in miniseed format
    [V] = ReadSacBinaryFile(Vfname); %vertical
    [NS] = ReadSacBinaryFile(NSfname); %North-south
    [EW] = ReadSacBinaryFile(EWfname); %East-West
end

if strcmp(extrafile, 'yes') == 1
    Vfname1 = strcat(path,'\',nm1,'_STATI_1_1.sac');
    NSfname1 = strcat(path,'\',nm1,'_STATI_1_2.sac');
    EWfname1 = strcat(path,'\',nm1,'_STATI_1_3.sac');
    [~,~,ext] = fileparts(Vfname1);
    %.sacBinary
    if strcmp(ext, '.sac') == 1 % use "rdmseed" if file is in miniseed format
        [V1] = ReadSacBinaryFile(Vfname1); %vertical
        [NS1] = ReadSacBinaryFile(NSfname1); %North-south
        [EW1] = ReadSacBinaryFile(EWfname1); %East-West
    end
    V = vertcat(V1,V);
    EW = vertcat(EW1,EW);
    NS = vertcat(NS1,NS);
end

if strcmp(Allplots, 'yes') == 1 || strcmp(Timeplot, 'yes') == 1
    timeseriesplot(NS,EW,V, fs,'sav','no','outpath',outpath)
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
    timeseriesplot(NS,EW,V, fs,'sav','yes','outpath',outpath)
end
%% Create structure
data.NS = NS;
data.EW = EW;
data.V = V;
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

if strcmp(wincut,'yes') == 1
    Vmatrix([which_win],:) = [];
    NSmatrix([which_win],:) = [];
    EWmatrix([which_win],:) = [];
    numwin = numwin - length(which_win);
end
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
% if strcmp(Allplots, 'yes') == 1 || strcmp(IUMagplot, 'yes') == 1
%     individmagrespplot(fax_HzN, Hmatrix2, Vmatrix2, fs, lowbound, outpath, sav)
% end

%% Average the un-smoothed magnitude responses
[ahatfhorz, sigmahorz, confinthighhorz, confintlowhorz] =  wavav(Hmatrix2);
[ahatfvert, sigmavert, confinthighvert, confintlowvert] =  wavav(Vmatrix2);

%% Plot averaged unfiltered magnitude responses (OUTPUT 3)
% if strcmp(Allplots, 'yes') == 1 || strcmp(AUMagplot, 'yes') == 1
%     averagedmagrespplot(fax_HzN, ahatfhorz, ahatfvert, fs,confinthighhorz, confintlowhorz, confinthighvert, confintlowvert, lowbound, outpath, sav)
% end

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
    individmagrespplot(fax_HzN, Hmatrix3, Vmatrix3, fs, 'Horizontal','Vertical', lowbound,upbound, outpath, sav, 'individ_mag_resp_plot')
end

%% Average the smoothed magnitude responses
[ahatfhorz, sigmahorz, confinthighhorz, confintlowhorz] =  wavav(Hmatrix3);
[ahatfvert, sigmavert, confinthighvert, confintlowvert] =  wavav(Vmatrix3);

%% Plot averaged, smoothed magnitude responses (OUTPUT 5)
if strcmp(Allplots, 'yes') == 1 || strcmp(AFMagplot, 'yes') == 1
    averagedmagrespplot(fax_HzN, ahatfhorz, ahatfvert, fs,confinthighhorz, confintlowhorz, confinthighvert, confintlowvert,'Horizontal','Vertical',lowbound,upbound, outpath, sav,'Averaged_mag_resps')
end

%% Compute the HVSR
H_V = zeros(numwin, N_2);
for iii = 1:numwin
    [H_V(iii,:)] = HV(Hmatrix3(iii,:),Vmatrix3(iii,:));
end

%% Plot all the HVSRs
if strcmp(Allplots, 'yes') == 1 || strcmp(HVSRplot, 'yes') == 1
    individplot(H_V, fax_HzN, statname, sav,lowbound,upbound,outpath,'individ_HVSR')
end
%% average the HVSR
[ahatf, sigma, confinthigh, confintlow] =  wavav(H_V);
data.ahatf = ahatf;
data.sigma = sigma;
data.confinthigh = confinthigh;
data.confintlow = confintlow;
data.freq = fax_HzN;
%% Plot the HVSR (OUTPUT 6)
if strcmp(Allplots, 'yes') == 1 || strcmp(HVSRplot, 'yes') == 1
    HVSRmicroplot(fax_HzN, ahatf, confinthigh, confintlow, statname, lowbound, upbound, outpath, sav, TTF)
end

%% Plot sigma
sig = figure;
plot(fax_HzN, sigma)
xlabel('Frequency (Hz)','FontSize', 18)
ylabel('\sigma','FontSize', 18)
title(strcat(statname,{' '},'\sigma'), 'FontSize', 18)
set(gca,'XScale','log','FontName', 'Times New Roman', 'FontSize', 18)
xlim([fax_HzN(lowbound) fax_HzN(upbound)])
%xlim([fax_HzN(1) 40])
ylim([0 1])
xticks([0.1 1 10])
xticklabels({'0.1','1','10', num2str(upbound)})

grid on 
box on

saveas(sig, strcat(outpath, '\', 'sigma.jpg'));
%% compute statistics on HVSR
lowbound1 = 4;
upbound1 = 39;
[~, upbound1] = min(abs(fax_HzN - upbound1));
[peakfreqs, peakamps, hpb, f1s, f2s, Areamat, proms, amps, peakind2, freqs, sigs, I1s, I2s] = peakiden(ahatf, fax_HzN, sigma, lowbound, upbound1);

[~, I] = max(amps);

datamat = vertcat(freqs,amps,proms,hpb,sigs);
datamatmax = datamat(:,I);
datamatmax = datamatmax';
datamat = datamat';
data.fn = datamatmax(1);
data.A = datamatmax(2);
data.prom = datamatmax(3);
data.hpb = datamatmax(4);
data.sigma = datamatmax(5);

%% now plot, make the figure and set the base
HVf = figure;
hold on
confidenceinterval=shadedplot(fax_HzN(1:length(fax_HzN)), confinthigh(1:length(fax_HzN)), confintlow(1:length(fax_HzN)),[.9,.9,.9],[1,1,1]);
hold on
ETF = plot(fax_HzN(1 :length(fax_HzN)), ahatf(1:length(fax_HzN)), 'Color', [0 0.30196 0.6588] , 'Linewidth', 1.5);
xlabel('Frequency (Hz)','FontSize', 18)
ylabel('Amplification','FontSize', 18)
title(strcat(statname), 'FontSize', 18)
xlim([fax_HzN(lowbound) fax_HzN(upbound)])
%xlim([fax_HzN(1) 40])
ylim([0.1 100])
xticks([0.1 1 10])
xticklabels({'0.1','1','10', num2str(upbound)})
yticks([0.1 1 10 100])
yticklabels({'0.1','1','10', '100'})
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
text(0.5 ,12,str, 'FontName', 'Times New Roman', 'FontSize', 18,'Color', 'black', 'HorizontalAlignment', 'right', 'EdgeColor','k','BackgroundColor', 'w')
saveas(HVf, strcat(outpath, '\', 'HVSR_w_data.jpg'));
%% Now save data
save(strcat(outpath,'data.mat'),'data')

%% Now get back into the right directory
cd(path)