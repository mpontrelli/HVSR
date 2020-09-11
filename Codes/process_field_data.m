%% Process data in the field all the way through HVSR processing steps 

close all
clear all

%% Inputs
fs = 100;
windowlen = 60;
numwin = 13;
windis = 25;
lowbound1 = 1.5;
upbound1 =  3;
LowCorner = 0.1;
HighCorner = fs/2 - 1;
Npoles = 4;
width = 0.5;
Filterplot = 'no'; % toggle off filter plot so it doesn't plot response three times
sampnum = windowlen*fs; 
windisnum = windis*fs;
sav = 'no';
outpath = 'no';
TTF = 'no';
statname = '3c';

%% Load data
load('C:\Users\Marshall\Box\Geohazards Research Group\2020_7_1\Charles_test\3c\3c');
V = data.V;
NS = data.NS;
EW = data.EW;

%% Now filter
[V] = Butter2(V, fs, 'LowCorner', LowCorner, 'HighCorner', HighCorner, 'Npoles', Npoles , 'Filterplot', Filterplot);
[NS] = Butter2(NS, fs, 'LowCorner', LowCorner, 'HighCorner', HighCorner, 'Npoles', Npoles , 'Filterplot', Filterplot);
[EW] = Butter2(EW, fs, 'LowCorner', LowCorner, 'HighCorner', HighCorner, 'Npoles', Npoles , 'Filterplot', Filterplot);

%% Create time series plot
timeseriesplot(NS,V,EW, fs)

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
    NSmatrix(iii,:) = 4*abs(fft(NSmatrix(iii,:)))/sampnum;
    EWmatrix(iii,:) = 4*abs(fft(EWmatrix(iii,:)))/sampnum;
end

%% compute geometric mean
Hmatrix = sqrt(EWmatrix.*NSmatrix);

%% Compute frequency axis
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
[~, lowbound] = min(abs(fax_HzN - lowbound1));
[~, upbound] = min(abs(fax_HzN - upbound1));

%% plot individual unfiltered magnitude responses (OUTPUT 2)
%individmagrespplot(fax_HzN, Hmatrix2, Vmatrix2, fs, lowbound, outpath, sav)


%% Average the un-smoothed magnitude responses
[ahatfhorz, sigmahorz, confinthighhorz, confintlowhorz] =  wavav(Hmatrix2);
[ahatfvert, sigmavert, confinthighvert, confintlowvert] =  wavav(Vmatrix2);

%% Plot averaged unfiltered magnitude responses (OUTPUT 3)
%averagedmagrespplot(fax_HzN, ahatfhorz, ahatfvert, fs,confinthighhorz, confintlowhorz, confinthighvert, confintlowvert, lowbound, outpath, sav)

%% compute smoothed magnitude responses
window = ceil((N/fs)*width); %width for smoothing filter in samples where 20 is the number of Hz on your x-axis
Hmatrix3 = zeros(numwin, N_2);
Vmatrix3 = zeros(numwin, N_2);
for iii = 1:numwin
    Vmatrix3(iii,:) = smooth(Vmatrix2(iii,:),window);
    Hmatrix3(iii,:) = smooth(Hmatrix2(iii,:),window);
end

%% Plot individual, smoothed magnitude responses (OUTPUT 4)
individmagrespplot(fax_HzN, Hmatrix3, Vmatrix3, fs, lowbound, outpath, sav)

%% Average the smoothed magnitude responses
[ahatfhorz, sigmahorz, confinthighhorz, confintlowhorz] =  wavav(Hmatrix3);
[ahatfvert, sigmavert, confinthighvert, confintlowvert] =  wavav(Vmatrix3);

%% Plot averaged, smoothed magnitude responses (OUTPUT 5)
averagedmagrespplot(fax_HzN, ahatfhorz, ahatfvert, fs,confinthighhorz, confintlowhorz, confinthighvert, confintlowvert, lowbound, outpath, sav)

%% Compute the HVSR
H_V = zeros(numwin, N_2);
for iii = 1:numwin
    [H_V(iii,:)] = HV(Hmatrix3(iii,:),Vmatrix3(iii,:));
end

%% plot individual HVSR
individplot(H_V, fax_HzN, '3c')
%% average the HVSR
[ahatf, sigma, confinthigh, confintlow] =  wavav(H_V);

%% Now quantify the shape of the peaks
[peakfreqs, peakamps, hpb, f1s, f2s, Areamat, proms, amps, peakind2, freqs, sigs, I1s, I2s] = peakiden(ahatf, fax_HzN, sigma, lowbound, upbound);
[~, I] = max(amps);

datamat = vertcat(freqs,amps,proms,hpb,sigs);
datamat = datamat(:,1);
datamat = datamat';

%% now plot, make the figure and set the base
HVSR = figure;
hold on
confidenceinterval=shadedplot(fax_HzN(1:length(fax_HzN)), confinthigh(1:length(fax_HzN)), confintlow(1:length(fax_HzN)),[.9,.9,.9],[1,1,1]);
hold on
ETF = plot(fax_HzN(1 :length(fax_HzN)), ahatf(1:length(fax_HzN)), 'Color', [0 0.30196 0.6588] , 'Linewidth', 1.5);
xlabel('Frequency (Hz)','FontSize', 18)
ylabel('Amplification','FontSize', 18)
title(strcat(statname), 'FontSize', 18)
set(gca,'YScale', 'log', 'XScale', 'log','FontName', 'Times New Roman', 'FontSize', 18)
xlim([fax_HzN(lowbound) fax_HzN(upbound)])
%xlim([fax_HzN(1) 40])
ylim([2 20])
xticks([2 3])
xticklabels({'2','3'})
yticks([0.1 2 10 100])
yticklabels({'0.1','2','10', '100'})
%makes figure full screen
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

grid on 
box on
hold on

%% now plot fundamental resonance
line([freqs(1),freqs(1)],[0.1, 100],'LineStyle', '--', 'color','k')
hold on


%% Now the hpb-prominence cross
hold on
plot([freqs(1),freqs(1)],[amps(1), amps(1) - proms(1)], 'c', 'Linewidth', 2)
hold on

freq_ar = peakfreqs{1};
amp_ar = peakamps{1};
ampd = amp_ar(I1s(1));
ampdd = amp_ar(I2s(1));
plot([f1s(1),f2s(1)],[ampd, ampdd], 'c', 'LineWidth', 2)
hold on

a = "Max peak freq =" + " "  + num2str(freqs,3);
b = strcat("Amp =" + " "  + num2str(amps,3));
c = strcat("HPB =" + " "  +num2str(hpb,3));
d = strcat("Prom =" + " "  +num2str(proms,3));
e = strcat("Area =" + " "  +num2str(Areamat,3));
f = strcat("\sigma =" + " "  +num2str(sigs,3));
str = {a, b, c, d, e, f};
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
q = find(ahatf == amps);

text(2.75,12,str, 'FontName', 'Times New Roman', 'FontSize', 18,'Color', 'black', 'HorizontalAlignment', 'right', 'EdgeColor','k','BackgroundColor', 'w')





