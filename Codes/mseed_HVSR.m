%% read DMC data
close all
clear all
%% necessary inputs
% Station
filename = 'C:\Users\mpontr01\Box\People\Marshall and Jeremy\NEHRP\DATA\weston_stations\NY_NJ\2020_001\bd_C30L_1577836800.msd';
% Filter
LowCorner = 0.1;
HighCorner = 10;
Npoles = 4;
Filterplot = 'no';

%%
x = rdmseed(filename);
[NS, EW, V] = openmseed(x);

%%
fs = x(1).SampleRate;
statname = x.StationIdentifierCode;


%% make all vectors the same length
a = length(NS);
b = length(EW);
c = length(V);
d = min([a,b,c]);
NS = NS(1:d);
EW = EW(1:d);
V = V(1:d);

%% Now filter
[NS] = Butter2(NS, fs, 'LowCorner', LowCorner, 'HighCorner', HighCorner, 'Npoles', Npoles , 'Filterplot', Filterplot);
[EW] = Butter2(EW, fs, 'LowCorner', LowCorner, 'HighCorner', HighCorner, 'Npoles', Npoles , 'Filterplot', Filterplot);
[V] = Butter2(V, fs, 'LowCorner', LowCorner, 'HighCorner', HighCorner, 'Npoles', Npoles , 'Filterplot', Filterplot);

%% plot time series
% find the maximum value to make bounds for plotting
a = max(abs(NS));
b = max(abs(EW));
c = max(abs(V));
cc = [a,b,c];
d = max(cc);

%% create a time vector and plot time series
time = ((1:length(NS))/fs)/3600;

figure;
    
% NS
subplot(3,1,1)
plot(time,NS)
title('NS')
ylabel('Counts')
xlim([0 length(NS)/(fs*3600)])
ylim([-d d])
grid on 
box on
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);

% EW
subplot(3,1,2)
plot(time,EW)
title('EW')
ylabel('Counts')
xlim([0 length(EW)/(fs*3600)])
ylim([-d d])
grid on 
box on
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

% V
subplot(3,1,3)
plot(time,V)
title('V')
ylabel('Counts')
xlabel('Time (Hours)')
xlim([0 length(V)/(fs*3600)])
ylim([-d d])
grid on 
box on
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

%% integrate for displacement
alpha = 0.25; %coefficients for Newmark-beta method
beta = 0.5;

% NS
dt = 1/fs;
dNS = zeros(1,length(NS)-1);
for i = 1:length(NS) - 1
    dNS(i+1) = dNS(i) + dt*((1-beta)*NS(i) + beta*NS(i+1));
end
dNS = dNS - mean(dNS);
[dNS] = Butter2(dNS, fs);
dNS = dNS*1000000; % to microns

% EW
dEW = zeros(1,length(EW)-1);
for i = 1:length(EW) - 1
    dEW(i+1) = dEW(i) + dt*((1-beta)*EW(i) + beta*EW(i+1));
end
dEW = dEW - mean(dEW);
[dEW] = Butter2(dEW, fs);
dEW = dEW*1000000; % to microns

% V
dV = zeros(1,length(V)-1);
for i = 1:length(V) - 1
    dV(i+1) = dV(i) + dt*((1-beta)*V(i) + beta*V(i+1));
end
dV = dV - mean(dV);
[dV] = Butter2(dV, fs);
dV = dV*1000000; % to microns

%% plot displacement
% find the maximum value to make bounds for plotting
a = max(abs(dNS));
b = max(abs(dEW));
c = max(abs(dV));
cc = [a,b, c];
d = max(cc);
figure;
    
% NS
subplot(3,1,1)
plot(time,dNS)
title('NS')
ylabel('Microns')
xlabel('Time (Hours)')
xlim([0 length(NS)/(fs*3600)])
ylim([-d d])
grid on 
box on
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);

% EW
subplot(3,1,2)
plot(time,dEW)
title('EW')
ylabel('Microns')
xlabel('Time (Hours)')
xlim([0 length(EW)/(fs*3600)])
ylim([-d d])
grid on 
box on
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

% V
subplot(3,1,3)
plot(time,dV)
title('V')
ylabel('Microns')
xlabel('Time (Hours)')
xlim([0 length(EW)/(fs*3600)])
ylim([-d d])
grid on 
box on
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

%% spectrogram inputs (edit this to make plots better)
% nfft = 2048/4;
% noverlap = round(nfft*0.50);
% win = hanning(nfft);
% 
% figure
% subplot(1,3,1)
% spectrogram(NS,win,noverlap,nfft,fs,'yaxis');
% title('NS')
% ylim([0.1 100])
% set(gca,'YScale', 'log','FontName', 'Times New Roman', 'FontSize', 14)
% 
% subplot(1,3,2)
% spectrogram(EW,win,noverlap,nfft,fs,'yaxis');
% title('EW')
% ylim([0.1 100])
% ylabel('')
% set(gca,'YScale', 'log','FontName', 'Times New Roman', 'FontSize', 14)
% 
% subplot(1,3,3)
% spectrogram(V,win,noverlap,nfft,fs,'yaxis');
% title('V')
% ylim([0.1 100])
% ylabel('')
% set(gca,'YScale', 'log','FontName', 'Times New Roman', 'FontSize', 14)
%set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
%% Now do magnitude responses and compare
% inputs
numwin = 100;
windowlen = 200;
windis = 20;
sampnum = windowlen*fs; 
windisnum = windis*fs;

k = [1,fs];
NSmatrix = zeros(numwin, sampnum);
EWmatrix = zeros(numwin, sampnum);
Vmatrix = zeros(numwin, sampnum);
for iii = 1:numwin
    NSmatrix(iii,:) = NS((k(1)):(k(2)*windowlen+k(1))-1);
    EWmatrix(iii,:) = EW((k(1)):(k(2)*windowlen+k(1))-1);
    Vmatrix(iii,:) = V((k(1)):(k(2)*windowlen+k(1))-1);
    k(1) = k(1)+ sampnum + windisnum;
end


%% window the data
win = hann(sampnum)';
for i = 1:numwin
    NSmatrix(i,:) = NSmatrix(i,:).*win;
    EWmatrix(i,:) = EWmatrix(i,:).*win;
    Vmatrix(i,:) = Vmatrix(i,:).*win;
end

%% Compute unfiltered magnitude responses
for iii = 1:numwin
    NSmatrix(iii,:) = 4*abs(fft(NSmatrix(iii,:)))/sampnum;
    EWmatrix(iii,:) = 4*abs(fft(EWmatrix(iii,:)))/sampnum;
    Vmatrix(iii,:) = 4*abs(fft(Vmatrix(iii,:)))/sampnum;
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
    mat1 = NSmatrix(iii,:);
    mat2 = EWmatrix(iii, :);
    mat3 = Vmatrix(iii, :);
    NSmatrix2(iii,:) = mat1(1 : N_2);
    EWmatrix2(iii,:) = mat2(1 : N_2);
    Vmatrix2(iii,:) = mat3(1 : N_2);
end

%% create upbound and lowbound in terms of sample number
upbound = HighCorner;
lowbound = LowCorner;
[~, lowbound] = min(abs(fax_HzN - lowbound));
[~, upbound] = min(abs(fax_HzN - upbound));
    
%% plot individual unfiltered magnitude responses (OUTPUT 2)
sav = 'no';
outpath = 'no';
individualunfiltered = figure;
subplot(1,3,1)
plot(fax_HzN, NSmatrix2, 'Color',  'k' , 'Linewidth', .5);
title('NS', 'FontSize', 20)
xlabel('Frequency (Hz)','FontSize', 18)
ylabel('Amplification','FontSize', 18)
set(gca,'YScale', 'log', 'XScale', 'log','FontName', 'Times New Roman', 'FontSize', 14)
xlim([0.1 10])
%ylim([10e-16 10e-5])
xticks([0.1 1 10])
xticklabels({'0.1','1','10'})
% yticks([0.1 1 10 100 1000])
% yticklabels({'0.1','1','10', '100', '1000'})
grid on 
box on

subplot(1,3,2)
plot(fax_HzN, EWmatrix2, 'Color',  'k' , 'Linewidth', .5);
title('EW', 'FontSize', 20)
xlabel('Frequency (Hz)','FontSize', 18)
ylabel('Amplification','FontSize', 18)
set(gca,'YScale', 'log', 'Xscale', 'log','FontName', 'Times New Roman', 'FontSize', 14)
xlim([0.1 10])
%ylim([10e-16 10e-5])
xticks([0.1 1 10])
xticklabels({'0.1','1','10'})
% yticks([0.1 1 10 100 1000])
% yticklabels({'0.1','1','10', '100', '1000'})
grid on 
box on

subplot(1,3,3)
plot(fax_HzN, Vmatrix2, 'Color',  'k' , 'Linewidth', .5);
title('V', 'FontSize', 20)
xlabel('Frequency (Hz)','FontSize', 18)
ylabel('Amplification','FontSize', 18)
set(gca,'YScale', 'log', 'Xscale', 'log','FontName', 'Times New Roman', 'FontSize', 14)
xlim([0.1 10])
%ylim([10e-16 10e-5])
xticks([0.1 1 10])
xticklabels({'0.1','1','10'})
% yticks([0.1 1 10 100 1000])
% yticklabels({'0.1','1','10', '100', '1000'})
grid on 
box on

%makes figure full screen
set(gcf, 'Units', 'Normalized', 'OuterPosition', [1, 0.5, 1, 0.5]);


%% Average the un-smoothed magnitude responses
[ahatfNS, sigmaNS, confinthighNS, confintlowNS] =  wavav(NSmatrix2);
[ahatfEW, sigmaEW, confinthighEW, confintlowEW] =  wavav(EWmatrix2);
[ahatfV, sigmaV, confinthighV, confintlowV] =  wavav(Vmatrix2);

%% now plot
averageunfiltered = figure;
subplot(1,3,1)
hold on
confidenceinterval=shadedplot(fax_HzN, confinthighNS, confintlowNS,[.9,.9,.9],'k');
hold on
plot(fax_HzN, ahatfNS, 'Color', [0 0.30196 0.6588] , 'Linewidth', 1.5);
title('NS')
xlabel('Frequency (Hz)')
ylabel('Amplification')
set(gca,'YScale', 'log', 'Xscale', 'log','FontName', 'Times New Roman', 'FontSize', 14)
xlim([0.1 10])
%ylim([10e-16 10e-5])
xticks([0.1 1 10])
xticklabels({'0.1','1','10'})
% yticks([0.1 1 10 100 1000])
% yticklabels({'0.1','1','10', '100', '1000'})
grid on 
box on

%makes figure full screen
set(gcf, 'Units', 'Normalized', 'OuterPosition', [1, 0.5, 1, 0.5]);


subplot(1,3,2)
hold on
confidenceinterval=shadedplot(fax_HzN, confinthighEW, confintlowEW,[.9,.9,.9],'k');
hold on
plot(fax_HzN, ahatfEW, 'Color', [0 0.30196 0.6588] , 'Linewidth', 1.5);
title('EW')
xlabel('Frequency (Hz)')
ylabel('Amplification')
set(gca,'YScale', 'log', 'Xscale', 'log','FontName', 'Times New Roman', 'FontSize', 14)
xlim([0.1 10])
%ylim([10e-16 10e-5])
xticks([0.1 1 10])
xticklabels({'0.1','1','10'})
% yticks([0.1 1 10 100 1000])
% yticklabels({'0.1','1','10', '100', '1000'})
grid on 
box on

subplot(1,3,3)
hold on
confidenceinterval=shadedplot(fax_HzN, confinthighV, confintlowV,[.9,.9,.9],'k');
hold on
plot(fax_HzN, ahatfV, 'Color', [0 0.30196 0.6588] , 'Linewidth', 1.5);
title('V')
xlabel('Frequency (Hz)')
ylabel('Amplification')
set(gca,'YScale', 'log', 'Xscale', 'log','FontName', 'Times New Roman', 'FontSize', 14)
xlim([0.1 10])
%ylim([10e-16 10e-5])
xticks([0.1 1 10])
xticklabels({'0.1','1','10'})
% yticks([0.1 1 10 100 1000])
% yticklabels({'0.1','1','10', '100', '1000'})
grid on 
box on

%makes figure full screen
set(gcf, 'Units', 'Normalized', 'OuterPosition', [1, 0.5, 1, 0.5]);

%% Now compute geometric mean
horz = sqrt(NSmatrix2.*EWmatrix2);
%% now compute HVSR
HV = horz ./Vmatrix2;
figure
plot(fax_HzN,HV, 'Color',  'k' , 'Linewidth', .5)
title('HVSR')
xlabel('Frequency (Hz)')
ylabel('Amplification')
set(gca,'YScale', 'log', 'Xscale', 'log','FontName', 'Times New Roman', 'FontSize', 14)
xlim([0.1 10])
ylim([0.1 100])
xticks([0.1 1 10])
xticklabels({'0.1','1','10'})
yticks([0.1 1 10 100])
yticklabels({'0.1','1','10', '100'})
grid on 
box on

%makes figure full screen
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

%% average HVSR
[ahatf, sigma, confinthigh, confintlow] =  wavav(HV);

%% Now plot
figure
hold on
confidenceinterval=shadedplot(fax_HzN, confinthigh, confintlow,[.9,.9,.9],'k');
hold on
plot(fax_HzN, ahatf, 'Color', [0 0.30196 0.6588] , 'Linewidth', 1.5);
title('HVSR')
xlabel('Frequency (Hz)')
ylabel('Amplification')
set(gca,'YScale', 'log', 'Xscale', 'log','FontName', 'Times New Roman', 'FontSize', 14)
xlim([0.1 10])
ylim([0.1 100])
xticks([0.1 1 10])
xticklabels({'0.1','1','10'})
yticks([0.1 1 10 100])
yticklabels({'0.1','1','10', '100'})
grid on 
box on

%% Now smooth
width = 0.15;
window = ceil((N/fs)*width); %width for smoothing filter in samples where 20 is the number of Hz on your x-axis
Vmatrix3 = zeros(numwin, N_2);
horz2 = zeros(numwin, N_2);
for iii = 1:numwin
    Vmatrix3(iii,:) = smooth(Vmatrix2(iii,:),window);
    horz2(iii,:) = smooth(horz(iii,:),window);
end

%% Now compute HVSR
HV = horz2 ./Vmatrix3;
figure
plot(fax_HzN,HV, 'Color',  'k' , 'Linewidth', .5)
title('HVSR')
xlabel('Frequency (Hz)')
ylabel('Amplification')
set(gca,'YScale', 'log', 'Xscale', 'log','FontName', 'Times New Roman', 'FontSize', 14)
xlim([0.1 10])
ylim([0.1 100])
xticks([0.1 1 10])
xticklabels({'0.1','1','10'})
yticks([0.1 1 10 100])
yticklabels({'0.1','1','10', '100'})
grid on 
box on

%makes figure full screen
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

%% average HVSR
[ahatf, sigma, confinthigh, confintlow] =  wavav(HV);

%% Now plot
figure
hold on
confidenceinterval=shadedplot(fax_HzN, confinthigh, confintlow,[.9,.9,.9],'k');
hold on
plot(fax_HzN, ahatf, 'Color', [0 0.30196 0.6588] , 'Linewidth', 1.5);
title(statname)
xlabel('Frequency (Hz)')
ylabel('Amplification')
set(gca,'YScale', 'log','FontName', 'Times New Roman', 'FontSize', 14)
%set(gca,'FontName', 'Times New Roman', 'FontSize', 14)
xlim([0.1 10])
ylim([0.1 100])
% xticks([0.1 1 10])
% xticklabels({'0.1','1','10'})
yticks([0.1 1 10 100])
yticklabels({'0.1','1','10', '100'})
grid on 
box on

%% Now plot sigma
figure
plot(fax_HzN, sigma)
title('Sigma')
xlabel('Frequency (Hz)')
ylabel('Sigma')
set(gca, 'Xscale', 'log','FontName', 'Times New Roman', 'FontSize', 14)
xlim([0.1 10])


%% Now interpolate onto a vector. This makes the hpb measurement better
% newfaxhz = linspace(0,10,100000);
% ahatf = interp1(fax_HzN, ahatf, newfaxhz);
% confinthigh = interp1(fax_HzN, confinthigh, newfaxhz);
% confintlow = interp1(fax_HzN, confintlow, newfaxhz);
% sigma = interp1(fax_HzN, sigma, newfaxhz);
% fax_HzN = newfaxhz;
%% Now quantify the shape of the peaks
[peakfreqs, peakamps, hpb, f1s, f2s, Areamat, proms, amps, peakind2, freqs, sigs, I1s, I2s] = peakiden(ahatf, fax_HzN, sigma, lowbound, upbound);


%% Now set some conditions to classify the peak
[~, freq_class] = sort(freqs, 'descend');
[~, amp_class] = sort(amps, 'descend');
[~, prom_class] = sort(proms, 'descend');
[~, hpb_class] = sort(hpb, 'descend');
[~, area_class] = sort(Areamat, 'descend');
[~, sig_class] = sort(sigs, 'descend');
classmatrix = vertcat(freq_class, amp_class, prom_class, hpb_class, area_class, sig_class);

%% Now save peak data
datamat = vertcat(freqs,amps,proms,hpb,Areamat,sigs);

%% now plot, make the figure and set the base
figure
xlim([0.1 10])
ylim([0.1 40])
xticks([.1 1 10])
xticklabels({'0.1', '1', '10'})
yticks([1 10 40])
yticklabels({ '1','10', '40'})
xlabel('Frequency (Hz)')
ylabel('Amplification')
title(statname)
set(gca,'YScale', 'log','XScale','log', 'FontName', 'Times New Roman', 'FontSize', 18)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
grid on
box on
hold on

%% start with the confidence interval
confidenceinterval=shadedplot(fax_HzN(lowbound:length(fax_HzN)), confinthigh(lowbound:length(confinthigh)), confintlow(lowbound:length(confintlow)),[.9,.9,.9],[1 1 1]);
hold on
    
%% Plot filled in peaks color coordinated based on their area
if length(amps) > 0
    for i = 1:length(amps)
        hold on
        freq_ar = peakfreqs{i};
        amp_ar = peakamps{i};
        if i == 1
            col = 'b';
        end
        if i == 2
            col = 'r';
        end
        if i == 3
            col = 'g';
        end
        if i > 3
            col = 'k'
        end
        fill(freq_ar, amp_ar, col, 'LineStyle','none','FaceAlpha',0.5)
        hold on
    end
end
if length(amps) > 0
    for i = 1:length(amps) 
        freq_ar = peakfreqs{i};
        amp_ar = peakamps{i};
        ampd = amp_ar(I1s(i));
        ampdd = amp_ar(I2s(i));
        plot([f1s(i),f2s(i)],[ampd, ampdd], 'c', 'LineWidth', 2)
        hold on
    end
end
hold on


%% now plot fundamental resonance
if length(amps) > 0
    for i = 1:length(f1s)
        line([freqs(i),freqs(i)],[0.1, 40],'LineStyle', '--', 'color','k')
        hold on
    end
end



%% Now the hpb-prominence cross
hold on
if length(amps) > 0
    for i = 1:length(amps)
        plot([freqs(i),freqs(i)],[amps(i), amps(i) - proms(i)], 'c', 'Linewidth', 2)
        hold on
    end
end

%% Now plot ahatf
plot(fax_HzN,ahatf, 'LineWidth', 2, 'Color', [0 0.5 0])
%% Now print some important info
if length(amps) > 0
    if length(amps) ==1
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
        if fax_HzN(q) < 1
            text(xlim(2)-5,ylim(2)-25,str, 'FontName', 'Times New Roman', 'FontSize', 18,'Color', 'black', 'HorizontalAlignment', 'right', 'EdgeColor','k','BackgroundColor', 'w')
        else
            text(xlim(1)+0.3,ylim(2)-25,str, 'FontName', 'Times New Roman', 'FontSize', 18,'Color', 'black', 'HorizontalAlignment', 'right', 'EdgeColor','k','BackgroundColor', 'w')
        end
        else
            [~, I] = max(amps);
        a = "Max peak freq =" + " "  + num2str(freqs(I),3);
        b = strcat("Amp =" + " "  + num2str(amps(I),3));
        c = strcat("HPB =" + " "  +num2str(hpb(I),3));
        d = strcat("Prom =" + " "  +num2str(proms(I),3));
        e = strcat("Area =" + " "  +num2str(Areamat(I),3));
        f = strcat("\sigma =" + " "  +num2str(sigs(I),3));
        str = {a, b, c, d, e, f};
        ylim=get(gca,'ylim');
        xlim=get(gca,'xlim');
        q = find(ahatf == amps(I));
        if fax_HzN(q) < 1
            text(xlim(2)-5,ylim(2)-25,str, 'FontName', 'Times New Roman', 'FontSize', 18,'Color', 'black', 'HorizontalAlignment', 'right', 'EdgeColor','k','BackgroundColor', 'w')
        else
            text(xlim(1)+0.3,ylim(2)-25,str, 'FontName', 'Times New Roman', 'FontSize', 18,'Color', 'black', 'HorizontalAlignment', 'right', 'EdgeColor','k','BackgroundColor', 'w')
        end

    end
else
    str = 'No peaks';
    text(0.3,15,str, 'FontName', 'Times New Roman', 'FontSize', 24,'Color', 'black', 'HorizontalAlignment', 'right', 'EdgeColor','k','BackgroundColor', 'w')

end
% str = {strcat('Peak 1 prom = ',{' '},num2str(proms(1))), strcat('Peak 2 prom = ',{' '},num2str(proms(2)))};
% stra = str{1}{1};
% strb = str{2}{1};
% strc = strcat('Peak 1 area = ',num2str(Areamat(1)));
% strd = strcat('Peak 1 hpb = ',num2str(hpbss(1)));
% stre = strcat('Peak 2 area = ',num2str(Areamat(2)));
% strf = strcat('Peak 2 hpb = ',num2str(hpbss(2)));
% strreal = {stra, strd, strc};
% text(.15,20,strreal, 'FontName', 'Times New Roman', 'FontSize', 12,'Color', 'red')
% strreal2 ={strb, strf,stre};
% text(.15,5,strreal2, 'FontName', 'Times New Roman', 'FontSize', 12, 'Color', 'blue')
% 
% % Now the arrow for freq and amp
% freq = strcat('fn = ',num2str(matrix(2,1)));
% amp = strcat('Amp = ',num2str(matrix(2,2)));
% strreal2 = {freq, amp};
% text(freqs(2)-5 ,amps(2) +1,strreal2,'FontName', 'Times New Roman', 'FontSize', 12, 'Color', 'blue' )
