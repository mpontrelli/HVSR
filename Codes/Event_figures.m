%% Event_figures
%% Author: Marshall Pontrelli
% Date: 4/27/2020

%%
close all
clear all

%%
CE32 = load('C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Processed_data2\CE32\CE3220170919181440.mat');
TP13 = load('C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Processed_data2\TP13\TP1320170919181440.mat');
%% inputs
fs = CE32.data.meta.instrument.fs;
time = CE32.data.processing.filtereddata.time_orig;
timelong = CE32.data.processing.filtereddata.time;
%% time series
% Acceleration CE32

NS = CE32.data.processing.filtereddata.acceleration.NS.waveform_orig;
EW = CE32.data.processing.filtereddata.acceleration.EW.waveform_orig;
V = CE32.data.processing.filtereddata.acceleration.EW.waveform_orig;
timeseriesplot(NS,EW, V, fs,'yaxislabel', 'm/s^2')

% Acceleration TP13
NS = TP13.data.processing.filtereddata.acceleration.NS.waveform_orig;
EW = TP13.data.processing.filtereddata.acceleration.EW.waveform_orig;
V = TP13.data.processing.filtereddata.acceleration.EW.waveform_orig;
timeseriesplot(NS,EW, V, fs,'yaxislabel', 'm/s^2')

%% time series with velocity and displacement CE32

NS = CE32.data.processing.filtereddata.acceleration.NS.waveform_orig;
[PGA, PGV, PGD, v, d] =  waveform_integrate(NS, fs);

% find the maximum value to make bounds for plotting
dd = max(abs(NS));

% Acceleration
subplot(3,1,1)
plot(time,NS)
ylabel('m/s^2')
xlim([0 length(NS)/fs])
ylim([-dd dd])
grid on 
box on
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);

% Velocity
subplot(3,1,2)
plot(time,v)
ylabel('m/s')
xlim([0 length(NS)/fs])
ylim([-dd dd])
grid on 
box on
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);

% displacement
subplot(3,1,3)
plot(time,d);
xlabel('Time (secs)')
ylabel('m')
xlim([0 length(NS)/fs])
ylim([-dd dd])
grid on 
box on
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);
    
%makes figure full screen and font in times
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
%% time series with velocity and displacement TP13
time = TP13.data.processing.filtereddata.time_orig;
NS = TP13.data.processing.filtereddata.acceleration.NS.waveform_orig;
[PGA, PGV, PGD, v, d] =  waveform_integrate(NS, fs);

% find the maximum value to make bounds for plotting
dd = max(abs(NS));

% Acceleration
subplot(3,1,1)
plot(time,NS)
ylabel('m/s^2')
xlim([0 length(NS)/fs])
ylim([-dd dd])
grid on 
box on
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);

% Velocity
subplot(3,1,2)
plot(time,v)
ylabel('m/s')
xlim([0 length(NS)/fs])
ylim([-dd dd])
grid on 
box on
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);

% displacement
subplot(3,1,3)
plot(time,d);
xlabel('Time (secs)')
ylabel('m')
xlim([0 length(NS)/fs])
ylim([-dd dd])
grid on 
box on
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);
    
%makes figure full screen and font in times

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
%% Arias intensity
% CE32
fs = CE32.data.meta.instrument.fs;
time = CE32.data.processing.filtereddata.time_orig;
NS_IaX = CE32.data.processing.filtereddata.acceleration.NS.arias.no_norm;
EW_IaX = CE32.data.processing.filtereddata.acceleration.EW.arias.no_norm;
V_IaX = CE32.data.processing.filtereddata.acceleration.V.arias.no_norm;
rot_IaX = CE32.data.processing.filtereddata.acceleration.rotated.arias.no_norm;
NS_arias_int = CE32.data.processing.filtereddata.acceleration.NS.arias.intensity;
NS_arias_D5 = CE32.data.processing.filtereddata.acceleration.NS.arias.D5;
NS_arias_D75 = CE32.data.processing.filtereddata.acceleration.NS.arias.D75;
NS_arias_D95 = CE32.data.processing.filtereddata.acceleration.NS.arias.D95;
NS_arias_D595 = CE32.data.processing.filtereddata.acceleration.NS.arias.D595;
NS_arias_D575 = CE32.data.processing.filtereddata.acceleration.NS.arias.D575;
NS_arias_rate = CE32.data.processing.filtereddata.acceleration.NS.arias.rate;

figure
NSplot = plot(time, NS_IaX, 'linewidth', 1.5);
hold on
EWplot = plot(time, EW_IaX,'linewidth', 1.5);
hold on
Vplot = plot(time, V_IaX,'linewidth', 1.5);
xlim([25 200])
ylim([0 0.5])
% 
% sigdurFN=line([NS_arias_D5 NS_arias_D5], [0 max(IaX)],'color','g');
% line([NS_arias_D95 NS_arias_D95], [0 max(IaX)],'color','g')
legend([NSplot, EWplot,Vplot],'NS', 'EW', 'V', 'location', 'east')
xlabel('Time (secs)')
ylabel('m/s')
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);
grid on

% TP13
fs = TP13.data.meta.instrument.fs;
time = TP13.data.processing.filtereddata.time_orig;
NS_IaX = TP13.data.processing.filtereddata.acceleration.NS.arias.no_norm;
EW_IaX = TP13.data.processing.filtereddata.acceleration.EW.arias.no_norm;
V_IaX = TP13.data.processing.filtereddata.acceleration.V.arias.no_norm;
rot_IaX = TP13.data.processing.filtereddata.acceleration.rotated.arias.no_norm;
NS_arias_int = TP13.data.processing.filtereddata.acceleration.NS.arias.intensity;
NS_arias_D5 = TP13.data.processing.filtereddata.acceleration.NS.arias.D5;
NS_arias_D75 = TP13.data.processing.filtereddata.acceleration.NS.arias.D75;
NS_arias_D95 = TP13.data.processing.filtereddata.acceleration.NS.arias.D95;
NS_arias_D595 = TP13.data.processing.filtereddata.acceleration.NS.arias.D595;
NS_arias_D575 = TP13.data.processing.filtereddata.acceleration.NS.arias.D575;
NS_arias_rate = TP13.data.processing.filtereddata.acceleration.NS.arias.rate;

figure
NSplot = plot(time, NS_IaX, 'linewidth', 1.5);
hold on
EWplot = plot(time, EW_IaX,'linewidth', 1.5);
hold on
Vplot = plot(time, V_IaX,'linewidth', 1.5);
xlim([25 200])
ylim([0 0.5])

% 
% sigdurFN=line([NS_arias_D5 NS_arias_D5], [0 max(IaX)],'color','g');
% line([NS_arias_D95 NS_arias_D95], [0 max(IaX)],'color','g')
legend([NSplot, EWplot,Vplot],'NS', 'EW', 'V', 'location', 'east')
xlabel('Time (secs)')
ylabel('m/s')
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);
grid on


%% Response spectra
period = CE32.data.processing.filtereddata.period_vec;
spectra_NS = CE32.data.processing.filtereddata.acceleration.NS.spectra.waveform;
spectra_EW = CE32.data.processing.filtereddata.acceleration.EW.spectra.waveform;
spectra_V = CE32.data.processing.filtereddata.acceleration.V.spectra.waveform;

subplot(3,1,1)
plot(period,spectra_NS)
ylabel('m/s^2')
grid on
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);
box on
xlim([period(2) 5])
ylim([0 4])

subplot(3,1,2)
plot(period,spectra_EW)
ylabel('m/s^2')
grid on
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);
box on
xlim([period(2) 5])
ylim([0 4])

subplot(3,1,3)
plot(period,spectra_V)
xlabel('Period (secs)')
ylabel('m/s^2')
grid on
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);
box on
xlim([period(2) 5])
ylim([0 4])
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

% TP13
period = TP13.data.processing.filtereddata.period_vec;
spectra_NS = TP13.data.processing.filtereddata.acceleration.NS.spectra.waveform;
spectra_EW = TP13.data.processing.filtereddata.acceleration.EW.spectra.waveform;
spectra_V = TP13.data.processing.filtereddata.acceleration.V.spectra.waveform;
figure
subplot(3,1,1)
plot(period,spectra_NS)
ylabel('m/s^2')
grid on
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);
box on
xlim([period(2) 5])
ylim([0 4])

subplot(3,1,2)
plot(period,spectra_EW)
ylabel('m/s^2')
grid on
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);
box on
xlim([period(2) 5])
ylim([0 4])

subplot(3,1,3)
plot(period,spectra_V)
xlabel('Period (secs)')
ylabel('m/s^2')
grid on
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);
box on
xlim([period(2) 5])
ylim([0 4])
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);


%% HVSR
%NS
freq = CE32.data.processing.filtereddata.freq_vec;
HV_NS = CE32.data.processing.filtereddata.acceleration.NS.HVSR.smooth.HV;
amps = CE32.data.processing.filtereddata.acceleration.NS.HVSR.smooth.freq_amp(:,2);
freqs = CE32.data.processing.filtereddata.acceleration.NS.HVSR.smooth.freq_amp(:,1);
f1s = CE32.data.processing.filtereddata.acceleration.NS.HVSR.smooth.hpb(:,2);
f2s = CE32.data.processing.filtereddata.acceleration.NS.HVSR.smooth.hpb(:,3);
proms = CE32.data.processing.filtereddata.acceleration.NS.HVSR.smooth.prom;
%%
areafreqs  = CE32.data.processing.filtereddata.acceleration.NS.HVSR.smooth.areafreqs;
areaamps  = CE32.data.processing.filtereddata.acceleration.NS.HVSR.smooth.areaamps;
%%

figure
% first plot the areas

for i = 1:length(areafreqs)
    freq_ar = areafreqs{i};
    amp_ar = areaamps{i};
    fill(freq_ar, amp_ar, 'r')
    hold on
end
hold on
plot(freq,HV_NS)
hold on

for i = 1:length(f1s)
    dddd = f1s(i);
    ampd = find(freq == dddd);
    ampd = HV_NS(ampd);
    ddddd = f2s(i);
    ampdd = find(freq == ddddd);
    ampdd = HV_NS(ampdd);
    plot([dddd,ddddd],[ampd, ampdd], 'k')
    hold on
end
hold on
for i = 1:length(amps)
    plot([freqs(i),freqs(i)],[amps(i), amps(i) - proms(i)], 'k')
    hold on
end
grid on
box on
xlim([0.1 5])
ylim([0.1 40])
xticks([.1 1 5])
xticklabels({'0.1', '1', '5'})
yticks([1 10 30])
yticklabels({ '1','10', '30'})
xlabel('Frequency (Hz)')
ylabel('Amplification')
title('NS')
set(gca,'YScale', 'log', 'XScale', 'log','FontName', 'Times New Roman', 'FontSize', 14)
%%

%EW
freq = CE32.data.processing.filtereddata.freq_vec;
HV_NS = CE32.data.processing.filtereddata.acceleration.EW.HVSR.smooth.HV;
amps = CE32.data.processing.filtereddata.acceleration.EW.HVSR.smooth.freq_amp(:,2);
freqs = CE32.data.processing.filtereddata.acceleration.EW.HVSR.smooth.freq_amp(:,1);
f1s = CE32.data.processing.filtereddata.acceleration.EW.HVSR.smooth.hpb(:,2);
f2s = CE32.data.processing.filtereddata.acceleration.EW.HVSR.smooth.hpb(:,3);
proms = CE32.data.processing.filtereddata.acceleration.EW.HVSR.smooth.prom;

figure
plot(freq,HV_NS)
hold on

for i = 1:length(f1s)
    dddd = f1s(i);
    ampd = find(freq == dddd);
    ampd = HV_NS(ampd);
    ddddd = f2s(i);
    ampdd = find(freq == ddddd);
    ampdd = HV_NS(ampdd);
    plot([dddd,ddddd],[ampd, ampdd], 'k')
    hold on
end
hold on
for i = 1:length(amps)
    plot([freqs(i),freqs(i)],[amps(i), amps(i) - proms(i)], 'k')
    hold on
end
grid on
box on
xlim([0.1 5])
ylim([0.1 40])
xticks([.1 1 5])
xticklabels({'0.1', '1', '5'})
yticks([1 10 30])
yticklabels({ '1','10', '30'})
xlabel('Frequency (Hz)')
ylabel('Amplification')
title('EW')
set(gca,'YScale', 'log', 'XScale', 'log','FontName', 'Times New Roman', 'FontSize', 14)

%complex
freq = CE32.data.processing.filtereddata.freq_vec;
HV_NS = CE32.data.processing.filtereddata.acceleration.complex.HVSR.smooth.HV;
amps = CE32.data.processing.filtereddata.acceleration.complex.HVSR.smooth.freq_amp(:,2);
freqs = CE32.data.processing.filtereddata.acceleration.complex.HVSR.smooth.freq_amp(:,1);
f1s = CE32.data.processing.filtereddata.acceleration.complex.HVSR.smooth.hpb(:,2);
f2s = CE32.data.processing.filtereddata.acceleration.complex.HVSR.smooth.hpb(:,3);
proms = CE32.data.processing.filtereddata.acceleration.complex.HVSR.smooth.prom;

figure
plot(freq,HV_NS)
hold on

for i = 1:length(f1s)
    dddd = f1s(i);
    ampd = find(freq == dddd);
    ampd = HV_NS(ampd);
    ddddd = f2s(i);
    ampdd = find(freq == ddddd);
    ampdd = HV_NS(ampdd);
    plot([dddd,ddddd],[ampd, ampdd], 'k')
    hold on
end
hold on
for i = 1:length(amps)
    plot([freqs(i),freqs(i)],[amps(i), amps(i) - proms(i)], 'k')
    hold on
end
grid on
box on
xlim([0.1 5])
ylim([0.1 40])
xticks([.1 1 5])
xticklabels({'0.1', '1', '5'})
yticks([1 10 35])
yticklabels({ '1','10', '30'})
xlabel('Frequency (Hz)')
ylabel('Amplification')
title('Complex')
set(gca,'YScale', 'log', 'XScale', 'log','FontName', 'Times New Roman', 'FontSize', 14)


%% SSR
% NS
CE32_NS = CE32.data.processing.filtereddata.acceleration.NS.mag_resps.smooth;
TP13_NS = TP13.data.processing.filtereddata.acceleration.NS.mag_resps.smooth;
NS_SSR = CE32_NS./TP13_NS;
%% create upbound and lowbound values
lowbound = 0.1;
upbound = 5;
[~, lowbound] = min(abs(freq - lowbound));
[~, upbound] = min(abs(freq - upbound));
[matrix, matrix1, peakind,ahatf1,newfaxhz, peakfreqs, peakamps, Areamat1] = peakiden(NS_SSR, freq, lowbound, upbound);
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
amps = matrix(:,2);
freqs = matrix(:,1);
f1s = hpb1(:,2);
f2s = hpb1(:,3);
proms = matrix1;
%%
figure
plot(freq,NS_SSR)
hold on
for i = 1:length(f1s)
    dddd = f1s(i);
    ampd = find(freq == dddd);
    ampd = NS_SSR(ampd);
    ddddd = f2s(i);
    ampdd = find(freq == ddddd);
    ampdd = NS_SSR(ampdd);
    plot([dddd,ddddd],[ampd, ampdd], 'k')
    hold on
end
hold on
for i = 1:length(amps)
    plot([freqs(i),freqs(i)],[amps(i), amps(i) - proms(i)], 'k')
    hold on
end
xlim([0.1 5])
ylim([0.1 40])
xticks([.1 1 5])
xticklabels({'0.1', '1', '5'})
yticks([1 10 35])
yticklabels({ '1','10', '30'})
xlabel('Frequency (Hz)')
ylabel('Amplification')
title('NS')
set(gca,'YScale', 'log', 'XScale', 'log','FontName', 'Times New Roman', 'FontSize', 14)
grid on
box on
hold on

% EW
CE32_NS = CE32.data.processing.filtereddata.acceleration.EW.mag_resps.smooth;
TP13_NS = TP13.data.processing.filtereddata.acceleration.EW.mag_resps.smooth;
NS_SSR = CE32_NS./TP13_NS;
% create upbound and lowbound values
lowbound = 0.1;
upbound = 5;
[~, lowbound] = min(abs(freq - lowbound));
[~, upbound] = min(abs(freq - upbound));
[matrix, matrix1, peakind,ahatf1,newfaxhz1, peakfreqs, peakamps, Areamat] = peakiden(NS_SSR, freq, lowbound, upbound);
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
amps = matrix(:,2);
freqs = matrix(:,1);
f1s = hpb1(:,2);
f2s = hpb1(:,3);
proms = matrix1;
figure
plot(freq,NS_SSR)
hold on
for i = 1:length(f1s)
    dddd = f1s(i);
    ampd = find(freq == dddd);
    ampd = NS_SSR(ampd);
    ddddd = f2s(i);
    ampdd = find(freq == ddddd);
    ampdd = NS_SSR(ampdd);
    plot([dddd,ddddd],[ampd, ampdd], 'k')
    hold on
end
hold on
for i = 1:length(amps)
    plot([freqs(i),freqs(i)],[amps(i), amps(i) - proms(i)], 'k')
    hold on
end
xlim([0.1 5])
ylim([0.1 40])
xticks([.1 1 5])
xticklabels({'0.1', '1', '5'})
yticks([1 10 35])
yticklabels({ '1','10', '30'})
xlabel('Frequency (Hz)')
ylabel('Amplification')
title('EW')
set(gca,'YScale', 'log', 'XScale', 'log','FontName', 'Times New Roman', 'FontSize', 14)
grid on
box on
hold on

% complex
CE32_NS = CE32.data.processing.filtereddata.acceleration.complex.mag_resps.smooth;
TP13_NS = TP13.data.processing.filtereddata.acceleration.complex.mag_resps.smooth;
NS_SSR = CE32_NS./TP13_NS;
% create upbound and lowbound values
lowbound = 0.1;
upbound = 5;
[~, lowbound] = min(abs(freq - lowbound));
[~, upbound] = min(abs(freq - upbound));
[matrix, matrix1, peakind,ahatf1,newfaxhz1, peakfreqs, peakamps, Areamat] = peakiden(NS_SSR, freq, lowbound, upbound);
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
amps = matrix(:,2);
freqs = matrix(:,1);
f1s = hpb1(:,2);
f2s = hpb1(:,3);
proms = matrix1;
figure
plot(freq,NS_SSR)
hold on
for i = 1:length(f1s)
    dddd = f1s(i);
    ampd = find(freq == dddd);
    ampd = NS_SSR(ampd);
    ddddd = f2s(i);
    ampdd = find(freq == ddddd);
    ampdd = NS_SSR(ampdd);
    plot([dddd,ddddd],[ampd, ampdd], 'k')
    hold on
end
hold on
for i = 1:length(amps)
    plot([freqs(i),freqs(i)],[amps(i), amps(i) - proms(i)], 'k')
    hold on
end
xlim([0.1 5])
ylim([0.1 40])
xticks([.1 1 5])
xticklabels({'0.1', '1', '5'})
yticks([1 10 35])
yticklabels({ '1','10', '30'})
xlabel('Frequency (Hz)')
ylabel('Amplification')
title('Complex')
set(gca,'YScale', 'log', 'XScale', 'log','FontName', 'Times New Roman', 'FontSize', 14)
grid on
box on
hold on

%% Now plot averages
% NS
clear hpb1
CE32 = load('C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Shape_statistics\CE32');
ahatf = CE32.shapedata.NS.ahatf;
confinthigh = CE32.shapedata.NS.confinthigh;
confintlow = CE32.shapedata.NS.confintlow;
areafreqs = CE32.shapedata.NS.area_freq;
areaamps = CE32.shapedata.NS.area_amp;
area = CE32.shapedata.NS.Areamat;
%%
lowbound = 0.1;
upbound = 5;
[~, lowbound] = min(abs(freq - lowbound));
[~, upbound] = min(abs(freq - upbound));
[matrix, matrix1, peakind,ahatf1,newfaxhz1, peakfreqs, peakamps, Areamat] = peakiden(ahatf', freq, lowbound, upbound);
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
amps = matrix(:,2);
freqs = matrix(:,1);
f1s = hpb1(:,2);
f2s = hpb1(:,3);
proms = matrix1;

figure
hold on
confidenceinterval=shadedplot(freq(lowbound:length(freq)), confinthigh(lowbound:length(confinthigh)), confintlow(lowbound:length(confintlow)),[.9,.9,.9],[1 1 1]);
hold on
for i = 1:length(areafreqs)
    freq_ar = areafreqs{i};
    amp_ar = areaamps{i};
    fill(freq_ar, amp_ar, 'r')
    alpha(.5)
    hold on
end
hold on


plot(freq,ahatf, 'LineWidth',2, 'Color',[0 0.5 0])
for i = 1:length(f1s)
    dddd = f1s(i);
    ampd = find(freq == dddd);
    ampd = ahatf(ampd);
    ddddd = f2s(i);
    ampdd = find(freq == ddddd);
    ampdd = ahatf(ampdd);
    plot([dddd,ddddd],[ampd, ampdd], 'k')
    hold on
end
hold on
for i = 1:length(amps)
    plot([freqs(i),freqs(i)],[amps(i), amps(i) - proms(i)], 'k')
    hold on
end
hold on
xlim([0.1 5])
ylim([0.1 40])
xticks([.1 1 5])
xticklabels({'0.1', '1', '5'})
yticks([1 10 100])
yticklabels({ '1','10', '100'})
xlabel('Frequency (Hz)')
ylabel('Amplification')
title('NS')
set(gca,'YScale', 'log', 'XScale', 'log','FontName', 'Times New Roman', 'FontSize', 14)
grid on
box on
hold on

%% EW
clear hpb1
CE32 = load('C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Shape_statistics\CE32');
ahatf = CE32.shapedata.EW.ahatf;
confinthigh = CE32.shapedata.EW.confinthigh;
confintlow = CE32.shapedata.EW.confintlow;
lowbound = 0.1;
upbound = 5;
[~, lowbound] = min(abs(freq - lowbound));
[~, upbound] = min(abs(freq - upbound));
[matrix, matrix1, peakind,ahatf1,newfaxhz1, peakfreqs, peakamps, Areamat] = peakiden(ahatf', freq, lowbound, upbound);
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
amps = matrix(:,2);
freqs = matrix(:,1);
f1s = hpb1(:,2);
f2s = hpb1(:,3);
proms = matrix1;

figure

hold on
confidenceinterval=shadedplot(freq(lowbound:length(freq)), confinthigh(lowbound:length(confinthigh)), confintlow(lowbound:length(confintlow)),[.9,.9,.9],[1 1 1]);
hold on
plot(freq,ahatf)
for i = 1:length(f1s)
    dddd = f1s(i);
    ampd = find(freq == dddd);
    ampd = ahatf(ampd);
    ddddd = f2s(i);
    ampdd = find(freq == ddddd);
    ampdd = ahatf(ampdd);
    plot([dddd,ddddd],[ampd, ampdd], 'k')
    hold on
end
hold on
for i = 1:length(amps)
    plot([freqs(i),freqs(i)],[amps(i), amps(i) - proms(i)], 'k')
    hold on
end
hold on
xlim([0.1 5])
ylim([0.1 40])
xticks([.1 1 5])
xticklabels({'0.1', '1', '5'})
yticks([1 10 100])
yticklabels({ '1','10', '100'})
xlabel('Frequency (Hz)')
ylabel('Amplification')
title('EW')
set(gca,'YScale', 'log', 'XScale', 'log','FontName', 'Times New Roman', 'FontSize', 14)
grid on
box on
hold on

% Complex
clear hpb1
CE32 = load('C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Shape_statistics\CE32');
ahatf = CE32.shapedata.complex.ahatf;
confinthigh = CE32.shapedata.complex.confinthigh;
confintlow = CE32.shapedata.complex.confintlow;
lowbound = 0.1;
upbound = 5;
[~, lowbound] = min(abs(freq - lowbound));
[~, upbound] = min(abs(freq - upbound));
[matrix, matrix1, peakind,ahatf1,newfaxhz1, peakfreqs, peakamps, Areamat] = peakiden(ahatf', freq, lowbound, upbound);
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
amps = matrix(:,2);
freqs = matrix(:,1);
f1s = hpb1(:,2);
f2s = hpb1(:,3);
proms = matrix1;

figure

hold on
confidenceinterval=shadedplot(freq(lowbound:length(freq)), confinthigh(lowbound:length(confinthigh)), confintlow(lowbound:length(confintlow)),[.9,.9,.9],[1 1 1]);
hold on
plot(freq,ahatf)
for i = 1:length(f1s)
    dddd = f1s(i);
    ampd = find(freq == dddd);
    ampd = ahatf(ampd);
    ddddd = f2s(i);
    ampdd = find(freq == ddddd);
    ampdd = ahatf(ampdd);
    plot([dddd,ddddd],[ampd, ampdd], 'k')
    hold on
end
hold on
for i = 1:length(amps)
    plot([freqs(i),freqs(i)],[amps(i), amps(i) - proms(i)], 'k')
    hold on
end
hold on
xlim([0.1 5])
ylim([0.1 100])
xticks([.1 1 5])
xticklabels({'0.1', '1', '5'})
yticks([1 10 100])
yticklabels({ '1','10', '100'})
xlabel('Frequency (Hz)')
ylabel('Amplification')
title('Complex')
set(gca,'YScale', 'log', 'XScale', 'log','FontName', 'Times New Roman', 'FontSize', 14)
grid on
box on
hold on

%% Now do SSR
codepath = 'C:\Users\mpontr01\Desktop\HVSR\Codes';
datapath = 'C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Processed_data2\CE32\';
datapath2 = 'C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Processed_data2\TP13\';
cd(datapath)
eventlist = dir;
eventlist = eventlist(3:length(eventlist));
cd(datapath2)
refeventlist = dir;
refeventlist = refeventlist(3:length(refeventlist));
counter = 0;
SR = [];
for i = 1 : length(eventlist)
    event = eventlist(i).name;
    filename = strcat(datapath,event);
    CE32 = load(filename);
    mag = CE32.data.processing.filtereddata.acceleration.NS.mag_resps.smooth;
    eventid = event(5:18);
    for ii = 1:length(refeventlist)
        refevent = refeventlist(ii).name;
        refid = refevent(5:18);
        if strcmp(eventid,refid) == 1
            counter = counter + 1;
            filename = strcat(datapath2,refevent);
            TP13 = load(filename);
            refmag = TP13.data.processing.filtereddata.acceleration.NS.mag_resps.smooth;
            if length(refmag) == 100000
                refmag = refmag(1:50000);
            end
            SR(counter,:) = mag./refmag;
        end
    end
end

freq = CE32.data.processing.filtereddata.freq_vec;
cd(codepath)
statname = 'CE32';
[ahatf, sigma, confinthigh, confintlow] =  wavav(SR);
lowbound = 0.1;
upbound = 5;
[~, lowbound] = min(abs(freq - lowbound));
[~, upbound] = min(abs(freq - upbound));
%%
[matrix, matrix1, peakind,ahatf1,newfaxhz1, peakfreqs, peakamps, Areamat] = peakiden(ahatf', freq, lowbound, upbound);
[taxstat, sigma1] = specratstat(peakind, matrix, matrix1, ahatf1, newfaxhz1, sigma, statname,lowbound, upbound);
%%
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
amps = matrix(:,2);
freqs = matrix(:,1);
f1s = hpb1(:,2);
f2s = hpb1(:,3);
proms = matrix1;
%%
figure

hold on
confidenceinterval=shadedplot(freq(lowbound:length(freq)), confinthigh(lowbound:length(confinthigh)), confintlow(lowbound:length(confintlow)),[.9,.9,.9],[1 1 1]);
hold on
plot(freq,ahatf)
for i = 1:length(f1s)
    dddd = f1s(i);
    ampd = find(freq == dddd);
    ampd = ahatf(ampd);
    ddddd = f2s(i);
    ampdd = find(freq == ddddd);
    ampdd = ahatf(ampdd);
    plot([dddd,ddddd],[ampd, ampdd], 'k')
    hold on
end
hold on
for i = 1:length(amps)
    plot([freqs(i),freqs(i)],[amps(i), amps(i) - proms(i)], 'k')
    hold on
end
hold on
xlim([0.1 5])
ylim([0.1 100])
xticks([.1 1 5])
xticklabels({'0.1', '1', '5'})
yticks([1 10 100])
yticklabels({ '1','10', '100'})
xlabel('Frequency (Hz)')
ylabel('Amplification')
title('NS')
set(gca,'YScale', 'log', 'XScale', 'log','FontName', 'Times New Roman', 'FontSize', 14)
grid on
box on
hold on

% EW
codepath = 'C:\Users\mpontr01\Desktop\HVSR\Codes';
datapath = 'C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Processed_data2\CE32\';
datapath2 = 'C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Processed_data2\TP13\';
cd(datapath)
eventlist = dir;
eventlist = eventlist(3:length(eventlist));
cd(datapath2)
refeventlist = dir;
refeventlist = refeventlist(3:length(refeventlist));
counter = 0;
SR = [];
for i = 1 : length(eventlist)
    event = eventlist(i).name;
    filename = strcat(datapath,event);
    CE32 = load(filename);
    mag = CE32.data.processing.filtereddata.acceleration.EW.mag_resps.smooth;
    eventid = event(5:18);
    for ii = 1:length(refeventlist)
        refevent = refeventlist(ii).name;
        refid = refevent(5:18);
        if strcmp(eventid,refid) == 1
            counter = counter + 1;
            filename = strcat(datapath2,refevent);
            TP13 = load(filename);
            refmag = TP13.data.processing.filtereddata.acceleration.EW.mag_resps.smooth;
            if length(refmag) == 100000
                refmag = refmag(1:50000);
            end
            SR(counter,:) = mag./refmag;
        end
    end
end

freq = CE32.data.processing.filtereddata.freq_vec;
cd(codepath)
[ahatf, sigma, confinthigh, confintlow] =  wavav(SR);
lowbound = 0.1;
upbound = 5;
[~, lowbound] = min(abs(freq - lowbound));
[~, upbound] = min(abs(freq - upbound));
[matrix, matrix1, peakind,ahatf1,newfaxhz1, peakfreqs, peakamps, Areamat] = peakiden(ahatf', freq, lowbound, upbound);
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
amps = matrix(:,2);
freqs = matrix(:,1);
f1s = hpb1(:,2);
f2s = hpb1(:,3);
proms = matrix1;

figure

hold on
confidenceinterval=shadedplot(freq(lowbound:length(freq)), confinthigh(lowbound:length(confinthigh)), confintlow(lowbound:length(confintlow)),[.9,.9,.9],[1 1 1]);
hold on
plot(freq,ahatf)
for i = 1:length(f1s)
    dddd = f1s(i);
    ampd = find(freq == dddd);
    ampd = ahatf(ampd);
    ddddd = f2s(i);
    ampdd = find(freq == ddddd);
    ampdd = ahatf(ampdd);
    plot([dddd,ddddd],[ampd, ampdd], 'k')
    hold on
end
hold on
for i = 1:length(amps)
    plot([freqs(i),freqs(i)],[amps(i), amps(i) - proms(i)], 'k')
    hold on
end
hold on
xlim([0.1 5])
ylim([0.1 100])
xticks([.1 1 5])
xticklabels({'0.1', '1', '5'})
yticks([1 10 100])
yticklabels({ '1','10', '100'})
xlabel('Frequency (Hz)')
ylabel('Amplification')
title('EW')
set(gca,'YScale', 'log', 'XScale', 'log','FontName', 'Times New Roman', 'FontSize', 14)
grid on
box on
hold on


% complex
codepath = 'C:\Users\mpontr01\Desktop\HVSR\Codes';
datapath = 'C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Processed_data2\CE32\';
datapath2 = 'C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Processed_data2\TP13\';
cd(datapath)
eventlist = dir;
eventlist = eventlist(3:length(eventlist));
cd(datapath2)
refeventlist = dir;
refeventlist = refeventlist(3:length(refeventlist));
counter = 0;
SR = [];
for i = 1 : length(eventlist)
    event = eventlist(i).name;
    filename = strcat(datapath,event);
    CE32 = load(filename);
    mag = CE32.data.processing.filtereddata.acceleration.complex.mag_resps.smooth;
    eventid = event(5:18);
    for ii = 1:length(refeventlist)
        refevent = refeventlist(ii).name;
        refid = refevent(5:18);
        if strcmp(eventid,refid) == 1
            counter = counter + 1;
            filename = strcat(datapath2,refevent);
            TP13 = load(filename);
            refmag = TP13.data.processing.filtereddata.acceleration.complex.mag_resps.smooth;
            if length(refmag) == 100000
                refmag = refmag(1:50000);
            end
            SR(counter,:) = mag./refmag;
        end
    end
end

freq = CE32.data.processing.filtereddata.freq_vec;
cd(codepath)
[ahatf, sigma, confinthigh, confintlow] =  wavav(SR);
lowbound = 0.1;
upbound = 5;
[~, lowbound] = min(abs(freq - lowbound));
[~, upbound] = min(abs(freq - upbound));
[matrix, matrix1, peakind,ahatf1,newfaxhz1, peakfreqs, peakamps, Areamat] = peakiden(ahatf', freq, lowbound, upbound);
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
amps = matrix(:,2);
freqs = matrix(:,1);
f1s = hpb1(:,2);
f2s = hpb1(:,3);
proms = matrix1;

figure

hold on
confidenceinterval=shadedplot(freq(lowbound:length(freq)), confinthigh(lowbound:length(confinthigh)), confintlow(lowbound:length(confintlow)),[.9,.9,.9],[1 1 1]);
hold on
plot(freq,ahatf)
for i = 1:length(f1s)
    dddd = f1s(i);
    ampd = find(freq == dddd);
    ampd = ahatf(ampd);
    ddddd = f2s(i);
    ampdd = find(freq == ddddd);
    ampdd = ahatf(ampdd);
    plot([dddd,ddddd],[ampd, ampdd], 'k')
    hold on
end
hold on
for i = 1:length(amps)
    plot([freqs(i),freqs(i)],[amps(i), amps(i) - proms(i)], 'k')
    hold on
end
hold on
xlim([0.1 5])
ylim([0.1 100])
xticks([.1 1 5])
xticklabels({'0.1', '1', '5'})
yticks([1 10 100])
yticklabels({ '1','10', '100'})
xlabel('Frequency (Hz)')
ylabel('Amplification')
title('Complex')
set(gca,'YScale', 'log', 'XScale', 'log','FontName', 'Times New Roman', 'FontSize', 14)
grid on
box on
hold on
