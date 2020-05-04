%% read DMC data
close all
clear all
%% necessary inputs
LowCorner = 0.1;
HighCorner = 15;
Npoles = 4;
Filterplot = 'no';

%% station 1
Station = 'DC23';
Network = 'ZN';
Component = 'EHZ';
time1 = '2015-06-07 06:07:00';
time2 = '2015-06-08 06:07:00';

[sampletimes1,trace1,stat1, fs1, sensitivity1, sensunits1] = getDMCData(time1,time2,Station,Network,Component);
[stat1] = Butter2(stat1, fs1, 'LowCorner', LowCorner, 'HighCorner', HighCorner, 'Npoles', Npoles , 'Filterplot', Filterplot);
stat1 = stat1/sensitivity1;
stat1 = stat1(1:length(stat1)-1);
%% station 2
Station = 'DC23';
Network = 'ZN';
Component = 'EHE';
time1 = '2015-06-07 06:07:00';
time2 = '2015-06-08 06:07:00';

[sampletimes2,trace2,stat2, fs2, sensitivity2, sensunits2] = getDMCData(time1,time2,Station,Network,Component);
[stat2] = Butter2(stat2, fs2, 'LowCorner', LowCorner, 'HighCorner', HighCorner, 'Npoles', Npoles , 'Filterplot', Filterplot);
stat2 = stat2/sensitivity2;
%stat2 = stat2(1:length(stat2)-1);
%% plot time series
% find the maximum value to make bounds for plotting
a = max(abs(stat1));
b = max(abs(stat2));
c = [a,b];
d = max(c);

%% create a time vector and plot time series
time1 = ((1:length(stat1))/fs1)/3600;
time2 = ((1:length(stat2))/fs2)/3600;

figure;
    
% station 1
subplot(2,1,1)
plot(time1,stat1)
title('Station 1')
ylabel(sensunits1)
xlabel('Time (Hours)')
xlim([0 length(stat1)/(fs1*3600)])
ylim([-d d])
grid on 
box on
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);

% station 2
subplot(2,1,2)
plot(time2,stat2)
title('Station 2')
ylabel(sensunits2)
xlabel('Time (Hours)')
xlim([0 length(stat2)/(fs2*3600)])
ylim([-d d])
grid on 
box on
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

%% integrate for displacement
alpha = 0.25; %coefficients for Newmark-beta method
beta = 0.5;

% station 1
dt = 1/fs1;
d1 = 0;
for i = 1:length(stat1) - 1
    d1(i+1) = d1(i) + dt*((1-beta)*stat1(i) + beta*stat1(i+1));
end
d1 = d1 - mean(d1);
[d1] = Butter2(d1, fs1);
d1 = d1*1000000; % to microns

% station 2
dt2 = 1/fs2;
d2 = 0;
for i = 1:length(stat2) - 1
    d2(i+1) = d2(i) + dt2*((1-beta)*stat2(i) + beta*stat2(i+1));
end
d2 = d2 - mean(d2);
[d2] = Butter2(d2, fs2);
d2 = d2*1000000; % to microns

%% plot displacement
% find the maximum value to make bounds for plotting
a = max(abs(d1));
b = max(abs(d2));
c = [a,b];
d = max(c);
figure;
    
% station 1
subplot(2,1,1)
plot(time1,d1)
title('Station 1')
ylabel('Microns')
xlabel('Time (Hours)')
xlim([0 length(stat1)/(fs1*3600)])
ylim([-d d])
grid on 
box on
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);

% station 2
subplot(2,1,2)
plot(time2,d2)
title('Station 2')
ylabel('Microns')
xlabel('Time (Hours)')
xlim([0 length(stat2)/(fs2*3600)])
ylim([-d d])
grid on 
box on
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

%% spectrogram inputs
nfft = 2048/4;
noverlap = round(nfft*0.50);
win = hanning(nfft);

figure
subplot(1,2,1)
spectrogram(stat1,win,noverlap,nfft,fs1,'yaxis');
title('station 1')
set(gca,'YScale', 'log','FontName', 'Times New Roman', 'FontSize', 14)

subplot(1,2,2)
spectrogram(stat2,win,noverlap,nfft,fs2,'yaxis');
title('station 2')
set(gca,'YScale', 'log','FontName', 'Times New Roman', 'FontSize', 14)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
%% Now do magnitude responses and compare
% inputs
numwin = 1000;
windowlen = 60;
windis = 20;

% station 1
sampnum = windowlen*fs1; 
windisnum = windis*fs1;

k = [1,fs1];
for iii = 1:numwin
    stat1matrix(iii,:) = stat1((k(1)):(k(2)*windowlen+k(1))-1);
    k(1) = k(1)+ sampnum + windisnum;
end

% station 2
sampnum = windowlen*fs2; 
windisnum = windis*fs2;

k = [1,fs2];
for iii = 1:numwin
    stat2matrix(iii,:) = stat2((k(1)):(k(2)*windowlen+k(1))-1);
    k(1) = k(1)+ sampnum + windisnum;
end

%% window the data
win = hann(sampnum)';
for i = 1:numwin
    stat1matrix(i,:) = stat1matrix(i,:).*win;
    stat2matrix(i,:) = stat2matrix(i,:).*win;
end

%% Compute unfiltered magnitude responses
for iii = 1:numwin
    stat1matrix(iii,:) = 4*abs(fft(stat1matrix(iii,:)))/sampnum;
    stat2matrix(iii,:) = 4*abs(fft(stat2matrix(iii,:)))/sampnum;
end

%Computing the frequency -axis
N = sampnum;
fax_binsN = (0 : N-1); %samples in NS component
fax_HzN1 = fax_binsN*fs1/N; %frequency axis NS (Hz)
N_2 = ceil(N/2); %half magnitude spectrum
fax_HzN = fax_HzN1(1 : N_2);
for iii = 1:numwin
    mat1 = stat1matrix(iii,:);
    mat2 = stat2matrix(iii, :);
    stat1matrix2(iii,:) = mat1(1 : N_2);
    stat2matrix2(iii,:) = mat2(1 : N_2);
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
subplot(1,2,1)
plot(fax_HzN, stat1matrix2, 'Color',  'k' , 'Linewidth', .5);
title('Station 1', 'FontSize', 20)
xlabel('Frequency (Hz)','FontSize', 18)
ylabel('Amplification','FontSize', 18)
set(gca,'YScale', 'log', 'XScale', 'log','FontName', 'Times New Roman', 'FontSize', 14)
xlim([0.1 10])
ylim([10e-16 10e-5])
xticks([0.1 1 10])
xticklabels({'0.1','1','10'})
% yticks([0.1 1 10 100 1000])
% yticklabels({'0.1','1','10', '100', '1000'})
grid on 
box on

subplot(1,2,2)
plot(fax_HzN, stat2matrix2, 'Color',  'k' , 'Linewidth', .5);
title('Station 2', 'FontSize', 20)
xlabel('Frequency (Hz)','FontSize', 18)
ylabel('Amplification','FontSize', 18)
set(gca,'YScale', 'log', 'Xscale', 'log','FontName', 'Times New Roman', 'FontSize', 14)
xlim([0.1 10])
ylim([10e-16 10e-5])
xticks([0.1 1 10])
xticklabels({'0.1','1','10'})
% yticks([0.1 1 10 100 1000])
% yticklabels({'0.1','1','10', '100', '1000'})
grid on 
box on

%makes figure full screen
set(gcf, 'Units', 'Normalized', 'OuterPosition', [1, 0.5, 1, 0.5]);


%% Average the un-smoothed magnitude responses
[ahatf1, sigma1, confinthigh1, confintlow1] =  wavav(stat1matrix2);
[ahatf2, sigma2, confinthigh2, confintlow2] =  wavav(stat2matrix2);

%% now plot
averageunfiltered = figure;
subplot(1,2,1)
hold on
confidenceinterval=shadedplot(fax_HzN, confinthigh1, confintlow1,[.9,.9,.9],'k');
hold on
plot(fax_HzN, ahatf1, 'Color', [0 0.30196 0.6588] , 'Linewidth', 1.5);
title('Station 1')
xlabel('Frequency (Hz)')
ylabel('Amplification')
set(gca,'YScale', 'log', 'Xscale', 'log','FontName', 'Times New Roman', 'FontSize', 14)
xlim([0.1 10])
ylim([10e-16 10e-5])
xticks([0.1 1 10])
xticklabels({'0.1','1','10'})
% yticks([0.1 1 10 100 1000])
% yticklabels({'0.1','1','10', '100', '1000'})
grid on 
box on

%makes figure full screen
set(gcf, 'Units', 'Normalized', 'OuterPosition', [1, 0.5, 1, 0.5]);


subplot(1,2,2)
hold on
confidenceinterval=shadedplot(fax_HzN, confinthigh2, confintlow2,[.9,.9,.9],'k');
hold on
plot(fax_HzN, ahatf2, 'Color', [0 0.30196 0.6588] , 'Linewidth', 1.5);
title('Station 2')
xlabel('Frequency (Hz)')
ylabel('Amplification')
set(gca,'YScale', 'log', 'Xscale', 'log','FontName', 'Times New Roman', 'FontSize', 14)
xlim([0.1 10])
ylim([10e-16 10e-5])
xticks([0.1 1 10])
xticklabels({'0.1','1','10'})
% yticks([0.1 1 10 100 1000])
% yticklabels({'0.1','1','10', '100', '1000'})
grid on 
box on

%makes figure full screen
set(gcf, 'Units', 'Normalized', 'OuterPosition', [1, 0.5, 1, 0.5]);

%% now take a ratio of the two responses with station 2 as numerator
rat = ahatf2 ./ahatf1;
figure
plot(fax_HzN,rat, 'LineWidth', 2)
title('Spectral ratio, stat 2 numerator')
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

%% now take a ratio of the two responses with station 1 as numerator
rat = ahatf1 ./ahatf2;
figure
plot(fax_HzN,rat, 'LineWidth', 2)
title('Spectral ratio, stat 1 numerator')
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

%% now compare the sigma vectors
sigma1 = smooth(sigma1);
sigma2 = smooth(sigma2);
figure
sig1 = plot(fax_HzN,sigma1, 'LineWidth', 2);
hold on
sig2 = plot(fax_HzN,sigma2, 'LineWidth', 2);
title('Sigma comparison')
xlabel('Frequency (Hz)')
ylabel('\sigma')
set(gca,'YScale', 'log', 'Xscale', 'log','FontName', 'Times New Roman', 'FontSize', 14)
xlim([0.1 10])
%ylim([0.1 100])
xticks([0.1 1 10])
xticklabels({'0.1','1','10'})
% yticks([0.1 1 10 100])
% yticklabels({'0.1','1','10', '100'})
legend([sig1, sig2], 'Station 1','Station 2')
grid on 
box on

%makes figure full screen
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

%% now take cross correlation
[c,lags] = xcorr(stat1, stat2,'normalized');

%%
figure
plot(lags/(3600*fs1),c)