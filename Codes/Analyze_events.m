%% Analyze events

% Author: Marshall Pontrelli
% Date: 5/15/2020

close all
clear all


%%

sitelist = {'CS78'};%, 'IM40', 'TE07','UI21','FJ74','TP13'};

%% CS78
filename = 'C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Processed_data2\CS78\CS7820170919181440';

CS78 = load(filename);
CS78ID = CS78.data.meta.event.ID;
CS78az = CS78.data.meta.event.azimuth;
CS78dist = CS78.data.meta.event.epi_dist;
CS78magrespNS = CS78.data.processing.filtereddata.acceleration.NS.mag_resps.smooth;
CS78wave = CS78.data.processing.filtereddata.acceleration.NS.waveform_orig;
CS78time = CS78.data.processing.filtereddata.time_orig;
% The next 3 lines normalize the waveform by the sqrt of its total energy
waveform1sq = CS78wave.^2;
waveform1norm = sqrt(sum(waveform1sq));
CS78wave = CS78wave/waveform1norm;

%% IM40
filename = 'C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Processed_data2\IM40\IM4020170919181440';

IM40 = load(filename);
IM40ID = IM40.data.meta.event.ID;
IM40az = IM40.data.meta.event.azimuth;
IM40dist = IM40.data.meta.event.epi_dist;
IM40magrespNS = IM40.data.processing.filtereddata.acceleration.NS.mag_resps.smooth;
IM40wave = IM40.data.processing.filtereddata.acceleration.NS.waveform_orig;
% The next 3 lines normalize the waveform by the sqrt of its total energy
waveform1sq = IM40wave.^2;
waveform1norm = sqrt(sum(waveform1sq));
IM40wave = IM40wave/waveform1norm;
%% TE07
% filename = 'C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Processed_data2\TE07\TE0720170919181440';
% 
% TE07 = load(filename);
% TE07ID = TE07.data.meta.event.ID;
% TE07az = TE07.data.meta.event.azimuth;
% TE07dist = TE07.data.meta.event.epi_dist;
% TE07magrespNS = TE07.data.processing.filtereddata.acceleration.NS.mag_resps.smooth;

%% UI21
filename = 'C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Processed_data2\UI21\UI2120170919181440';

UI21 = load(filename);
UI21ID = UI21.data.meta.event.ID;
UI21az = UI21.data.meta.event.azimuth;
UI21dist = UI21.data.meta.event.epi_dist;
UI21magrespNS = UI21.data.processing.filtereddata.acceleration.NS.mag_resps.smooth;
UI21wave = UI21.data.processing.filtereddata.acceleration.NS.waveform_orig;
UI21time = UI21.data.processing.filtereddata.time_orig;
% The next 3 lines normalize the waveform by the sqrt of its total energy
waveform1sq = UI21wave.^2;
waveform1norm = sqrt(sum(waveform1sq));
UI21wave = UI21wave/waveform1norm;

%% TP13
filename = 'C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Processed_data2\TP13\TP1320170919181440';

TP13 = load(filename);
TP13ID = TP13.data.meta.event.ID;
TP13az = TP13.data.meta.event.azimuth;
TP13dist = TP13.data.meta.event.epi_dist;
TP13magrespNS = TP13.data.processing.filtereddata.acceleration.NS.mag_resps.smooth;
TP13wave = TP13.data.processing.filtereddata.acceleration.NS.waveform_orig;
TP13time = TP13.data.processing.filtereddata.time_orig;
% The next 3 lines normalize the waveform by the sqrt of its total energy
waveform1sq = TP13wave.^2;
waveform1norm = sqrt(sum(waveform1sq));
TP13wave = TP13wave/waveform1norm;

%% FJ74
filename = 'C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Processed_data2\FJ74\FJ7420170919181440';

FJ74 = load(filename);
FJ74ID = FJ74.data.meta.event.ID;
FJ74az = FJ74.data.meta.event.azimuth;
FJ74dist = FJ74.data.meta.event.epi_dist;
FJ74magrespNS = FJ74.data.processing.filtereddata.acceleration.NS.mag_resps.smooth;
FJ74wave = FJ74.data.processing.filtereddata.acceleration.NS.waveform_orig;
FJ74time = FJ74.data.processing.filtereddata.time_orig;
% The next 3 lines normalize the waveform by the sqrt of its total energy
waveform1sq = FJ74wave.^2;
waveform1norm = sqrt(sum(waveform1sq));
FJ74wave = FJ74wave/waveform1norm;

%% Now load frequency vector
freq = FJ74.data.processing.filtereddata.freq_vec;

%% Now plot

fig = figure;
xlim([0.1 10])
%ylim([0.1 100])
xticks([.1 1 10])
xticklabels({'0.1', '1', '10'})
%yticks([1 10 100])
%yticklabels({ '1','10', '100'})
xlabel('Frequency (Hz)')
ylabel('Amplification')
%title(statname)
set(gca,'YScale', 'log','XScale','log', 'FontName', 'Times New Roman', 'FontSize', 18)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
grid on
box on
hold on
CS78 = plot(freq, CS78magrespNS,'Color', 'm');
hold on
IM40 = plot(freq, IM40magrespNS(1:50000),'Color', 'r');
hold on
UI21 = plot(freq, UI21magrespNS,'Color', 'k');
hold on
TP13 = plot(freq, TP13magrespNS,'Color', 'b');
hold on
FJ74 = plot(freq, FJ74magrespNS,'Color', 'c');

legend([CS78, FJ74, IM40, TP13, UI21],'CS78', 'FJ74', 'IM40', 'TP13', 'UI21')


%% Now correlate
[c12,lags12]=xcorr(UI21wave,TP13wave);
[maxc12,indc12]=max(c12);
%lagstime1=lags12(indc12)*delta;
figure
subplot(3,1,1)
plot(TP13time,TP13wave)
ylim([-0.1, 0.1])
subplot(3,1,2)
plot(UI21time,UI21wave)
ylim([-0.1, 0.1])
subplot(3,1,3)
plot(lags12,c12)