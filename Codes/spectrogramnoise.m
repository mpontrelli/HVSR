%% spectrogram
% based off Brian Traceys code shown in class 10/31/2019

%fnm = 'ghost01.wav';  cax = [-120 -30];
%close all
close all
clear all

outpath = 'C:\Users\mpontr01\Desktop\boston_site_response\Reports\Noise_comparison_908DP_Tufts_7a_11_1_2019';
%1a
% Vfname = 'C:\Users\mpontr01\Desktop\boston_site_response\field_deployments\Tufts_Campus\Stations\Line_A\1a\Data\2019303120000006_1a____1_1.sac';
% NSfname = 'C:\Users\mpontr01\Desktop\boston_site_response\field_deployments\Tufts_Campus\Stations\Line_A\1a\Data\2019303120000006_1a____1_2.sac';
% EWfname = 'C:\Users\mpontr01\Desktop\boston_site_response\field_deployments\Tufts_Campus\Stations\Line_A\1a\Data\2019303120000006_1a____1_3.sac';

% 2a
% Vfname = 'C:\Users\mpontr01\Desktop\boston_site_response\field_deployments\Tufts_Campus\Stations\Line_A\2a\Data\2019303\9B27\1\2019303124614006_2a____1_1.sac';
% NSfname = 'C:\Users\mpontr01\Desktop\boston_site_response\field_deployments\Tufts_Campus\Stations\Line_A\2a\Data\2019303\9B27\1\2019303124614006_2a____1_2.sac';
% EWfname = 'C:\Users\mpontr01\Desktop\boston_site_response\field_deployments\Tufts_Campus\Stations\Line_A\2a\Data\2019303\9B27\1\2019303124614006_2a____1_3.sac';

%Danahey park (Yilar's data)
% Vfname = 'B:\Erkan Yilar\Ambient Noise Data\Files\Ambient Noise Data Set\08.15.2014\Dan\Danehy Park Cambridge\Ground\Original Files\2014227134758005_T4260_1_1.sac';
% NSfname =  'B:\Erkan Yilar\Ambient Noise Data\Files\Ambient Noise Data Set\08.15.2014\Dan\Danehy Park Cambridge\Ground\Original Files\2014227134758005_T4260_1_2.sac';
% EWfname =  'B:\Erkan Yilar\Ambient Noise Data\Files\Ambient Noise Data Set\08.15.2014\Dan\Danehy Park Cambridge\Ground\Original Files\2014227134758005_T4260_1_3.sac';


%Danahey park (Our data)
% Vfname = 'C:\Users\mpontr01\Box\People\Marshall_and_Ashkan\Inversion_project\Data\908DP\Microtremor\Data\2000004171152005_908DP_1_1.sac';
% NSfname =  'C:\Users\mpontr01\Box\People\Marshall_and_Ashkan\Inversion_project\Data\908DP\Microtremor\Data\2000004171152005_908DP_1_2.sac';
% EWfname =  'C:\Users\mpontr01\Box\People\Marshall_and_Ashkan\Inversion_project\Data\908DP\Microtremor\Data\2000004171152005_908DP_1_3.sac';

% 7a
Vfname = 'C:\Users\mpontr01\Desktop\boston_site_response\field_deployments\Tufts_Campus\Stations\Line_A\7a\Data\2019292111241005_09B27_1_1.sac';
NSfname = 'C:\Users\mpontr01\Desktop\boston_site_response\field_deployments\Tufts_Campus\Stations\Line_A\7a\Data\2019292111241005_09B27_1_2.sac';
EWfname = 'C:\Users\mpontr01\Desktop\boston_site_response\field_deployments\Tufts_Campus\Stations\Line_A\7a\Data\2019292111241005_09B27_1_3.sac';

[xV] = ReadSacBinaryFile(Vfname); %vertical
[xNS] = ReadSacBinaryFile(NSfname); %North-south
[xEW] = ReadSacBinaryFile(EWfname); %East-West
% CE32
%load('C:\Users\mpontr01\Box\People\Marshall and Jeremy\HVSR_info\CE32.mat');

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
[xH] =  complex_time(xNS, xEW);
soundsc(xV/2, 4000)
audiowrite('C:\Users\mpontr01\Desktop\7a.wav',xV,4000);

fs = 100;

%% spectrogram inputs
nfft = 2048/4;
noverlap = round(nfft*0.01);
win = boxcar(nfft);
%win = blackman(nfft);
%win = hanning(nfft);

%% Vertical
sig = xV;
vert = figure;
subplot(5,1,1:4)
spectrogram(sig,win,noverlap,nfft,fs,'yaxis');
colorbar('horiz')
%caxis(cax)
title('Vertical')

subplot(5,1,5)
t = (1:length(sig))./fs;
plot(t,sig)
xlim([0 t(end)])  % try to get time axes of two plots to match...
%saveas(vert, strcat(outpath, '\', '7avertspectro.jpg'));
%% North south
sig = xNS;
NS = figure;
subplot(5,1,1:4)
spectrogram(sig,win,noverlap,nfft,fs,'yaxis');
colorbar('horiz')
%caxis(cax)
title('North south')

subplot(5,1,5)
t = (1:length(sig))./fs;
plot(t,sig)
xlim([0 t(end)])  % try to get time axes of two plots to match...
%saveas(NS, strcat(outpath, '\', '7aNSspectro.jpg'));
%% East west
sig = xEW;
EW = figure;
subplot(5,1,1:4)
spectrogram(sig,win,noverlap,nfft,fs,'yaxis');
colorbar('horiz')
%caxis(cax)
title('East west')

subplot(5,1,5)
t = (1:length(sig))./fs;
plot(t,sig)
xlim([0 t(end)])  % try to get time axes of two plots to match...

%saveas(EW, strcat(outpath, '\', '7aEWspectro.jpg'));

%% complex
sig = xH;
xomp = figure;
subplot(5,1,1:4)
spectrogram(sig,win,noverlap,nfft,fs,'yaxis');
colorbar('horiz')
%caxis(cax)
title('East west')

subplot(5,1,5)
t = (1:length(sig))./fs;
plot(t,sig)
xlim([0 t(end)])  % try to get time axes of two plots to match...