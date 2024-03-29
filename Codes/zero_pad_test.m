%% Test zero padding on HVSR
% Conclusion, zero pad as much as you want just make sure you sero pad on
% both sides and you divede byt the number of samples of the origional
% signal, not the new sample length after zero padding.
%close all
clear all

filename1 = 'C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Data\CE32\CE3220170919181440';
filename2 = 'C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Data\CE32\CE3220170908044918';
fs = 100;
[xNS,xV,xEW] = readfile1(filename1);
[xNS2,xV2,xEW2] = readfile1(filename2);

a = length(xV);
b = length(xV2);
c = abs(b-a);

%% zero pad the smaller vector
xNS1 = vertcat(zeros(1,1000)', xNS, zeros(1,1000)');
xV1 = vertcat(zeros(1,1000)', xV, zeros(1,1000)');
xEW1 = vertcat(zeros(1,1000)', xEW, zeros(1,1000)');
b = length(xEW1);

% %% window the data
win1 = hann(a);
win = hann(b);

xNS_win = xNS.*win1;
xEW_win = xEW.*win1;
xV_win = xV.*win1;

xNS1_win = xNS1.*win;
xEW1_win = xEW1.*win;
xV1_win = xV1.*win;

% xNS2_win = xNS2.*win;
% xEW2_win = xEW2.*win;
% xV2_win = xV2.*win;

%% Now do HVSRs
XNS = 4*abs(fft(xNS_win))/a;
XEW = 4*abs(fft(xEW_win))/a;
V = 4*abs(fft(xV_win))/a;
XNS1 = 4*abs(fft(xNS1_win))/a;
XEW1 = 4*abs(fft(xEW1_win))/a;
V1 = 4*abs(fft(xV1_win))/a;
% XNS2 = 4*abs(fft(xNS2_win))/b;
% XEW2 = 4*abs(fft(xEW2_win))/b;
% V2 = 4*abs(fft(xV2_win))/b;

%% freq vec for zero padded signals
N = b;
fax_binsN = (0 : N-1); %samples in NS component
fax_HzN1 = fax_binsN*fs/N; %frequency axis NS (Hz)
N_2 = ceil(N/2); %half magnitude spectrum
fax_HzN = fax_HzN1(1 : N_2);

%% freq vec for zero padded signals
N = a;
fax_binsN = (0 : N-1); %samples in NS component
fax_HzN1 = fax_binsN*fs/N; %frequency axis NS (Hz)
N_3 = ceil(N/2); %half magnitude spectrum
fax_HzN2 = fax_HzN1(1 : N_3);

%%
XNS_fin = XNS(1 : N_3);
V_fin = V(1 : N_3);

XNS1_fin = XNS1(1 : N_2);
V1_fin = V1(1 : N_2);

% XNS2_fin = XNS2(1 : N_2);
% V2_fin = V2(1 : N_2);

HV1 = XNS1_fin./V1_fin;
% HV2 = XNS2_fin./V2_fin;
HV3 = XNS_fin./V_fin;

%% This plots the HVSR from the orogional record and the zero-padded HVSR
% In the noisier parts of the signal, the comparison is off a little bit,
% but for the most part, I think the two match up well
figure
plot(fax_HzN, HV1)
hold on
% plot(fax_HzN, HV2)
% hold on
plot(fax_HzN2, HV3)

%% Here is the same with magnitude response
figure
plot(fax_HzN, XNS1_fin)
hold on
% plot(fax_HzN, HV2)
% hold on
plot(fax_HzN2, XNS_fin)
