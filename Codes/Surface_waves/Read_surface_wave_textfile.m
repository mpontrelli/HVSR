%% Read_surface_wave_textfile

% Reads and plots the text file output from the RAS-24 system to be used
% for surface wave analysis

%% Author: Marshall Pontrelli
% Date: 9/11/2020
close all
clear all

%% inputs
fs = 500

%% import the data
filename = 'C:\Users\mpontr01\Box\Projects\Surface_wave\9_11\data1.txt';
A = importdata(filename);

%% figure out how big A is
a = size(A);
a = a(2);

%% now plot
figure
for i = 1:a
    sta = A(:,i);
    sta = sta(1:100);
    subplot(1,a,i)
    plot(sta)
    set(gca,'xtick',[], 'ytick', [])
    view([90, 90])
    ti = num2str(i);
    title(ti)
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

%% now make a dispersion curve
d1 = 5.94;
d2 = 7.47;

%% pick right stations
sta1 = A(:, 1);
sta2 = A(:, 6);

%% window
win = hann(length(sta1))';
sta1 = sta1.*win';
sta2 = sta2.*win';

%% Take ffts
STA1 = angle(fft(sta1)/length(sta1));
STA2 = angle(fft(sta2)/length(sta2));

%% compute frequency vector
N = length(sta1);
fax_binsN = (0 : N-1); %samples in NS component
fax_HzN1 = fax_binsN*fs/N; %frequency axis NS (Hz)
N_2 = ceil(N/2); %half magnitude spectrum
fax_HzN = fax_HzN1(1 : N_2);
STA1_2 = STA1(1: N_2);
STA2_2 = STA2(1: N_2);

%% now plot
figure
plot(fax_HzN, STA1_2)
hold on
plot(fax_HzN, STA2_2)

%% now do dispersion curve (Kramer 1996 page 204 eqs 6.27 - 6.29
phi_f = STA2_2 - STA1_2;
del_t_f = phi_f'./ (2*pi*fax_HzN);
del_d = d2 - d1;
vR_f = (1/del_d)*del_t_f;
lamR_f = vR_f./fax_HzN;

%% now plot
figure
plot(vR_f, lamR_f)

