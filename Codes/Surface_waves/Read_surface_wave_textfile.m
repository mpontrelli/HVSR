%% Read_surface_wave_textfile

% Reads and plots the text file output from the RAS-24 system to be used
% for surface wave analysis

%% Author: Marshall Pontrelli
% Date: 9/11/2020
close all
clear all

%% inputs
fs = 500;
datapath = 'C:\Users\mpontr01\Box\Projects\Surface_wave\9_14\';
%% import the data
% ACCESSING THE DATA
% go into the data folder and get a list of stations
path = pwd;
cd(datapath)
tracelist = dir;
tracelist = tracelist(3:length(tracelist));
% start the for loop that goes through all the station folders
B = zeros(501, 24);
for eee = 1:length(tracelist)
    trace = tracelist(eee); % read the folder info
    tracename = trace.name;
    filename = strcat(datapath,tracename);
    A = importdata(filename);
    B = A + B;
    q{eee} = importdata(filename);
end
B = B/length(tracelist);

%% figure out how big B is
a = size(B);
a = a(2);



%% now plot
figure
for i = 1:24
    sta = B(:,i);
    sta = sta(1:100);
    subplot(1,a,i)
    plot(sta)
    set(gca,'xtick',[], 'ytick', [])
    view([90, 90])
    ti = num2str(i);
    title(ti)
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

%% Now plot each trace
a = size(q);
a = a(2);

%% now plot all traces to be stacked
figure
for i = 1:a
    tr = q{i};
    for j = 1:24
        sta = tr(:,j);
        sta = sta(1:100);
        subplot(1,24,j)
        plot(sta)
        set(gca,'xtick',[], 'ytick', [])
        view([90, 90])
        ti = num2str(j);
        title(ti)
        hold on
    end
end
%% now make a dispersion curve
Sx = B(:,2);
N = length(Sx);
dref = 3; % station 2, 3 meters from the trigger

Sx = fft(Sx)/ N;
Sx = conj(Sx);

fax_binsN = (0 : N-1); %samples in NS component
omega = fax_binsN*2*pi/N; %frequency axis NS (radians)
Vrf_plot = [];
lam_theta_plot = [];
for i = 1:10
    ds = dref + i;
    Sy = B(:,i+2);
    Sy = fft(Sy)/N;
    Sy = real(Sy);
    Gxy = Sx.*Sy;
    thetaxy = atan(imag(Gxy)./real(Gxy));
    tf = thetaxy./omega;
    lam_theta = thetaxy.^-1*((ds - dref));
    Vrf = omega'.*lam_theta;
    ind = find(lam_theta /3 < ds - dref & lam_theta > 2*lam_theta);
    Vrf = Vrf(ind);
    lam_theta = lam_theta(ind);
    Vrf_plot = vertcat(Vrf_plot,Vrf);
    lam_theta_plot = vertcat(lam_theta_plot,lam_theta);
end

%% Now plot
figure
plot(omega(ind).^-1,abs(Vrf), '*')

% 
% %% pick right stations
% sta1 = A(:, 1);
% sta2 = A(:, 6);
% 
% %% window
% win = hann(length(sta1))';
% sta1 = sta1.*win';
% sta2 = sta2.*win';
% 
% %% Take ffts
% STA1 = angle(fft(sta1)/length(sta1));
% STA2 = angle(fft(sta2)/length(sta2));
% 
% %% compute frequency vector
% N = length(sta1);
% fax_binsN = (0 : N-1); %samples in NS component
% fax_HzN1 = fax_binsN*fs/N; %frequency axis NS (Hz)
% N_2 = ceil(N/2); %half magnitude spectrum
% fax_HzN = fax_HzN1(1 : N_2);
% STA1_2 = STA1(1: N_2);
% STA2_2 = STA2(1: N_2);
% 
% %% now plot
% figure
% plot(fax_HzN, STA1_2)
% hold on
% plot(fax_HzN, STA2_2)
% 
% %% now do dispersion curve (Kramer 1996 page 204 eqs 6.27 - 6.29
% phi_f = STA2_2 - STA1_2;
% del_t_f = phi_f'./ (2*pi*fax_HzN);
% del_d = d2 - d1;
% vR_f = (1/del_d)*del_t_f;
% lamR_f = vR_f./fax_HzN;
% 
% %% now plot
% figure
% plot(vR_f, lamR_f)

