%% Surface_wave_simulate

% Reads and plots the text file output from the RAS-24 system to be used
% for surface wave analysis

%% Author: Marshall Pontrelli
% Date: 10/21/2020
%close all
clear all

%% begin
%Create a sin wave of many frequencies
%inputs
A = 1; %Amplitude
lam = 0.9; % wavelength
fs = 100;
f = 5;
t_vec = 0:1/fs:1;
x_vec = 10:2:56;
time = zeros(1,1000);
%mat = zeros(length(t_vec),length(x_vec));
%xxx = zeros(length(t_vec));
figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

for i = 1:length(x_vec)
    x = x_vec(i);
    for j = 1:length(t_vec)
        t = t_vec(j);
        xxx(j) = A*cos(2*pi*f*t-(2*pi/lam)*x);
        
    end
    mat(:,i) = xxx;
    % Extract positive and negative part
     yp = (xxx + abs(xxx))/2;
     yn = (xxx - abs(xxx))/2;
    subplot(1,24,i)
    
    hold on
    
    sta
    hold on 
    if i ==1
        area(t_vec,yp,'Facecolor','r')
    else
        area(t_vec,yp,'Facecolor','k')
    end
    view([90, 90])
    hold off
    disp(i)
    
end

%% plot horizontals
%figure
%plot(x_vec,(mat(20,:)))
% for i = 1:length(t_vec)
%     t = t_vec(i);
%     for j = 1:length(x_vec)
%         x = x_vec(j);
%         xxx(j) = A*cos(2*pi*f*t-(2*pi/lam)*x);
%         xxxx(j) = A*cos(2*pi*f2*t-(2*pi/lam2)*x);
%     end
%     new_x = xxx;% + xxxx;
% %     plot(x_vec,new_x)
% %     drawnow
%     mat(i,:) = new_x;
% end
% 
% 
% %% Now plot
% % figure
% % plot(time,xx)
% 
% 
% %% Now plot
% figure
% plot(x_vec,mat(1,:))
% 
% figure
% plot(t_vec,mat(:,1))
% 
% %% Now fft2
% %Y = fft2(mat);
% figure
% imagesc(mat)
% %% Now compute phase response
% 
% % theta = 0.9; %phase angle (Hz)
% % fs = 100; %sampling frequency in Hz
% % lam = 1; % wavelength (meters)
% % x = 1; % position (m)
% % t = 0:1/fs:10; %time vector between 0 and 100 seconds sampled at fs
% % x = A*cos(2*pi*f*t-(2*pi/lam)*x);
% % subplot(1,3,1)
% % Fig1 = plot(t,x);
% % %view([90, 90])
% % set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% % 
% % %% Decomposing signals using the fft
%% space
% X2 = fft(mat(:,1)); %taking the fft
% N = length(X2); %we need this for our frequency axis
% X2_mag = abs(X2)/(N); %magnitude response
% X2_phase = angle(X2)/(2*pi); %phase response
% bins = 0 : N-1; %samples in NS component
% freq = bins*fs/N; %frequency axis NS (Hz)
% figure
% subplot(2,1,1)
% plot(freq,X2_mag)
% ylim([0 20])
% title('Magnitude response')
% xlabel('Frequency (Hz)')
% ylabel('Amplitude (samples)')
% grid on
% subplot(2,1,2)
% Fig3 = plot(freq, X2_phase);
% title('Phase response')
% xlabel('Frequency (Hz)')
% ylabel('Phase (2\pi radians)')
% grid on
% 
% %% time
% X2 = fft(mat(:,1)); %taking the fft
% N = length(X2); %we need this for our frequency axis
% X2_mag = abs(X2)/(N); %magnitude response
% X2_phase_2 = angle(X2)/(2*pi); %phase response
% bins = 0 : N-1; %samples in NS component
% freq = bins*fs/N; %frequency axis NS (Hz)
% figure
% subplot(2,1,1)
% plot(freq,X2_mag)
% %ylim([0 20])
% title('Magnitude response')
% xlabel('Frequency (Hz)')
% ylabel('Amplitude (samples)')
% grid on
% subplot(2,1,2)
% Fig3 = plot(freq, X2_phase_2);
% title('Phase response')
% xlabel('Frequency (Hz)')
% ylabel('Phase (2\pi radians)')
% grid on
% 
% %% 
% z = X2_phase_2.* X2_phase;
% 
% figure
% plot(z)