close all
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Problem 1


%Inputs
alpha = 0.25; %coefficients for Newmark-beta method
beta = 0.5;
width = 0.1; % width of smoothing filter in Hz
bracketcut = 0.05*9.8; %This is the cutoff for the bolt bracketed duration in units of acceleration (not gs)
statname = 'CO56';
%output = 'C:\Users\mpontr01\Box Sync\Box Sync\tspring_2019\Earthquake Engineering\Assignments\Assignment 2\Figures';

%% Read the data and make time vector

[xNS,xV,xEW, fs] = readfile1(strcat('C:\Users\mpontr01\Box\2020_1_spring\SSA\Time domain\Puebla_data\', statname));
[xNS, xV, xEW] = Butter(xNS, xV, xEW, fs); %filter the data
% xNS = xNS - mean(xNS);
% [xNS] = Butter2(xNS);
n = length(xNS);
win = hann(n);
dt = 1/fs;
time = 0:dt:(n - 1)*dt;

%% find PGA
[PGANS,PGAV,PGAEW] = PGA(xNS,xV,xEW);
a = max([PGANS PGAV PGAEW]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot North-south component acceleration
%acceleration
figure
%sgtitle('North-South','FontName', 'Times New Roman', 'FontSize', 18, 'FontWeight', 'bold', 'color','k');
subplot(3,2,1)
plot(time, xNS);set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, 'XColor', 'k');
ylabel('(cm/s^2)','color','k'); 
grid on; box on; xlim([0 time(end)]); %ylim([(-a -5) (a +5) ]);

%% now do magnitude response
xNSw = xNS.*win;
XNS = fft(xNSw); %taking the fft
N = length(XNS); %we need this for our frequency axis
XNS_mag = 4*abs(XNS)/(N); %magnitude response
window = ceil((N/fs)*width);
XNS_mag = smooth(XNS_mag, window);
bins = 0 : N-1; %samples in NS component
freq = bins*fs/N; %frequency axis NS (Hz)
subplot(3,2,2)
freqvert = plot(freq, XNS_mag);set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, 'Yscale','log', 'Xscale', 'log', 'XColor', 'k');
grid on; box on;
xlim([0.03 10]); ylim([1e-4 10]);
xticks([0.1 1 10])
xticklabels({'0.1','1','10'})
yticks([0.001 0.01 0.1 1 10])
yticklabels({'0.001','0.01','0.1','1', '10'})

%% PGA = max(abs(FN));
%disp(strcat('PGA FN = ', num2str(PGA)));
%bracketed duration
%q = (FN > bracketcut) | (FN < -1* bracketcut);
%qq = find(q == 1);
%bdfirst = time(qq(1));
%bdlast = time(qq(length(qq)));
%boltbrackdur = bdlast-bdfirst;
%line([0 time(length(time))], [bracketcut bracketcut],'color','k')
%line([0 time(length(time))], [-bracketcut -bracketcut],'color','k')
%line([bdfirst bdfirst], [-PGA PGA],'color','r')
%line([bdlast bdlast], [-PGA PGA],'color','r')
%disp(strcat('Bracketed duration FN = ', num2str(boltbrackdur)));
%hold off
%% integrate for velocity (Newmark-Beta method)
v = 0;
for i = 1:length(xNS) - 1
    v(i+1) = v(i) + dt*((1-beta)*xNS(i) + beta*xNS(i+1));
end
v = v - mean(v);
[vNS] = Butter2(v);
subplot(3,2,3)
plot(time, vNS);set(gca, 'FontName', 'Times New Roman', 'FontSize', 14,'XColor', 'k');
ylabel('(cm/s)', 'color','k');
grid on; box on; xlim([0 time(end)]); 


%% now do magnitude response of velocity
vw = vNS';
vw = vw.*win;
V = fft(vw); %taking the fft
N = length(V); %we need this for our frequency axis
V_mag = 4*abs(V)/(N); %magnitude response
window = ceil((N/fs)*width);
V_mag = smooth(V_mag, window);
bins = 0 : N-1; %samples in NS component
freq = bins*fs/N; %frequency axis NS (Hz)
subplot(3,2,4)
freqvert = plot(freq, V_mag);set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, 'Yscale','log', 'Xscale', 'log','XColor', 'k');
grid on; box on;
xlim([0.03 10]); ylim([1e-4 10]);
xticks([0.1 1 10])
xticklabels({'0.1','1','10'})
yticks([0.001 0.01 0.1 1 10])
yticklabels({'0.001','0.01','0.1','1', '10'})
%PGV = max(abs(v));
%disp(strcat('PGV FN = ', num2str(PGV)));
%% Integrate for displacement 
dN = 0;
for i = 1:length(xEW) - 1
    dN(i+1) = dN(i) + (dt)*v(i) + (dt)^2 * ((0.5-alpha)*xNS(i) + alpha*xNS(i+1));
end
dN = dN - mean(dN);
[dNS] = Butter2(dN);
subplot(3,2,5)
fig1 = plot(time, dNS);set(gca, 'FontName', 'Times New Roman', 'FontSize', 14,'XColor', 'k');
xlabel('time (secs)','color','k'); ylabel('(cm)', 'color','k');
grid on; box on; xlim([0 time(end)]); 

%% now do magnitude response of displacement
dNS = dNS';
dNSw = dNS.*win;
DNS = fft(dNSw); %taking the fft
N = length(DNS); %we need this for our frequency axis
DNS_mag = 4*abs(DNS)/(N); %magnitude response
window = ceil((N/fs)*width);
DNS_mag = smooth(DNS_mag, window);
bins = 0 : N-1; %samples in NS component
freq = bins*fs/N; %frequency axis NS (Hz)
subplot(3,2,6)
freqvert = plot(freq, DNS_mag);set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, 'Yscale','log', 'Xscale', 'log', 'XColor','k');
xlabel('Frequency (Hz)', 'color','k'); 
grid on; box on;
xlim([0.03 10]); ylim([1e-4 10]);
xticks([0.1 1 10])
xticklabels({'0.1','1','10'})
yticks([0.001 0.01 0.1 1 10])
yticklabels({'0.001','0.01','0.1','1', '10'})
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

%% Plot East-west component acceleration
%acceleration
figure
%sgtitle('East-West','FontName', 'Times New Roman', 'FontSize', 18, 'FontWeight', 'bold','color','k');
subplot(3,2,1)
plot(time, xEW);set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, 'XColor','k');
ylabel('(cm/s^2)', 'color', 'k'); 
grid on; box on; xlim([0 time(end)]); %ylim([(-a -5) (a +5) ]);

%% now do magnitude response
xEWw = xEW.*win;
XEW = fft(xEWw); %taking the fft
N = length(XEW); %we need this for our frequency axis
XEW_mag = 4*abs(XEW)/(N); %magnitude response
window = ceil((N/fs)*width);
XEW_mag = smooth(XEW_mag, window);
bins = 0 : N-1; %samples in NS component
freq = bins*fs/N; %frequency axis NS (Hz)
subplot(3,2,2)
freqvert = plot(freq, XEW_mag);set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, 'Yscale','log', 'Xscale', 'log', 'XColor','k');
grid on; box on;
xlim([0.03 10]); ylim([1e-4 10]);
xticks([0.1 1 10])
xticklabels({'0.1','1','10'})
yticks([0.001 0.01 0.1 1 10])
yticklabels({'0.001','0.01','0.1','1', '10'})

%% integrate for velocity (Newmark-Beta method)
vEW = 0;
for i = 1:length(xNS) - 1
    vEW(i+1) = vEW(i) + dt*((1-beta)*xEW(i) + beta*xEW(i+1));
end
vEW = vEW - mean(vEW);
[vEW] = Butter2(vEW);
subplot(3,2,3)
plot(time, vEW);set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, 'XColor','k');
ylabel('(cm/s)', 'color', 'k');
grid on; box on; xlim([0 time(end)]); 


%% now do magnitude response of velocity
vEWw = vEW';
vEWw = vEWw.*win;
VEW = fft(vEWw); %taking the fft
N = length(VEW); %we need this for our frequency axis
VEW_mag = 4*abs(VEW)/(N); %magnitude response
window = ceil((N/fs)*width);
VEW_mag = smooth(VEW_mag, window);
bins = 0 : N-1; %samples in NS component
freq = bins*fs/N; %frequency axis NS (Hz)
subplot(3,2,4)
freqvert = plot(freq, VEW_mag);set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, 'Yscale','log', 'Xscale', 'log', 'XColor','k');
grid on; box on;
xlim([0.03 10]); ylim([1e-4 10]);
xticks([0.1 1 10])
xticklabels({'0.1','1','10'})
yticks([0.001 0.01 0.1 1 10])
yticklabels({'0.001','0.01','0.1','1', '10'})
%PGV = max(abs(v));
%disp(strcat('PGV FN = ', num2str(PGV)));
%% Integrate for displacement 
dEW = 0;
for i = 1:length(xEW) - 1
    dEW(i+1) = dEW(i) + (dt)*vEW(i) + (dt)^2 * ((0.5-alpha)*xEW(i) + alpha*xEW(i+1));
end
dEW = dEW - mean(dEW);
[dEW] = Butter2(dEW);
subplot(3,2,5)
fig1 = plot(time, dEW);set(gca, 'FontName', 'Times New Roman', 'FontSize', 14,'XColor','k');
xlabel('Time (secs)', 'color','k'); ylabel('(cm)', 'color','k');
grid on; box on; xlim([0 time(end)]); 

%% now do magnitude response of displacement
dEW = dEW';
dEWw = dEW.*win;
DEW = fft(dEWw); %taking the fft
N = length(DEW); %we need this for our frequency axis
DEW_mag = 4*abs(DEW)/(N); %magnitude response
window = ceil((N/fs)*width);
DEW_mag = smooth(DEW_mag, window);
bins = 0 : N-1; %samples in NS component
freq = bins*fs/N; %frequency axis NS (Hz)
subplot(3,2,6)
freqvert = plot(freq, DEW_mag);set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, 'Yscale','log', 'Xscale', 'log', 'XColor','k');
xlabel('Frequency (Hz)'), 'color', 'k';
grid on; box on;
xlim([0.03 10]); ylim([1e-4 10]);
xticks([0.1 1 10])
xticklabels({'0.1','1','10'})
yticks([0.001 0.01 0.1 1 10])
yticklabels({'0.001','0.01','0.1','1', '10'})
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

%% Plot Vertical component acceleration
%acceleration
figure
%sgtitle('Vertical','FontName', 'Times New Roman', 'FontSize', 18, 'FontWeight', 'bold', 'color','k');
subplot(3,2,1)
plot(time, xV);set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, 'XColor','k');
ylabel('(cm/s^2)', 'color','k'); 
grid on; box on; xlim([0 time(end)]); %ylim([(-a -5) (a +5) ]);

%% now do magnitude response
xVw = xV.*win;
XV = fft(xVw); %taking the fft
N = length(XV); %we need this for our frequency axis
XV_mag = 4*abs(XV)/(N); %magnitude response
window = ceil((N/fs)*width);
XV_mag = smooth(XV_mag, window);
bins = 0 : N-1; %samples in NS component
freq = bins*fs/N; %frequency axis NS (Hz)
subplot(3,2,2)
freqvert = plot(freq, XV_mag);set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, 'Yscale','log', 'Xscale', 'log', 'XColor','k')
grid on; box on;
xlim([0.03 10]); ylim([1e-4 10]);
xticks([0.1 1 10])
xticklabels({'0.1','1','10'})
yticks([0.001 0.01 0.1 1 10])
yticklabels({'0.001','0.01','0.1','1', '10'})

%% integrate for velocity (Newmark-Beta method)
vV = 0;
for i = 1:length(xV) - 1
    vV(i+1) = vV(i) + dt*((1-beta)*xV(i) + beta*xV(i+1));
end
vV = vV - mean(vV);
[vV] = Butter2(vV);
subplot(3,2,3)
plot(time, vV);set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, 'XColor','k');
ylabel('(cm/s)', 'color','k');
grid on; box on; xlim([0 time(end)]); 


%% now do magnitude response of velocity
vVw = vV';
vVw = vVw.*win;
VV = fft(vVw); %taking the fft
N = length(VV); %we need this for our frequency axis
VV_mag = 4*abs(VV)/(N); %magnitude response
window = ceil((N/fs)*width);
VV_mag = smooth(VV_mag, window);
bins = 0 : N-1; %samples in NS component
freq = bins*fs/N; %frequency axis NS (Hz)
subplot(3,2,4)
freqvert = plot(freq, VV_mag);set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, 'Yscale','log', 'Xscale', 'log', 'XColor','k');
grid on; box on;
xlim([0.03 10]); ylim([1e-4 10]);
xticks([0.1 1 10])
xticklabels({'0.1','1','10'})
yticks([0.001 0.01 0.1 1 10])
yticklabels({'0.001','0.01','0.1','1', '10'})
%PGV = max(abs(v));
%disp(strcat('PGV FN = ', num2str(PGV)));
%% Integrate for displacement 
dV = 0;
for i = 1:length(xV) - 1
    dV(i+1) = dV(i) + (dt)*vV(i) + (dt)^2 * ((0.5-alpha)*xV(i) + alpha*xV(i+1));
end
dV = dV - mean(dV);
[dV] = Butter2(dV);
subplot(3,2,5)
fig1 = plot(time, dV);set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, 'XColor','k');
xlabel('Time (secs)', 'color','k'); ylabel('(cm)', 'color','k');
grid on; box on; xlim([0 time(end)]); 

%% now do magnitude response of displacement
dV = dV';
dVw = dV.*win;
DV = fft(dVw); %taking the fft
N = length(DV); %we need this for our frequency axis
DV_mag = 4*abs(DV)/(N); %magnitude response
window = ceil((N/fs)*width);
DV_mag = smooth(DV_mag, window);
bins = 0 : N-1; %samples in NS component
freq = bins*fs/N; %frequency axis NS (Hz)
subplot(3,2,6)
freqvert = plot(freq, DV_mag);set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, 'Yscale','log', 'Xscale', 'log', 'XColor','k');
xlabel('Frequency (Hz)', 'color','k'); 
grid on; box on;
xlim([0.03 10]); ylim([1e-4 10]);
xticks([0.1 1 10])
xticklabels({'0.1','1','10'})
yticks([0.001 0.01 0.1 1 10])
yticklabels({'0.001','0.01','0.1','1', '10'})
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
vEW = vEW';
vNS = vNS';
vV = vV';
save(strcat('C:\Users\mpontr01\Box\2020_1_spring\SSA\Time domain\Data_mat_files\',statname,'.mat'),'xNS', 'vNS', 'dNS', 'xEW', 'vEW', 'dEW', 'xV', 'vV', 'dV')
% %% plot particle motion
% figure
% title('Particle motion')
% xlabel('North-south')
% ylabel('East-west')
% grid on
% plot(dEW,dNS, 'o','markeredgecolor', 'k', 'markerFacecolor', 'k','markersize', 0.2)
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
% xlabel('east-west disp'); ylabel('north-south disp'); title('particle motion')
% hold on; grid on
% % for i = 2:length(dNS)
% %     plot(dEW(i),dNS(i), 'o','markeredgecolor', 'k', 'markerFacecolor', 'k','markersize', 0.2)
% %     hold on
% % end
% 
% 
% %% do displacement HVSR
% HV_disp = DEW_mag./DV_mag;
% figure 
% freqvert = plot(freq, HV_disp, 'linewidth',1.2);set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, 'Yscale','log', 'Xscale', 'log');
% xlabel('Frequency (Hz)'); ylabel('Magnitude response'); title('HVSR displacement');
% grid on; box on;
% xlim([0.03 10]); ylim([1e-2 100]);
% xticks([0.1 1 10])
% xticklabels({'0.1','1','10'})
% yticks([0.001 0.01 0.1 1 10])
% yticklabels({'0.001','0.01','0.1','1', '10'})
% %% plot particle motion in 3d
% figure
% title('Particle motion')
% subplot(2,1,1)
% plot3(dEW(1),dNS(1),dV(1), 'o','markeredgecolor', 'k', 'markerFacecolor', 'k','markersize', 0.4)
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
% xlabel('east-west disp'); ylabel('north-south disp');zlabel('vertical disp'); title('particle motion')
% hold on; grid on
% xlim([-20 20]);ylim([-20 20]);zlim([-3 3]);
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% subplot(2,1,2)
% plot(time(1),dEW(1), 'o','markeredgecolor', 'k', 'markerFacecolor', 'k','markersize', 0.4)
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
% xlabel('time (secs)'); ylabel('displacement (cm)');
% hold on; grid on
% xlim([0 time(end)]); ylim([(-a -5) (a +5) ]);
% 
% for i = 6000:length(dNS)
%     subplot(2,1,1)
%     plot3(dEW(i),dNS(i),dV(i), 'o','markeredgecolor', 'k', 'markerFacecolor', 'k','markersize', 0.4)
%     hold on
% end


%% do displacement HVSR
HV_disp = DEW_mag./DV_mag;
figure 
freqvert = plot(freq, HV_disp, 'linewidth',1.2);set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, 'Yscale','log', 'Xscale', 'log');
xlabel('Frequency (Hz)'); ylabel('Magnitude response'); title('HVSR displacement');
grid on; box on;
xlim([0.03 10]); ylim([1e-2 100]);
xticks([0.1 1 10])
xticklabels({'0.1','1','10'})
yticks([0.001 0.01 0.1 1 10])
yticklabels({'0.001','0.01','0.1','1', '10'})
%% plot particle motion in 3d
figure
title('Particle motion')
plot3(dEW(1),dNS(1),dV(1), 'o','markeredgecolor', 'k', 'markerFacecolor', 'k','markersize', 0.4)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
xlabel('east-west disp'); ylabel('north-south disp');zlabel('vertical disp'); title('particle motion')
hold on; grid on
xlim([-20 20]);ylim([-20 20]);zlim([-3 3]);
for i = 7000:length(dNS)
    plot3(dEW(i),dNS(i),dV(i), 'o','markeredgecolor', 'k', 'markerFacecolor', 'k','markersize', 0.4)
    hold on
    pause(0.00001)
    xlim([-20 20]);ylim([-20 20]);zlim([-3 3]);
    disp(i)
end

%     
%     xlim([-20 20]);ylim([-20 20]);zlim([-3 3]);
%     subplot(2,1,2)
%     plot(time(i),dEW(i), 'o','markeredgecolor', 'k', 'markerFacecolor', 'k','markersize', 0.4)
%     set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
%     xlabel('time (secs)'); ylabel('displacement (cm)');
%     hold on; grid on
%     xlim([0 time(end)]); ylim([-15 15 ]);
%     pause(0.001)
% end
