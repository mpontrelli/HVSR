
filename='test_nrattle_02mar11.nrattle_amps4plot.out';
M = dlmread(filename,'',19,0);
freq=M(:,1);
amps=M(:,3);

% plot(freq,amps)
% title('Theoretical Transfer Function')
% xlabel('Frequency (Hz)','FontSize', 18)
% ylabel('Amplification','FontSize', 18)
% set(gca,'FontSize',20,'YScale', 'log')
% xlim([0 40])
% ylim([0.1 100])
% grid on

[maxamp,I]=max(amps);
disp(strcat('maximum amplitude = ', num2str(maxamp)))
FSF=freq(I);
disp(strcat('fundamental frequency = ', num2str(FSF)))
figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [1, 0.5, 1, 0.5]);
subplot(1,2,1)
TTF = plot(freq,amps);
hold on
ETF = plot(fax_HzN, HV , 'Linewidth', 1.5);
xlim([0.2 10])
xlabel('Frequency (Hz)','FontSize', 18)
ylabel('Amplification','FontSize', 18)

ylim([0.1 10])
xticks([.1 1 10])
xticklabels({'0.1', '1', '10'})
yticks([0.1 1 10 100])
yticklabels({'0.1', '1','10', '100'})
title('NP 8040 TTF-HVSR comparison')
legend([TTF,ETF], 'TTF', 'HVSR', 'FontName', 'Times New Roman', 'FontSize', 18, 'location','southeast')

set(gca,'YScale', 'log', 'XScale', 'log','FontName', 'Times New Roman', 'FontSize', 14)
grid on 
box on


h = [50,10, 50];
v = [250, 600, 1000];
h_new = 0;
for i = 1: length(h) % Sum up the thicknesses
    h_new(i+1) = h_new(i) + h(i);
end
h_new = repelem(h_new,2);
h_new(1) = []; 
h_new(end) = [];
v_new = repelem(v,2); 
subplot(1,2,2)
vel = plot(v_new, h_new, 'linewidth', 2, 'color','k');

xlabel('Velocity (m/s)')
ylabel('Depth (m)')
title('TTF velocity profile')
xlim([0 v_new(end) + v_new(1)])
ylim([0, 100])
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, 'Ydir','reverse');
grid on
box on
name = strcat(filepath,'TTF_HVSR.jpg');
saveas(vel, name, 'jpg');

%% calculate thomps
[pks,locs] = findpeaks(amps);
a = locs(1);
b = locs(3);
amps_cor = amps(a:b);
HV_cor = HV(a:b);
corrcoef(amps_cor, HV_cor)