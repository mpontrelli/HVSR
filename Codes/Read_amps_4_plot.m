
filename=strcat(pwd, '\', 'NRATTLE', '\', 'test_nrattle_02mar11.nrattle_amps4plot.out');
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