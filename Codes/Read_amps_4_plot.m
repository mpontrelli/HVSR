
filename='test_nrattle_02mar11.nrattle_amps4plot.out';
M = dlmread(filename,'',19,0);
freq=M(:,1);
amps=M(:,3);

% plot(freq,amps)
% xlabel('Frequency (hz)')
% ylabel('Amplitude')
% title('Theoretical Transfer Function')
% grid on
%ylim([0 15])

[maxamp,I]=max(amps);
disp(strcat('maximum amplitude = ', num2str(maxamp)))
FSF=freq(I);
disp(strcat('fundamental frequency = ', num2str(FSF)))