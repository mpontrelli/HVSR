function HVSRplot(ahatf, newfaxhz, confinthigh, confintlow, station)

figure
hold on
confidenceinterval=shadedplot(newfaxhz(10:length(newfaxhz)-1), confinthigh(10:length(confinthigh)-1), confintlow(10:length(confintlow)-1),[.9,.9,.9],'k');
hold on
ETF = plot(newfaxhz(10:length(newfaxhz)-1), ahatf(10:length(ahatf)-1), 'Color', [0.149 0.45098 0] , 'Linewidth', 1.5);
title(strcat(station), 'FontSize', 20)
xlabel('Frequency (Hz)','FontSize', 18)
ylabel('Amplification','FontSize', 18)
set(gca,'FontSize',20,'YScale', 'log')
xlim([0 5])
ylim([0.1 100])
grid on    
end