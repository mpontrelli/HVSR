function HVSRplot(ahatf, newfaxhz, sigma, confinthigh, confintlow, station)

newfaxhz=newfaxhz(10:length(newfaxhz)-1);
confinthigh=confinthigh(10:length(confinthigh)-1);
confintlow=confintlow(10:length(confintlow)-1);

figure
hold on
confidenceinterval=shadedplot(newfaxhz, confinthigh, confintlow,[.9,.9,.9],'k');
hold on
ahatf=ahatf(10:length(ahatf)-1);
ETF = plot(newfaxhz, ahatf, 'Color', [0.149 0.45098 0] , 'Linewidth', 1.5);
title(strcat(station), 'FontSize', 20)
xlabel('Frequency (Hz)','FontSize', 18)
ylabel('Amplification','FontSize', 18)
set(gca,'FontSize',20,'YScale', 'log')
xlim([0 5])
ylim([0.1 100])
grid on    
end