function HVSRplot(ahatf, newfaxhz, confinthigh, confintlow, lowbound, statname)
figure
hold on
confidenceinterval=shadedplot(newfaxhz(lowbound:length(newfaxhz)), confinthigh(lowbound:length(confinthigh)), confintlow(lowbound:length(confintlow)),[.9,.9,.9],'k');
hold on
ETF = plot(newfaxhz(lowbound:length(newfaxhz)), ahatf(lowbound:length(ahatf)), 'Color', [0.149 0.45098 0] , 'Linewidth', 1.5);
title(strcat(statname, ' HVSR'), 'FontSize', 20)
xlabel('Frequency (Hz)','FontSize', 18)
ylabel('Ampli fication','FontSize', 18)
set(gca,'FontSize',20,'YScale', 'log')
xlim([newfaxhz(lowbound) 49])
%ylim([0.1 100])
grid on    
hold off
end