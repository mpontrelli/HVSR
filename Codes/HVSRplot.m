function HVSRplot(ahatf, newfaxhz, confinthigh, confintlow, lowbound, statname)
figure
hold on
confidenceinterval=shadedplot(newfaxhz(lowbound:length(newfaxhz)), confinthigh(lowbound:length(confinthigh)), confintlow(lowbound:length(confintlow)),[.9,.9,.9],'k');
hold on
ETF = plot(newfaxhz(lowbound:length(newfaxhz)), ahatf(lowbound:length(ahatf)), 'Color', [0 0.30196 0.6588] , 'Linewidth', 1.5);
title(statname)
xlabel('Frequency (Hz)','FontSize', 18)
ylabel('Amplification','FontSize', 18)
set(gca,'FontSize',20,'YScale', 'log', 'XScale', 'log')
xlim([newfaxhz(lowbound) 20])
ylim([0.1 100])
xticks([.1 1 10])
xticklabels({'0.1', '1', '10'})
yticks([0.1 1 10 100])
yticklabels({'0.1', '1','10', '100'})
grid on 
box on
end