function wavplot(ahatf, newfaxhz, confinthigh, confintlow, statname)
figure
hold on
confidenceinterval=shadedplot(newfaxhz, confinthigh, confintlow,[.9,.9,.9],'k');
hold on
ETF = plot(newfaxhz, ahatf, 'Color', [0.149 0.45098 0] , 'Linewidth', 1.5);
title(strcat(statname), 'FontSize', 20)
xlabel('Frequency (Hz)','FontSize', 18)
ylabel('Amplification','FontSize', 18)
set(gca,'FontSize',20,'YScale', 'log')
%xlim([0 5])
%ylim([0.1 100])
grid on    
end