function magrespplot(ahatf, newfaxhz, confinthigh, confintlow, statname, xbounds, ybounds, comp)
figure
hold on
confidenceinterval=shadedplot(newfaxhz, confinthigh, confintlow,[.9,.9,.9],'k');
hold on
ETF = plot(newfaxhz, ahatf, 'Color', [0.149 0.45098 0] , 'Linewidth', 1.5);
title(strcat(statname, ' Magnitude response', ' ', comp), 'FontSize', 20)
xlabel('Frequency (Hz)','FontSize', 18)
ylabel('Mag resp','FontSize', 18)
set(gca,'FontSize',20,'YScale', 'log')
xlim(xbounds)
ylim(ybounds)
grid on    
end