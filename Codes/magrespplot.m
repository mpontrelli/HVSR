function magrespplot(ahatf, newfaxhz, confinthigh, confintlow, statname, xbounds, ybounds, comp, taxstat, plotcolor,lowbound, plotnum)
freq = taxstat{1,2};
num = find(newfaxhz == freq);
max1 = max(ahatf(num));
disp(strcat('max=',num2str(max1)));
subplot(1,2,plotnum)
hold on
confidenceinterval=shadedplot(newfaxhz, confinthigh, confintlow,[.9,.9,.9],'k');
hold on
ETF = plot(newfaxhz, ahatf, 'Color', plotcolor , 'Linewidth', 1.5);
hold on
%plot(newfaxhz(num), max1, 'o', 'markeredgecolor','k', 'markerfacecolor', 'k', 'markersize', 12)
title(strcat(statname, ' ', comp), 'FontSize', 20)
xlabel('Frequency (Hz)','FontSize', 18)
ylabel('Mag resp','FontSize', 18)
set(gca,'FontSize',20,'YScale', 'log', 'XScale','log')
xlim([newfaxhz(lowbound) 20])
xticks([.1 1 10])
xticklabels({'0.1', '1', '10'})
ylim([10e-5 10])
grid on    
%makes figure full screen
set(gcf, 'Units', 'Normalized', 'OuterPosition', [1, 0.5, 1, 0.5]);
end