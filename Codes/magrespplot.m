function magrespplot(ahatf, newfaxhz, confinthigh, confintlow, statname, xbounds, ybounds, comp, taxstat)
freq = taxstat{1,2};
num = find(newfaxhz == freq);
max1 = max(ahatf(num));
disp(strcat('max=',num2str(max1)));
figure
hold on
confidenceinterval=shadedplot(newfaxhz, confinthigh, confintlow,[.9,.9,.9],'k');
hold on
ETF = plot(newfaxhz, ahatf, 'Color', [0.5176 0 0.6588] , 'Linewidth', 1.5);
hold on
plot(newfaxhz(num), max1, 'o', 'markeredgecolor','k', 'markerfacecolor', 'k', 'markersize', 12)
title(strcat(statname, ' ', comp), 'FontSize', 20)
xlabel('Frequency (Hz)','FontSize', 18)
ylabel('Mag resp','FontSize', 18)
set(gca,'FontSize',20,'YScale', 'log')
xlim([0 10])
ylim([10e-5 10])
grid on    
end