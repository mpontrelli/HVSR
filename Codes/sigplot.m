
figure
freq = taxstat{1,2};
num = find(newfaxhz == freq);

sigma1 = smooth(sigma, (length(sigma)/fsmin)/2);
plot(newfaxhz(lowbound:length(newfaxhz)), sigma1(lowbound:length(sigma1)), 'linewidth',2, 'color', [0.5176, 0, 0.6588])
hold on
plot(newfaxhz(num), sigma1(num), 'o', 'markeredgecolor','k', 'markerfacecolor', 'k', 'markersize', 12)
title('AO24 sigma')
xlabel('Frequency (Hz)','FontSize', 18)
ylabel('Sigma','FontSize', 18)
set(gca,'FontSize',20)
xlim([0 10])
ylim([0.15 0.65])
grid on