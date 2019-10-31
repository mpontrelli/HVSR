%% averagedmagrespplot

%% Author: Marshall Pontrelli
% Date: 10/30/2019  

%% Start
function averagedmagrespplot(x, h, v, fs, confinthighhorz, confintlowhorz, confinthighvert, confintlowvert, lowbound, outpath, sav)
averageunfiltered = figure;
subplot(1,2,1)
hold on
confidenceinterval=shadedplot(x, confinthighhorz, confintlowhorz,[.9,.9,.9],'k');
hold on
plot(x, h, 'Color', [0 0.30196 0.6588] , 'Linewidth', 1.5);
title('Horizontal', 'FontSize', 20)
xlabel('Frequency (Hz)','FontSize', 18)
ylabel('Amplification','FontSize', 18)
set(gca,'FontSize',20,'YScale', 'log', 'Xscale', 'log')
xlim([x(1) (fs/2)])
ylim([0.1 1000])
xticks([0.1 1 10 fs/2])
xticklabels({'0.1','1','10', num2str(fs/2)})
yticks([0.1 1 10 100 1000])
yticklabels({'0.1','1','10', '100', '1000'})
grid on 
box on


subplot(1,2,2)
hold on
confidenceinterval=shadedplot(x, confinthighvert, confintlowvert,[.9,.9,.9],'k');
hold on
plot(x, v, 'Color', [0 0.30196 0.6588] , 'Linewidth', 1.5);
title('Vertical', 'FontSize', 20)
xlabel('Frequency (Hz)','FontSize', 18)
ylabel('Amplification','FontSize', 18)
set(gca,'FontSize',20,'YScale', 'log', 'Xscale', 'log')
xlim([x(1) (fs/2)])
ylim([0.1 1000])
xticks([0.1 1 10 fs/2])
xticklabels({'0.1','1','10', num2str(fs/2)})
yticks([0.1 1 10 100 1000])
yticklabels({'0.1','1','10', '100', '1000'})
grid on 
box on

%makes figure full screen
set(gcf, 'Units', 'Normalized', 'OuterPosition', [1, 0.5, 1, 0.5]);

%save
if strcmp(sav, 'yes') == 1
    saveas(averageunfiltered, strcat(outpath, '\', 'averageunfiltered.jpg'));
    saveas(averageunfiltered, strcat(outpath, '\', 'averageunfiltered.fig'));
end