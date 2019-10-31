%% individmagrespplot
% plot the individual magnitude responses


%% Author: Marshall Pontrelli
% Date: 10/30/2019 

%% Start
function individmagrespplot(x, h, v, fs, lowbound, outpath, sav)
individualunfiltered = figure;
subplot(1,2,1)
plot(x, h, 'Color',  'k' , 'Linewidth', .5);
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
plot(x, v, 'Color',  'k' , 'Linewidth', .5);
title('Vertical', 'FontSize', 20)
xlabel('Frequency (Hz)','FontSize', 18)
ylabel('Amplification','FontSize', 18)
set(gca,'FontSize',20,'YScale', 'log', 'Xscale', 'log')
xlim([x(1) fs/2])
ylim([0.1 1000])
xticks([.1 1 10 fs/2])
xticklabels({'.1', '1','10', num2str(fs/2)})
yticks([0.1 1 10 100 1000])
yticklabels({'0.1','1','10', '100', '1000'})
grid on 
box on

%makes figure full screen
set(gcf, 'Units', 'Normalized', 'OuterPosition', [1, 0.5, 1, 0.5]);

%save
if strcmp(sav, 'yes') == 1
    saveas(individualunfiltered, strcat(outpath, '\', 'individualunfiltered.jpg'));
    saveas(individualunfiltered, strcat(outpath, '\', 'individualunfiltered.fig'));
end