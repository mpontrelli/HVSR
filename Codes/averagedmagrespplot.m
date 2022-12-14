%% averagedmagrespplot

%% Author: Marshall Pontrelli
% Date: 10/30/2019  

%% Start
function averagedmagrespplot(fax_HzN, h, v, fs, confinthighhorz, confintlowhorz, confinthighvert, confintlowvert, left_lab,right_lab, lowbound,upbound, outpath, sav, out_name)
averageunfiltered = figure;
subplot(1,2,1)
hold on
confidenceinterval=shadedplot(fax_HzN, confinthighhorz, confintlowhorz,[.9,.9,.9],[1,1,1]);
hold on
plot(fax_HzN, h, 'Color', [0 0.30196 0.6588] , 'Linewidth', 1.5);
title('', 'FontSize', 14)
xlabel('Frequency (Hz)','FontSize', 14)
ylabel('Amplification','FontSize', 14)
set(gca,'YScale', 'log', 'XScale', 'log','FontName', 'Times New Roman', 'FontSize', 14)
xlim([0.1 49])
ylim([0.1 1000])
xticks([0.1 1 10 fs/2])
xticklabels({'0.1','1','10', num2str(fs/2)})
yticks([0.1 1 10 100 1000])
yticklabels({'0.1','1','10', '100', '1000'})
grid on 
box on


subplot(1,2,2)
hold on
confidenceinterval=shadedplot(fax_HzN, confinthighvert, confintlowvert,[.9,.9,.9],[1,1,1]);
hold on
plot(fax_HzN, v, 'Color', [0 0.30196 0.6588] , 'Linewidth', 1.5);
title('', 'FontSize', 14)
xlabel('Frequency (Hz)','FontSize', 14)
ylabel('Amplification','FontSize', 14)
set(gca,'YScale', 'log', 'XScale', 'log','FontName', 'Times New Roman', 'FontSize', 14)
xlim([0.1 49])
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
%if strcmp(sav, 'yes') == 1
    %saveas(averageunfiltered, strcat(outpath, '\', strcat(out_name,'averageunfiltered.jpg')));
    %saveas(averageunfiltered, strcat(outpath, '\', strcat(out_name,'averageunfiltered.fig')));
%end