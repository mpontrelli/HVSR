%% individmagrespplot
% plot the individual magnitude responses


%% Author: Marshall Pontrelli
% Date: 10/30/2019 

%% Start
function individmagrespplot(fax_HzN, h, v, fs, left_lab,right_lab, lowbound, upbound, outpath, sav, out_name)
individualunfiltered = figure;
subplot(1,2,1)
plot(fax_HzN, h, 'Color',  'k' , 'Linewidth', .5);
title(left_lab, 'FontSize', 20)
xlabel('Frequency (Hz)','FontSize', 18)
ylabel('Amplification','FontSize', 18)
set(gca,'YScale', 'log', 'XScale', 'log','FontName', 'Times New Roman', 'FontSize', 14)
xlim([fax_HzN(lowbound) fax_HzN(upbound)])
ylim([0.1 1000])
xticks([0.1 1 10 fs/2])
xticklabels({'0.1','1','10', num2str(fs/2)})
yticks([0.1 1 10 100 1000])
yticklabels({'0.1','1','10', '100', '1000'})
grid on 
box on

subplot(1,2,2)
plot(fax_HzN, v, 'Color',  'k' , 'Linewidth', .5);
title(right_lab, 'FontSize', 20)
xlabel('Frequency (Hz)','FontSize', 18)
ylabel('Amplification','FontSize', 18)
set(gca,'YScale', 'log', 'XScale', 'log','FontName', 'Times New Roman', 'FontSize', 14)
xlim([fax_HzN(lowbound) fax_HzN(upbound)])
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
    saveas(individualunfiltered, strcat(outpath, '\', strcat(out_name,'individualunfiltered.jpg')));
    saveas(individualunfiltered, strcat(outpath, '\', strcat(out_name,'individualunfiltered.fig')));
end