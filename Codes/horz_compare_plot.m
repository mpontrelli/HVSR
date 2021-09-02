%% horz_compare_plot

%% Author: Marshall Pontrelli
% Date: 7/30/2021

%% Start

function horz_compare_plot(fax_HzN, NS_horz, EW_horz, statname, lowbound, upbound, outpath, fs, sav)
HVSR = figure;
hold on

NS = plot(fax_HzN(1 :length(fax_HzN)), NS_horz(1:length(fax_HzN)), 'Color', [0 0.30196 0.6588] , 'Linewidth', 1.5);
hold on
EW = plot(fax_HzN(1 :length(fax_HzN)), EW_horz(1:length(fax_HzN)), 'Color', [0.6588 0.3412 0] , 'Linewidth', 1.5);
xlabel('Frequency (Hz)','FontSize', 14)
ylabel('Amplification','FontSize', 14)
title(strcat(statname), 'FontSize', 14)
title('Horizontal comparison', 'FontSize', 14)
xlabel('Frequency (Hz)','FontSize', 14)
ylabel('Amplification','FontSize', 14)
set(gca,'YScale', 'log', 'XScale', 'log','FontName', 'Times New Roman', 'FontSize', 14)
xlim([fax_HzN(lowbound) fax_HzN(upbound)])
ylim([0.1 1000])
xticks([0.1 1 10 fs/2])
xticklabels({'0.1','1','10', num2str(fs/2)})
yticks([0.1 1 10 100 1000])
yticklabels({'0.1','1','10', '100', '1000'})
grid on 
box on
legend('NS','EW')

grid on 
box on
hold on

if strcmp(sav, 'yes') == 1
    saveas(EW, strcat(outpath, '\', 'horz_compare.jpg'));
    saveas(HVSR, strcat(outpath, '\', 'horz_compare.fig'));
end
end