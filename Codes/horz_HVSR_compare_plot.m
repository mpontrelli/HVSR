%% horz_HVSR_compare_plot

%% Author: Marshall Pontrelli
% Date: 7/30/2021

%% Start
function horz_HVSR_compare_plot(fax_HzN, NS,EW, lowbound, upbound, outpath, sav)
HVSR = figure;
hold on
NS = plot(fax_HzN(1 :length(fax_HzN)), NS(1:length(fax_HzN)), 'Color', [0 0.30196 0.6588] , 'Linewidth', 1.5);
hold on
EW = plot(fax_HzN(1 :length(fax_HzN)), EW(1:length(fax_HzN)), 'Color', [0.6588 0.3412 0] , 'Linewidth', 1.5);
xlabel('Frequency (Hz)','FontSize', 14)
ylabel('Amplification','FontSize', 14)
title('HVSR horizontal comparison','FontSize', 14)
set(gca,'YScale', 'log','XScale','log','FontName', 'Times New Roman', 'FontSize', 14)
xlim([fax_HzN(lowbound) fax_HzN(upbound)])
%xlim([fax_HzN(1) 40])
ylim([0.1 100])
xticks([0.1 1 10])
xticklabels({'0.1','1','10', num2str(upbound)})
yticks([0.1 1 10 100])
yticklabels({'0.1','1','10', '100'})
legend('NS','EW')

grid on 
box on
hold on




%makes figure full screen
%set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

%save
if strcmp(sav, 'yes') == 1
    saveas(HVSR, strcat(outpath, '\', 'HVSR_horz_compare.jpg'));
    saveas(HVSR, strcat(outpath, '\', 'HVSR_horz_compare.fig'));
end
end