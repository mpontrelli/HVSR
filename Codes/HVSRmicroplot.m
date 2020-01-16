%% HVSRmicropplot

%% Author: Marshall Pontrelli
% Date: 10/30/2019  

%% Start
function HVSRmicroplot(fax_HzN, ahatf, confinthigh, confintlow, statname, lowbound, upbound, outpath, sav, TTF)
HVSR = figure;
hold on
confidenceinterval=shadedplot(fax_HzN(1:length(fax_HzN)), confinthigh(1:length(fax_HzN)), confintlow(1:length(fax_HzN)),[.9,.9,.9],'k');
hold on
ETF = plot(fax_HzN(1 :length(fax_HzN)), ahatf(1:length(fax_HzN)), 'Color', [0 0.30196 0.6588] , 'Linewidth', 1.5);
xlabel('Frequency (Hz)','FontSize', 18)
ylabel('Amplification','FontSize', 18)
title(strcat(statname), 'FontSize', 18)
set(gca,'FontSize',20,'YScale', 'log', 'XScale','log')
xlim([fax_HzN(lowbound) fax_HzN(upbound)])
set(gca,'FontSize',20,'YScale', 'log')
%xlim([fax_HzN(1) 40])
ylim([0.1 100])
% xticks([lowbound 0.1 1 10 upbound])
% xticklabels({num2str(lowbound),'1','10', num2str(upbound)})
yticks([0.1 1 10 100])
yticklabels({'0.1','1','10', '100'})

grid on 
box on
hold on



%% if you want to plot TTF from NRATTLE
if strcmp(TTF, 'yes') == 1
    Read_amps_4_plot
    TTF = plot(freq, amps, 'Color', 'k', 'linewidth', 2);
    legend([ETF, TTF], {'HVSR', 'TTF'})
end

%makes figure full screen
%set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

%save
if strcmp(sav, 'yes') == 1
    saveas(HVSR, strcat(outpath, '\', 'HVSR.jpg'));
    saveas(HVSR, strcat(outpath, '\', 'HVSR.fig'));
end
end