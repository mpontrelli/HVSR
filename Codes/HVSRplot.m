%% HVSRplot
% plot the HVSR lognormal average with upper and lower confidence intervals

    % INPUTS
    
    % ahatf - lognormal average spectral ratio
    
    % newfaxhz - frequency vector
    
    % confinthigh - upper 95% confidence interval
    
    % confintlow - lower 95% confidence interval
    
    % lowbound - low resolvable frequency
    
    % statname - station name (title of plot)
    
    % OUTPUTS
    
    % A plot of the HVSR
    
%% Author: Marshall Pontrelli
% Date: developed between September, 2017 and August, 2019
%% Start  
function ETF = HVSRplot(ahatf, newfaxhz, confinthigh, confintlow, lowbound, statname, plotcolor)
figure
hold on
confidenceinterval=shadedplot(newfaxhz(lowbound:length(newfaxhz)), confinthigh(lowbound:length(confinthigh)), confintlow(lowbound:length(confintlow)),[.9,.9,.9],[1 1 1]);
hold on
ETF = plot(newfaxhz(lowbound:length(newfaxhz)), ahatf(lowbound:length(ahatf)), 'Color', plotcolor , 'Linewidth', 1.5);
title(statname)
xlabel('Frequency (Hz)','FontSize', 18)
ylabel('Amplification','FontSize', 18)
set(gca,'YScale', 'log', 'XScale', 'log','FontName', 'Times New Roman', 'FontSize', 14)
xlim([newfaxhz(lowbound) 20])
ylim([0.1 100])
xticks([.1 1 10])
xticklabels({'0.1', '1', '10'})
yticks([0.1 1 10 100])
yticklabels({'0.1', '1','10', '100'})
grid on 
box on
end