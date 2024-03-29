%% Plot HVSR - SSR together Mexico_city
close all
clear all


sitelist = {'AE02','AL01','AO24', 'AP68','AR14','AU11','AU46','BA49','BL45','BO39',...
    'CA20', 'CA59', 'CB43','CC55','CE18', 'CE23','CE32','CH84','CI05','CJ03',...
    'CJ04','CO47', 'CO56','CP28','CS78','CT64', 'CU80', 'DM12','DR16', 'DX37','EO30','ES57','EX08','EX09','EX12','FJ74','GA62',...
    'GC38','GR27','HA41','HJ72', 'IB22','IM40','JA43','JC54','LI33', 'LI58', 'LV17','ME52',...
    'MI15', 'MY19', 'NZ20', 'NZ31','PA34', 'PD42','PE10', 'RI76', 'RM48',...
    'SI53', 'SP51', 'TE07', 'TH35', 'TL08', 'TL55', 'TP13', 'UC44', 'VG09', 'VM29',...
    'XP06'};


% Load SSR
load(strcat('C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\SSR_Shape_statistics\',...
    statname, '\',statname, '-TP13'))

ahatfSSR = data.complex.ahatf';
sigmaSSR = data.complex.sigma';
confinthighSSR = data.complex.confinthigh';
confintlowSSR  = data.complex.confintlow';

% load HVSR
load(strcat('C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Shape_statistics\',statname));
ahatfHV = shapedata.complex.ahatf';
sigmaHV = shapedata.complex.sigma';
confinthighHV = shapedata.complex.confinthigh';
confintlowHV  = shapedata.complex.confintlow';

%%

% load frequency vector
load('C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Processed_data2\AE02\AE0219900511234349')
fax_HzN = data.processing.filtereddata.freq_vec;


% Now some inputs
upbound = 10;
lowbound = 0.1;
[~, lowbound] = min(abs(fax_HzN - lowbound));
[~, upbound] = min(abs(fax_HzN - upbound));

%% Now find fundamental resonance
[peakfreqsSSR, peakampsSSR, hpb, f1s, f2s, Areamat, proms, ampsSSR, peakind2, freqsSSR, sigs, I1s, I2s] = peakiden(ahatfSSR, fax_HzN, sigmaSSR, lowbound, upbound);
[peakfreqsHVSR, peakampsHVSR, hpb, f1s, f2s, Areamat, proms, ampsHVSR, peakind2, freqsHVSR, sigs, I1s, I2s] = peakiden(ahatfHV, fax_HzN, sigmaHV, lowbound, upbound);


%% now plot, make the figure and set the base
figure
xlim([0.1 10])
ylim([0.1 100])
xticks([.1 1 10])
xticklabels({'0.1', '1', '10'})
yticks([1 10 100])
yticklabels({ '1','10', '100'})
xlabel('Frequency (Hz)')
ylabel('Amplification')
title(statname)
set(gca,'YScale', 'log','XScale','log', 'FontName', 'Times New Roman', 'FontSize', 18)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
grid on
box on
hold on
%% start with the SSR confidence interval
fr = fax_HzN(lowbound:length(fax_HzN))';
cohr = confinthighSSR(lowbound:length(fax_HzN))';
colr = confintlowSSR(lowbound:length(fax_HzN))';
x_plot =[fr, fliplr(fr)];
y_plot = [cohr, fliplr(colr)];

SSRconf = fill(x_plot, y_plot, 1,'facecolor', [0 0 0.5], 'edgecolor', 'none', 'facealpha', 0.4);

% confidenceinterval=shadedplot(fr, cohr, colr,[0 0 0.5],[1 1 1], 'FaceAlpha',0.5);
hold on

%% Now HVSR confidence interval
fr = fax_HzN(lowbound:length(fax_HzN))';
cohr = confinthighHV(lowbound:length(fax_HzN))';
colr = confintlowHV(lowbound:length(fax_HzN))';
x_plot =[fr, fliplr(fr)];
y_plot = [cohr, fliplr(colr)];

HVconf =fill(x_plot, y_plot, 1,'facecolor', [0 0.5 0], 'edgecolor', 'none', 'facealpha', 0.2);

% confidenceinterval=shadedplot(fr, cohr, colr,[.9,.9,.9],[1 1 1], 'FaceAlpha',0.5);
hold on

%% Now plot HVSR
HV_plot = plot(fax_HzN,ahatfHV, 'LineWidth', 2, 'Color', [0 0.5 0]);
SSR_plot = plot(fax_HzN,ahatfSSR, 'LineWidth', 2, 'Color', [0 0 0.5]);



%% Now add the fundamental resonance

% SSR
if length(ampsSSR) > 0  
    [SSRmax, Issr] = max(ampsSSR);
    SSRfn = line([freqsSSR(Issr),freqsSSR(Issr)],[0.1, 100],'LineStyle', '--', 'color',[0 0 0.5]);
end

% SSR
if length(ampsHVSR) > 0  
    [HVSRmax, Issr] = max(ampsHVSR);
    HVfn = line([freqsHVSR(Issr),freqsHVSR(Issr)],[0.1, 100],'LineStyle', '--', 'color',[0 0.5 0]);
end
hold on

%% Now do cross correlation
ahatfHV2 = ahatfHV - mean(ahatfHV);
ahatfSSR2 = ahatfSSR - mean(ahatfSSR);
[c,lags] = xcorr(ahatfHV2, ahatfSSR2,'normalized');
% figure
% plot(lags/5000,c)

[maxr, I] = max(c);
lagmax = lags(I)/5000;

r = c(50000);

datamat = [r,maxr,lagmax];

%% now paste info
a = "r =" + " "  + num2str(r,3);
b = "max r =" + " "  + num2str(maxr,3);
c = strcat("lag =" + " "  +num2str(lagmax,3) + " " + "hz");
str = {a, b, c};
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
       
text(2.25,30,str, 'FontName', 'Times New Roman', 'FontSize', 18,'Color', 'black', 'HorizontalAlignment', 'left', 'EdgeColor','k','BackgroundColor', 'w')
%% legend

legend([HV_plot, HVconf, HVfn, SSR_plot, SSRconf, SSRfn], 'HV', 'HV_{conf}', 'HV_{fn}', 'SSR', 'SSR_{conf}','SSR_{fn}')


