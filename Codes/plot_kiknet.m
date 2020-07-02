close all
clear all


load('C:\Users\mpontr01\Box\Projects\Kik-Net\Processed_data\TKCH08\TKCH080012080614.mat')
fax_HzN = data.processing.filtereddata.freq_vec;
load('C:\Users\mpontr01\Box\Projects\Kik-Net\Shape_data\TKCH08')
statname = 'TKCH08';
upbound = 20;
lowbound = 0.1;
[~, lowbound] = min(abs(fax_HzN - lowbound));
[~, upbound] = min(abs(fax_HzN - upbound));

%% HV
ahatf = shapedata.HV.geo2.ahatf';
sigma = shapedata.HV.geo2.sigma';
confinthigh = shapedata.HV.geo2.confinthigh';
confintlow = shapedata.HV.geo2.confintlow';

peakfreqs = shapedata.HV.geo2.shapes.peak_freqs;
peakamps = shapedata.HV.geo2.shapes.peak_amps;
hpb = shapedata.HV.geo2.shapes.hpb;
f1s = shapedata.HV.geo2.shapes.f1s;
f2s = shapedata.HV.geo2.shapes.f2s;
Areamat = shapedata.HV.geo2.shapes.Area;
proms = shapedata.HV.geo2.shapes.proms;
amps = shapedata.HV.geo2.shapes.amps;
peakind2 = shapedata.HV.geo2.shapes.peakind;
freqs = shapedata.HV.geo2.shapes.freqs;
sigs = shapedata.HV.geo2.shapes.sigs;
I1s = shapedata.HV.geo2.shapes.I1s;
I2s = shapedata.HV.geo2.shapes.I2s;

%% Now set some conditions to classify the peak
[~, freq_class] = sort(freqs, 'descend');
[~, amp_class] = sort(amps, 'descend');
[~, prom_class] = sort(proms, 'descend');
[~, hpb_class] = sort(hpb, 'descend');
[~, area_class] = sort(Areamat, 'descend');
[~, sig_class] = sort(sigs, 'descend');
classmatrix = vertcat(freq_class, amp_class, prom_class, hpb_class, area_class, sig_class);

%% Now save peak data
datamat = vertcat(freqs,amps,proms,hpb,Areamat,sigs);

%% now plot, make the figure and set the base
fin_fig = figure;
xlim([0.1 20])
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

%% start with the confidence interval
fr = fax_HzN(lowbound:length(fax_HzN))';
cohr = confinthigh(lowbound:length(fax_HzN))';
colr = confintlow(lowbound:length(fax_HzN))';
confidenceinterval=shadedplot(fr, cohr, colr,[.9,.9,.9],[1 1 1]);
hold on
    
%% Plot filled in peaks color coordinated based on their area
if length(amps) > 0
    for i = 1:length(amps)
        hold on
        freq_ar = peakfreqs{i};
        amp_ar = peakamps{i};
        if i == 1
            col = 'b';
        end
        if i == 2
            col = 'r';
        end
        if i == 3
            col = 'g';
        end
        if i > 3
            col = 'k';
        end
        fill(freq_ar, amp_ar, col, 'LineStyle','none','FaceAlpha',0.5)
        hold on
    end
end
if length(amps) > 0
    for i = 1:length(amps) 
        freq_ar = peakfreqs{i};
        amp_ar = peakamps{i};
        ampd = amp_ar(I1s(i));
        ampdd = amp_ar(I2s(i));
        plot([f1s(i),f2s(i)],[ampd, ampdd], 'c', 'LineWidth', 2)
        hold on
    end
end
hold on


%% now plot fundamental resonance
if length(amps) > 0
    for i = 1:length(f1s)
        line([freqs(i),freqs(i)],[0.1, 100],'LineStyle', '--', 'color','k')
        hold on
    end
end



%% Now the hpb-prominence cross
hold on
if length(amps) > 0
    for i = 1:length(amps)
        plot([freqs(i),freqs(i)],[amps(i), amps(i) - proms(i)], 'c', 'Linewidth', 2)
        hold on
    end
end

%% Now plot ahatf
plot(fax_HzN,ahatf, 'LineWidth', 2, 'Color', [0 0.5 0])
%% Now print some important info
if length(amps) > 0
    if length(amps) ==1
        a = "Max peak freq =" + " "  + num2str(freqs,3);
        b = strcat("Amp =" + " "  + num2str(amps,3));
        c = strcat("HPB =" + " "  +num2str(hpb,3));
        d = strcat("Prom =" + " "  +num2str(proms,3));
        e = strcat("Area =" + " "  +num2str(Areamat,3));
        f = strcat("\sigma =" + " "  +num2str(sigs,3));
        str = {a, b, c, d, e, f};
        ylim=get(gca,'ylim');
        xlim=get(gca,'xlim');
        q = find(ahatf == amps);
        if fax_HzN(q) < 1
             text(xlim(2)-5,ylim(2)-60,str, 'FontName', 'Times New Roman', 'FontSize', 18,'Color', 'black', 'HorizontalAlignment', 'right', 'EdgeColor','k','BackgroundColor', 'w')
        else
             text(xlim(1)+0.3,ylim(2)-60,str, 'FontName', 'Times New Roman', 'FontSize', 18,'Color', 'black', 'HorizontalAlignment', 'right', 'EdgeColor','k','BackgroundColor', 'w')
        end
        else
            [~, I] = max(amps);
            a = "Max peak freq =" + " "  + num2str(freqs(I),3);
            b = strcat("Amp =" + " "  + num2str(amps(I),3));
            c = strcat("HPB =" + " "  +num2str(hpb(I),3));
            d = strcat("Prom =" + " "  +num2str(proms(I),3));
            e = strcat("Area =" + " "  +num2str(Areamat(I),3));
            f = strcat("\sigma =" + " "  +num2str(sigs(I),3));
            str = {a, b, c, d, e, f};
            ylim=get(gca,'ylim');
            xlim=get(gca,'xlim');
            q = find(ahatf == amps(I));
            if fax_HzN(q) < 1
                text(xlim(2)-5,ylim(2)-60,str, 'FontName', 'Times New Roman', 'FontSize', 18,'Color', 'black', 'HorizontalAlignment', 'right', 'EdgeColor','k','BackgroundColor', 'w')
            else
                text(xlim(1)+0.3,ylim(2)-60,str, 'FontName', 'Times New Roman', 'FontSize', 18,'Color', 'black', 'HorizontalAlignment', 'right', 'EdgeColor','k','BackgroundColor', 'w')
            end

    end
    else
       str = 'No peaks';
       ext(0.3,15,str, 'FontName', 'Times New Roman', 'FontSize', 24,'Color', 'black', 'HorizontalAlignment', 'right', 'EdgeColor','k','BackgroundColor', 'w')
    
end
clear xlim
clear ylim

%% SSR
ahatf = shapedata.BSR.geo.ahatf';
sigma = shapedata.BSR.geo.sigma';
confinthigh = shapedata.BSR.geo.confinthigh';
confintlow = shapedata.BSR.geo.confintlow';

peakfreqs = shapedata.BSR.geo.shapes.peak_freqs;
peakamps = shapedata.BSR.geo.shapes.peak_amps;
hpb = shapedata.BSR.geo.shapes.hpb;
f1s = shapedata.BSR.geo.shapes.f1s;
f2s = shapedata.BSR.geo.shapes.f2s;
Areamat = shapedata.BSR.geo.shapes.Area;
proms = shapedata.BSR.geo.shapes.proms;
amps = shapedata.BSR.geo.shapes.amps;
peakind2 = shapedata.BSR.geo.shapes.peakind;
freqs = shapedata.BSR.geo.shapes.freqs;
sigs = shapedata.BSR.geo.shapes.sigs;
I1s = shapedata.BSR.geo.shapes.I1s;
I2s = shapedata.BSR.geo.shapes.I2s;

%% Now set some conditions to classify the peak
[~, freq_class] = sort(freqs, 'descend');
[~, amp_class] = sort(amps, 'descend');
[~, prom_class] = sort(proms, 'descend');
[~, hpb_class] = sort(hpb, 'descend');
[~, area_class] = sort(Areamat, 'descend');
[~, sig_class] = sort(sigs, 'descend');
classmatrix = vertcat(freq_class, amp_class, prom_class, hpb_class, area_class, sig_class);

%% Now save peak data
datamat = vertcat(freqs,amps,proms,hpb,Areamat,sigs);

%% now plot, make the figure and set the base
fin_fig = figure;
xlim([0.1 20])
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

%% start with the confidence interval
fr = fax_HzN(lowbound:length(fax_HzN))';
cohr = confinthigh(lowbound:length(fax_HzN))';
colr = confintlow(lowbound:length(fax_HzN))';
confidenceinterval=shadedplot(fr, cohr, colr,[.9,.9,.9],[1 1 1]);
hold on
    
%% Plot filled in peaks color coordinated based on their area
if length(amps) > 0
    for i = 1:length(amps)
        hold on
        freq_ar = peakfreqs{i};
        amp_ar = peakamps{i};
        if i == 1
            col = 'b';
        end
        if i == 2
            col = 'r';
        end
        if i == 3
            col = 'g';
        end
        if i > 3
            col = 'k';
        end
        fill(freq_ar, amp_ar, col, 'LineStyle','none','FaceAlpha',0.5)
        hold on
    end
end
if length(amps) > 0
    for i = 1:length(amps) 
        freq_ar = peakfreqs{i};
        amp_ar = peakamps{i};
        ampd = amp_ar(I1s(i));
        ampdd = amp_ar(I2s(i));
        plot([f1s(i),f2s(i)],[ampd, ampdd], 'c', 'LineWidth', 2)
        hold on
    end
end
hold on


%% now plot fundamental resonance
if length(amps) > 0
    for i = 1:length(f1s)
        line([freqs(i),freqs(i)],[0.1, 100],'LineStyle', '--', 'color','k')
        hold on
    end
end



%% Now the hpb-prominence cross
hold on
if length(amps) > 0
    for i = 1:length(amps)
        plot([freqs(i),freqs(i)],[amps(i), amps(i) - proms(i)], 'c', 'Linewidth', 2)
        hold on
    end
end

%% Now plot ahatf
plot(fax_HzN,ahatf, 'LineWidth', 2, 'Color', [0 0.5 0])
%% Now print some important info
if length(amps) > 0
    if length(amps) ==1
        a = "Max peak freq =" + " "  + num2str(freqs,3);
        b = strcat("Amp =" + " "  + num2str(amps,3));
        c = strcat("HPB =" + " "  +num2str(hpb,3));
        d = strcat("Prom =" + " "  +num2str(proms,3));
        e = strcat("Area =" + " "  +num2str(Areamat,3));
        f = strcat("\sigma =" + " "  +num2str(sigs,3));
        str = {a, b, c, d, e, f};
        ylim=get(gca,'ylim');
        xlim=get(gca,'xlim');
        q = find(ahatf == amps);
        if fax_HzN(q) < 1
             text(xlim(2)-5,ylim(2)-60,str, 'FontName', 'Times New Roman', 'FontSize', 18,'Color', 'black', 'HorizontalAlignment', 'right', 'EdgeColor','k','BackgroundColor', 'w')
        else
             text(xlim(1)+0.3,ylim(2)-60,str, 'FontName', 'Times New Roman', 'FontSize', 18,'Color', 'black', 'HorizontalAlignment', 'right', 'EdgeColor','k','BackgroundColor', 'w')
        end
        else
            [~, I] = max(amps);
            a = "Max peak freq =" + " "  + num2str(freqs(I),3);
            b = strcat("Amp =" + " "  + num2str(amps(I),3));
            c = strcat("HPB =" + " "  +num2str(hpb(I),3));
            d = strcat("Prom =" + " "  +num2str(proms(I),3));
            e = strcat("Area =" + " "  +num2str(Areamat(I),3));
            f = strcat("\sigma =" + " "  +num2str(sigs(I),3));
            str = {a, b, c, d, e, f};
            ylim=get(gca,'ylim');
            xlim=get(gca,'xlim');
            q = find(ahatf == amps(I));
            if fax_HzN(q) < 1
                text(xlim(2)-5,ylim(2)-60,str, 'FontName', 'Times New Roman', 'FontSize', 18,'Color', 'black', 'HorizontalAlignment', 'right', 'EdgeColor','k','BackgroundColor', 'w')
            else
                text(xlim(1)+0.3,ylim(2)-60,str, 'FontName', 'Times New Roman', 'FontSize', 18,'Color', 'black', 'HorizontalAlignment', 'right', 'EdgeColor','k','BackgroundColor', 'w')
            end

    end
    else
       str = 'No peaks';
       ext(0.3,15,str, 'FontName', 'Times New Roman', 'FontSize', 24,'Color', 'black', 'HorizontalAlignment', 'right', 'EdgeColor','k','BackgroundColor', 'w')
    
end
clear xlim
clear ylim
%% Now plot all together
%% Now load TTF
freq = xlsread('C:\Users\mpontr01\Box\2020_1_spring\Research\Proposals\HV_classification\TFs', 1, 'A1:A1000');
TF = xlsread('C:\Users\mpontr01\Box\2020_1_spring\Research\Proposals\HV_classification\TFs', 1, 'B1:B1000');
TF = interp1(freq, TF, fax_HzN);

%% loadSSR (actually BSR)
ahatfSSR = shapedata.BSR.geo.ahatf';
sigmaSSR = shapedata.BSR.geo.sigma';
confinthighSSR = shapedata.BSR.geo.confinthigh';
confintlowSSR = shapedata.BSR.geo.confintlow';
ampsSSR = shapedata.BSR.geo.shapes.amps;
freqsSSR = shapedata.BSR.geo.shapes.freqs;


%% HV
ahatfHV = shapedata.HV.geo2.ahatf';
sigmaHV = shapedata.HV.geo2.sigma';
confinthighHV = shapedata.HV.geo2.confinthigh';
confintlowHV = shapedata.HV.geo2.confintlow';
ampsHV = shapedata.HV.geo2.shapes.amps;
freqsHV = shapedata.HV.geo2.shapes.freqs;


% now plot
figure
xlim([0.1 20])
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
TTF = plot(fax_HzN,TF, 'LineWidth', 2, 'Color', 'k');


%% Now add the fundamental resonance

% SSR
if length(ampsSSR) > 0  
    [SSRmax, Issr] = max(ampsSSR);
    SSRfn = line([freqsSSR(Issr),freqsSSR(Issr)],[0.1, 100],'LineStyle', '--', 'color',[0 0 0.5]);
end

% SSR
if length(ampsHV) > 0  
    [HVSRmax, Issr] = max(ampsHV);
    HVfn = line([freqsHV(Issr),freqsHV(Issr)],[0.1, 100],'LineStyle', '--', 'color',[0 0.5 0]);
end
hold on

%% Now do cross correlation
[pks,locs] = findpeaks(TF);
I1 = locs(1);
I2 = locs(4);
ahatfHV = ahatfHV(I1:I2);
ahatfSSR = ahatfSSR(I1:I2);
TTF2 = TF(I1:I2);
ahatfHV2 = ahatfHV - mean(ahatfHV);
ahatfSSR2 = ahatfSSR - mean(ahatfSSR);
TTF2 = TTF2 - mean(TTF2);

% Between HVSR and BSR
r = corrcoef(ahatfHV2, ahatfSSR2);
r = r(2,1);


% Between HVSR and TTF
r1 = corrcoef(ahatfHV2, TTF2);
r1 = r1(2,1);
% Between BSR and TTF
r2 = corrcoef(ahatfSSR2, TTF2);
r2 = r2(2,1);
rs = [r r1 r2];
%% now paste info
a = "HVSR-BSR =" + " "  + num2str(r,3);
b = "HVSR-TTF =" + " "  + num2str(r1,3);
c = strcat("BSR-TTF =" + " "  +num2str(r2,3));
str = {a, b, c};
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
       
text(0.3,30,str, 'FontName', 'Times New Roman', 'FontSize', 18,'Color', 'black', 'HorizontalAlignment', 'left', 'EdgeColor','k','BackgroundColor', 'w')
%% legend
legend([HV_plot, HVconf, HVfn, SSR_plot, SSRconf, SSRfn, TTF], 'HV', 'HV_{conf}', 'HV_{fn}', 'SSR', 'SSR_{conf}','SSR_{fn}', 'TTF', 'location', 'northwest')


