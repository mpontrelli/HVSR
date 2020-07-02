%% Try azimuths
% Author: Marshall Pontrelli
% Date: 5/26/2020

% For HVSR taxonomy USGS proposal 2020
close all 
clear all


% go into the data folder and get a list of stations
codepath = 'C:\Users\mpontr01\Desktop\HVSR\Codes';
datapath = 'C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Processed_data2\CE32';
figpath = 'C:\Users\mpontr01\Box\2020_1_spring\Research\Proposals\HV_classification\figures\HVs\';
% 'AE02','AL01','AO24', 'AP68','AR14','AU11','AU46','BA49','BL45','BO39',...
%     'CA20', 'CA59', 'CB43','CC55', 'CE23','CE32','CH84','CI05','CJ03',...
%     'CJ04','CO47', 'CO56', 'CU80', 'DM12','DR16', 'DX37','EO30','ES57','EX08','EX09','EX12','GA62',...
%     'GC38','GR27','HA41','HJ72', 'IB22','JA43','JC54','LI33', 'LI58', 'LV17','ME52',...
%    'MI15', 'MY19', 'NZ20', 'NZ31', 'PD42','PE10', 'RI76', 'RM48',...
%    'SI53', 'SP51', 'TH35', 'TL08', 'TL55', 'UC44', 'VG09', 'VM29',...
%    'XP06'};
sitelist = {'ME52'};
compare_mat = [];
mat6 = [];
mat7 = [];
for qqqq = 1:length(sitelist)
    station_now = sitelist{qqqq};
datapath = strcat('C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Processed_data2\', station_now);
cd(datapath)
eventlist = dir;
eventlist = eventlist(3:length(eventlist));
%     event_num{i,1} = station;
%     event_num{i,2} = length(eventlist);
lat = zeros(length(eventlist), 1);
lon = zeros(length(eventlist), 1);
az = zeros(length(eventlist), 1);
mag = zeros(length(eventlist), 1);
ID = {};
loc = zeros(length(eventlist), 1);
for j = 1: length(eventlist)
    filename = eventlist(j);
    filename = strcat(filename.folder, '\', filename.name);
    load(filename)
    lat(j) = data.meta.event.lat;
    lon(j) = data.meta.event.lon;
    az(j) = data.meta.event.azimuth;
    mag(j) = data.meta.event.mag;
    ID{j} = data.meta.event.ID;
    station3{j} = data.meta.station.name;
    loc(j) = j;
    clear data
    disp(j)
end

%% 

mat = [];
mat = horzcat(lat,lon,mag,az, loc);

%%
max_az = max(az);
min_az = min(az);

mat = mat(mat(:,4) >= 80, :);
mat = mat(mat(:,4) <= 270, :);
mat6 = vertcat(mat6,mat);
figure
histogram(mat(:,4))

%% 
a = (max(mat(:,4)) - min(mat(:,4))) / 3;
cut1 = min(mat(:,4));
cut2 = cut1 + a;
cut3 = cut2 + a;
cut4 = cut3 + a;
compare_mat(qqqq,31) = cut1;
compare_mat(qqqq,32) = cut2;
compare_mat(qqqq,33) = cut3;
compare_mat(qqqq,34) = cut4;
counter1 = 0;
counter2 = 0;
counter3 = 0;
az1HV = [];
az1HVEW = [];
az2HV = [];
az2HVEW = [];
az3HV = [];
az3HVEW = [];
HV_all = [];
HVEW_all = [];
az1 = [];
az2 = [];
az3 = [];

for i = 1:length(mat)
    filename = eventlist(i);
    filename = strcat(filename.folder, '\', filename.name);
    load(filename)
    HV = data.processing.filtereddata.acceleration.NS.HVSR.smooth.HV;
        if length(HV) == 100000
            HV = HV(1:50000);
        end
    HV_all(i,:) = HV;
    HV = data.processing.filtereddata.acceleration.EW.HVSR.smooth.HV;
        if length(HV) == 100000
            HV = HV(1:50000);
        end
    HVNS_all(i,:) = HV;
    if mat(i,4) < cut2 && mat(i,4) >= cut1
        counter1 = counter1+1;
        az1(counter1,:) = mat(i,:);
        idnow = mat(i,5);
        filename = eventlist(idnow);
        filename = strcat(filename.folder, '\', filename.name);
        load(filename)
        HV = data.processing.filtereddata.acceleration.NS.HVSR.smooth.HV;
        if length(HV) == 100000
            HV = HV(1:50000);
        end
        az1HV(counter1,:) = HV;
        HV = data.processing.filtereddata.acceleration.EW.HVSR.smooth.HV;
        if length(HV) == 100000
            HV = HV(1:50000);
        end
        az1HVEW(counter1,:) = HV;

    end
    if mat(i,4) < cut3 && mat(i,4) >= cut2
        counter2 = counter2+1;
        az2(counter2,:) = mat(i,:);
        idnow = mat(i,5);
        filename = eventlist(idnow);
        filename = strcat(filename.folder, '\', filename.name);
        load(filename)
        HV = data.processing.filtereddata.acceleration.NS.HVSR.smooth.HV;
        if length(HV) == 100000
            HV = HV(1:50000);
        end
        az2HV(counter2,:) = HV;
        HV = data.processing.filtereddata.acceleration.EW.HVSR.smooth.HV;
        if length(HV) == 100000
            HV = HV(1:50000);
        end
        az2HVEW(counter2,:) = HV;
    end
    if mat(i,4) < cut4 && mat(i,4) >= cut3
        counter3 = counter3+1;
        az3(counter3,:) = mat(i,:);
        idnow = mat(i,5);
        filename = eventlist(idnow);
        filename = strcat(filename.folder, '\', filename.name);
        load(filename)
        HV = data.processing.filtereddata.acceleration.NS.HVSR.smooth.HV;
        if length(HV) == 100000
            HV = HV(1:50000);
        end
        az3HV(counter3,:) = HV;
        HV = data.processing.filtereddata.acceleration.EW.HVSR.smooth.HV;
        if length(HV) == 100000
            HV = HV(1:50000);
        end
        az3HVEW(counter3,:) = HV;
    end
end

statname = 'Southeast azimuth NS';
fax_HzN = data.processing.filtereddata.freq_vec;
upbound = 10;
lowbound = 0.1;
[~, lowbound] = min(abs(fax_HzN - lowbound));
[~, upbound] = min(abs(fax_HzN - upbound));

%% Now plat az1
clear xlim
clear ylim

cd(codepath)
[ahatf, sigma, confinthigh, confintlow] =  wavav(az1HV);
ahatf = ahatf(1:50000);
sigma = sigma(1:50000);
confinthigh = confinthigh(1:50000);
confintlow = confintlow(1:50000);

[peakfreqs, peakamps, hpb, f1s, f2s, Areamat, proms, amps, peakind2, freqs, sigs, I1s, I2s] = peakiden(ahatf, fax_HzN, sigma, lowbound, upbound);
[~, I] = max(amps);
compare_mat(qqqq,6) = freqs(I);
compare_mat(qqqq,7) = amps(I);
compare_mat(qqqq,8) = hpb(I);
compare_mat(qqqq,9) = proms(I);
compare_mat(qqqq,10) = sigs(I);

[~, freq_class] = sort(freqs, 'descend');
[~, amp_class] = sort(amps, 'descend');
[~, prom_class] = sort(proms, 'descend');
[~, hpb_class] = sort(hpb, 'descend');
[~, area_class] = sort(Areamat, 'descend');
[~, sig_class] = sort(sigs, 'descend');
classmatrix = vertcat(freq_class, amp_class, prom_class, hpb_class, area_class, sig_class);


datamat = vertcat(freqs,amps,proms,hpb,Areamat,sigs);


fin_fig1 = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

xlim([0.1 10])
ylim([0.1 100])
xticks([.1 1 10])
xticklabels({'0.1', '1', '10'})
yticks([1 10 100])
yticklabels({ '1','10', '100'})
ylabel('Amplification')
title(statname)
set(gca,'YScale', 'log','XScale','log', 'FontName', 'Times New Roman', 'FontSize', 24)

grid on
box on
hold on

cd(codepath)
[ahatf_all, sigma_all, confinthigh_all, confintlow_all] =  wavav(HV_all);
ahatf_all = ahatf_all(1:50000);
sigma_all = sigma_all(1:50000);
confinthigh_all = confinthigh_all(1:50000);
confintlow_all = confintlow_all(1:50000);

[peakfreqs_all, peakamps_all, hpb_all, f1s_all, f2s_all, Areamat_all, proms_all, amps_all, peakind2_all, freqs_all, sigs_all, I1s_all, I2s_all] = peakiden(ahatf_all, fax_HzN, sigma_all, lowbound, upbound);
[~, I] = max(amps_all);
compare_mat(qqqq,1) = freqs_all(I);
compare_mat(qqqq,2) = amps_all(I);
compare_mat(qqqq,3) = hpb_all(I);
compare_mat(qqqq,4) = proms_all(I);
compare_mat(qqqq,5) = sigs_all(I);
[~, I] = max(amps_all);
sig_all = sigs_all(I);
mat7(qqqq,1) = sig_all;


fr = fax_HzN(lowbound:length(fax_HzN))';
cohr = confinthigh_all(lowbound:length(fax_HzN));
colr = confintlow_all(lowbound:length(fax_HzN));
x_plot =[fr, fliplr(fr)];
y_plot = [cohr, fliplr(colr)];




statname = 'West azimuth EW';
fax_HzN = data.processing.filtereddata.freq_vec;
upbound = 10;
lowbound = 0.1;
[~, lowbound] = min(abs(fax_HzN - lowbound));
[~, upbound] = min(abs(fax_HzN - upbound));


cd(codepath)
[ahatf, sigma, confinthigh, confintlow] =  wavav(az3HV);
ahatf = ahatf(1:50000);
sigma = sigma(1:50000);
confinthigh = confinthigh(1:50000);
confintlow = confintlow(1:50000);

[peakfreqs, peakamps, hpb, f1s, f2s, Areamat, proms, amps, peakind2, freqs, sigs, I1s, I2s] = peakiden(ahatf, fax_HzN, sigma, lowbound, upbound);
[~, I] = max(amps);
compare_mat(qqqq,31) = freqs(I);
compare_mat(qqqq,32) = amps(I);
compare_mat(qqqq,33) = hpb(I);
compare_mat(qqqq,34) = proms(I);
compare_mat(qqqq,35) = sigs(I);

[~, freq_class] = sort(freqs, 'descend');
[~, amp_class] = sort(amps, 'descend');
[~, prom_class] = sort(proms, 'descend');
[~, hpb_class] = sort(hpb, 'descend');
[~, area_class] = sort(Areamat, 'descend');
[~, sig_class] = sort(sigs, 'descend');
classmatrix = vertcat(freq_class, amp_class, prom_class, hpb_class, area_class, sig_class);


datamat = vertcat(freqs,amps,proms,hpb,Areamat,sigs);


fin_fig3 = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

xlim([0.92 1.5])
ylim([2 30])
% xticks([.1 1 10])
% xticklabels({'0.1', '1', '10'})
yticks([1 10 100])
yticklabels({ '1','10', '100'})
ylabel('Amplification')
title(statname)
set(gca,'YScale', 'log','XScale','log', 'FontName', 'Times New Roman', 'FontSize', 24)

grid on
box on
hold on


fr = fax_HzN(lowbound:length(fax_HzN))';
cohr = confinthigh_all(lowbound:length(fax_HzN));
colr = confintlow_all(lowbound:length(fax_HzN));
x_plot =[fr, fliplr(fr)];
y_plot = [cohr, fliplr(colr)];

all_conf = fill(x_plot, y_plot, 1,'facecolor', [0 0.5 0], 'edgecolor', 'none', 'facealpha', 0.2);

% confidenceinterval=shadedplot(fr, cohr, colr,[.9,.9,.9],[1 1 1], 'FaceAlpha',0.5);
hold on

all_curve = plot(fax_HzN,ahatf_all, 'LineWidth', 2, 'Color', [0 0.5 0]);

fr = fax_HzN(lowbound:length(fax_HzN))';
cohr = confinthigh(lowbound:length(fax_HzN));
colr = confintlow(lowbound:length(fax_HzN));
x_plot =[fr, fliplr(fr)];
y_plot = [cohr, fliplr(colr)];

az_conf = fill(x_plot, y_plot, 1,'facecolor', [0 0 0.5], 'edgecolor', 'none', 'facealpha', 0.4);

% confidenceinterval=shadedplot(fr, cohr, colr,[0 0 0.5],[1 1 1], 'FaceAlpha',0.5);
hold on
%     








az_curve = plot(fax_HzN,ahatf, 'LineWidth', 2, 'Color', [0 0 0.5]);




if length(amps) > 0
    if length(amps) ==1
        f = strcat("Azimuth \sigma =" + " "  +num2str(sigs,3));
        mat7(qqqq,7) = sigs;
        g = strcat("All event \sigma =" + " "  +num2str(sig_all,3));
        str = {f,g};
        ylim=get(gca,'ylim');
        xlim=get(gca,'xlim');
        q = find(ahatf == amps);
        if fax_HzN(q) < 1
             text(xlim(2)-6,ylim(2)-60,str, 'FontName', 'Times New Roman', 'FontSize', 24,'Color', 'black', 'HorizontalAlignment', 'right', 'EdgeColor','k','BackgroundColor', 'w')
        else
             text(xlim(1)+0.3,ylim(2)-60,str, 'FontName', 'Times New Roman', 'FontSize', 24,'Color', 'black', 'HorizontalAlignment', 'right', 'EdgeColor','k','BackgroundColor', 'w')
        end
        else
            [~, I] = max(amps);
            f = strcat("Azimuth \sigma =" + " "  +num2str(sigs(I),3));
            mat7(qqqq,7) = sigs(I);
            g = strcat("All event \sigma =" + " "  +num2str(sig_all,3));
            str = {f, g};
            ylim=get(gca,'ylim');
            xlim=get(gca,'xlim');
            q = find(ahatf == amps(I));
            if fax_HzN(q) < 1
                text(xlim(2)-6,ylim(2)-60,str, 'FontName', 'Times New Roman', 'FontSize', 24,'Color', 'black', 'HorizontalAlignment', 'right', 'EdgeColor','k','BackgroundColor', 'w')
            else
                text(xlim(1)+0.3,ylim(2)-60,str, 'FontName', 'Times New Roman', 'FontSize', 24,'Color', 'black', 'HorizontalAlignment', 'right', 'EdgeColor','k','BackgroundColor', 'w')
            end

    end
    else
       str = 'No peaks';
       ext(0.3,15,str, 'FontName', 'Times New Roman', 'FontSize', 8,'Color', 'black', 'HorizontalAlignment', 'right', 'EdgeColor','k','BackgroundColor', 'w')
    
end
legend([all_curve, all_conf, az_curve, az_conf], 'HV all', 'HV all_{conf}', 'HV az', 'Hv az_{conf}')

% clear xlim
% clear ylim

ahatfaz3 = ahatf;

% saveas(fin_fig1, strcat(figpath,station_now,'southeast_EW','.jpg'))
% saveas(fin_fig2, strcat(figpath,station_now,'southwest_EW','.jpg'))
% saveas(fin_fig3, strcat(figpath,station_now,'west_EW','.jpg'))


end