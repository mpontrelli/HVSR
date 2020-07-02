%% Geo_HVSR
% Perform HVSR for geometric means on Mexico City

%% process Mexico City RACM data using code Rosetta and put it into a .mat file
% for each event at each station in folder "Processed_data2"
close all
clear all
%% ACCESS THE DATA
% go into the data folder and get a list of stations
codepath = 'C:\Users\mpontr01\Desktop\HVSR\Codes';
datapath = 'C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\geomean';
outpath = 'C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\geomean_HVSR';
cd(datapath)
stationlist = dir;
stationlist = stationlist(3:length(stationlist));
event_num = {};
for i = 20%1 : length(stationlist)
    HV_mat = [];
    disp(i)
    station = stationlist(i).name;
    cd(strcat(datapath, '\', station))
    eventlist = dir;
    eventlist = eventlist(3:length(eventlist));
    tic
    for j = 1: length(eventlist)
        filename = eventlist(j);
        event_name = filename.name;
        load(event_name)
        HV_mat(j,:) = dat.HV;

    end
%     cd(codepath)
%     [ahatf, sigma, confinthigh, confintlow] =  wavav(HV_mat);
%     data.ahatf = ahatf;
%     data.sigma = sigma;
%     data.confinthigh = confinthigh;
%     data.confintlow = confintlow;
%     fname = strcat(outpath,'\',station);
%     save(fname, 'data')
    
end

%% plot 
load('C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Processed_data2\AE02\AE0219900511234349')
filename = 'C:\Users\mpontr01\Box\Pontrelli_et_al_2020\Final_tables.xlsx';
newfaxhz = data.processing.filtereddata.freq_vec;
statname = 'CE32';
cd(codepath)
        fin_fig = figure;
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
plot(newfaxhz, HV_mat, 'Color',  'k' , 'Linewidth', .5);
