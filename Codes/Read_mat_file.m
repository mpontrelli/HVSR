%% process .mat file created in creat_mat_file and compute average HVSR and Thompson statistics

close all
clear all
plotcolor = 'k';
lowbound = 100;
%% ACCESS THE DATA
% go into the data folder and get a list of stations
codepath = 'C:\Users\mpontr01\Desktop\HVSR\Codes';
datapath = 'C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Processed_data2';
figpath = 'C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Figures\';
shapepath = 'C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Shape_statistics\';

cd(datapath)
stationlist = dir;
stationlist = stationlist(3:length(stationlist));
event_num = {};

for i = 1: length(stationlist)
    HV_EW_mat = [];
    HV_NS_mat = [];
    HV_comp_mat = [];
    station = stationlist(i).name;
    cd(strcat(datapath, '\', station))
    eventlist = dir;
    eventlist = eventlist(3:length(eventlist));
%     event_num{i,1} = station;
%     event_num{i,2} = length(eventlist);
    for j = 1: length(eventlist)
        disp(j)
        filename = eventlist(j);
        filename = strcat(filename.folder, '\', filename.name);
        load(filename)
        cd(codepath)
        if length(data.processing.filtereddata.acceleration.EW.HVSR.smooth.HV) == 50000
            HV_EW_mat(j,:) = data.processing.filtereddata.acceleration.EW.HVSR.smooth.HV;
            HV_NS_mat(j,:) = data.processing.filtereddata.acceleration.NS.HVSR.smooth.HV;
            HV_comp_mat(j,:) = data.processing.filtereddata.acceleration.complex.HVSR.smooth.HV;
        else
            hv_short = data.processing.filtereddata.acceleration.EW.HVSR.smooth.HV;
            hv_shortEW = hv_short(1:50000);
            HV_EW_mat(j,:) = hv_shortEW;
            hv_short = data.processing.filtereddata.acceleration.NS.HVSR.smooth.HV;
            hv_shortNS = hv_short(1:50000);
            HV_NS_mat(j,:) = hv_shortNS;
            hv_short = data.processing.filtereddata.acceleration.complex.HVSR.smooth.HV;
            hv_shortcomp = hv_short(1:50000);
            HV_comp_mat(j,:) = hv_shortcomp;
            clear data
        end
        
    end
    load(filename)
    statname = data.meta.station.name;
    if length(data.processing.filtereddata.freq_vec) == 100000
        freq_short = data.processing.filtereddata.freq_vec;
        freq = freq_short(1:50000);
    else
        freq = data.processing.filtereddata.freq_vec;
    end
    % remove infinite rows
    HV_EW_mat(any(isinf(HV_EW_mat),2),:) = [];
    HV_NS_mat(any(isinf(HV_NS_mat),2),:) = [];
    HV_comp_mat(any(isinf(HV_comp_mat),2),:) = [];

    [ahatfEW, sigmaEW, confinthighEW, confintlowEW] =  wavav(HV_EW_mat);
    shapedata.EW.ahatf = ahatfEW;
    shapedata.EW.sigma = sigmaEW;
    shapedata.EW.confinthigh = confinthighEW;
    shapedata.EW.confintlow = confintlowEW;
    [ahatfNS, sigmaNS, confinthighNS, confintlowNS] =  wavav(HV_NS_mat);
    shapedata.NS.ahatf = ahatfNS;
    shapedata.NS.sigma = sigmaNS;
    shapedata.NS.confinthigh = confinthighNS;
    shapedata.NS.confintlow = confintlowNS;
    [ahatfcomp, sigmacomp, confinthighcomp, confintlowcomp] =  wavav(HV_comp_mat);
    shapedata.complex.ahatf = ahatfcomp;
    shapedata.complex.sigma = sigmacomp;
    shapedata.complex.confinthigh = confinthighcomp;
    shapedata.complex.confintlow = confintlowcomp;
    %% EW
    title = strcat(statname, {' '}, 'EW');
    EW = HVSRplot(ahatfEW, freq, confinthighEW, confintlowEW, lowbound, title, plotcolor);
    individplot(HV_EW_mat, freq, statname)
    upbound = 20000;
    [matrix, matrix1, peakind,ahatf1,newfaxhz1] = peakiden(ahatfEW, freq', lowbound, upbound);
    [taxstatEW, sigma1EW] = specratstat(peakind, matrix, matrix1, ahatf1, newfaxhz1, sigmaEW, statname,lowbound, upbound);
    saveas(EW, strcat(figpath, statname, '\', 'EWHVSR.jpg'));
    shapedata.EW.shapes = taxstatEW;
    shapedata.EW.sig_vec = sigma1EW;
    %% NS
    title = strcat(statname, {' '}, 'NS');
    NS = HVSRplot(ahatfNS, freq, confinthighNS, confintlowNS, lowbound, title, plotcolor);
    individplot(HV_NS_mat, freq, statname)
    upbound = 20000;
    [matrix, matrix1, peakind,ahatf1,newfaxhz1] = peakiden(ahatfNS, freq', lowbound, upbound);
    [taxstatNS, sigma1NS] = specratstat(peakind, matrix, matrix1, ahatf1, newfaxhz1, sigmaNS, statname,lowbound, upbound);
    saveas(NS, strcat(figpath, statname, '\', 'NSHVSR.jpg'));
    shapedata.NS.shapes = taxstatNS;
    shapedata.NS.sig_vec = sigma1NS;
    %% complex
    title = strcat(statname, {' '}, 'Complex');
    comp = HVSRplot(ahatfcomp, freq, confinthighcomp, confintlowcomp, lowbound, title, plotcolor);
    individplot(HV_comp_mat, freq, statname)
    upbound = 20000;
    [matrix, matrix1, peakind,ahatf1,newfaxhz1] = peakiden(ahatfcomp, freq', lowbound, upbound);
    [taxstatcomp, sigma1comp] = specratstat(peakind, matrix, matrix1, ahatf1, newfaxhz1, sigmacomp, statname,lowbound, upbound);
    saveas(comp, strcat(figpath, statname, '\', 'compHVSR.jpg'));
    shapedata.complex.shapes = taxstatcomp;
    shapedata.complex.sig_vec = sigma1comp;
    shapedata.sigma_freq = newfaxhz1;
    %% now save the file
    save(strcat(shapepath, statname, '.mat'), 'shapedata')
    %close all
end

