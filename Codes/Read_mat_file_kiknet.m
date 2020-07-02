%% read_mat_file_kiknet
% process .mat file created in creat_mat_file and compute average HVSR and Thompson statistics

close all
clear all
plotcolor = 'k';
lowbound = 100;
%% ACCESS THE DATA
% go into the data folder and get a list of stations
codepath = 'C:\Users\mpontr01\Desktop\HVSR\Codes';
datapath = 'C:\Users\mpontr01\Box\Projects\Kik-Net\Processed_data';
figpath = 'C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Figures\';
shapepath = 'C:\Users\mpontr01\Box\Projects\Kik-Net\Shape_data\';

cd(datapath)
stationlist = dir;
stationlist = stationlist(3:length(stationlist));
event_num = {};
%%
for i = 1: length(stationlist)
    HV_EW1 = [];
    HV_EW2 = [];
    HV_NS1 = [];
    HV_NS2 = [];    
    HV_comp1 = [];
    HV_comp2 = [];
    HV_geo1 = [];
    HV_geo2 = [];
    HV_rot1 = [];
    HV_rot2 = [];
    BR_NS = [];
    BR_EW = [];
    BR_comp = [];
    BR_geo = [];
    BR_rot = [];

    station = stationlist(i).name;
    cd(strcat(datapath, '\', station))
    eventlist = dir;
    eventlist = eventlist(3:length(eventlist));
    for j = 1:100% length(eventlist)
        disp(j)
        filename = eventlist(j);
        filename = strcat(filename.folder, '\', filename.name);
        load(filename)
        cd(codepath)

        HV_EW1(j,:) = data.processing.filtereddata.acceleration.EW1.HVSR.smooth.HV;
        HV_EW2(j,:) = data.processing.filtereddata.acceleration.EW2.HVSR.smooth.HV;
        HV_NS1(j,:) = data.processing.filtereddata.acceleration.NS1.HVSR.smooth.HV;
        HV_NS2(j,:) = data.processing.filtereddata.acceleration.NS2.HVSR.smooth.HV;
        HV_comp1(j,:) = data.processing.filtereddata.acceleration.complex1.HVSR.smooth.HV;
        HV_comp2(j,:) = data.processing.filtereddata.acceleration.complex2.HVSR.smooth.HV;
        HV_geo1(j,:) = data.processing.filtereddata.acceleration.geo_mean1.HVSR.smooth.HV;
        HV_geo2(j,:) = data.processing.filtereddata.acceleration.geo_mean2.HVSR.smooth.HV;
        BR_NS(j,:) = data.processing.filtereddata.acceleration.BSR.NS;
        BR_EW(j,:) = data.processing.filtereddata.acceleration.BSR.EW;
        BR_comp(j,:) = data.processing.filtereddata.acceleration.BSR.complex;
        BR_geo(j,:) = data.processing.filtereddata.acceleration.BSR.geo_mean;
        BR_rot(j,:) = data.processing.filtereddata.acceleration.BSR.rot;        
        clear data
    end
    
    %%
    load(filename)
    statname = data.meta.station.name;
    freq = data.processing.filtereddata.freq_vec;
    % Now some inputs
    upbound = 20;
    lowbound = 0.1;
    [~, lowbound] = min(abs(freq - lowbound));
    [~, upbound] = min(abs(freq - upbound));



    %% EW1
    [ahatfEW1, sigmaEW1, confinthighEW1, confintlowEW1] =  wavav(HV_EW1);
    [peakfreqs, peakamps, hpb, f1s, f2s, Areamat, proms, amps, peakind2, freqs, sigs, I1s, I2s] = peakiden(ahatfEW1, freq, sigmaEW1, lowbound, upbound);
    shapedata.HV.EW1.shapes.peak_freqs = peakfreqs;
    shapedata.HV.EW1.shapes.peak_amps = peakamps;
    shapedata.HV.EW1.shapes.hpb = hpb;
    shapedata.HV.EW1.shapes.f1s = f1s;
    shapedata.HV.EW1.shapes.f2s = f2s;
    shapedata.HV.EW1.shapes.Area = Areamat;
    shapedata.HV.EW1.shapes.proms = proms;
    shapedata.HV.EW1.shapes.amps = amps;
    shapedata.HV.EW1.shapes.peakind = peakind2;
    shapedata.HV.EW1.shapes.freqs = freqs;
    shapedata.HV.EW1.shapes.sigs = sigs;
    shapedata.HV.EW1.shapes.I1s = I1s;
    shapedata.HV.EW1.shapes.I2s = I2s;
    shapedata.HV.EW1.ahatf = ahatfEW1;
    shapedata.HV.EW1.sigma = sigmaEW1;
    shapedata.HV.EW1.confinthigh = confinthighEW1;
    shapedata.HV.EW1.confintlow = confintlowEW1;
    
    %% EW2
    [ahatfEW2, sigmaEW2, confinthighEW2, confintlowEW2] =  wavav(HV_EW2);
    [peakfreqs, peakamps, hpb, f1s, f2s, Areamat, proms, amps, peakind2, freqs, sigs, I1s, I2s] = peakiden(ahatfEW2, freq, sigmaEW2, lowbound, upbound);
    shapedata.HV.EW2.shapes.peak_freqs = peakfreqs;
    shapedata.HV.EW2.shapes.peak_amps = peakamps;
    shapedata.HV.EW2.shapes.hpb = hpb;
    shapedata.HV.EW2.shapes.f1s = f1s;
    shapedata.HV.EW2.shapes.f2s = f2s;
    shapedata.HV.EW2.shapes.Area = Areamat;
    shapedata.HV.EW2.shapes.proms = proms;
    shapedata.HV.EW2.shapes.amps = amps;
    shapedata.HV.EW2.shapes.peakind = peakind2;
    shapedata.HV.EW2.shapes.freqs = freqs;
    shapedata.HV.EW2.shapes.sigs = sigs;
    shapedata.HV.EW2.shapes.I1s = I1s;
    shapedata.HV.EW2.shapes.I2s = I2s;
    shapedata.HV.EW2.ahatf = ahatfEW2;
    shapedata.HV.EW2.sigma = sigmaEW2;
    shapedata.HV.EW2.confinthigh = confinthighEW2;
    shapedata.HV.EW2.confintlow = confintlowEW2;
    
    %% NS1
    [ahatfNS1, sigmaNS1, confinthighNS1, confintlowNS1] =  wavav(HV_NS1);
    [peakfreqs, peakamps, hpb, f1s, f2s, Areamat, proms, amps, peakind2, freqs, sigs, I1s, I2s] = peakiden(ahatfNS1, freq, sigmaNS1, lowbound, upbound);
    shapedata.HV.NS1.shapes.peak_freqs = peakfreqs;
    shapedata.HV.NS1.shapes.peak_amps = peakamps;
    shapedata.HV.NS1.shapes.hpb = hpb;
    shapedata.HV.NS1.shapes.f1s = f1s;
    shapedata.HV.NS1.shapes.f2s = f2s;
    shapedata.HV.NS1.shapes.Area = Areamat;
    shapedata.HV.NS1.shapes.proms = proms;
    shapedata.HV.NS1.shapes.amps = amps;
    shapedata.HV.NS1.shapes.peakind = peakind2;
    shapedata.HV.NS1.shapes.freqs = freqs;
    shapedata.HV.NS1.shapes.sigs = sigs;
    shapedata.HV.NS1.shapes.I1s = I1s;
    shapedata.HV.NS1.shapes.I2s = I2s;
    shapedata.HV.NS1.ahatf = ahatfNS1;
    shapedata.HV.NS1.sigma = sigmaNS1;
    shapedata.HV.NS1.confinthigh = confinthighNS1;
    shapedata.HV.NS1.confintlow = confintlowNS1;

    %% NS2
    [ahatfNS2, sigmaNS2, confinthighNS2, confintlowNS2] =  wavav(HV_NS2);
    [peakfreqs, peakamps, hpb, f1s, f2s, Areamat, proms, amps, peakind2, freqs, sigs, I1s, I2s] = peakiden(ahatfNS2, freq, sigmaNS2, lowbound, upbound);
    shapedata.HV.NS2.shapes.peak_freqs = peakfreqs;
    shapedata.HV.NS2.shapes.peak_amps = peakamps;
    shapedata.HV.NS2.shapes.hpb = hpb;
    shapedata.HV.NS2.shapes.f1s = f1s;
    shapedata.HV.NS2.shapes.f2s = f2s;
    shapedata.HV.NS2.shapes.Area = Areamat;
    shapedata.HV.NS2.shapes.proms = proms;
    shapedata.HV.NS2.shapes.amps = amps;
    shapedata.HV.NS2.shapes.peakind = peakind2;
    shapedata.HV.NS2.shapes.freqs = freqs;
    shapedata.HV.NS2.shapes.sigs = sigs;
    shapedata.HV.NS2.shapes.I1s = I1s;
    shapedata.HV.NS2.shapes.I2s = I2s;
    shapedata.HV.NS2.ahatf = ahatfNS2;
    shapedata.HV.NS2.sigma = sigmaNS2;
    shapedata.HV.NS2.confinthigh = confinthighNS2;
    shapedata.HV.NS2.confintlow = confintlowNS2;
    
    %% Complex 1
    [ahatfcomp1, sigmacomp1, confinthighcomp1, confintlowcomp1] =  wavav(HV_comp1);
    [peakfreqs, peakamps, hpb, f1s, f2s, Areamat, proms, amps, peakind2, freqs, sigs, I1s, I2s] = peakiden(ahatfcomp1, freq, sigmacomp1, lowbound, upbound);
    shapedata.HV.comp1.shapes.peak_freqs = peakfreqs;
    shapedata.HV.comp1.shapes.peak_amps = peakamps;
    shapedata.HV.comp1.shapes.hpb = hpb;
    shapedata.HV.comp1.shapes.f1s = f1s;
    shapedata.HV.comp1.shapes.f2s = f2s;
    shapedata.HV.comp1.shapes.Area = Areamat;
    shapedata.HV.comp1.shapes.proms = proms;
    shapedata.HV.comp1.shapes.amps = amps;
    shapedata.HV.comp1.shapes.peakind = peakind2;
    shapedata.HV.comp1.shapes.freqs = freqs;
    shapedata.HV.comp1.shapes.sigs = sigs;
    shapedata.HV.comp1.shapes.I1s = I1s;
    shapedata.HV.comp1.shapes.I2s = I2s;
    shapedata.HV.comp1.ahatf = ahatfcomp1;
    shapedata.HV.comp1.sigma = sigmacomp1;
    shapedata.HV.comp1.confinthigh = confinthighcomp1;
    shapedata.HV.comp1.confintlow = confintlowcomp1;    
    
    %% Complex2
    [ahatfcomp2, sigmacomp2, confinthighcomp2, confintlowcomp2] =  wavav(HV_comp2);
    [peakfreqs, peakamps, hpb, f1s, f2s, Areamat, proms, amps, peakind2, freqs, sigs, I1s, I2s] = peakiden(ahatfcomp2, freq, sigmacomp2, lowbound, upbound);
    shapedata.HV.comp2.shapes.peak_freqs = peakfreqs;
    shapedata.HV.comp2.shapes.peak_amps = peakamps;
    shapedata.HV.comp2.shapes.hpb = hpb;
    shapedata.HV.comp2.shapes.f1s = f1s;
    shapedata.HV.comp2.shapes.f2s = f2s;
    shapedata.HV.comp2.shapes.Area = Areamat;
    shapedata.HV.comp2.shapes.proms = proms;
    shapedata.HV.comp2.shapes.amps = amps;
    shapedata.HV.comp2.shapes.peakind = peakind2;
    shapedata.HV.comp2.shapes.freqs = freqs;
    shapedata.HV.comp2.shapes.sigs = sigs;
    shapedata.HV.comp2.shapes.I1s = I1s;
    shapedata.HV.comp2.shapes.I2s = I2s;
    shapedata.HV.comp2.ahatf = ahatfcomp2;
    shapedata.HV.comp2.sigma = sigmacomp2;
    shapedata.HV.comp2.confinthigh = confinthighcomp2;
    shapedata.HV.comp2.confintlow = confintlowcomp2;   
    
    %% Geo1
    [ahatfgeo1, sigmageo1, confinthighgeo1, confintlowgeo1] =  wavav(HV_geo1);
    [peakfreqs, peakamps, hpb, f1s, f2s, Areamat, proms, amps, peakind2, freqs, sigs, I1s, I2s] = peakiden(ahatfgeo1, freq, sigmageo1, lowbound, upbound);
    shapedata.HV.geo1.shapes.peak_freqs = peakfreqs;
    shapedata.HV.geo1.shapes.peak_amps = peakamps;
    shapedata.HV.geo1.shapes.hpb = hpb;
    shapedata.HV.geo1.shapes.f1s = f1s;
    shapedata.HV.geo1.shapes.f2s = f2s;
    shapedata.HV.geo1.shapes.Area = Areamat;
    shapedata.HV.geo1.shapes.proms = proms;
    shapedata.HV.geo1.shapes.amps = amps;
    shapedata.HV.geo1.shapes.peakind = peakind2;
    shapedata.HV.geo1.shapes.freqs = freqs;
    shapedata.HV.geo1.shapes.sigs = sigs;
    shapedata.HV.geo1.shapes.I1s = I1s;
    shapedata.HV.geo1.shapes.I2s = I2s;
    shapedata.HV.geo1.ahatf = ahatfgeo1;
    shapedata.HV.geo1.sigma = sigmageo1;
    shapedata.HV.geo1.confinthigh = confinthighgeo1;
    shapedata.HV.geo1.confintlow = confintlowgeo1;   
    
    %% Geo2
    [ahatfgeo2, sigmageo2, confinthighgeo2, confintlowgeo2] =  wavav(HV_geo2);
    [peakfreqs, peakamps, hpb, f1s, f2s, Areamat, proms, amps, peakind2, freqs, sigs, I1s, I2s] = peakiden(ahatfgeo2, freq, sigmageo2, lowbound, upbound);
    shapedata.HV.geo2.shapes.peak_freqs = peakfreqs;
    shapedata.HV.geo2.shapes.peak_amps = peakamps;
    shapedata.HV.geo2.shapes.hpb = hpb;
    shapedata.HV.geo2.shapes.f1s = f1s;
    shapedata.HV.geo2.shapes.f2s = f2s;
    shapedata.HV.geo2.shapes.Area = Areamat;
    shapedata.HV.geo2.shapes.proms = proms;
    shapedata.HV.geo2.shapes.amps = amps;
    shapedata.HV.geo2.shapes.peakind = peakind2;
    shapedata.HV.geo2.shapes.freqs = freqs;
    shapedata.HV.geo2.shapes.sigs = sigs;
    shapedata.HV.geo2.shapes.I1s = I1s;
    shapedata.HV.geo2.shapes.I2s = I2s;
    shapedata.HV.geo2.ahatf = ahatfgeo2;
    shapedata.HV.geo2.sigma = sigmageo2;
    shapedata.HV.geo2.confinthigh = confinthighgeo2;
    shapedata.HV.geo2.confintlow = confintlowgeo2; 
    
    %% Now BSRs
    %% NS
    [ahatfNS, sigmaNS, confinthighNS, confintlowNS] =  wavav(BR_NS);
    [peakfreqs, peakamps, hpb, f1s, f2s, Areamat, proms, amps, peakind2, freqs, sigs, I1s, I2s] = peakiden(ahatfNS, freq, sigmaNS, lowbound, upbound);
    shapedata.BSR.NS.shapes.peak_freqs = peakfreqs;
    shapedata.BSR.NS.shapes.peak_amps = peakamps;
    shapedata.BSR.NS.shapes.hpb = hpb;
    shapedata.BSR.NS.shapes.f1s = f1s;
    shapedata.BSR.NS.shapes.f2s = f2s;
    shapedata.BSR.NS.shapes.Area = Areamat;
    shapedata.BSR.NS.shapes.proms = proms;
    shapedata.BSR.NS.shapes.amps = amps;
    shapedata.BSR.NS.shapes.peakind = peakind2;
    shapedata.BSR.NS.shapes.freqs = freqs;
    shapedata.BSR.NS.shapes.sigs = sigs;
    shapedata.BSR.NS.shapes.I1s = I1s;
    shapedata.BSR.NS.shapes.I2s = I2s;
    shapedata.BSR.NS.ahatf = ahatfNS;
    shapedata.BSR.NS.sigma = sigmaNS;
    shapedata.BSR.NS.confinthigh = confinthighNS;
    shapedata.BSR.NS.confintlow = confintlowNS;     
    
    %% EW
    [ahatfEW, sigmaEW, confinthighEW, confintlowEW] =  wavav(BR_EW);
    [peakfreqs, peakamps, hpb, f1s, f2s, Areamat, proms, amps, peakind2, freqs, sigs, I1s, I2s] = peakiden(ahatfEW, freq, sigmaEW, lowbound, upbound);
    shapedata.BSR.EW.shapes.peak_freqs = peakfreqs;
    shapedata.BSR.EW.shapes.peak_amps = peakamps;
    shapedata.BSR.EW.shapes.hpb = hpb;
    shapedata.BSR.EW.shapes.f1s = f1s;
    shapedata.BSR.EW.shapes.f2s = f2s;
    shapedata.BSR.EW.shapes.Area = Areamat;
    shapedata.BSR.EW.shapes.proms = proms;
    shapedata.BSR.EW.shapes.amps = amps;
    shapedata.BSR.EW.shapes.peakind = peakind2;
    shapedata.BSR.EW.shapes.freqs = freqs;
    shapedata.BSR.EW.shapes.sigs = sigs;
    shapedata.BSR.EW.shapes.I1s = I1s;
    shapedata.BSR.EW.shapes.I2s = I2s;
    shapedata.BSR.EW.ahatf = ahatfEW;
    shapedata.BSR.EW.sigma = sigmaEW;
    shapedata.BSR.EW.confinthigh = confinthighEW;
    shapedata.BSR.EW.confintlow = confintlowEW;   
    
    %% complex
    [ahatfcomp, sigmacomp, confinthighcomp, confintlowcomp] =  wavav(BR_comp);
    [peakfreqs, peakamps, hpb, f1s, f2s, Areamat, proms, amps, peakind2, freqs, sigs, I1s, I2s] = peakiden(ahatfcomp, freq, sigmacomp, lowbound, upbound);
    shapedata.BSR.comp.shapes.peak_freqs = peakfreqs;
    shapedata.BSR.comp.shapes.peak_amps = peakamps;
    shapedata.BSR.comp.shapes.hpb = hpb;
    shapedata.BSR.comp.shapes.f1s = f1s;
    shapedata.BSR.comp.shapes.f2s = f2s;
    shapedata.BSR.comp.shapes.Area = Areamat;
    shapedata.BSR.comp.shapes.proms = proms;
    shapedata.BSR.comp.shapes.amps = amps;
    shapedata.BSR.comp.shapes.peakind = peakind2;
    shapedata.BSR.comp.shapes.freqs = freqs;
    shapedata.BSR.comp.shapes.sigs = sigs;
    shapedata.BSR.comp.shapes.I1s = I1s;
    shapedata.BSR.comp.shapes.I2s = I2s;
    shapedata.BSR.comp.ahatf = ahatfcomp;
    shapedata.BSR.comp.sigma = sigmacomp;
    shapedata.BSR.comp.confinthigh = confinthighcomp;
    shapedata.BSR.comp.confintlow = confintlowcomp;   
    
    %% geo
    [ahatfgeo, sigmageo, confinthighgeo, confintlowgeo] =  wavav(BR_geo);
    [peakfreqs, peakamps, hpb, f1s, f2s, Areamat, proms, amps, peakind2, freqs, sigs, I1s, I2s] = peakiden(ahatfgeo, freq, sigmageo, lowbound, upbound);
    shapedata.BSR.geo.shapes.peak_freqs = peakfreqs;
    shapedata.BSR.geo.shapes.peak_amps = peakamps;
    shapedata.BSR.geo.shapes.hpb = hpb;
    shapedata.BSR.geo.shapes.f1s = f1s;
    shapedata.BSR.geo.shapes.f2s = f2s;
    shapedata.BSR.geo.shapes.Area = Areamat;
    shapedata.BSR.geo.shapes.proms = proms;
    shapedata.BSR.geo.shapes.amps = amps;
    shapedata.BSR.geo.shapes.peakind = peakind2;
    shapedata.BSR.geo.shapes.freqs = freqs;
    shapedata.BSR.geo.shapes.sigs = sigs;
    shapedata.BSR.geo.shapes.I1s = I1s;
    shapedata.BSR.geo.shapes.I2s = I2s;
    shapedata.BSR.geo.ahatf = ahatfgeo;
    shapedata.BSR.geo.sigma = sigmageo;
    shapedata.BSR.geo.confinthigh = confinthighgeo;
    shapedata.BSR.geo.confintlow = confintlowgeo;  
    
    %% Rotated
    [ahatfrot, sigmarot, confinthighrot, confintlowrot] =  wavav(BR_rot);
    [peakfreqs, peakamps, hpb, f1s, f2s, Areamat, proms, amps, peakind2, freqs, sigs, I1s, I2s] = peakiden(ahatfrot, freq, sigmarot, lowbound, upbound);
    shapedata.BSR.rot.shapes.peak_freqs = peakfreqs;
    shapedata.BSR.rot.shapes.peak_amps = peakamps;
    shapedata.BSR.rot.shapes.hpb = hpb;
    shapedata.BSR.rot.shapes.f1s = f1s;
    shapedata.BSR.rot.shapes.f2s = f2s;
    shapedata.BSR.rot.shapes.Area = Areamat;
    shapedata.BSR.rot.shapes.proms = proms;
    shapedata.BSR.rot.shapes.amps = amps;
    shapedata.BSR.rot.shapes.peakind = peakind2;
    shapedata.BSR.rot.shapes.freqs = freqs;
    shapedata.BSR.rot.shapes.sigs = sigs;
    shapedata.BSR.rot.shapes.I1s = I1s;
    shapedata.BSR.rot.shapes.I2s = I2s;
    shapedata.BSR.rot.ahatf = ahatfrot;
    shapedata.BSR.rot.sigma = sigmarot;
    shapedata.BSR.rot.confinthigh = confinthighrot;
    shapedata.BSR.rot.confintlow = confintlowrot;    
    %% now save the file
    save(strcat(shapepath, statname, '.mat'), 'shapedata')
end