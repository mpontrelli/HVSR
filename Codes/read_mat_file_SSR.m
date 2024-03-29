%% process .mat file created in creat_mat_file and compute average HVSR and Thompson statistics
% read_matfile_SSR
close all
clear all
plotcolor = 'k';
lowbound = 100;
%% ACCESS THE DATA
% go into the data folder and get a list of stations

sitelist = {'AE02','AL01','AO24', 'AP68','AR14','AU11','AU46','BA49','BL45','BO39',...
    'CA20', 'CA59', 'CB43','CC55','CE18', 'CE23','CE32','CH84','CI05','CJ03',...
    'CJ04','CO47', 'CO56','CP28','CS78','CT64', 'CU80', 'DM12','DR16', 'DX37','EO30','ES57','EX08','EX09','EX12','FJ74','GA62',...
    'GC38','GR27','HA41','HJ72', 'IB22','IM40','JA43','JC54','LI33', 'LI58', 'LV17','ME52',...
    'MI15', 'MY19', 'NZ20', 'NZ31','PA34', 'PD42','PE10', 'RI76',... 
    'RM48', 'SI53', 'SP51', 'TE07', 'TH35', 'TL08', 'TL55', 'TP13', 'UC44','UI21', 'VG09', 'VM29',...
    'XP06'};


reflist = {'TP13'};
codepath = 'C:\Users\Marshall\Desktop\HVSR\Codes';
outpath = 'C:\Users\Marshall\Box\Data\Ground motion\Mexico CIty\SSR_Shape_statistics\';
for a = 1:length(sitelist)
    currsite = sitelist{a};
    disp(currsite)
    datapath = strcat('C:\Users\Marshall\Box\Data\Ground motion\Mexico CIty\Processed_data2\', currsite, '\');
    figpath = strcat('C:\Users\Marshall\Box\Data\Ground motion\Mexico CIty\SSR_figures\',currsite);
    cd(datapath)
    eventlist = dir;
    eventlist = eventlist(3:length(eventlist));
    for b = 1:length(reflist)
        currref = reflist{b};
        if strcmp(currsite,currref) ==0
            datapath2 = strcat('C:\Users\Marshall\Box\Data\Ground motion\Mexico CIty\Processed_data2\',currref,'\');
            cd(datapath2)
            refeventlist = dir;
            refeventlist = refeventlist(3:length(refeventlist));
            counter = 0;
            SR_NS = [];
            SR_EW = [];
            SR_complex = [];
            for i = 1 : length(eventlist)
                event = eventlist(i).name;
                filename = strcat(datapath,event);
                soft = load(filename);
                NSmag = soft.data.processing.filtereddata.acceleration.NS.mag_resps.smooth;
                if length(NSmag) == 100000
                    NSmag = NSmag(1:50000);
                end
                EWmag = soft.data.processing.filtereddata.acceleration.EW.mag_resps.smooth;
                if length(EWmag) == 100000
                    EWmag = EWmag(1:50000);
                end
                compmag = soft.data.processing.filtereddata.acceleration.complex.mag_resps.smooth;
                if length(compmag) == 100000
                    compmag = compmag(1:50000);
                end
                eventid = event(5:18);
                for ii = 1:length(refeventlist)
                    refevent = refeventlist(ii).name;
                    refid = refevent(5:18);
                    if strcmp(eventid,refid) == 1
                        counter = counter + 1;
                        filename = strcat(datapath2,refevent);
                        hard = load(filename);
                        NSrefmag = hard.data.processing.filtereddata.acceleration.NS.mag_resps.smooth;
                        if length(NSrefmag) == 100000
                            NSrefmag = NSrefmag(1:50000);
                        end
                        EWrefmag = hard.data.processing.filtereddata.acceleration.EW.mag_resps.smooth;
                        if length(EWrefmag) == 100000
                            EWrefmag = EWrefmag(1:50000);
                        end
                        comprefmag = hard.data.processing.filtereddata.acceleration.complex.mag_resps.smooth;
                        if length(comprefmag) == 100000
                            comprefmag = comprefmag(1:50000);
                        end
                    SR_NS(counter,:) = NSmag./NSrefmag;
                    SR_EW(counter,:) = EWmag./EWrefmag;
                    SR_complex(counter,:) = compmag./comprefmag;
                    end
                end
      
            end
        statname =  strcat(currsite,'-',currref);
        if length(soft.data.processing.filtereddata.freq_vec) == 100000
            freq_short = soft.data.processing.filtereddata.freq_vec;
            freq = freq_short(1:50000);
        else
            freq = soft.data.processing.filtereddata.freq_vec;
        end
        % remove infinite rows
        SR_NS(any(isinf(SR_NS),2),:) = [];
        SR_EW(any(isinf(SR_EW),2),:) = [];
        SR_complex(any(isinf(SR_complex),2),:) = [];

         cd(codepath)
        [ahatfEW, sigmaEW, confinthighEW, confintlowEW] =  wavav(SR_EW);
        shapedata.EW.ahatf = ahatfEW;
        shapedata.EW.sigma = sigmaEW;
        shapedata.EW.confinthigh = confinthighEW;
        shapedata.EW.confintlow = confintlowEW;
        [ahatfNS, sigmaNS, confinthighNS, confintlowNS] =  wavav(SR_NS);
        shapedata.NS.ahatf = ahatfNS;
        shapedata.NS.sigma = sigmaNS;
        shapedata.NS.confinthigh = confinthighNS;
        shapedata.NS.confintlow = confintlowNS;
        [ahatfcomp, sigmacomp, confinthighcomp, confintlowcomp] =  wavav(SR_complex);
        shapedata.complex.ahatf = ahatfcomp;
        shapedata.complex.sigma = sigmacomp;
        shapedata.complex.confinthigh = confinthighcomp;
        shapedata.complex.confintlow = confintlowcomp;
%         %% EW
%         title = strcat(statname, {' '}, 'EW');
%         title = title{1};
%         EW = HVSRplot(ahatfEW, freq, confinthighEW, confintlowEW, lowbound, title, plotcolor);
%         individplot(SR_EW, freq, statname)
%         upbound = 20000;
%         [matrix, matrix1, peakind,ahatf1,newfaxhz1, peakfreqs, peakamps, Areamat] = peakiden(ahatfEW', freq, lowbound, upbound);
%         [taxstatEW, sigma1EW] = specratstat(peakind, matrix, matrix1, ahatf1, newfaxhz1, sigmaEW, statname,lowbound, upbound);
%         saveas(EW, strcat(figpath, '\',title,'.jpg'));
%         shapedata.EW.shapes = taxstatEW;
%         shapedata.EW.sig_vec = sigma1EW;
%         shapedata.EW.area_freq = peakfreqs;
%         shapedata.EW.area_amp = peakamps;
%         shapedata.EW.Areamat = Areamat;
%         %% NS
%         title = strcat(statname, {' '}, 'NS');
%         title = title{1};
%         NS = HVSRplot(ahatfNS, freq, confinthighNS, confintlowNS, lowbound, title, plotcolor);
%         individplot(SR_NS, freq, statname)
%         upbound = 20000;
%         [matrix, matrix1, peakind,ahatf1,newfaxhz1, peakfreqs, peakamps, Areamat] = peakiden(ahatfNS', freq, lowbound, upbound);
%         [taxstatNS, sigma1NS] = specratstat(peakind, matrix, matrix1, ahatf1, newfaxhz1, sigmaNS, statname,lowbound, upbound);
%         saveas(NS, strcat(figpath, '\',title,'.jpg'));
%         shapedata.NS.shapes = taxstatNS;
%         shapedata.NS.sig_vec = sigma1NS;
%         shapedata.NS.area_freq = peakfreqs;
%         shapedata.NS.area_amp = peakamps;
%         shapedata.NS.Areamat = Areamat;
%         %% complex
%         title = strcat(statname, {' '}, 'Complex');
%         title = title{1};
%         comp = HVSRplot(ahatfcomp, freq, confinthighcomp, confintlowcomp, lowbound, title, plotcolor);
%         individplot(SR_complex, freq, statname)
%         upbound = 20000;
%         [matrix, matrix1, peakind,ahatf1,newfaxhz1, peakfreqs, peakamps, Areamat] = peakiden(ahatfcomp', freq, lowbound, upbound);
%         [taxstatcomp, sigma1comp] = specratstat(peakind, matrix, matrix1, ahatf1, newfaxhz1, sigmacomp, statname,lowbound, upbound);
%         saveas(comp, strcat(figpath, '\',title,'.jpg'));
%         shapedata.complex.shapes = taxstatcomp;
%         shapedata.complex.sig_vec = sigma1comp;
%         shapedata.complex.area_freq = peakfreqs;
%         shapedata.complex.area_amp = peakamps;
%         shapedata.complex.Areamat = Areamat; 
%         shapedata.sigma_freq = newfaxhz1;
        
        %% now compute distance and azimuth between stations
        softlat = soft.data.meta.station.lat;
        softlon = soft.data.meta.station.lon;
        reflat = hard.data.meta.station.lat;
        reflon = hard.data.meta.station.lon;
        
        az = azimuth(softlat,softlon,reflat, reflon);
        shapedata.azimuth = az;
        
        dist = deg2km(distance(softlat,softlon,reflat, reflon));
        shapedata.dist = dist;
        
        shapedata.numstations = counter;
        %% now save the file
        fname = strcat(outpath,currsite,'\',statname, '.mat');
        parsave(fname, shapedata)
        clear SR_NS 
        clear SR_EW
        clear SR_complex
        close all   
        
        
        end
    end
    
end

