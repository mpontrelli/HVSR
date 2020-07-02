%% Geo_HVSR
% Perform HVSR for geometric means on Mexico City

%% process Mexico City RACM data using code Rosetta and put it into a .mat file
% for each event at each station in folder "Processed_data2"
close all
clear all
%% ACCESS THE DATA
% go into the data folder and get a list of stations
codepath = 'C:\Users\mpontr01\Desktop\HVSR\Codes';
datapath = 'C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\geomean\';
refpath = 'C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\geomean\TP13';
outpath = 'C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\geomean_SSR';
cd(datapath)
stationlist = { 'AE02','AL01','AO24', 'AP68','AR14','AU11','AU46','BA49','BL45','BO39',...
     'CA20', 'CA59', 'CB43','CC55', 'CE23','CE32','CH84','CI05','CJ03',...
     'CJ04','CO47', 'CO56', 'CU80', 'DM12','DR16', 'DX37','EO30','ES57','EX08','EX09','EX12','GA62',...
     'GC38','GR27','HA41','HJ72', 'IB22','JA43','JC54','LI33', 'LI58', 'LV17','ME52',...
    'MI15', 'MY19', 'NZ20', 'NZ31', 'PD42','PE10', 'RI76', 'RM48',...
    'SI53', 'SP51', 'TH35', 'TL08', 'TL55', 'UC44', 'VG09', 'VM29',...
    'XP06'};
cd(refpath)
reflist = dir;
reflist = reflist(3:length(reflist));
event_num = {};
for i = 1 : length(stationlist)
    SSR_mat = [];
    disp(i)
    station = stationlist{i};
    cd(strcat(datapath, station))
    eventlist = dir;
    eventlist = eventlist(3:length(eventlist));
    tic
    counter = 0;
    for j = 1: length(eventlist)
        filename = eventlist(j);
        event_name = filename.name;
        ID = event_name(5:18);
        
        for jj = 1:length(reflist)
            refID = reflist(jj);
            refID1 = refID.name;
            refID = refID1(5:18);
            if strcmp(ID,refID) == 1
                counter = counter +1;
                load(strcat(datapath,station,'\',event_name))
                numerator = dat.geomean;
                load(strcat(refpath,'\',refID1))
                denominator = dat.geomean;
                SR = numerator./denominator;
                SSR_mat(counter,:) = SR;
            end
        end
    end
    cd(codepath)
    [ahatf, sigma, confinthigh, confintlow] =  wavav(SSR_mat);
    data.ahatf = ahatf;
    data.sigma = sigma;
    data.confinthigh = confinthigh;
    data.confintlow = confintlow;
    data.num = counter;
    fname = strcat(outpath,'\',station);
    save(fname, 'data')
    
end