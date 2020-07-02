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
outpath_SE = 'C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\geomean_HVSR_SE';
outpath_SW = 'C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\geomean_HVSR_SW';
outpath_W = 'C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\geomean_HVSR_W';
cd(datapath)
stationlist = dir;
stationlist = stationlist(3:length(stationlist));
event_num = {};
for i = 1 : length(stationlist)
    HV_mat_SE = [];
    HV_mat_SW = [];
    HV_mat_W = [];
    counter1 = 0;
    counter2 = 0;
    counter3 = 0;
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
        az = dat.az;
        if az >=111 && az <= 164
            counter1 = counter1 + 1;
            HV_mat_SE(counter1,:) = dat.HV;
        end
        if az >164 && az <= 217
            counter2 = counter2 + 1;
            HV_mat_SW(counter2,:) = dat.HV;
        end
        if az > 217 && az <= 270
            counter3 = counter3 + 1;
            HV_mat_W(counter3,:) = dat.HV;
        end
    end
    cd(codepath)
    [ahatf, sigma, confinthigh, confintlow] =  wavav(HV_mat_SE);
    data.ahatf = ahatf;
    data.sigma = sigma;
    data.confinthigh = confinthigh;
    data.confintlow = confintlow;
    data.num = counter1;
    fname = strcat(outpath_SE,'\',station);
    save(fname, 'data')
    
    clear data
    [ahatf, sigma, confinthigh, confintlow] =  wavav(HV_mat_SW);
    data.ahatf = ahatf;
    data.sigma = sigma;
    data.confinthigh = confinthigh;
    data.confintlow = confintlow;
    data.num = counter2;
    fname = strcat(outpath_SW,'\',station);
    save(fname, 'data')
    
    clear data
    [ahatf, sigma, confinthigh, confintlow] =  wavav(HV_mat_W);
    data.ahatf = ahatf;
    data.sigma = sigma;
    data.confinthigh = confinthigh;
    data.confintlow = confintlow;
    data.num = counter3;
    fname = strcat(outpath_W,'\',station);
    save(fname, 'data')
    
    
end