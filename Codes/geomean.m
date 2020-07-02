%% Compute geometric means

%% process Mexico City RACM data using code Rosetta and put it into a .mat file
% for each event at each station in folder "Processed_data2"
close all
clear all
%% ACCESS THE DATA
% go into the data folder and get a list of stations
codepath = 'C:\Users\mpontr01\Desktop\HVSR\Codes';
datapath = 'C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Processed_data2';
outpath = 'C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\geomean';
cd(datapath)
stationlist = dir;
stationlist = stationlist(3:length(stationlist));
event_num = {};
for i = 58 : length(stationlist)
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
        EW = data.processing.filtereddata.acceleration.EW.mag_resps.smooth;
        NS = data.processing.filtereddata.acceleration.NS.mag_resps.smooth;
        V = data.processing.filtereddata.acceleration.V.mag_resps.smooth;
        geo_mean = sqrt(EW.*NS);
        HV = geo_mean./V;
        if length(HV) == 100000
            HV = HV(1:50000);
        end
        dat.HV = geo_mean./V;
        if length(V) == 100000
            V = V(1:50000);
        end
        dat.geomean = geo_mean;
        dat.az = data.meta.event.azimuth;
        dat.dist = data.meta.event.epi_dist;
        dat.mag = data.meta.event.mag;
        dat.lat = data.meta.event.lat;
        dat.lon = data.meta.event.lon;
        fname = strcat(outpath,'\',station,'\',event_name);
        save(fname, 'dat')
    end
    toc
end