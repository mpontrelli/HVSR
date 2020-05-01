%% process Mexico City RACM data using code Rosetta and put it into a .mat file
% for each event at each station in folder "Processed_data2"
close all
clear all
%% ACCESS THE DATA
% go into the data folder and get a list of stations
codepath = 'C:\Users\mpontr01\Desktop\HVSR\Codes';
datapath = 'C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Data';
outpath = 'C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Processed_data2';
cd(datapath)
stationlist = dir;
stationlist = stationlist(3:length(stationlist));
event_num = {};
for i = 1% : length(stationlist)
    station = stationlist(i).name;
    cd(strcat(datapath, '\', station))
    eventlist = dir;
    eventlist = eventlist(3:length(eventlist));
    event_num{i,1} = station;
    event_num{ i,2} = length(eventlist);
    tic
    for j = 1%1: length(eventlist)
        %disp(j)
        filename = eventlist(j);
        event_name = filename.name;
        filename = strcat(filename.folder, '\', filename.name);
        cd(codepath)
        data = Rosetta(filename);
        fname = strcat(outpath,'\',station,'\',event_name, '.mat');
        parsave(fname, data)
    end
    toc
end