close all
clear all

%% Open up the folder with the data
folder = 'C:\Users\mpontr01\Box\Projects\Kik-Net\Data\TKCH08\';
outpath = 'C:\Users\mpontr01\Box\Projects\Kik-Net\Processed_data\';
codepath = 'C:\Users\mpontr01\Desktop\HVSR\Codes';
station = 'TKCH08\';
cd(folder)
files = dir;
files = files(3:length(files));
name_mat = [];
counter = 0;
for file = files'
    a = strcat(folder,file.name);
    [filepath,name,ext] = fileparts(a);
    if strcmp(ext, '.EW1') ==1
        counter = counter + 1;
        name_mat{counter} = name;
    end
end

%% now loop through events
parfor iii = 1:length(name_mat)
    %% Now save
    event = name_mat{iii};
    cd(codepath)
    data = Rosetta_Kik_net(folder,event)
    fname = strcat(outpath,station,event,'.mat');
    parsave(fname, data);
    fclose('all')
end