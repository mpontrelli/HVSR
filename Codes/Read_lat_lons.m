%% Read lat-lons
close all
clear all


sitelist = {'AE02','AL01','AO24', 'AP68','AR14','AU11','AU46','BA49','BL45','BO39',...
    'CA20', 'CA59', 'CB43','CC55','CE18', 'CE23','CE32','CH84','CI05','CJ03',...
    'CJ04','CO47', 'CO56','CP28','CS78','CT64', 'CU80', 'DM12','DR16', 'DX37','EO30','ES57','EX08','EX09','EX12','FJ74','GA62',...
    'GC38','GR27','HA41','HJ72', 'IB22','IM40','JA43','JC54','LI33', 'LI58', 'LV17','ME52',...
    'MI15', 'MY19', 'NZ20', 'NZ31','PA34', 'PD42','PE10', 'RI76', 'RM48',...
    'SI53', 'SP51', 'TE07','TH35', 'TL08', 'TL55','TP13', 'UC44', 'UI21', 'VG09', 'VM29',...
    'XP06'};

filepath = 'C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Processed_data2\';
mat = [];

for ii = 1:length(sitelist)
    statname = sitelist{ii};
    station = strcat(filepath,statname);
    cd(station)
    files = dir;
    files = files(3:length(files));
    file = files(end);
    filename = strcat(station,'\',file.name);
    load(filename)
    lat = data.meta.station.lat;
    lon = data.meta.station.lon;
    mat(ii,1) = lat;
    mat(ii,2) = lon;
    clear data
    disp(ii)
end