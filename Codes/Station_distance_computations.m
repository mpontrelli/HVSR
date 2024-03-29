% This code takes a network of stations and computes the distance to each
% station from each station as well as the distance north (or south) and
% east (or west)

% Author - Marshall Pontrelli
% Date - 1/13/2020

%% start
close all
clear all
% If you have an event, plug in the epicenter and compute the distance from
% the epicenter to the station
Epicenter_lat = 18.40;
Epicenter_lon = -98.72;
event_depth = 57; %km

Station_lon = -99.1353;
Station_lat = 19.4198;

dist = deg2km(distance(Station_lat, Station_lon,Epicenter_lat, Epicenter_lon))
az = azimuth(Station_lat, Station_lon,Epicenter_lat, Epicenter_lon)

% station distance
% import all the station latitudes and longitudes
T = readtable('C:\Users\mpontr01\Box\2020_1_spring\SSA\Time domain\Tables\Stations.xlsx','Range','B1:C62');
A = table2array(T);
A(61,:) = [];
Asize = size(A);
len = Asize(1);
matrix = zeros(len,4);
for i = 1:len
    Stat2_lon = A(i,1);
    Stat2_lat = A(i,2);

    % [dist, azimuth] =  hypdist(Station_lat,Stat2_lat, Station_lon, Stat2_lon)
    matrix(i,1) = deg2km(distance(Station_lat, Station_lon,Stat2_lat, Stat2_lon));
    matrix(i,2) = azimuth(Station_lat, Station_lon,Stat2_lat, Stat2_lon);

    % NS, compute NS distance
    Stat2_lon2 = Station_lon;
    Stat2_lat2 = Stat2_lat;
    matrix(i,3) = deg2km(distance(Station_lat, Station_lon,Stat2_lat2, Stat2_lon2));
    % EW, compute EW distance
    Stat2_lon3 = Stat2_lon;
    Stat2_lat3 = Station_lat;
    matrix(i,4) = deg2km(distance(Station_lat, Station_lon,Stat2_lat3, Stat2_lon3));
end