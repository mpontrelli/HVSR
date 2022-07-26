%% f0 = Vs/4d calculation
% This code takes an f0, an overburden Vs, and a bedrock Vs and calcualtes
% Vs30
% Marshall Pontrelli
% 6/9/2022

%% 
close all
clear all

%% Start
f0 = 3.18;
v1 = 220;
d = v1/(4*f0);

cd(strcat('C:\Users\',getenv('username'),'\Desktop\HVSR\Codes'));
%% Vs = 180 proglacial fine
d1 = d;
d2 = 200-d1;
d_vec = [d1,d2];
v = [v1,2500];
[vs_30, class] =  Vs_30(d_vec,v);

%% Braganza equation
%vs30_Braganza = 30/((1/(4*f0))*(1-(v1/2500))+30/2500)

