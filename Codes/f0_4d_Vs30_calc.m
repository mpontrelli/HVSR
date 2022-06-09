%% f0 = Vs/4d calculation
% This code takes an f0, an overburden Vs, and a bedrock Vs and calcualtes
% Vs30
% Marshall Pontrelli
% 6/9/2022

%% 
close all
clear all

%% Start
f0 = 2.7;
Vs = 180;
d = Vs/(4*f0);

cd(strcat('C:\Users\',getenv('username'),'\Desktop\HVSR\Codes'));
%% Vs = 180 proglacial fine
v1 = 180;
for i = 1:200
    d1 = i;
    d2 = 200-d1;
    d_vec = [d1,d2];
    v = [v1,2500];
    [vs_30, class] =  Vs_30(d_vec,v);
    f0(i) = v1/(4*d1);
    Vs30_vec(i) = vs_30;
end