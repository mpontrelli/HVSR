%% Lysmer_input

% Inputs soil profile into raylee_lysmer

%% Author: Marshall Pontrelli
% Date: 11/6/2020

close all
clear all

%% inputs
Nn = 2; % number of elements in solid part of model
Nnf = 0; % number of elements in fluid part of model
hv = [5,5]; % vector of grid spacings for solid (meters)
hfv = 0; % vector of grid spacings for fluid (meters)
f =60; % frequency (Hz)
modn = 1; % which mode (1=fundamental, 2=first overtone, etc)
vsv = [150, 170]; % shear velocity model in solid, a vector (m/s)
vpv = vsv*1.87; % compressional velocity model in solid, a vector (m/s)
rhov = [1600, 1700]; % density model in solid, a vector (kg/m^3)
vpfv = 0; % compressional velocity model in fluid, a vector (m/s)
rhofv = 0; % density model in fluid, a vector (kg/m^3)

[kk, vpk, vgk, ev] = ... 
    raylee_lysmer(Nn,vsv,vpv,rhov,f,hv,modn,Nnf,vpfv,rhofv,hfv)