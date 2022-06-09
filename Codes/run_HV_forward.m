% Run HV_forward
close all
clear all
n = 2;
t = [10;10];
Vp = [500;600];
Vs = [200;210];
rho = [2100;2150];
Vp_base = [3000];
Vs_base = 2000;
rho_base = 2700;
fmin1 = 0.1;
fmax1 = 20;
NR = 10;
NL = 10;
dk = 500;
log = 0;
[HV,f] = HV_forward(n,t, Vp, Vs, rho, Vp_base, Vs_base, rho_base, fmin1, fmax1, NR, NL, dk, log);