%% solution for geomechanics final
% Marshall Pontrelli
% 5/8/2020
close all
clear all

basement_shear = 800;
basement_rho = 2.7;
basement_z = 0.05;
basement_h = 50;

cs = [70, basement_shear];
rho = [1.5, basement_rho];
zeta = [0.05,basement_z];
h = [50, basement_h];

velocity_profile(h, cs)
%% TEST
h_cur = 50;
cs_bot = 800;
cs_top = 70;
zeta_bot = 0.05;
zeta_top = 0.05;
rho_bot = 2.7;
rho_top = 1.5;
omega = linspace(0, 100, 1000)/(2*pi);
for j = 1:length(omega)
    o = omega(j);
    amp(j) = 1/(abs(cos(2*pi*h_cur*o/(cs_top*(1+1i*zeta_top))) + (rho_top*(cs_top*(1+1i*zeta_top)))/(rho_bot*(cs_bot*(1+1i*zeta_bot)))*1i*sin(2*pi*h_cur*o/(cs_top*(1+1i*zeta_top))))); 
end
figure
plot(omega,amp)
