%% Try to model waveforms
close all
clear all
% Author: Marshall Pontrelli
% Date: 5/15/2020
cs = 100;
f = 1;
w = f*2*pi;
k = w/cs;
t = 2;


x = linspace(0,100,20);
y = linspace(0,100,20);
z = linspace(0,100,20);
[X,Y,Z] = meshgrid(x,y,z);
U0 = real(exp(1i*(w*t - k*X - k*Y - k*Z)))+real(exp(1i*(w*t - k*X + k*Y - k*z)));
colorbar
surf(X,Y,Z, U0)
F(ii) = getframe(gcf);
close all

%%
for ii = 1:100
    t = (ii - 1)*0.01;
    x = linspace(0,100,20);
    y = linspace(0,100,20);
    [X,Y] = meshgrid(x,y);
    U0 = real(exp(1i*(w*t - k*X - k*Y)))+real(exp(1i*(w*t - k*X + k*Y)));
    colorbar
    surf(X,Y,U0)
    F(ii) = getframe(gcf);
    close all
end

%%
movie(F)