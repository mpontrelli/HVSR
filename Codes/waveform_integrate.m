%% waveform_integrate
% Integrate for velocity and displacement given a acceleration time history

%% Author: Marshall Pontrelli
% Date: 2/18/2020 based on code from Spring 2019 for Earthquake Engineering 
% CEE 247 which was edited in summer 2019 and fall 2019
% edited - 2/18/2020 - formatting and editing to make compatable with
% updated code.


function [PGA, PGV, PGD, v, d] =  waveform_integrate(x, fs)

    PGA = max(abs(x));
    %Inputs
    alpha = 0.25; %coefficients for Newmark-beta method
    beta = 0.5;
    LowCorner = 0.1;
    HighCorner = 49;
    Npoles = 4;
    Filterplot = 'no';
    dt = 1/fs;
    %% integrate for velocity (Newmark-Beta method)
    v = 0;
    for i = 1:length(x) - 1
        v(i+1) = v(i) + dt*((1-beta)*x(i) + beta*x(i+1));
    end
    v = v - mean(v);
    [v] = Butter2(v, fs);
    PGV = max(abs(v));
    %% Integrate for displacement 
    d = 0;
    for i = 1:length(x) - 1
        d(i+1) = d(i) + (dt)*v(i) + (dt)^2 * ((0.5-alpha)*x(i) + alpha*x(i+1));
    end
    d = d - mean(d);
    [d] = Butter2(d, fs);
    PGD = max(abs(d));
end 