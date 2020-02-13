%% Response_Spectra

% Compute response spectra following the procedure used in section 6.6.5
% in Chopra "Dynamics of Structures" using scaling factors for deformation
% spectra to acceleration spectra in section 6.6.3 of the same book.
% Integrate the acceleration time series using the Crank-Nicholson
% approximation. The acceleration is the pseudo-response spectra.

% Theory

% Crank-Nicholson approximation

% Pseudo acceleration response spectrum calculation

    % INPUTS
    
    % x - Input ground motion
    
    % fs - sampling frequency of ground motion (used to solve for dt by
    % doing 1/fs)
    
    % zeta - damping ratio in percent (this gets converted to a decimal)

    % OUTPUTS
    
    % y - a vector containing the acceleration response spectra
    

%% Author: Marshall Pontrelli
% Date: 2/11/2020

% Edits 2/12 - removed plots and computed values at important periods
%% start 
function [y, displacement, velocity, period, max1, per_max, resp_02_secs, resp_05_secs, resp_1_secs, resp_2_secs, resp_5_secs    ] = Response_Spectra(x, fs, zeta)
    alpha = 0.5; %value for Crank-Nicholson approximation
    dt = 1/fs;
    zeta = zeta /100;

    %% Loop over SDOF systems of differing period between 0 and 10 seconds
    for ii = 1:1000
        Tn = (ii-1)*0.01;
        % inputs for the stucture
        wn = 2*pi / Tn; %circular frequency

        % coefficients for the crank nicholson approximation
        a1 = 1 + 0.5*dt*2*zeta*wn + 0.25*dt^2*wn^2;
        a2 = 1 - 0.5*dt*2*zeta*wn - 0.25*dt^2*wn^2;
        a3 = -dt*wn^2;
    
        % velocity and displacement using the Crank-Nicholson approximation
        v = 0;
        u = 0;
        for i = 1:length(x)-1
            v(i+1) = (1/a1)*(a2*v(i)+ a3*u(i) + dt*(alpha*x(i+1)+(1-alpha)*x(i)));
            u(i+1) = u(i) + (alpha*v(i+1) + (1 - alpha) *v(i))*dt;
        end

        displacement(ii) = max(abs(u)); % displacement response spectra
        velocity(ii) = max(abs(v)); % velocity response spectra
        y(ii) = max(abs(u)) * (2*pi/Tn)^2; % pseudo-acceleration response spectra using equaiton 6.6.3 from Chopra
        period(ii) = Tn; % vector of periods
    end
%     y_freq = flip(y);
%     freq = 1/period;
    
    %% now find some important periods
    
    % maximum
    [max1, I] = max(y);
    per_max = period(I);
    
    % 0.2 secs
    [~, I] = min(abs(period - 0.2));
    resp_02_secs = y(I);
    
    % 0.5 secs
    [~, I] = min(abs(period - 0.5));
    resp_05_secs = y(I);    
    % 1 sec
    [~, I] = min(abs(period - 1));
    resp_1_secs = y(I);
    % 2 secs
    [~, I] = min(abs(period - 2));
    resp_2_secs = y(I);    
    % 5 secs
    [~, I] = min(abs(period - 5));
    resp_5_secs = y(I);
end
