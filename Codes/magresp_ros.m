%% magresp_ros

% Compute the magnitude response of a time series and the corresponding
% frequency vector. Only outputs half the magnitude response below the
% frequency vector and multiplies by 2 to account for cos(x) = 1/2 (e^jw +
% e^0jw) so we recover the amplitude of the sinusoid at a single frequency
% instead of 2.

    % INPUTS
    
    % x - times series input
    
    % fs - sampling frequency of the record
    
    % ff - positive value, one over which is the percent of the magnitude
    % response you want out.
    
    % OUTPUTS
    
    % X - magnitude response
    
    % fax_HzN - frequency vector of X
    

%% Author: Marshall Pontrelli
% Date: developed between September, 2017 and August, 2019, altered
% 2/20/2020 from magresp to work with Rosetta
%% Start   
function [X_mag, X_mag_smooth]=  magresp_ros(x, fs, len)
    width = 0.3;
    N = length(x);
    win = hann(N);
    x_win = x.*win;
    X_mag = 4*abs(fft(x_win))/len;
    window = ceil((N/fs)*width); %width for smoothing filter in samples where 20 is the number of Hz on your x-axis
    X_mag_smooth = smooth(X_mag,window);
    N_2 = floor(N/2); %half magnitude spectrum
    X_mag = X_mag(1 : N_2);
    X_mag_smooth = X_mag_smooth(1 : N_2);
end