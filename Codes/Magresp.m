%% Magresp
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
% Date: developed between September, 2017 and August, 2019
%% Start   
function [X_mag, fax_HzN]=  Magresp(x, fs, ff)
N = length(x); %length of signal
X = fft(x); %fft of signal
X_mag = 4*abs(X)/(N); %normalized magnitude spectra (multiply by 4 is to account for hanning window and one-sided response
fax_binsN = [0 : N-1]; % number of frequency samples in the record
fax_HzN1 = fax_binsN*fs/N; %frequency axis NS (Hz)
N_2 = ceil(N/ff) - fs; %half magnitude spectrum - 1 hz
fax_HzN = fax_HzN1(1 : N_2);
X_mag = X_mag(1: N_2);
end