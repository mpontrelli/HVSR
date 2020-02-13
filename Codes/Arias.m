%% Arias
% compute the Arias intensity, D_5-95, D_5-75 and rate of increase of Arias
% intensity from a time vector, an input strong ground motion and a
% sampling frequency. The input ground motion must be in g's. 

    % INPUTS
    % time - time vector to index for duration calculations
    
    % X - a strong ground motion in g's
    
    % fs - sampling frequency of the record
    
    % OUTPUTS
    
    % Iaval - Arias intensity
    
    % D5D95 - significant duration a normalized arias intensities between 5
    % and 95 
    
    % D5D75 - same as above but for 75
    
    % rate_arias - rate of increase defined as Ia_90 / D_5-95
    
    % Ianorm - vector of the normailed Arias Intensity plot if you want to
    % plot it later. 
    
    
%% Author: Marshall Pontrelli
% Date: Spring 2019 for Earthquake Engineering CEE 247, summer 2019
% edited - 2/12/2020 - formatting and editing to make compatable with
% updated code.


function [Iaval, D5D95, D5D75, rate_arias, Ianorm] =  Arias(time, X, fs)
    dt = 1/fs;
    IaX = 0;
    for i = 1:length(X)-1
        IaX(i+1) = IaX(i) + (X(i)^2+X(i+1)^2)*dt/2;   
    end
    IaX = (IaX * pi/(2*9.8));
    Iaval = max(IaX);
    Ianorm = IaX / Iaval;
    [~,I5] = min(abs((Ianorm - 0.05)));
    D5 = time(I5);
    [~,I95] = min(abs((Ianorm - 0.95)));
    D95 = time(I95);
    D5D95 = D95 - D5;
    [~,I75] = min(abs((Ianorm - 0.75)));
    D75 = time(I75);
    D5D75 = D75 - D5;
    [~,I90] = min(abs((Ianorm - 0.90)));
    rate_arias = IaX(I90)/D5D95;
end