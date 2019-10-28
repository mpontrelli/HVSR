%% waveformcut
% waveform cut cuts all the records from the dataset to the length of the
% record with the lowest number of samples. It came about after questions
% about the validity of upsampling using linear interpolation. It now is an
% option the user can input. It finds the peak value in the record, which
% is likely the peak value of the shear wave, and cuts to 3/4 the number of
% samples after the peak and 1/4 the number of samples before the peak.

    % INPUTS
    
    % xNS - north south component of the ground motion record

    % xV - vertical component of the ground motion record
        
    % xEW - east west component of the ground motion record
    
    % recmin - minimum number of samples in the records. This is an output
    % of statrecinfo
    
    % OUTPUTS
    
    % xNS - cut north south component of the ground motion record

    % xV - cut vertical component of the ground motion record
        
    % xEW - cut east west component of the ground motion record
    
%% Author: Marshall Pontrelli
% Date: developed between September, 2017 and August, 2019
%% Start   
function [xNS,xV,xEW] = waveformcut(xNS,xV,xEW, recmin)
    if recmin == length(xNS) % if the record being tested is the shortest record, skip
        xNS = xNS;
        xV = xV;
        xEW = xEW;
        return
    end
    
    % find the peaks
    [~, a]= max(xNS);
    [~, b] = max(xV);
    [~, c] = max(xEW);
    % number of samples to cut to
    precuta = ceil(recmin/4);
    postcuta = recmin - precuta;
    precutb = precuta;
    postcutb = postcuta;
    precutc = precuta;
    postcutc = postcuta;
    % This set of conditional statements assures that the cut waveform does
    % not fall below 0 or above the length of the waveform. It cuts around
    % the maximum, preferable with a ratio of 1/4 to the left of the max and
    % 3/4 to the right of the max, but given that these ratios are less than
    % zero or greater than the max number of samples respectively, the
    % difference between these bounds and the cut value is added or
    % subtracted to the ideal cut value. See photo
    if a - precuta < 0
        q = abs(a - precuta) + 1;
        precuta = precuta - q;
        postcuta = postcuta + q;
    end
    if a + postcuta > length(xNS)
        q = a + postcuta - length(xNS);
        postcuta = postcuta - q;
        precuta = precuta + q;
    end
    xNS = xNS(a - precuta: a + postcuta - 1);
    if b - precutb < 0
        q = abs(b - precutb) + 1;
        precutb = precutb - q;
        postcutb = postcutb + q;
    end
    if b + postcutb > length(xV)
        q = b + postcutb - length(xV);
        postcutb = postcutb - q;
        precutb = precutb + q;
    end
    xV = xV(b - precutb: b + postcutb -1);
    if c - precutc < 0
        q = abs(c - precutc) + 1;
        precutc = precutc - q;
        postcutc = postcutc + q;
    end
    if c + postcutc > length(xEW)
        q = c + postcutc - length(xEW);
        postcutc = postcutc - q;
        precutc = precutc + q;
    end
    xEW = xEW(c - precutc: c + postcutc - 1);      
end