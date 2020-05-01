%% peakiden
% indentify significant peaks above a certain prominence threshold.
% Prominence is defined in the matlab "findpeaks" documentation.
% https://www.mathworks.com/help/signal/ref/findpeaks.html#bufbbs1-p

    % INPUTS
    
    % ahatf - a spectral ratio, can be averaged from microtremor or just a
    % normal transfer function
    
    % newfaxhz - vector of frequencies the same length as ahatf
    
    % lowbound - the value of the lowest frequency of interest
    
    % OUTPUTS
    
    % matrix - column 1 is the frequency of the significant peak, column 2
    % is the amplitude of the significant peak
    
    % matrix1 - column 1 is the prominence of the significant peak, column
    % 2 is the width of the significant peak
    
    % peakind - peak index
    
    % ahatf1 - the clipped transfer function vector to take the lowbound and
    % upbound values into account 
    
    % newfaxhz1 - the clipped frequency vector to take the lowbound and
    % upbound values into account
    

%% Author: Marshall Pontrelli
% Date: developed during summer 2019 
% Update - 1/16/2020 added comments while I was coming back to debug this
% and specratstat - Marshall

% Update - 4/29/2020 added peak areas

function [matrix, matrix1, peakind,ahatf1,newfaxhz1, peakfreqs, peakamps, Areamat] = peakiden(ahatf, newfaxhz, lowbound, upbound)
% find the values in the frequency vector closest to the lowbound and
% upbound values

newfaxhz1 = newfaxhz(lowbound:upbound);
ahatf1 = ahatf(lowbound: upbound);

%% Determine if peak is a peak
[peaks,locs,w,p] = findpeaks(ahatf1);
counter = 0;
% loop through the prominence vector
for ii = 1:length(p)
    if peaks(ii) - p(ii) < peaks(ii)/sqrt(2) % if the peak minus the prominence is less than the peak/sqrt(2), this is a significant peak
        counter = counter+1;
        o(counter) = ii; %this peak is not a sig peak, so minus 1 to get rid of it
    end 
end 
topprom = p(o)';
%topwidth = w(o)';
A = peaks(o)'; peakind = locs(o)'; fn = newfaxhz1(peakind)';
matrix = [fn',A'];
matrix1 = topprom';%, topwidth];


%% now compute area under the peak and output the polygon for area
AA = A - topprom;
peakfreqs = cell(length(peakind),1);
peakamps = cell(length(peakind),1);
for f = 1:length(peakind)
    loc = peakind(f);
    cur_amp = AA(f);
    %move down signal to the right
    for i = 1:length(newfaxhz1)
        ii = loc + i;
        k = ahatf1(ii);
        if 0.01 - (k - cur_amp) > 0
            I1 = ii;
            break
        end
    end
    %move down signal to the left
    for i = 1:length(newfaxhz1)
        ii = loc - i;
        k = ahatf1(ii);
        if 0.01 - (k - cur_amp) > 0
            I2 = ii;
            break
        end
    end
    % Now make the vectors
    peak_freqs = newfaxhz1(I2:I1)';
    peak_amps = ahatf1(I2:I1)';
    
    % Now compute areas using trapezoids
    for i = 1:length(peak_freqs) - 1
        if i == 1
            Area(i) = ((peak_freqs(i+1) - peak_freqs(i)) * peak_amps(i+1))/2;
        else
            Area(i) = ((peak_freqs(i+1) - peak_freqs(i)) * peak_amps(i)) + ((peak_freqs(i+1) - peak_freqs(i)) * (peak_amps(i) - peak_amps(i +1)))/2;
        end          
    end
    Areamat(f) = sum(Area);
    
    % Now save
    peakamps{f} = peak_amps;
    peakfreqs{f} = peak_freqs;
    
    clear peak_freqs
    clear peak_amps
    clear Area
end
end