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

% Update 5/3/2020 - made peak definition more rigorous and changed outputs
% to reflect new classification scheme

function [peakfreqs, peakamps, hpb, f1s, f2s, Areamat, proms, amps, peakind2, freqs, sigs, I1s, I2s] = peakiden(ahatf, newfaxhz, sigma, lowbound, upbound)
% find the values in the frequency vector closest to the lowbound and
% upbound values
newfaxhz1 = newfaxhz(lowbound:upbound);
ahatf1 = ahatf(lowbound: upbound);
sigma1 = sigma(lowbound: upbound);

%% Step 1 in determining if peak is a peak, does it have a halfpower bandwidth?
[peaks,locs,w,p] = findpeaks(ahatf1);
counter = 0;
% loop through the prominence vector
o =0;
for ii = 1:length(p)
    if peaks(ii) - p(ii) < peaks(ii)/sqrt(2) % if the peak minus the prominence is less than the peak/sqrt(2), this is a significant peak
        counter = counter+1;
        o(counter) = ii; %this peak is not a sig peak, so minus 1 to get rid of it
    end
end 
if o == 0
    I2s = [];
    I1s = [];
    sigs = [];
    hpb = [];
    f1s = [];
    f2s = [];
    Areamat = [];
    proms = [];
    amps = [];
    peakind2 = [];
    freqs = [];
    peakamps = [];
    peakfreqs = []; 
    disp('No peaks')
else
topprom = p(o)';
A = peaks(o)'; peakind = locs(o)'; fn = newfaxhz1(peakind)';



%% Step 2: Find the interval minima and make the peak shape, then figure out how many peaks 
% there are across that peak
AA = A - topprom;
peakfreqs = cell(length(peakind),1);
peakamps = cell(length(peakind),1);
counter = 0;
for f = 1:length(peakind)
    loc = peakind(f);
    cur_amp = AA(f);
    %move down signal to the right
    for i = 1:length(newfaxhz1)
        ii = loc + i;
        k = ahatf1(ii);
        if 0.001 - (k - cur_amp) > 0
            I1 = ii;
            break
        end
    end
    %move down signal to the left
    for i = 1:length(newfaxhz1)
        ii = loc - i;
        k = ahatf1(ii);
        if 0.001 - (k - cur_amp) > 0
            I2 = ii;
            break
        end
    end
    % Now make the vectors
    peak_freqs = newfaxhz1(I2:I1)';
    peak_amps = ahatf1(I2:I1)';
    peak_sigs = sigma1(I2: I1)';
   
    %% Now figure out the shape of the peak
    [peaks,locs,~,p] = findpeaks(peak_amps);
    if min(peaks./p) < 2 % Then it's a real peak
        counter = counter + 1;      
        %% Now compute areas using trapezoids
        for i = 1:length(peak_freqs) - 1
            if i == 1
                Area(i) = ((peak_freqs(i+1) - peak_freqs(i)) * abs((peak_amps(i+1)-peak_amps(i))))/2;
            else
                Area(i) = (((peak_freqs(i+1) - peak_freqs(i)) * abs((peak_amps(i)-peak_amps(1)))) + ((peak_freqs(i+1) - peak_freqs(i)) * abs((peak_amps(i+1) - peak_amps(i -1))))/2);
            end          
        end
        
        %% now compute halfpower bandwidth
        amp = A(f)/sqrt(2);
        Ind =find(peak_amps == A(f));
    
        %move down signal to the right
        for i = 1:length(peak_freqs)
            ii = Ind + i;
            k = peak_amps(ii);
            if k < amp
                I2 = ii;
                f2 = peak_freqs(I2);
                break
            end
        end
        %move down signal to the left
        for i = 1:length(peak_freqs)
            ii = Ind - i;
            k = peak_amps(ii);
            if k < amp
                I1 = ii;
                f1 = peak_freqs(I1);
                break
            end
        end
        
        %% Now sigma
        I2s(counter) = I2;
        I1s(counter) = I1;
        sigs(counter) = median(peak_sigs);
        hpb(counter) = f2 - f1;
        f1s(counter) = f1;
        f2s(counter) = f2;
        %% and save all of the peak stats
        Areamat(counter) = sum(Area);
        proms(counter) = topprom(f);
        amps(counter) = A(f);
        peakind2(counter) = peakind(f);
        freqs(counter) = fn(f);
        peakamps{counter} = peak_amps;
        peakfreqs{counter} = peak_freqs; 
        clear peak_freqs
        clear peak_amps
        clear Area
    end
end
if counter == 0
    I2s = [];
    I1s = [];
    sigs = [];
    hpb = [];
    f1s = [];
    f2s = [];
    Areamat = [];
    proms = [];
    amps = [];
    peakind2 = [];
    freqs = [];
    peakamps = [];
    peakfreqs = []; 
    disp('No peaks')
else
% now compute slopes
% % for i = 1:length(ahatf1) - 1
% %     slope(i) = (ahatf1(i+1) - ahatf1(i))/ (newfaxhz1(i+1) - newfaxhz1(i));
% %     newnewfreq(i) = (newfaxhz1(i+1) + newfaxhz1(i))/2;
% % end
% % slope = smooth(slope);
% % slope = smooth(slope);
% % 
% % % Now compute slopes of slope
% % for i = 1:length(slope) - 1
% %     slope1(i) = (slope(i+1) - slope(i))/ (newnewfreq(i+1) - newnewfreq(i));
% %     nnnfreq(i) = (newnewfreq(i+1) + newnewfreq(i))/2;
% % end
end
end
