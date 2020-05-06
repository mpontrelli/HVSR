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

%% Step 1 in determining if peak is a peak using the whole signal, does it have a halfpower bandwidth?
[peaks,locs,~,p] = findpeaks(ahatf1);
counter = 0;
% loop through the prominence vector and find all the peaks that have hpbs
o =0;
for ii = 1:length(p)
    if peaks(ii) - p(ii) < peaks(ii)/sqrt(2) % if the peak minus the prominence is less than the peak/sqrt(2), this is a significant peak
        counter = counter+1;
        o(counter) = ii; %this peak is a sig peak, keep it
    end
end 
if o == 0 % If a station has no peaks with hpbs, then it has no peaks
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
else % else there are some hpb peaks, let's take a look
    topprom = p(o)'; % save the prominences, Amps, peakindexes and fns of all peaks with hpbs
    A = peaks(o)'; peakind = locs(o)'; fn = newfaxhz1(peakind)';


    %% Step 2: Now find the interval minima and make the peak shape, 
    % then 
    counter = 0;
    for f = 1:length(peakind) % loop across all peaks with hpbs
        loc = peakind(f); % peak the peak we're on
        cur_amp = A(f); % and the current amp
        % Now find the right interval minima
        for i = 1:(length(newfaxhz1) - loc - 1)
            ii = loc + i;
            k = ahatf1(ii);
            if k >= cur_amp
                I1 = ii;
                break
            else
                I1 = length(newfaxhz1);
            end
        end
        right_int = ahatf1(loc:I1);
        [r_min, ~] = min(right_int);
        r_min = find(ahatf1 == r_min);
        %move down signal to the left
        % Now find the left interval minima
        for i = 1:(loc - 1)
            ii = loc - i;
            k = ahatf1(ii);
            if k >= cur_amp
                I1 = ii;
                break
            else
                I1 = 1;
            end
        end
        left_int = ahatf1(I1:loc);
        [l_min, ~] = min(left_int);
        l_min = find(ahatf1 == l_min);
        
        % now find the highest minimum 
        [~, I] = max([ahatf1(l_min), ahatf1(r_min)]);
        if I == 1 % that means the left interval min is the highest
            for i = 1: length(ahatf1)
                q = loc + i +5;
                gg = ahatf1(q);
                if gg <= ahatf1(l_min)
                    r_min = q;
                    break
                end
            end
        end
        if I == 2 % that means the right interval min is the highest
            for i = 1: length(ahatf1)
                q = loc - i - 5;
                gg = ahatf1(q);
                if gg <= ahatf1(r_min)
                    l_min = q;
                    break
                end
            end
        end
        
        % Now make the vectors
        peak_freqs = newfaxhz1(l_min:r_min)';
        peak_amps = ahatf1(l_min:r_min)';
        peak_sigs = sigma1(l_min:r_min)';
   
        %% Now take the peak that is constrained between the left and right
        % intervals and see if it has any peaks within it that make the hpb
        % criteria
%         [peaks2,locs2,~,p2] = findpeaks(peak_amps);
%         o =0;
%         counter2 = 0;
%         for ii = 1:length(p2)
%             if peaks2(ii) - p2(ii) < peaks2(ii)/sqrt(2) % if the peak minus the prominence is less than the peak/sqrt(2), this is a significant peak
%                 counter2 = counter2+1;
%                 o(counter2) = ii; %this peak is not a sig peak, so minus 1 to get rid of it
%             end
%         end 
%         
%         % now find all the statistics from the peaks within the peak
%         topprom2 = p2(o)';
%         A2 = peaks2(o)'; peakind2 = locs2(o)'; fn2 = peak_freqs(peakind2)';
% 
%         %% Step 2: Find the interval minima and make the peak shape, then figure out how many peaks 
%         % there are across that peak
%         for ff = 1:length(peakind2)
%             loc = peakind2(ff); % peak the peak we're on
%             cur_amp = A2(ff); % and the current amp
%             % Now find the right interval minima
%             for i = 1:(length(peak_freqs) - loc - 5)
%             ii = loc + i +5;
%             k = peak_amps(ii);
%             if k >= cur_amp
%                 I1 = ii;
%                 break
%             else
%                 I1 = length(peak_freqs);
%             end
%         end
%         right_int = peak_amps(loc:I1);
%         [r_min, ~] = min(right_int);
%         r_min = find(peak_amps == r_min);
%         %move down signal to the left
%         % Now find the left interval minima
%         for i = 1:(loc - 5)
%             ii = loc - i + 5;
%             k = peak_amps(ii);
%             if k >= cur_amp
%                 I1 = ii;
%                 break
%             else
%                 I1 = 1;
%             end
%         end
%         left_int = peak_amps(I1:loc);
%         [l_min, ~] = min(left_int);
%         l_min = find(peak_amps == l_min);
%         
%         % now find the highest minimum 
%         [~, I] = max([peak_amps(l_min), peak_amps(r_min)]);
%         if I == 1 % that means the left interval min is the highest
%             for i = 1: length(peak_amps)
%                 q = loc + i +5;
%                 gg = peak_amps(q);
%                 if gg <= peak_amps(l_min)
%                     r_min = q;
%                     break
%                 end
%             end
%         end
%         if I == 2 % that means the right interval min is the highest
%             for i = 1: length(peak_amps)
%                 q = loc - i - 5;
%                 gg = peak_amps(q);
%                 if gg <= peak_amps(r_min)
%                     l_min = q;
%                     break
%                 end
%             end
%         end
%         % Now make the vectors
%         peak_freqs2 = peak_freqs(l_min:r_min)';
%         peak_amps2 = peak_amps(l_min:r_min)';
%         peak_sigs2 = peak_sigs(l_min:r_min)';
    
        %% Now compute areas using trapezoids
            for i = 1:length(peak_freqs) - 1
                if i == 1
                    Area(i) = ((peak_freqs(i+1) - peak_freqs(i)) * abs((peak_amps(i+1)-peak_amps(i))))/2;
                else
                    Area(i) = (((peak_freqs(i+1) - peak_freqs(i)) * abs((peak_amps(i)-peak_amps(1)))) + ((peak_freqs(i+1) - peak_freqs(i)) * abs((peak_amps(i+1) - peak_amps(i -1))))/2);
                end          
            end
        
            %% now compute halfpower bandwidth
            amp2 = cur_amp/sqrt(2);
            Ind2 =find(peak_amps == A(f));

            %move down signal to the right
            for i = 1:length(peak_freqs) - 1
                ii = Ind2 + i;
                k = peak_amps(ii);
                if k < amp2
                    I2 = ii;
                    f2 = peak_freqs(I2);
                    break
                end
            end

            %move down signal to the left
            for i = 1:length(peak_freqs) - 1
                ii = Ind2 - i;
                k = peak_amps(ii);
                if k <= amp2
                    I1 = ii;
                    f1 = peak_freqs(I1);
                    break
                end
            end
            counter = counter +1;
        
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
%           clear peak_freqs
%           clear peak_amps
%           clear Area
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
end % loop across all peaks with hpbs

