function [matrix, matrix1, peakind,ahatf1,newfaxhz1] = peakiden(ahatf, newfaxhz, lowbound, upbound, fsmin)
% N = length(ahatf)-lowbound; %length North_South_Component
% width = 2; %width for triangle moving average filter in hz
% q = ceil((fsmin/2-1)*width); %width for triangle moving average filter in samples
% e = smooth(ahatf, q, 'moving');
% e = e(lowbound: length(e));
newfaxhz1 = newfaxhz(lowbound: end);
ahatf1 = ahatf(lowbound:end);

%% Determine if peak is a peak
[peaks,locs,w,p] = findpeaks(ahatf1);
counter = 0;
for ii = 1:length(p)
    if peaks(ii) - p(ii) < peaks(ii)/sqrt(2)
        counter = counter+1;
        o(counter) = ii; %this peak is not a sig peak, so minus 1 to get rid of it
    end 
end 
topprom = p(o)';
topwidth = w(o)';
A = peaks(o)'; peakind = locs(o)'; fn = newfaxhz1(peakind)';
matrix = [fn,A];
matrix1 = [topprom, topwidth];

%plot(fn, A, 'o', 'markeredgecolor', 'g', 'markerfacecolor', 'g')
% for hh = 1:length(amps)
%     freqpeak(hh) = newfaxhz1(amplocs(hh));
% end
% [trough,troughloc] = findpeaks(-1*e);
% for hh = 1:length(trough)
%     freqtrough(hh) = newfaxhz1(troughloc(hh));
% end
% if freqtrough(1) > freqpeak(1)
%     Y = 1;
%     X = 1;
%     trough = horzcat(Y,trough');
%     troughloc = horzcat(X,troughloc');
% end
% if freqpeak(length(freqpeak)) > freqtrough(length(freqtrough))
%     Y = -1;
%     X = newfaxhz1(length(newfaxhz1));
%     trough = horzcat(trough, Y);
%     troughloc = horzcat(troughloc,X);
% end
% for hh = 1:length(trough)
%     freqtrough(hh) = newfaxhz1(troughloc(hh));
% end
% for k = 1:length(amps)
%     height = amps(k)/sqrt(2);
%     lefttrough = -1 * trough(k);
%     righttrough = -1 * trough(k + 1);
%     if lefttrough < height && righttrough < height
%         count  = count + 1;
%         peakfreq(count) = newfaxhz1(amplocs(k));
%         peakamp(count) = amps(k); 
%         amplocs2(count) = amplocs(k);
%     end
% end
% figure 
% plot(newfaxhz1,ahatf1)
% hold on
% plot(fn,A,'go')

end