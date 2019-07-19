function [matrix,peakind,ahatf1,newfaxhz1] = peakiden(ahatf, newfaxhz, lowbound, fsmin)
% N = length(ahatf)-lowbound; %length North_South_Component
% width = 2; %width for triangle moving average filter in hz
% q = ceil((fsmin/2-1)*width); %width for triangle moving average filter in samples
% e = smooth(ahatf, q, 'moving');
% e = e(lowbound: length(e));
newfaxhz1 = newfaxhz(lowbound: length(newfaxhz));
ahatf1 = ahatf(lowbound:length(ahatf));

peakfreq = [];
peakamp = []; 
amplocs2 = [];
%% Determine if peak is a peak
count = 0;
[peaks,locs,w,p] = findpeaks(ahatf1,newfaxhz1);
[B,I] = sort(p);
for ii = 1:length(B)
    if B(ii) > 1.5
        o = ii;
        break
    end 
end 
B=B(o:end);
w = w(I);
topprom = B;
topwidth = w(o:end);
ind = I(o:end);
A = peaks(ind); fn = locs(ind);
[fn,I2] = sort(fn);
A = A(I2); topprom = topprom(I2); topwidth = topwidth(I2);
A=A'; topprom = topprom'; topwidth = topwidth';fn = fn';
matrix = [fn,A,topprom, topwidth];
for j = 1:length(fn)
    peakind(j)=find(newfaxhz1==fn(j));
end
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