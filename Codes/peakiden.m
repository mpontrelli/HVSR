function [peakamp, peakfreq, amplocs2] = peakiden(ahatf, newfaxhz)
N = length(ahatf); %length North_South_Component
width = .1; %width for triangle moving average filter in hz
q = ceil((N/20)*width); %width for triangle moving average filter in samples
e = smooth(ahatf, q, 'moving');
%% Determine if peak is a peak
count = 0;
[amps,amplocs] = findpeaks(e);
for hh = 1:length(amps)
    freqpeak(hh) = newfaxhz(amplocs(hh));
end
[trough,troughloc] = findpeaks(-1*e);
for hh = 1:length(trough)
    freqtrough(hh) = newfaxhz(troughloc(hh));
end
if freqtrough(1) > freqpeak(1)
    Y = 1;
    X = 1;
    trough = horzcat(Y,trough');
    troughloc = horzcat(X,troughloc');
end
if freqpeak(length(freqpeak)) > freqtrough(length(freqtrough))
    Y = -1;
    X = 19991;
    trough = horzcat(trough, Y);
    troughloc = horzcat(troughloc,X);
end
for hh = 1:length(trough)
    freqtrough(hh) = newfaxhz(troughloc(hh));
end
for k = 1:length(amps)
    height = amps(k)/sqrt(2);
    lefttrough = -1 * trough(k);
    righttrough = -1 * trough(k + 1);
    if lefttrough < height && righttrough < height
        count  = count + 1;
        peakfreq(count) = newfaxhz(amplocs(k));
        peakamp(count) = amps(k); 
        amplocs2(count) = amplocs(k);
    end
end
figure 
end