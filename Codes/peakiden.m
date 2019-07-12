function [peakamp, peakfreq, amplocs2] = peakiden(ahatf, newfaxhz, lowbound)
N = length(ahatf)-lowbound; %length North_South_Component
width = .1; %width for triangle moving average filter in hz
q = ceil((N/20)*width); %width for triangle moving average filter in samples
e = smooth(ahatf, q, 'moving');
e = e(lowbound: length(e));
newfaxhz1 = newfaxhz(lowbound: length(newfaxhz));
disp(length(newfaxhz1))
%% Determine if peak is a peak
count = 0;
[amps,amplocs] = findpeaks(e);
for hh = 1:length(amps)
    freqpeak(hh) = newfaxhz1(amplocs(hh));
end
[trough,troughloc] = findpeaks(-1*e);
for hh = 1:length(trough)
    freqtrough(hh) = newfaxhz1(troughloc(hh));
end
if freqtrough(1) > freqpeak(1)
    Y = 1;
    X = 1;
    trough = horzcat(Y,trough');
    troughloc = horzcat(X,troughloc');
end
if freqpeak(length(freqpeak)) > freqtrough(length(freqtrough))
    Y = -1;
    X = newfaxhz1(length(newfaxhz1));
    trough = horzcat(trough, Y);
    troughloc = horzcat(troughloc,X);
end
for hh = 1:length(trough)
    freqtrough(hh) = newfaxhz1(troughloc(hh));
end
for k = 1:length(amps)
    height = amps(k)/sqrt(2);
    lefttrough = -1 * trough(k);
    righttrough = -1 * trough(k + 1);
    if lefttrough < height && righttrough < height
        count  = count + 1;
        peakfreq(count) = newfaxhz1(amplocs(k));
        peakamp(count) = amps(k); 
        amplocs2(count) = amplocs(k);
    end
end
figure 
end