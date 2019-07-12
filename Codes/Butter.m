function [xNS, xV, xEW] = Butter(xNS, xV, xEW, fs)
% Here are the filter parameters
HighCorner = (fs / 2) - 1.01;
Npoles = 4;  % Corner for 1 pass of the two-pass filter
%Filter Design
fN = (fs / 2) - 1;
Highcut = HighCorner / fN;
[bb1, aa1] = butter(Npoles, Highcut);
%Filter soft site (Filename 1)
xNS = filtfilt(bb1,aa1,xNS - mean(xNS));
xV = filtfilt(bb1,aa1,xV - mean(xV));
xEW=filtfilt(bb1,aa1,xEW - mean(xEW));
end