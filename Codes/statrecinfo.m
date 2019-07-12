function [fsmin, recmax] = statrecinfo(files, station)
counter1 = 0;
lengthvec = [];
fsvec = [];
for file = files'
    counter1 = counter1 + 1;
    filename = strcat(station,'\',file.name);
    [xNS,xV,xEW, fs] = readfile1(filename);
    lengthvec(counter1, :) = length(xNS);
    fsvec(counter1, :) = fs;
end
fsmin = min(fsvec);
recmax = max(lengthvec);
end