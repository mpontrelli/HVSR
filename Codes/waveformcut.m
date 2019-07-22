function [xNS,xV,xEW] = waveformcut(xNS,xV,xEW, recmin)
    if recmin == length(xNS)
        xNS = xNS;
        xV = xV;
        xEW = xEW;
        return
    end
    [~, a]= max(xNS);
    [~, b] = max(xV);
    [~, c] = max(xEW);
    precuta = ceil(recmin/4);
    postcuta = recmin - precuta;
    precutb = precuta;
    postcutb = postcuta;
    precutc = precuta;
    postcutc = postcuta;
    %This set of conditional statements assures that the cut waveform does
    %not fall below 0 or above the length of the waveform. It cuts around
    %the maximum, preferable with a ratio of 1/4 to the left of the max and
    %3/4 to the right of the max, but given that these ratios are less than
    %zero or greater than the max number of samples respectively, the
    %difference between these bounds and the cut value is added or
    %subtracted to the ideal cut value. See photo
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