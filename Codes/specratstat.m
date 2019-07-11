function [taxstat] = specratstat(peakamp, peakfreq, amplocs2, ahatf, newfaxhz, sigma)
Taxstat = [];
for f = 1:length(peakamp)
    taxstat(f,1) = f;
    taxstat(f,2) = peakfreq(f); 
    A = peakamp(f);
    taxstat(f,3) = A;
    amploc2 = amplocs2(f);
    [I1, I2, f1, f2, hpb] =  HalfPowerBand2(A, amploc2, newfaxhz, ahatf); 
    taxstat(f,4) = hpb;
    a = sigma(I1:I2);
    sigmai = median(a);
    taxstat(f,5) = sigmai;
end
end