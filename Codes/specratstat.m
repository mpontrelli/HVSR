function [taxstat] = specratstat(peakamp, peakfreq, amplocs2, ahatf, newfaxhz, sigma, statname)
taxstat = {};
for f = 1:length(peakamp)
    taxstat{f,1} = statname;
    taxstat{f,2} = num2str(f);
    taxstat{f,3} = num2str(peakfreq(f)); 
    A = peakamp(f);
    taxstat{f,4} = num2str(A);
    amploc2 = amplocs2(f);
    [I1, I2, f1, f2, hpb] =  HalfPowerBand2(A, amploc2, newfaxhz, ahatf); 
    taxstat{f,5} = num2str(hpb);
    a = sigma(I1:I2);
    sigmai = median(a);
    taxstat{f,6} = num2str(sigmai);
end
end