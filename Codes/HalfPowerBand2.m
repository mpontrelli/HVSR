function [I1, I2, f1, f2, hpb]=  HalfPowerBand2(A, amploc2, newfaxhz, ahatf)
    %Half power bandwidth
    
    amp = A/sqrt(2);
    
    %move down signal to the right
    for i = 1:length(newfaxhz)
        ii = amploc2 + i;
        k = ahatf(ii);
        if k < amp
            I2 = ii;
            f2 = newfaxhz(I2);
            break
        end
    end
    %move down signal to the left
    for i = 1:length(newfaxhz)
        ii = amploc2 - i;
        k = ahatf(ii);
        if k < amp
            I1 = ii;
            f1 = newfaxhz(I1);
            break
        end
    end
    hpb = f2-f1;
end