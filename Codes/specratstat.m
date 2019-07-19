function [taxstat] = specratstat(peakind,matrix, ahatf1, newfaxhz1, sigma, statname,lowbound)
[m] = length(peakind);
sigma1 = sigma(lowbound:length(sigma)); 
sigvec = []; 
for f = 1:m
    loc = peakind(f);
    disp(loc)
    A = matrix(f,2);
    disp(A)
    [I1, I, f1, f2, hpb] =  HalfPowerBand2(A, loc, newfaxhz1, ahatf1);
    taxstat{f,1} = statname;
    taxstat{f,2} = num2str(hpb);
    a = sigma1(I1:I);
    sigmai = median(a);
    taxstat{f,3} = num2str(sigmai);
end
matrix = num2cell(matrix);
taxstat=[taxstat,matrix];
end