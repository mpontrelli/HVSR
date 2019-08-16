function [taxstat] = specratstat(peakind, matrix, matrix1, ahatf1, newfaxhz1, sigma, statname,lowbound)
[m] = length(peakind);
sigma1 = sigma(lowbound:length(sigma));  
for f = 1:m
    loc = peakind(f);
    A = matrix(f,2);
    [I1, I, f1, f2, hpb] =  HalfPowerBand2(A, loc, newfaxhz1, ahatf1);
    taxstat{f,1} = statname;
    taxstat1(f,1) = hpb;
    a = sigma1(I1:I);
    sigmai = median(a);
    taxstat1(f,2) = sigmai;
end
matrix = num2cell(matrix);
matrix1 = num2cell(matrix1);
taxstat1 = num2cell(taxstat1);
taxstat=[taxstat, matrix, taxstat1, matrix1];
end