%% peakiden
% indentify significant peaks above a certain prominence threshold.
% Prominence is defined in the matlab "findpeaks" documentation.
% https://www.mathworks.com/help/signal/ref/findpeaks.html#bufbbs1-p

    % INPUTS
    
    % peakind - peak index, comes from peakiden
    
    % matrix - matrix with frequencies and amplitudes of significant peaks,
    % from peakiden
    
    % matrix1 - matrix with prominence and widths of significant peaks from
    % peakiden
    
    % ahatf1 - spectral ratio of interest, likely from peakiden
    
    % newfaxhz1 - vector of frequencies corresponding to ahatf1, likely
    % from peakiden
    
    % sigma - Vector of standard deviations, likely from wavav
    
    % statname - station name, this gets put in the first column
    
    % lowbound - lowest frequency of interest
    
    % OUTPUTS
    
    % taxstat - a matrix containing information about all the significant peaks 
    % with: column 1) name of the station, column 2) the frequency of the
    % significant peak, column 3) the amplitude of the significant peak,
    % column 4) the halfpower bandwidth of the significant peak, 5) the
    % interevent variability of each peak across the halfpower bandwidth,
    % 6) the prominence of the significant peak, anf 7) the width of the
    % significant peak. The last two have definitions in the matlab
    % documentation for findpeaks. 
    

%% Author: Marshall Pontrelli
% Date: developed during summer 2019 
% Update - 1/16/2020 added comments while I was coming back to debug this
% and peakiden - Marshall

function [taxstat] = specratstat(peakind, matrix, matrix1, ahatf1, newfaxhz, newfaxhz1, sigma, statname,lowbound, upbound)
[m] = length(peakind);
sigma1 = sigma(find(newfaxhz == lowbound): find(newfaxhz == upbound)); 
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