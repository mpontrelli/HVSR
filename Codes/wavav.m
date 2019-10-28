%% wavav
% Compute the lognormal average, standard deviation and 95% confidence interval
% of a spectral ratio of magnitude responses at each frequency
% using equations 2-4 from Thompson at al. 2012 which come from Steidl et
% al. 1994. 

    % INPUTS
    
    % SR_matrix - Matix containing a spectral ratio in each row. Each
    % column must have the amplification value correponding to the same
    % frequency at each spectral ratio.
    
    % OUTPUTS
    
    % ahatf - lognormal average spectral ratio
    
    % sigma - standard deviation
    
    % confinthigh - upper 95% confidence interval
    
    %confintlow - lower 95% confidence interval

%% Author: Marshall Pontrelli
% Date: developed between September, 2017 and August, 2019   

%% Start
function [ahatf, sigma, confinthigh, confintlow] =  wavav(SR_matrix)
size1 = size(SR_final_matrix);
len = size1(1);
q = log(SR_final_matrix);
ahatf = exp(nansum(q)/len);

for i = 1:len
    q(i,:) = (log(HV_final_matrix(i,:))- log(ahatf)).^2;
end
sigma = sqrt(nansum(q)/len);

confinthigh = exp(log(ahatf)+1.96*sigma);
confintlow = exp(log(ahatf)-1.96*sigma);
end