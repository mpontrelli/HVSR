%Thompson et al. 2012

function [ahatf, sigma, confinthigh, confintlow] =  wavav(HV_final_matrix)
size1 = size(HV_final_matrix);
len = size1(1);
q = log(HV_final_matrix);
ahatf = exp(nansum(q)/len);

for i = 1:len
    q(i,:) = (log(HV_final_matrix(i,:))- log(ahatf)).^2;
end
sigma = sqrt(nansum(q)/len);

confinthigh = exp(log(ahatf)+1.96*sigma);
confintlow = exp(log(ahatf)-1.96*sigma);
end