function [interevent_variability, goodness_of_fit] = Thompstat(ETF,TTF, sigma)
[~, locs] = findpeaks(TTF);
loc1 = locs(1);
loc2 = locs(4);
interevent_variability = sigma(loc1:loc2);
interevent_variability = median(interevent_variability);
ETFcut = ETF(loc1:loc2);
TTFcut = TTF(loc1:loc2);
goodness_of_fit = corrcoef(ETFcut, TTFcut);
goodness_of_fit = goodness_of_fit(2);
end