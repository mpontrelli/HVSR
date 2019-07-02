function [PGANS,PGAV,PGAEW] = PGA(xNS,xV,xEW)
PGANS = max(abs(xNS))/981;
PGAV = max(abs(xV))/981;
PGAEW = max(abs(xEW)/981);
end 