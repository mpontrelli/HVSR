% Compute the distance and azimuth between two points on earth given their latitude and
% longitude using the haversine formula, found on the website 
% https://www.omnicalculator.com/other/azimuth#how-to-calculate-the-azimuth-from-latitude-and-longitude
% save for future, but recognize that commands: dist = deg2km(distance(lat1,
% lon1, lat2, lon2) and az = azimuth(lat1, lon1, lat2, lon2) do the same
% exact thing, This formulation of azimuth is not yet complete, distance
% works. 

% Update - 1/13/2020 - The azimuth calculation is now correct and the whole
% code is operational. I used formulations in degrees and radians and don't
% feel like trouble shooting to do it on one or the other, but if you're
% ever bored, come back and make the code more efficient by commiting to
% radians or degrees. The results are identical to azimuth and deg2km
% matlab commands

% Author - Marshall Pontrelli
% Date - 5/30/2019
function [dist, azimuth] =  hypdist(lat1, lat2, lon1, lon2)

% First do azimuth (computation in degrees)
dlat = lat2-lat1;
dlon = lon2 - lon1;
azimuth = atan2d((sind((dlon))*cosd(lat2)), (cosd(lat1) * sind(lat2) - sind(lat1) *cosd(lat2) * cosd((dlon))));
   
% Now do distance in km (computation in radians)
lat1 = lat1 .* pi / 180;
lat2 = lat2 .* pi / 180;
lon1 = lon1 .* pi / 180;
lon2 = lon2 .* pi / 180;
dlat = lat2-lat1;
dlon = lon2 - lon1;
R = 6371; %radius of the earth in km
a = (sin(dlat./2)).^2 + cos(lat1) .* cos(lat2) .* (sin(dlon./2)).^2;
c = 2 .* asin(sqrt(a));
dist = R * c;
end