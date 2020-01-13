%Compute the distance and azimuth between two points on earth given their latitude and
%longitude using the haversine formula, found on the website 
%https://www.omnicalculator.com/other/azimuth#how-to-calculate-the-azimuth-from-latitude-and-longitude
%save for future, but recognize that commands: dist = deg2km(distance(lat1,
%lon1, lat2, lon2) and az = azimuth(lat1, lon1, lat2, lon2) do the same
%exact thing, This formulation of azimuth is not yet complete, distance
%works.
function [dist, azimuth] =  hypdist(lat1, lat2, lon1, lon2)
    lat1 = lat1 .* pi / 180;
    lat2 = lat2 .* pi / 180;
    lon1 = lon1 .* pi / 180;
    lon2 = lon2 .* pi / 180;
    dlat = lat2-lat1;
    dlon = lon2 - lon1;
    R = 6371; %radius of the earth in km
    a = (sin(dlat./2)).^2 + cos(lat1) .* cos(lat2) .* (sin(dlon./2)).^2;
    c = 2 .* asin(sqrt(a));
%     a = (sind((lat2-lat1))).^2 + cosd(lat1) * cosd(lat2) * (sind((lon2-lon1))).^2;
%     c = 2 * atan2d(sqrt(a),sqrt(1-a));
    dist = R * c;
%     
    azimuth = atan2d((sind((lon2-lon1))*cosd(lat2)), (cosd(lat1 * sind(lat2) - sind(lat1) *cosd(lat2) * cosd((lon2-lon1)))));
    

end