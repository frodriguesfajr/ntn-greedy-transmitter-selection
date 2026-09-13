function [lat_deg, lon_deg, h_m] = ecef_to_geodetic_wgs84(rECEF_m)
%ECEF_TO_GEODETIC_WGS84 Convert ECEF to WGS-84 geodetic coordinates.

a = 6378137.0;
f = 1 / 298.257223563;
e2 = f * (2 - f);

x = rECEF_m(1);
y = rECEF_m(2);
z = rECEF_m(3);

lon = atan2(y,x);
p = hypot(x,y);

lat = atan2(z,p*(1-e2));

for k = 1:10
    N = a / sqrt(1-e2*sin(lat)^2);
    h = p/cos(lat) - N;

    latNew = atan2(z, ...
        p*(1 - e2*N/(N+h)));

    if abs(latNew-lat) < 1e-13
        lat = latNew;
        break;
    end

    lat = latNew;
end

N = a / sqrt(1-e2*sin(lat)^2);
h_m = p/cos(lat) - N;

lat_deg = rad2deg(lat);
lon_deg = rad2deg(lon);

end