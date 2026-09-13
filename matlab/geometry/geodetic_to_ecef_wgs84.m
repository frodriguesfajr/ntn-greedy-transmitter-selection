function rECEF_m = geodetic_to_ecef_wgs84(lat_deg, lon_deg, h_m)
%GEODETIC_TO_ECEF_WGS84 Convert geodetic coordinates to WGS-84 ECEF.

a = 6378137.0;                 % semi-major axis [m]
f = 1 / 298.257223563;
e2 = f * (2 - f);

lat = deg2rad(lat_deg);
lon = deg2rad(lon_deg);

N = a / sqrt(1 - e2*sin(lat)^2);

x = (N + h_m) * cos(lat) * cos(lon);
y = (N + h_m) * cos(lat) * sin(lon);
z = (N*(1 - e2) + h_m) * sin(lat);

rECEF_m = [x; y; z];
end