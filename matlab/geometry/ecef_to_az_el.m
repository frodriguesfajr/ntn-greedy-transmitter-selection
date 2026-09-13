function [az_deg, el_deg, range_m] = ecef_to_az_el( ...
    rSatECEF_m, lat_deg, lon_deg, h_m)
%ECEF_TO_AZ_EL Compute azimuth, elevation, and range.

rRxECEF_m = geodetic_to_ecef_wgs84(lat_deg, lon_deg, h_m);

dr = rSatECEF_m(:) - rRxECEF_m(:);

lat = deg2rad(lat_deg);
lon = deg2rad(lon_deg);

R = [ ...
    -sin(lon),              cos(lon),             0;
    -sin(lat)*cos(lon),    -sin(lat)*sin(lon),   cos(lat);
     cos(lat)*cos(lon),     cos(lat)*sin(lon),   sin(lat)];

enu = R * dr;

east  = enu(1);
north = enu(2);
up    = enu(3);

az_deg = mod(atan2d(east,north),360);
el_deg = atan2d(up,hypot(east,north));

range_m = norm(dr);
end