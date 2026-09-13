function rHAPS_ECEF_m = haps_from_az_el( ...
    latRx_deg, lonRx_deg, hRx_m, ...
    az_deg, el_deg, hHAPS_m)
%HAPS_FROM_AZ_EL Generate HAPS ECEF position from receiver az/el.
%
% The LOS is extended until the WGS-84 geodetic height equals hHAPS_m.

rRx = geodetic_to_ecef_wgs84( ...
    latRx_deg,lonRx_deg,hRx_m);

lat = deg2rad(latRx_deg);
lon = deg2rad(lonRx_deg);
az  = deg2rad(az_deg);
el  = deg2rad(el_deg);

%% LOS unit vector in ENU
uENU = [ ...
    sin(az)*cos(el);
    cos(az)*cos(el);
    sin(el)];

%% ECEF -> ENU rotation
R = [ ...
    -sin(lon),              cos(lon),             0;
    -sin(lat)*cos(lon),    -sin(lat)*sin(lon),   cos(lat);
     cos(lat)*cos(lon),     cos(lat)*sin(lon),   sin(lat)];

%% ENU -> ECEF
uECEF = R.' * uENU;
uECEF = uECEF / norm(uECEF);

%% Initial flat-Earth slant-range estimate
s0 = (hHAPS_m-hRx_m)/sin(el);

heightError = @(s) localHeight( ...
    rRx + s*uECEF) - hHAPS_m;

%% Find an upper bracket
sUpper = 2*s0;

while heightError(sUpper) < 0
    sUpper = 2*sUpper;
end

%% Solve for the LOS range giving the requested HAPS height
s = fzero(heightError,[0 sUpper]);

rHAPS_ECEF_m = rRx + s*uECEF;

end


function h = localHeight(rECEF)
[~,~,h] = ecef_to_geodetic_wgs84(rECEF);
end