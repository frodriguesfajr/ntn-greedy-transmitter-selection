clear;
clc;

thisFile = mfilename('fullpath');
geometryDir = fileparts(thisFile);

addpath(geometryDir);

%% Receiver
latRx = -22.8596582;
lonRx = -43.2303236;
hRx   = 10;

%% Synthetic HAPS geometry
azHAPS = [45 135 225 315];
elHAPS = [80 30 30 25];

hHAPS = 20000;  % m

fprintf('\n');
fprintf('============================================================\n');
fprintf('SYNTHETIC HAPS GEOMETRY\n');
fprintf('============================================================\n');

for k = 1:4

    rHAPS = haps_from_az_el( ...
        latRx,lonRx,hRx, ...
        azHAPS(k),elHAPS(k),hHAPS);

    [az,el,range_m] = ecef_to_az_el( ...
        rHAPS,latRx,lonRx,hRx);

    [lat,lon,h] = ecef_to_geodetic_wgs84(rHAPS);

    fprintf(['HAPS%d  Az=%7.3f deg  El=%7.3f deg  ' ...
             'Range=%8.3f km  Height=%8.3f km\n'], ...
        k,az,el,range_m/1000,h/1000);

end

fprintf('============================================================\n');