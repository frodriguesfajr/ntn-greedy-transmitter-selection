clear;
clc;

%% Repository paths
thisFile = mfilename('fullpath');
geometryDir = fileparts(thisFile);
matlabDir = fileparts(geometryDir);
repoRoot = fileparts(matlabDir);

addpath(fullfile(repoRoot,'external','sgp4'));
addpath(fullfile(repoRoot,'matlab','propagation'));
addpath(geometryDir);

%% Receiver
latRx_deg = -22.8596582;
lonRx_deg = -43.2303236;
hRx_m     = 10;

%% Reference epoch
epochUTC = datetime(2026,8,1,12,0,0, ...
    'TimeZone','UTC');

%% Nominal orbital pool
leo = [ ...
    65175 66864 65348 59340 ...
    67736 59960 65279 48482 ];

meo = [ ...
    62362 40081 44115 43231 ...
    43234 40080 40079 ];

geo = [ ...
    41589 41904 42692 43562 ...
    43175 43228 38087 ];

idsByType = {leo, meo, geo};
types = {'LEO','MEO','GEO'};

mask_deg = 5;

fprintf('\n');
fprintf('============================================================\n');
fprintf('NOMINAL ORBITAL VISIBILITY\n');
fprintf('Epoch: 01-Aug-2026 12:00:00 UTC\n');
fprintf('Receiver: %.7f deg, %.7f deg, %.1f m\n', ...
    latRx_deg,lonRx_deg,hRx_m);
fprintf('Elevation mask: %.1f deg\n',mask_deg);
fprintf('============================================================\n');

nVisible = 0;
nTotal   = 0;

for g = 1:numel(types)

    ids = idsByType{g};

    fprintf('\n%s\n',types{g});
    fprintf('------------------------------------------------------------\n');

    for k = 1:numel(ids)

        norad = ids(k);

        tleFile = fullfile( ...
            repoRoot,'tle',lower(types{g}), ...
            sprintf('%d',norad), ...
            sprintf('%d.tle',norad));

        out = propagate_tle_to_ecef(tleFile,epochUTC);

        [az,el,range_m] = ecef_to_az_el( ...
            1000*out.rECEF_km, ...
            latRx_deg,lonRx_deg,hRx_m);

        visible = el >= mask_deg;

        fprintf('%5d  Az=%7.2f deg  El=%7.2f deg  R=%9.1f km  %s\n', ...
            norad,az,el,range_m/1000, ...
            string(visible));

        nVisible = nVisible + visible;
        nTotal = nTotal + 1;
    end
end

fprintf('\n============================================================\n');
fprintf('Visible satellites: %d / %d\n',nVisible,nTotal);
fprintf('============================================================\n');