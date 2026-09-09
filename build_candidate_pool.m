clear;
clc;

%% ============================================================
% BUILD CORRECTED 26-TRANSMITTER CANDIDATE POOL
%
% Composition is kept identical to the original paper pool:
%   4 HAPS
%   8 LEO
%   7 MEO
%   7 GEO
%
% The only change is the orbital state source:
%   OLD: MATLAB satelliteScenario SGP4
%   NEW: Vallado SGP4 / WGS-72
%
% HAPS positions are kept unchanged.
% =============================================================

rootDir = fileparts(mfilename('fullpath'));

if isempty(rootDir)
    rootDir = pwd;
end

cd(rootDir);
setup_paths;

fprintf('\n============================================================\n');
fprintf('BUILD 26-TX CANDIDATE POOL - VALLADO CORRECTION\n');
fprintf('============================================================\n');

%% ============================================================
% Old candidate pool
%
% Used here only to preserve:
%   - exact 26-Tx composition
%   - HAPS positions
%   - labels and ordering
%
% This dependency will be removed from the final standalone
% reproducibility repository after validation.
% =============================================================

oldFile = fullfile( ...
    'C:\Repository\master-thesis-ntn-positioning', ...
    'matlab_code', ...
    'results_scenario4_candidate_pool', ...
    'scenario4_candidate_pool_26.mat');

assert(isfile(oldFile), ...
    'Original candidate pool not found:\n%s',oldFile);

Old = load(oldFile);

Told = Old.Tpool;

assert(height(Told) == 26, ...
    'Original pool should contain 26 transmitters.');

%% ============================================================
% Corrected Vallado orbital catalog
% =============================================================

catalogFile = fullfile( ...
    rootDir, ...
    'results', ...
    'orbital_catalog.mat');

assert(isfile(catalogFile), ...
    'Corrected Vallado catalog not found:\n%s',catalogFile);

New = load(catalogFile);

Tcatalog = New.TOrbitalCatalog;

analysisTimeUTC = New.analysisTimeUTC;

UserPositionECEF = New.UserPositionECEF(:);

userLat_deg = New.userLat_deg;
userLon_deg = New.userLon_deg;
userH_m     = New.userH_m;

%% ============================================================
% Start from the exact original composition/order
% =============================================================

Tpool = Told;

% Normalize string variables
Tpool.Label        = string(Tpool.Label);
Tpool.Architecture = string(Tpool.Architecture);
Tpool.NORAD_CAT_ID = string(Tpool.NORAD_CAT_ID);
Tpool.ObjectName   = string(Tpool.ObjectName);

Tcatalog.Architecture = string(Tcatalog.Architecture);
Tcatalog.NORAD_CAT_ID = string(Tcatalog.NORAD_CAT_ID);

%% ============================================================
% Replace ONLY orbital ECEF states
% =============================================================

orbitalArchitectures = ["LEO","MEO","GEO"];

for k = 1:height(Tpool)

    arch = Tpool.Architecture(k);

    if any(arch == orbitalArchitectures)

        norad = Tpool.NORAD_CAT_ID(k);

        idx = ...
            Tcatalog.Architecture == arch & ...
            Tcatalog.NORAD_CAT_ID == norad;

        assert(sum(idx) == 1, ...
            'Could not uniquely match %s NORAD %s.', ...
            arch,norad);

        Tpool.X_ECEF_m(k) = Tcatalog.X_ECEF_m(idx);
        Tpool.Y_ECEF_m(k) = Tcatalog.Y_ECEF_m(idx);
        Tpool.Z_ECEF_m(k) = Tcatalog.Z_ECEF_m(idx);

    end

end

%% ============================================================
% Recompute geometry for all 26 transmitters
% =============================================================

Npool = height(Tpool);

PositionECEF = zeros(Npool,3);
uLOS         = zeros(Npool,3);

Azimuth_deg   = zeros(Npool,1);
Elevation_deg = zeros(Npool,1);
SlantRange_m  = zeros(Npool,1);
Visible       = false(Npool,1);

for k = 1:Npool

    rTx = [ ...
        Tpool.X_ECEF_m(k);
        Tpool.Y_ECEF_m(k);
        Tpool.Z_ECEF_m(k)];

    PositionECEF(k,:) = rTx.';

    [az,el,range_m,u] = lookAnglesECEF( ...
        rTx, ...
        UserPositionECEF, ...
        userLat_deg, ...
        userLon_deg);

    Azimuth_deg(k)   = az;
    Elevation_deg(k) = el;
    SlantRange_m(k)  = range_m;

    uLOS(k,:) = u.';

    Visible(k) = el >= Tpool.Mask_deg(k);

end

%% Update table geometry

Tpool.Azimuth_deg   = Azimuth_deg;
Tpool.Elevation_deg = Elevation_deg;
Tpool.SlantRange_km = SlantRange_m/1000;
Tpool.Visible       = Visible;

%% ============================================================
% Compatibility variables
% =============================================================

PoolIndex    = Tpool.PoolIndex;
Label        = Tpool.Label;
Architecture = Tpool.Architecture;
NORAD        = Tpool.NORAD_CAT_ID;
ObjectName   = Tpool.ObjectName;
Mask_deg     = Tpool.Mask_deg;

%% ============================================================
% Composition checks
% =============================================================

assert(sum(Architecture=="HAPS") == 4);
assert(sum(Architecture=="LEO")  == 8);
assert(sum(Architecture=="MEO")  == 7);
assert(sum(Architecture=="GEO")  == 7);

assert(all(Visible), ...
    'At least one transmitter in the corrected nominal pool is not visible.');

%% ============================================================
% Compare against old 26-Tx pool
% =============================================================

rOld = [ ...
    Told.X_ECEF_m, ...
    Told.Y_ECEF_m, ...
    Told.Z_ECEF_m];

rNew = PositionECEF;

positionDifference_m = sqrt(sum((rNew-rOld).^2,2));

elevationDifference_deg = ...
    Tpool.Elevation_deg - Told.Elevation_deg;

azimuthDifference_deg = mod( ...
    Tpool.Azimuth_deg - Told.Azimuth_deg + 180, ...
    360) - 180;

rangeDifference_m = ...
    1000*(Tpool.SlantRange_km - Told.SlantRange_km);

visibilityChanged = ...
    Tpool.Visible ~= Told.Visible;

%% ============================================================
% Summary
% =============================================================

fprintf('\nPool composition:\n');
fprintf('HAPS : %d\n',sum(Architecture=="HAPS"));
fprintf('LEO  : %d\n',sum(Architecture=="LEO"));
fprintf('MEO  : %d\n',sum(Architecture=="MEO"));
fprintf('GEO  : %d\n',sum(Architecture=="GEO"));

architectures = ["HAPS","LEO","MEO","GEO"];

fprintf('\n============================================================\n');
fprintf('OLD POOL vs CORRECTED VALLADO POOL\n');
fprintf('============================================================\n');

for a = 1:numel(architectures)

    arch = architectures(a);
    idx = Architecture == arch;

    fprintf('\n%s (%d)\n',arch,sum(idx));

    fprintf('Max position difference  = %.3f m\n', ...
        max(positionDifference_m(idx)));

    fprintf('Max |azimuth difference| = %.6f deg\n', ...
        max(abs(azimuthDifference_deg(idx))));

    fprintf('Max |elevation diff.|     = %.6f deg\n', ...
        max(abs(elevationDifference_deg(idx))));

    fprintf('Max |range difference|    = %.3f m\n', ...
        max(abs(rangeDifference_m(idx))));

    fprintf('Visibility changes        = %d\n', ...
        sum(visibilityChanged(idx)));

end

%% ============================================================
% Print corrected pool
% =============================================================

fprintf('\n============================================================\n');
fprintf('CORRECTED 26-TX POOL\n');
fprintf('============================================================\n');

disp(Tpool(:,{ ...
    'PoolIndex', ...
    'Label', ...
    'Architecture', ...
    'NORAD_CAT_ID', ...
    'ObjectName', ...
    'Azimuth_deg', ...
    'Elevation_deg', ...
    'SlantRange_km', ...
    'Visible'}));

%% ============================================================
% Save
% =============================================================

resultsDir = fullfile(rootDir,'results');

if ~exist(resultsDir,'dir')
    mkdir(resultsDir);
end

matFile = fullfile( ...
    resultsDir, ...
    'candidate_pool_26.mat');

csvFile = fullfile( ...
    resultsDir, ...
    'candidate_pool_26.csv');

save(matFile, ...
    'Tpool', ...
    'PoolIndex', ...
    'Label', ...
    'Architecture', ...
    'NORAD', ...
    'ObjectName', ...
    'PositionECEF', ...
    'Azimuth_deg', ...
    'Elevation_deg', ...
    'SlantRange_m', ...
    'Mask_deg', ...
    'Visible', ...
    'uLOS', ...
    'analysisTimeUTC', ...
    'UserPositionECEF', ...
    'userLat_deg', ...
    'userLon_deg', ...
    'userH_m', ...
    'positionDifference_m', ...
    'azimuthDifference_deg', ...
    'elevationDifference_deg', ...
    'rangeDifference_m', ...
    'visibilityChanged');

writetable(Tpool,csvFile);

fprintf('\nMAT file:\n%s\n',matFile);
fprintf('\nCSV file:\n%s\n',csvFile);

fprintf('\n============================================================\n');
fprintf('CORRECTED 26-TX CANDIDATE POOL COMPLETED\n');
fprintf('============================================================\n');


%% ============================================================
% Local function: ECEF -> azimuth/elevation/range/LOS
% =============================================================

function [az_deg,el_deg,range_m,uLOS] = ...
    lookAnglesECEF( ...
    rTx_m, ...
    rUser_m, ...
    userLat_deg, ...
    userLon_deg)

    rTx_m   = rTx_m(:);
    rUser_m = rUser_m(:);

    d = rTx_m-rUser_m;

    range_m = norm(d);

    if range_m <= 0
        error('Invalid transmitter/user geometry.');
    end

    uLOS = d/range_m;

    lat = userLat_deg*pi/180;
    lon = userLon_deg*pi/180;

    east = ...
        -sin(lon)*d(1) + ...
         cos(lon)*d(2);

    north = ...
        -sin(lat)*cos(lon)*d(1) ...
        -sin(lat)*sin(lon)*d(2) ...
        +cos(lat)*d(3);

    up = ...
         cos(lat)*cos(lon)*d(1) ...
        +cos(lat)*sin(lon)*d(2) ...
        +sin(lat)*d(3);

    az_deg = mod( ...
        atan2(east,north)*180/pi, ...
        360);

    el_deg = atan2( ...
        up, ...
        hypot(east,north))*180/pi;

end