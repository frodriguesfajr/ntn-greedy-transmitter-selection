clear;
clc;

%% ============================================================
%  BUILD ORBITAL CATALOG - VALLADO SGP4
%
%  Corrected orbital propagation for the paper
%
%  Epoch:
%    01-Aug-2026 12:00:00 UTC
%
%  Orbital catalog:
%    20 LEO
%     7 MEO
%     7 GEO
%
%  Propagation:
%    Vallado SGP4
%    WGS-72
%
%  TEME -> ECEF:
%    Vallado transformation with fixed IERS EOP
%
%  No Aerospace Toolbox or Satellite Communications Toolbox
%  required.
% =============================================================

%% Repository root

rootDir = fileparts(mfilename('fullpath'));

if isempty(rootDir)
    rootDir = pwd;
end

cd(rootDir);

fprintf('\n============================================================\n');
fprintf('BUILD ORBITAL CATALOG - VALLADO SGP4\n');
fprintf('============================================================\n');
fprintf('Repository:\n%s\n',rootDir);

%% Configure paths

setup_paths;

assert(exist('propagate_tle_vallado','file') == 2, ...
    'propagate_tle_vallado.m was not found.');

assert(exist('paper_epoch_eop','file') == 2, ...
    'paper_epoch_eop.m was not found.');

assert(exist('sgp4','file') == 2, ...
    'Vallado sgp4.m was not found.');

fprintf('\nSGP4 implementation:\n%s\n',which('sgp4'));

%% ============================================================
% Analysis epoch
% =============================================================

analysisTimeUTC = datetime( ...
    2026,8,1,12,0,0, ...
    'TimeZone','UTC');

fprintf('\nAnalysis epoch : %s\n',string(analysisTimeUTC));

%% ============================================================
% Earth Orientation Parameters
% =============================================================

eop = paper_epoch_eop();

fprintf('\nEarth orientation parameters\n');
fprintf('xp       = %.7f arcsec\n',eop.xp_arcsec);
fprintf('yp       = %.7f arcsec\n',eop.yp_arcsec);
fprintf('UT1-UTC  = %.7f s\n',eop.dut1_s);
fprintf('LOD      = %.8f s\n',eop.lod_s);

%% ============================================================
% User position
% =============================================================

userLat_deg = -22.8596582;
userLon_deg = -43.2303236;
userH_m     = 10.0;

UserPositionECEF = geodeticToECEFWGS84( ...
    userLat_deg, ...
    userLon_deg, ...
    userH_m);

fprintf('\nUser position\n');
fprintf('Latitude  = %.7f deg\n',userLat_deg);
fprintf('Longitude = %.7f deg\n',userLon_deg);
fprintf('Height    = %.3f m\n',userH_m);

fprintf('\nUser ECEF [m]\n');
fprintf('X = %.3f\n',UserPositionECEF(1));
fprintf('Y = %.3f\n',UserPositionECEF(2));
fprintf('Z = %.3f\n',UserPositionECEF(3));

%% ============================================================
% Catalog configuration
% =============================================================

architectures = ["LEO","MEO","GEO"];

% Orbital transmitter elevation mask used in the paper
orbitalMask_deg = 5.0;

expectedCounts.LEO = 20;
expectedCounts.MEO = 7;
expectedCounts.GEO = 7;

%% ============================================================
% Output arrays
% =============================================================

Architecture = strings(0,1);
NORAD_CAT_ID = strings(0,1);
ObjectName   = strings(0,1);
TLEFile      = strings(0,1);

X_TEME_km = zeros(0,1);
Y_TEME_km = zeros(0,1);
Z_TEME_km = zeros(0,1);

VX_TEME_kms = zeros(0,1);
VY_TEME_kms = zeros(0,1);
VZ_TEME_kms = zeros(0,1);

X_ECEF_m = zeros(0,1);
Y_ECEF_m = zeros(0,1);
Z_ECEF_m = zeros(0,1);

VX_ECEF_mps = zeros(0,1);
VY_ECEF_mps = zeros(0,1);
VZ_ECEF_mps = zeros(0,1);

Azimuth_deg   = zeros(0,1);
Elevation_deg = zeros(0,1);
SlantRange_km = zeros(0,1);

Mask_deg = zeros(0,1);
Visible  = false(0,1);

uX = zeros(0,1);
uY = zeros(0,1);
uZ = zeros(0,1);

TLEAge_h = zeros(0,1);
SGP4Error = zeros(0,1);

%% ============================================================
% Propagate all orbital transmitters
% =============================================================

for a = 1:numel(architectures)

    arch = architectures(a);

    archDir = fullfile( ...
        rootDir, ...
        'tle', ...
        char(arch));

    assert(isfolder(archDir), ...
        'TLE directory not found: %s',archDir);

    satDirs = dir(archDir);
    satDirs = satDirs([satDirs.isdir]);

    folderNames = string({satDirs.name});

    % Keep only folders whose names are numeric NORAD IDs
    noradNumeric = str2double(folderNames);

    valid = ...
        ~ismember(folderNames,[".",".."]) & ...
        ~isnan(noradNumeric);

    satDirs = satDirs(valid);
    noradNumeric = noradNumeric(valid);

    % Deterministic ordering
    [~,order] = sort(noradNumeric);
    satDirs = satDirs(order);

    fprintf('\n------------------------------------------------------------\n');
    fprintf('%s\n',arch);
    fprintf('------------------------------------------------------------\n');
    fprintf('TLE files found: %d\n',numel(satDirs));

    switch arch

        case "LEO"
            assert(numel(satDirs) == expectedCounts.LEO, ...
                'Expected %d LEO satellites, found %d.', ...
                expectedCounts.LEO,numel(satDirs));

        case "MEO"
            assert(numel(satDirs) == expectedCounts.MEO, ...
                'Expected %d MEO satellites, found %d.', ...
                expectedCounts.MEO,numel(satDirs));

        case "GEO"
            assert(numel(satDirs) == expectedCounts.GEO, ...
                'Expected %d GEO satellites, found %d.', ...
                expectedCounts.GEO,numel(satDirs));

    end

    for k = 1:numel(satDirs)

        noradText = string(satDirs(k).name);

        tleFile = fullfile( ...
            archDir, ...
            char(noradText), ...
            [char(noradText) '.tle']);

        assert(isfile(tleFile), ...
            'TLE file not found: %s',tleFile);

        %% ----------------------------------------------------
        % Vallado SGP4
        % -----------------------------------------------------

        V = propagate_tle_vallado( ...
            tleFile, ...
            analysisTimeUTC, ...
            eop);

        %% ----------------------------------------------------
        % Geometry relative to user
        % -----------------------------------------------------

        [az_deg,el_deg,range_m,uLOS] = ...
            lookAnglesECEF( ...
            V.rECEF_m, ...
            UserPositionECEF, ...
            userLat_deg, ...
            userLon_deg);

        %% ----------------------------------------------------
        % Store result
        % -----------------------------------------------------

        Architecture(end+1,1) = arch;
        NORAD_CAT_ID(end+1,1) = noradText;
        ObjectName(end+1,1)   = string(V.name);

        % Save relative path so repository is portable
        TLEFile(end+1,1) = string(fullfile( ...
            'tle', ...
            char(arch), ...
            char(noradText), ...
            [char(noradText) '.tle']));

        % TEME
        X_TEME_km(end+1,1) = V.rTEME_km(1);
        Y_TEME_km(end+1,1) = V.rTEME_km(2);
        Z_TEME_km(end+1,1) = V.rTEME_km(3);

        VX_TEME_kms(end+1,1) = V.vTEME_kms(1);
        VY_TEME_kms(end+1,1) = V.vTEME_kms(2);
        VZ_TEME_kms(end+1,1) = V.vTEME_kms(3);

        % ECEF
        X_ECEF_m(end+1,1) = V.rECEF_m(1);
        Y_ECEF_m(end+1,1) = V.rECEF_m(2);
        Z_ECEF_m(end+1,1) = V.rECEF_m(3);

        VX_ECEF_mps(end+1,1) = V.vECEF_mps(1);
        VY_ECEF_mps(end+1,1) = V.vECEF_mps(2);
        VZ_ECEF_mps(end+1,1) = V.vECEF_mps(3);

        % User geometry
        Azimuth_deg(end+1,1)   = az_deg;
        Elevation_deg(end+1,1) = el_deg;
        SlantRange_km(end+1,1) = range_m/1000;

        Mask_deg(end+1,1) = orbitalMask_deg;
        Visible(end+1,1)  = el_deg >= orbitalMask_deg;

        uX(end+1,1) = uLOS(1);
        uY(end+1,1) = uLOS(2);
        uZ(end+1,1) = uLOS(3);

        TLEAge_h(end+1,1) = V.tsince_min/60;
        SGP4Error(end+1,1) = V.sgp4_error;

        fprintf( ...
            'NORAD %-6s  %-24s  El = %+8.3f deg  Visible = %d\n', ...
            char(noradText), ...
            char(string(V.name)), ...
            el_deg, ...
            el_deg >= orbitalMask_deg);

    end
end

%% ============================================================
% Complete orbital catalog
% =============================================================

CatalogIndex = (1:numel(Architecture)).';

TOrbitalCatalog = table( ...
    CatalogIndex, ...
    Architecture, ...
    NORAD_CAT_ID, ...
    ObjectName, ...
    TLEFile, ...
    X_TEME_km, ...
    Y_TEME_km, ...
    Z_TEME_km, ...
    VX_TEME_kms, ...
    VY_TEME_kms, ...
    VZ_TEME_kms, ...
    X_ECEF_m, ...
    Y_ECEF_m, ...
    Z_ECEF_m, ...
    VX_ECEF_mps, ...
    VY_ECEF_mps, ...
    VZ_ECEF_mps, ...
    Azimuth_deg, ...
    Elevation_deg, ...
    SlantRange_km, ...
    Mask_deg, ...
    Visible, ...
    uX,uY,uZ, ...
    TLEAge_h, ...
    SGP4Error);

%% Architecture-specific tables

TLEO20 = TOrbitalCatalog( ...
    TOrbitalCatalog.Architecture=="LEO",:);

TMEO7 = TOrbitalCatalog( ...
    TOrbitalCatalog.Architecture=="MEO",:);

TGEO7 = TOrbitalCatalog( ...
    TOrbitalCatalog.Architecture=="GEO",:);

%% ============================================================
% Sanity checks
% =============================================================

assert(height(TOrbitalCatalog) == 34, ...
    'Expected 34 orbital objects, found %d.', ...
    height(TOrbitalCatalog));

assert(height(TLEO20) == 20, ...
    'Expected 20 LEO objects.');

assert(height(TMEO7) == 7, ...
    'Expected 7 MEO objects.');

assert(height(TGEO7) == 7, ...
    'Expected 7 GEO objects.');

assert(all(TOrbitalCatalog.SGP4Error == 0), ...
    'At least one SGP4 propagation returned an error.');

assert(all(isfinite(TOrbitalCatalog.X_ECEF_m)));
assert(all(isfinite(TOrbitalCatalog.Y_ECEF_m)));
assert(all(isfinite(TOrbitalCatalog.Z_ECEF_m)));

%% ============================================================
% Save results
% =============================================================

resultsDir = fullfile(rootDir,'results');

if ~exist(resultsDir,'dir')
    mkdir(resultsDir);
end

matFile = fullfile( ...
    resultsDir, ...
    'orbital_catalog.mat');

csvFile = fullfile( ...
    resultsDir, ...
    'orbital_catalog.csv');

save(matFile, ...
    'TOrbitalCatalog', ...
    'TLEO20', ...
    'TMEO7', ...
    'TGEO7', ...
    'analysisTimeUTC', ...
    'UserPositionECEF', ...
    'userLat_deg', ...
    'userLon_deg', ...
    'userH_m', ...
    'orbitalMask_deg', ...
    'eop');

writetable(TOrbitalCatalog,csvFile);

%% ============================================================
% Summary
% =============================================================

fprintf('\n============================================================\n');
fprintf('CATALOG SUMMARY\n');
fprintf('============================================================\n');

for a = 1:numel(architectures)

    arch = architectures(a);

    idx = TOrbitalCatalog.Architecture == arch;

    fprintf('\n%s\n',arch);

    fprintf('Total   : %d\n',sum(idx));

    fprintf('Visible : %d\n', ...
        sum(TOrbitalCatalog.Visible(idx)));

    fprintf('El min  : %.3f deg\n', ...
        min(TOrbitalCatalog.Elevation_deg(idx)));

    fprintf('El max  : %.3f deg\n', ...
        max(TOrbitalCatalog.Elevation_deg(idx)));

    fprintf('TLE age min/max : %.3f / %.3f h\n', ...
        min(TOrbitalCatalog.TLEAge_h(idx)), ...
        max(TOrbitalCatalog.TLEAge_h(idx)));

end

fprintf('\nTOTAL ORBITAL OBJECTS : %d\n', ...
    height(TOrbitalCatalog));

fprintf('TOTAL VISIBLE         : %d\n', ...
    sum(TOrbitalCatalog.Visible));

fprintf('\nMAT file:\n%s\n',matFile);

fprintf('\nCSV file:\n%s\n',csvFile);

fprintf('\n============================================================\n');
fprintf('CORRECTED ORBITAL CATALOG COMPLETED\n');
fprintf('============================================================\n');


%% ============================================================
% Local function: geodetic -> ECEF, WGS-84
% =============================================================

function rECEF_m = ...
    geodeticToECEFWGS84(lat_deg,lon_deg,h_m)

    a = 6378137.0;
    f = 1/298.257223563;

    e2 = f*(2-f);

    lat = lat_deg*pi/180;
    lon = lon_deg*pi/180;

    sinLat = sin(lat);
    cosLat = cos(lat);

    N = a/sqrt(1-e2*sinLat^2);

    x = (N+h_m)*cosLat*cos(lon);
    y = (N+h_m)*cosLat*sin(lon);
    z = (N*(1-e2)+h_m)*sinLat;

    rECEF_m = [x;y;z];

end


%% ============================================================
% Local function: ECEF geometry
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