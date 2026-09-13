clear;
clc;
format long;

%% ============================================================
% BUILD ORBITAL CATALOG - VALLADO SGP4
%
% Canonical location:
%   matlab/setup/build_orbital_catalog.m
%
% Output:
%   results/orbital_catalog.mat
%   results/orbital_catalog.csv
%
% Required repository structure:
%   tle/LEO/<NORAD>/<NORAD>.tle
%   tle/MEO/<NORAD>/<NORAD>.tle
%   tle/GEO/<NORAD>/<NORAD>.tle
%
% Propagation:
%   Vallado SGP4 / WGS-72
%
% No Aerospace Toolbox or Satellite Communications Toolbox required.
% ============================================================

%% Repository root

thisFile = mfilename('fullpath');

if isempty(thisFile)
    error('This script must be executed from a saved .m file.');
end

currentDir = fileparts(thisFile);
repoRoot = currentDir;

while ~isfolder(fullfile(repoRoot,'.git'))

    parentDir = fileparts(repoRoot);

    if strcmp(parentDir,repoRoot)
        error('Could not locate Git repository root from:\n%s',thisFile);
    end

    repoRoot = parentDir;
end

%% Configure paths

setupFile = fullfile( ...
    repoRoot,'matlab','setup','setup_paths.m');

assert(isfile(setupFile), ...
    'setup_paths.m not found:\n%s',setupFile);

addpath(fullfile(repoRoot,'matlab','setup'),'-begin');
setup_paths;

assert(exist('propagate_tle_vallado','file')==2, ...
    'propagate_tle_vallado.m was not found.');

assert(exist('paper_epoch_eop','file')==2, ...
    'paper_epoch_eop.m was not found.');

assert(exist('sgp4','file')==2, ...
    'Vallado sgp4.m was not found.');

fprintf('\n============================================================\n');
fprintf('BUILD ORBITAL CATALOG\n');
fprintf('============================================================\n');
fprintf('Repository:\n%s\n',repoRoot);
fprintf('\nSGP4 implementation:\n%s\n',which('sgp4'));

%% Analysis epoch

analysisTimeUTC = datetime( ...
    2026,8,1,12,0,0, ...
    'TimeZone','UTC');

fprintf('\nAnalysis epoch : %s\n',string(analysisTimeUTC));

%% Earth orientation parameters

eop = paper_epoch_eop();

fprintf('\nEarth orientation parameters\n');
fprintf('xp       = %.7f arcsec\n',eop.xp_arcsec);
fprintf('yp       = %.7f arcsec\n',eop.yp_arcsec);
fprintf('UT1-UTC  = %.7f s\n',eop.dut1_s);
fprintf('LOD      = %.8f s\n',eop.lod_s);

%% User position

userLat_deg = -22.8596582;
userLon_deg = -43.2303236;
userH_m = 10.0;

UserPositionECEF = ...
    geodeticToECEFWGS84( ...
    userLat_deg,userLon_deg,userH_m);

fprintf('\nUser position\n');
fprintf('Latitude  = %.7f deg\n',userLat_deg);
fprintf('Longitude = %.7f deg\n',userLon_deg);
fprintf('Height    = %.3f m\n',userH_m);

fprintf('\nUser ECEF [m]\n');
fprintf('X = %.3f\n',UserPositionECEF(1));
fprintf('Y = %.3f\n',UserPositionECEF(2));
fprintf('Z = %.3f\n',UserPositionECEF(3));

%% Catalog configuration

architectures = ["LEO","MEO","GEO"];

expectedCounts.LEO = 20;
expectedCounts.MEO = 7;
expectedCounts.GEO = 7;

orbitalMask_deg = 5.0;

tleRoot = fullfile(repoRoot,'tle');

assert(isfolder(tleRoot), ...
    'TLE root not found:\n%s',tleRoot);

%% Output arrays

Architecture = strings(0,1);
NORAD_CAT_ID = strings(0,1);
ObjectName = strings(0,1);
TLEFile = strings(0,1);

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

Azimuth_deg = zeros(0,1);
Elevation_deg = zeros(0,1);
SlantRange_km = zeros(0,1);

Mask_deg = zeros(0,1);
Visible = false(0,1);

uX = zeros(0,1);
uY = zeros(0,1);
uZ = zeros(0,1);

TLEAge_h = zeros(0,1);
SGP4Error = zeros(0,1);

%% Propagate all orbital transmitters

for a = 1:numel(architectures)

    arch = architectures(a);

    archDir = fullfile( ...
        tleRoot,char(arch));

    assert(isfolder(archDir), ...
        'TLE directory not found: %s',archDir);

    satDirs = dir(archDir);
    satDirs = satDirs([satDirs.isdir]);

    folderNames = string({satDirs.name});

    noradNumeric = str2double(folderNames);

    valid = ...
        ~ismember(folderNames,[".",".."]) & ...
        ~isnan(noradNumeric);

    satDirs = satDirs(valid);
    noradNumeric = noradNumeric(valid);

    [~,order] = sort(noradNumeric);
    satDirs = satDirs(order);

    switch arch
        case "LEO"
            expectedCount = expectedCounts.LEO;
        case "MEO"
            expectedCount = expectedCounts.MEO;
        case "GEO"
            expectedCount = expectedCounts.GEO;
        otherwise
            error('Unsupported architecture: %s',arch);
    end

    assert(numel(satDirs)==expectedCount, ...
        'Expected %d %s satellites, found %d.', ...
        expectedCount,arch,numel(satDirs));

    fprintf('\n------------------------------------------------------------\n');
    fprintf('%s\n',arch);
    fprintf('------------------------------------------------------------\n');
    fprintf('TLE files found: %d\n',numel(satDirs));

    for k = 1:numel(satDirs)

        noradText = string(satDirs(k).name);

        tleFile = fullfile( ...
            archDir, ...
            char(noradText), ...
            [char(noradText) '.tle']);

        assert(isfile(tleFile), ...
            'TLE file not found: %s',tleFile);

        V = propagate_tle_vallado( ...
            tleFile,analysisTimeUTC,eop);

        assert(isfield(V,'rECEF_m'), ...
            'Propagation result does not contain rECEF_m.');

        assert(isfield(V,'rTEME_km'), ...
            'Propagation result does not contain rTEME_km.');

        if isfield(V,'sgp4_error')
            sgp4Error = double(V.sgp4_error);
        else
            sgp4Error = 0;
        end

        assert(sgp4Error==0, ...
            'SGP4 error %d for NORAD %s.', ...
            sgp4Error,noradText);

        [az_deg,el_deg,range_m,uLOS] = ...
            lookAnglesECEF( ...
            V.rECEF_m, ...
            UserPositionECEF, ...
            userLat_deg, ...
            userLon_deg);

        Architecture(end+1,1) = arch; %#ok<SAGROW>
        NORAD_CAT_ID(end+1,1) = noradText; %#ok<SAGROW>

        if isfield(V,'name')
            ObjectName(end+1,1) = string(V.name); %#ok<SAGROW>
        else
            ObjectName(end+1,1) = ""; %#ok<SAGROW>
        end

        TLEFile(end+1,1) = string(fullfile( ...
            'tle',char(arch),char(noradText), ...
            [char(noradText) '.tle'])); %#ok<SAGROW>

        X_TEME_km(end+1,1) = V.rTEME_km(1); %#ok<SAGROW>
        Y_TEME_km(end+1,1) = V.rTEME_km(2); %#ok<SAGROW>
        Z_TEME_km(end+1,1) = V.rTEME_km(3); %#ok<SAGROW>

        if isfield(V,'vTEME_kms')
            VX_TEME_kms(end+1,1) = V.vTEME_kms(1); %#ok<SAGROW>
            VY_TEME_kms(end+1,1) = V.vTEME_kms(2); %#ok<SAGROW>
            VZ_TEME_kms(end+1,1) = V.vTEME_kms(3); %#ok<SAGROW>
        else
            VX_TEME_kms(end+1,1) = NaN; %#ok<SAGROW>
            VY_TEME_kms(end+1,1) = NaN; %#ok<SAGROW>
            VZ_TEME_kms(end+1,1) = NaN; %#ok<SAGROW>
        end

        X_ECEF_m(end+1,1) = V.rECEF_m(1); %#ok<SAGROW>
        Y_ECEF_m(end+1,1) = V.rECEF_m(2); %#ok<SAGROW>
        Z_ECEF_m(end+1,1) = V.rECEF_m(3); %#ok<SAGROW>

        if isfield(V,'vECEF_mps')
            VX_ECEF_mps(end+1,1) = V.vECEF_mps(1); %#ok<SAGROW>
            VY_ECEF_mps(end+1,1) = V.vECEF_mps(2); %#ok<SAGROW>
            VZ_ECEF_mps(end+1,1) = V.vECEF_mps(3); %#ok<SAGROW>
        else
            VX_ECEF_mps(end+1,1) = NaN; %#ok<SAGROW>
            VY_ECEF_mps(end+1,1) = NaN; %#ok<SAGROW>
            VZ_ECEF_mps(end+1,1) = NaN; %#ok<SAGROW>
        end

        Azimuth_deg(end+1,1) = az_deg; %#ok<SAGROW>
        Elevation_deg(end+1,1) = el_deg; %#ok<SAGROW>
        SlantRange_km(end+1,1) = range_m/1000; %#ok<SAGROW>

        Mask_deg(end+1,1) = orbitalMask_deg; %#ok<SAGROW>
        Visible(end+1,1) = el_deg>=orbitalMask_deg; %#ok<SAGROW>

        uX(end+1,1) = uLOS(1); %#ok<SAGROW>
        uY(end+1,1) = uLOS(2); %#ok<SAGROW>
        uZ(end+1,1) = uLOS(3); %#ok<SAGROW>

        if isfield(V,'tsince_min')
            TLEAge_h(end+1,1) = V.tsince_min/60; %#ok<SAGROW>
        else
            TLEAge_h(end+1,1) = NaN; %#ok<SAGROW>
        end

        SGP4Error(end+1,1) = sgp4Error; %#ok<SAGROW>

        fprintf( ...
            'NORAD %-6s  %-24s  El = %+8.3f deg  Visible = %d\n', ...
            char(noradText), ...
            char(ObjectName(end)), ...
            el_deg, ...
            Visible(end));
    end
end

%% Build table

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

TLEO20 = ...
    TOrbitalCatalog(TOrbitalCatalog.Architecture=="LEO",:);

TMEO7 = ...
    TOrbitalCatalog(TOrbitalCatalog.Architecture=="MEO",:);

TGEO7 = ...
    TOrbitalCatalog(TOrbitalCatalog.Architecture=="GEO",:);

%% Sanity checks

assert(height(TOrbitalCatalog)==34, ...
    'Expected 34 orbital objects, found %d.', ...
    height(TOrbitalCatalog));

assert(height(TLEO20)==20);
assert(height(TMEO7)==7);
assert(height(TGEO7)==7);

assert(all(TOrbitalCatalog.SGP4Error==0));

assert(all(isfinite(TOrbitalCatalog.X_ECEF_m)));
assert(all(isfinite(TOrbitalCatalog.Y_ECEF_m)));
assert(all(isfinite(TOrbitalCatalog.Z_ECEF_m)));

%% Save

resultsDir = fullfile(repoRoot,'results');

if ~isfolder(resultsDir)
    mkdir(resultsDir);
end

matFile = fullfile( ...
    resultsDir,'orbital_catalog.mat');

csvFile = fullfile( ...
    resultsDir,'orbital_catalog.csv');

save( ...
    matFile, ...
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

%% Summary

fprintf('\n============================================================\n');
fprintf('CATALOG SUMMARY\n');
fprintf('============================================================\n');

for a = 1:numel(architectures)

    arch = architectures(a);
    idx = TOrbitalCatalog.Architecture==arch;

    fprintf('\n%s\n',arch);
    fprintf('Total   : %d\n',sum(idx));
    fprintf('Visible : %d\n',sum(TOrbitalCatalog.Visible(idx)));
    fprintf('El min  : %.3f deg\n', ...
        min(TOrbitalCatalog.Elevation_deg(idx)));
    fprintf('El max  : %.3f deg\n', ...
        max(TOrbitalCatalog.Elevation_deg(idx)));
end

fprintf('\nTOTAL ORBITAL OBJECTS : %d\n', ...
    height(TOrbitalCatalog));

fprintf('TOTAL VISIBLE         : %d\n', ...
    sum(TOrbitalCatalog.Visible));

fprintf('\nMAT file:\n%s\n',matFile);
fprintf('\nCSV file:\n%s\n',csvFile);

fprintf('\n============================================================\n');
fprintf('ORBITAL CATALOG COMPLETED\n');
fprintf('============================================================\n');

%% ============================================================
% Local functions
% =============================================================

function rECEF_m = ...
    geodeticToECEFWGS84(lat_deg,lon_deg,h_m)

    a = 6378137.0;
    f = 1/298.257223563;
    e2 = f*(2-f);

    lat = deg2rad(lat_deg);
    lon = deg2rad(lon_deg);

    N = a/sqrt(1-e2*sin(lat)^2);

    x = (N+h_m)*cos(lat)*cos(lon);
    y = (N+h_m)*cos(lat)*sin(lon);
    z = (N*(1-e2)+h_m)*sin(lat);

    rECEF_m = [x;y;z];
end

function [az_deg,el_deg,range_m,uLOS] = ...
    lookAnglesECEF( ...
        rTx_m, ...
        rUser_m, ...
        userLat_deg, ...
        userLon_deg)

    rTx_m = rTx_m(:);
    rUser_m = rUser_m(:);

    d = rTx_m-rUser_m;
    range_m = norm(d);

    assert(range_m>0, ...
        'Invalid transmitter/user geometry.');

    uLOS = d/range_m;

    lat = deg2rad(userLat_deg);
    lon = deg2rad(userLon_deg);

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

    az_deg = mod(atan2d(east,north),360);

    el_deg = ...
        atan2d(up,hypot(east,north));
end
