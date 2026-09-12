%% scenario4_task1_build_temporal_geometry.m
% TASK CODE 1 - STEP 1
% Build and validate four temporal geometry sets for the paper.
%
% Paper repository:
%   /MATLAB Drive/ntn-greedy-transmitter-selection
%
% Epochs:
%   t = 0, 30, 60, 90 s
% relative to:
%   01-Aug-2026 12:00:00 UTC
%
% This script DOES NOT run the link budget, FTS, or BTS yet.
% It only:
%   1) loads the nominal 26-candidate pool;
%   2) propagates LEO/MEO/GEO candidates with the archived TLEs and
%      Vallado SGP4;
%   3) keeps the four synthetic HAPS fixed in ECEF;
%   4) computes ECEF position, slant range, azimuth, elevation,
%      visibility, and geometry matrix H at each epoch;
%   5) validates t = 0 against the stored nominal candidate pool;
%   6) saves the four geometry sets into one MAT file.
%
% Main output:
%   results/temporal_geometry_4epochs/scenario4_geometry_4epochs.mat
%
% Additional inspection outputs:
%   results/temporal_geometry_4epochs/geometry_summary.csv
%   results/temporal_geometry_4epochs/geometry_t000s.csv
%   results/temporal_geometry_4epochs/geometry_t030s.csv
%   results/temporal_geometry_4epochs/geometry_t060s.csv
%   results/temporal_geometry_4epochs/geometry_t090s.csv
%
% Required existing project files:
%   setup_paths.m
%   paper_epoch_eop.m
%   propagate_tle_vallado.m
%   results/candidate_pool_26.csv
%   tle/<Architecture>/<NORAD>/<NORAD>.tle
%
% The resulting MAT file is intended to be the common geometry input for
% the next steps:
%   Step 2: temporal link budget -> R(t)
%   Step 3: No Selection / FTS / BTS comparison
%
% -------------------------------------------------------------------------

close all;
clear;
clc;
format long;

%% ========================================================================
% Repository
% =========================================================================

paperRoot = '/MATLAB Drive/ntn-greedy-transmitter-selection';

if ~isfolder(paperRoot)
    error('Paper repository not found: %s', paperRoot);
end

cd(paperRoot);

if ~isfile(fullfile(paperRoot,'setup_paths.m'))
    error('setup_paths.m not found in paper repository.');
end

setup_paths;

candidateCsv = fullfile(paperRoot,'results','candidate_pool_26.csv');

outDir = fullfile( ...
    paperRoot,'results','temporal_geometry_4epochs');

if ~isfile(candidateCsv)
    error('Candidate-pool file not found: %s', candidateCsv);
end

if ~isfolder(outDir)
    mkdir(outDir);
end

if exist('paper_epoch_eop','file') == 0
    error('paper_epoch_eop.m is not available on the MATLAB path.');
end

if exist('propagate_tle_vallado','file') == 0
    error('propagate_tle_vallado.m is not available on the MATLAB path.');
end

%% ========================================================================
% Fixed experiment configuration
% =========================================================================

analysisTimeUTC = datetime( ...
    2026,8,1,12,0,0, ...
    'TimeZone','UTC');

timeOffset_s = [0;30;60;90];
timeUTC = analysisTimeUTC + seconds(timeOffset_s);
timeUTC.Format = 'yyyy-MM-dd HH:mm:ss';

user.lat_deg = -22.8596582;
user.lon_deg = -43.2303236;
user.h_m = 10;

mask.HAPS_deg = 15;
mask.LEO_deg  = 5;
mask.MEO_deg  = 5;
mask.GEO_deg  = 5;

tleRoot = fullfile(paperRoot,'tle');

%% ========================================================================
% Load nominal candidate pool
% =========================================================================

Tpool = readtable(candidateCsv,'TextType','string');

if ~ismember("NORAD_CAT_ID",string(Tpool.Properties.VariableNames))
    if ismember("NORAD",string(Tpool.Properties.VariableNames))
        Tpool.NORAD_CAT_ID = string(Tpool.NORAD);
    else
        error('Candidate pool has no NORAD identifier column.');
    end
end

requiredColumns = [ ...
    "PoolIndex","Label","Architecture","NORAD_CAT_ID", ...
    "X_ECEF_m","Y_ECEF_m","Z_ECEF_m"];

availableColumns = string(Tpool.Properties.VariableNames);
missingColumns = setdiff(requiredColumns,availableColumns);

if ~isempty(missingColumns)
    error('Missing candidate-pool column(s): %s', ...
        strjoin(missingColumns,", "));
end

Tpool = sortrows(Tpool,'PoolIndex');

N = height(Tpool);

if N ~= 26
    error('Expected 26 candidates; found %d.',N);
end

PoolIndex = double(Tpool.PoolIndex(:));
Label = string(Tpool.Label(:));
Architecture = upper(string(Tpool.Architecture(:)));
NORAD_CAT_ID = string(Tpool.NORAD_CAT_ID(:));

nominalPositionECEF_m = [ ...
    double(Tpool.X_ECEF_m), ...
    double(Tpool.Y_ECEF_m), ...
    double(Tpool.Z_ECEF_m)];

validArchitectures = ["HAPS","LEO","MEO","GEO"];

if any(~ismember(Architecture,validArchitectures))
    bad = unique(Architecture(~ismember(Architecture,validArchitectures)));
    error('Unsupported architecture(s): %s',strjoin(bad,", "));
end

%% ========================================================================
% Receiver ECEF
% =========================================================================

userPositionECEF_m = geodeticToECEFWGS84( ...
    user.lat_deg,user.lon_deg,user.h_m);

%% ========================================================================
% Propagate the same 26 candidate identities at four epochs
% =========================================================================

Ntime = numel(timeOffset_s);
PositionECEF_m = nan(N,3,Ntime);

eop = paper_epoch_eop();

fprintf('\n============================================================\n');
fprintf('TASK CODE 1 - STEP 1: TEMPORAL GEOMETRY\n');
fprintf('============================================================\n');
fprintf('Repository     : %s\n',paperRoot);
fprintf('Reference UTC  : %s\n',char(string(analysisTimeUTC)));
fprintf('Epoch offsets  : %s s\n',mat2str(timeOffset_s.'));
fprintf('Candidates     : %d\n',N);
fprintf('============================================================\n\n');

for i = 1:N

    if Architecture(i) == "HAPS"

        % Synthetic HAPS are stationary in this short-term experiment.
        positionHistory_m = repmat( ...
            nominalPositionECEF_m(i,:),Ntime,1);

        fprintf('Fixed       %-7s HAPS\n',Label(i));

    else

        norad = strip(NORAD_CAT_ID(i));

        if ismissing(norad) || strlength(norad) == 0
            error('Missing NORAD ID for %s.',Label(i));
        end

        tleFile = fullfile( ...
            tleRoot,Architecture(i),norad,norad+".tle");

        if ~isfile(tleFile)
            error('Archived TLE missing for %s: %s', ...
                Label(i),tleFile);
        end

        positionHistory_m = nan(Ntime,3);

        fprintf('Propagating %-7s %-3s NORAD %s\n', ...
            Label(i),Architecture(i),norad);

        for k = 1:Ntime

            state = propagate_tle_vallado( ...
                tleFile,timeUTC(k),eop);

            if isfield(state,'sgp4_error') && ...
                    state.sgp4_error ~= 0
                error( ...
                    'SGP4 error %d for %s at %s.', ...
                    state.sgp4_error, ...
                    Label(i), ...
                    char(string(timeUTC(k))));
            end

            if ~isfield(state,'rECEF_m')
                error('propagate_tle_vallado did not return rECEF_m.');
            end

            positionHistory_m(k,:) = state.rECEF_m(:).';
        end
    end

    PositionECEF_m(i,:,:) = permute( ...
        positionHistory_m,[3 2 1]);
end

%% ========================================================================
% Validate propagated t = 0 positions
% =========================================================================

positionDifferenceT0_m = vecnorm( ...
    PositionECEF_m(:,:,1)-nominalPositionECEF_m,2,2);

maxPositionDifferenceT0_m = max(positionDifferenceT0_m);

fprintf('\nMaximum propagated-vs-nominal ECEF difference at t=0: %.6f m\n', ...
    maxPositionDifferenceT0_m);

% Same tolerance already used in the current paper-figure script.
if maxPositionDifferenceT0_m > 0.05
    error([ ...
        'The propagated geometry does not reproduce the nominal candidate ' ...
        'pool at t=0. Maximum difference = %.6f m.'], ...
        maxPositionDifferenceT0_m);
end

%% ========================================================================
% Build four geometry sets
% =========================================================================

epochs = repmat(struct( ...
    'TimeOffset_s',[], ...
    'TimeUTC',[], ...
    'PositionECEF_m',[], ...
    'SlantRange_m',[], ...
    'Azimuth_deg',[], ...
    'Elevation_deg',[], ...
    'ElevationMask_deg',[], ...
    'Visible',[], ...
    'VisibleIndices',[], ...
    'VisibleCount',[], ...
    'H',[], ...
    'VisibleH',[], ...
    'VisibleRank',[], ...
    'GeometryTable',[]),Ntime,1);

VisibleCount = zeros(Ntime,1);
VisibleRank = zeros(Ntime,1);
VisibleLabels = strings(Ntime,1);

for k = 1:Ntime

    positionNow_m = PositionECEF_m(:,:,k);

    [H,slantRange_m,azimuth_deg,elevation_deg] = ...
        geometryFromECEF( ...
            positionNow_m, ...
            userPositionECEF_m, ...
            user.lat_deg, ...
            user.lon_deg);

    elevationMask_deg = nan(N,1);

    elevationMask_deg(Architecture=="HAPS") = mask.HAPS_deg;
    elevationMask_deg(Architecture=="LEO")  = mask.LEO_deg;
    elevationMask_deg(Architecture=="MEO")  = mask.MEO_deg;
    elevationMask_deg(Architecture=="GEO")  = mask.GEO_deg;

    visible = elevation_deg >= elevationMask_deg;
    visibleIdx = find(visible);

    Hvisible = H(visibleIdx,:);
    rankVisible = rank(Hvisible);

    VisibleCount(k) = numel(visibleIdx);
    VisibleRank(k) = rankVisible;
    VisibleLabels(k) = strjoin(Label(visibleIdx),";");

    geometryTable = table( ...
        PoolIndex, ...
        Label, ...
        Architecture, ...
        NORAD_CAT_ID, ...
        positionNow_m(:,1), ...
        positionNow_m(:,2), ...
        positionNow_m(:,3), ...
        slantRange_m, ...
        azimuth_deg, ...
        elevation_deg, ...
        elevationMask_deg, ...
        visible, ...
        H(:,1), ...
        H(:,2), ...
        H(:,3), ...
        H(:,4), ...
        'VariableNames',{ ...
        'PoolIndex', ...
        'Label', ...
        'Architecture', ...
        'NORAD_CAT_ID', ...
        'X_ECEF_m', ...
        'Y_ECEF_m', ...
        'Z_ECEF_m', ...
        'SlantRange_m', ...
        'Azimuth_deg', ...
        'Elevation_deg', ...
        'ElevationMask_deg', ...
        'Visible', ...
        'H1', ...
        'H2', ...
        'H3', ...
        'H4'});

    epochs(k).TimeOffset_s = timeOffset_s(k);
    epochs(k).TimeUTC = timeUTC(k);
    epochs(k).PositionECEF_m = positionNow_m;
    epochs(k).SlantRange_m = slantRange_m;
    epochs(k).Azimuth_deg = azimuth_deg;
    epochs(k).Elevation_deg = elevation_deg;
    epochs(k).ElevationMask_deg = elevationMask_deg;
    epochs(k).Visible = visible;
    epochs(k).VisibleIndices = visibleIdx;
    epochs(k).VisibleCount = numel(visibleIdx);
    epochs(k).H = H;
    epochs(k).VisibleH = Hvisible;
    epochs(k).VisibleRank = rankVisible;
    epochs(k).GeometryTable = geometryTable;

    epochFile = fullfile( ...
        outDir, ...
        sprintf('geometry_t%03ds.csv',timeOffset_s(k)));

    writetable(geometryTable,epochFile);

    fprintf('\n--- Epoch t = %3d s ---\n',timeOffset_s(k));
    fprintf('UTC           : %s\n',char(string(timeUTC(k))));
    fprintf('Visible       : %d / %d\n',numel(visibleIdx),N);
    fprintf('rank(Hvisible): %d\n',rankVisible);

    if rankVisible < 4
        warning('Visible geometry is rank deficient at t=%d s.', ...
            timeOffset_s(k));
    end
end

%% ========================================================================
% Save summary and MAT package
% =========================================================================

Tsummary = table( ...
    timeOffset_s, ...
    timeUTC, ...
    VisibleCount, ...
    VisibleRank, ...
    VisibleLabels, ...
    'VariableNames',{ ...
    'TimeOffset_s', ...
    'TimeUTC', ...
    'VisibleCount', ...
    'VisibleRank', ...
    'VisibleLabels'});

summaryCsv = fullfile(outDir,'geometry_summary.csv');
writetable(Tsummary,summaryCsv);

metadata = struct;
metadata.Repository = paperRoot;
metadata.AnalysisTimeUTC = analysisTimeUTC;
metadata.TimeOffset_s = timeOffset_s;
metadata.UserLatitude_deg = user.lat_deg;
metadata.UserLongitude_deg = user.lon_deg;
metadata.UserHeight_m = user.h_m;
metadata.UserPositionECEF_m = userPositionECEF_m;
metadata.ElevationMasks_deg = mask;
metadata.CandidateCount = N;
metadata.CandidateLabels = Label;
metadata.CandidateArchitectures = Architecture;
metadata.CandidateNORAD = NORAD_CAT_ID;
metadata.MaxPositionDifferenceT0_m = maxPositionDifferenceT0_m;
metadata.SGP4Implementation = 'Vallado implementation included in repository';
metadata.HAPSModel = 'Fixed ECEF synthetic HAPS';
metadata.Purpose = [ ...
    'Temporal geometry input for subsequent link-budget, ', ...
    'FTS, and BTS evaluation.'];

matFile = fullfile( ...
    outDir,'scenario4_geometry_4epochs.mat');

save( ...
    matFile, ...
    'epochs', ...
    'metadata', ...
    'Tsummary', ...
    'PoolIndex', ...
    'Label', ...
    'Architecture', ...
    'NORAD_CAT_ID', ...
    'PositionECEF_m', ...
    'nominalPositionECEF_m', ...
    'userPositionECEF_m', ...
    'timeOffset_s', ...
    'timeUTC', ...
    '-v7.3');

%% ========================================================================
% Final console report
% =========================================================================

fprintf('\n============================================================\n');
fprintf('TEMPORAL GEOMETRY COMPLETED\n');
fprintf('============================================================\n');
disp(Tsummary(:,["TimeOffset_s","TimeUTC","VisibleCount","VisibleRank"]));

fprintf('MAT file:\n  %s\n',matFile);
fprintf('Summary CSV:\n  %s\n',summaryCsv);
fprintf('Maximum t=0 position difference: %.6f m\n', ...
    maxPositionDifferenceT0_m);
fprintf('============================================================\n');

%% ========================================================================
% Local functions
% =========================================================================

function positionECEF_m = ...
    geodeticToECEFWGS84(lat_deg,lon_deg,h_m)

    a_m = 6378137.0;
    flattening = 1/298.257223563;
    eccentricitySquared = flattening*(2-flattening);

    lat = deg2rad(lat_deg);
    lon = deg2rad(lon_deg);

    primeVerticalRadius_m = ...
        a_m/sqrt( ...
        1-eccentricitySquared*sin(lat)^2);

    positionECEF_m = [ ...
        (primeVerticalRadius_m+h_m)*cos(lat)*cos(lon), ...
        (primeVerticalRadius_m+h_m)*cos(lat)*sin(lon), ...
        (primeVerticalRadius_m*(1-eccentricitySquared)+h_m)*sin(lat)];
end

function [H,slantRange_m,azimuth_deg,elevation_deg] = ...
    geometryFromECEF( ...
        transmitterECEF_m, ...
        userECEF_m, ...
        userLat_deg, ...
        userLon_deg)

    userECEF_m = userECEF_m(:).';

    delta_m = transmitterECEF_m-userECEF_m;

    slantRange_m = vecnorm(delta_m,2,2);

    if any(~isfinite(slantRange_m)) || ...
       any(slantRange_m <= eps)
        error('Invalid transmitter-to-user range.');
    end

    % Unit vector from receiver to transmitter.
    unitReceiverToTransmitter = ...
        delta_m./slantRange_m;

    % The pseudorange derivative with respect to receiver position is the
    % opposite direction: (p_u - p_i)/||p_i-p_u||.
    H = [ ...
        -unitReceiverToTransmitter, ...
        ones(size(transmitterECEF_m,1),1)];

    [azimuth_deg,elevation_deg] = ...
        ecefAzEl( ...
            transmitterECEF_m, ...
            userECEF_m, ...
            userLat_deg, ...
            userLon_deg);
end

function [azimuth_deg,elevation_deg] = ...
    ecefAzEl( ...
        transmitterECEF_m, ...
        userECEF_m, ...
        userLat_deg, ...
        userLon_deg)

    lat = deg2rad(userLat_deg);
    lon = deg2rad(userLon_deg);

    rotationECEFToENU = [ ...
        -sin(lon),             cos(lon),            0
        -sin(lat)*cos(lon),   -sin(lat)*sin(lon),   cos(lat)
         cos(lat)*cos(lon),    cos(lat)*sin(lon),   sin(lat)];

    delta_m = ...
        transmitterECEF_m-userECEF_m(:).';

    enu_m = ...
        (rotationECEFToENU*delta_m.').';

    azimuth_deg = ...
        mod(atan2d(enu_m(:,1),enu_m(:,2)),360);

    elevation_deg = ...
        atan2d( ...
            enu_m(:,3), ...
            hypot(enu_m(:,1),enu_m(:,2)));
end
