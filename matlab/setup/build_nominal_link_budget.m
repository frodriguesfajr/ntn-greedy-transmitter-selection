clear;
clc;

%% ============================================================
% BUILD NOMINAL LINK BUDGET - PORTABLE P.2145 + P.676-13
%
% Geometry:
%   Corrected Vallado SGP4 / WGS-72 orbital states
%
% Link-budget model:
%   Same RF assumptions used in the paper pipeline.
%
% Atmospheric loss:
%   Portable ITU-R P.2145 + ITU-R P.676-13 implementation.
%   No Satellite Communications Toolbox is required.
%
% Notes:
%   - Annual-median (50th percentile) P.2145 meteorology is used.
%   - Historical MATLAB/P.676-12 comparisons belong only to the
%     validation script diagnose_gas_model.m and are intentionally
%     excluded from this public simulation pipeline.
% =============================================================

thisFile = mfilename('fullpath');

if isempty(thisFile)
    error('This script must be executed from a saved .m file.');
end

setupDir = fileparts(thisFile);
matlabDir = fileparts(setupDir);
repoRoot = fileparts(matlabDir);

assert(isfolder(repoRoot), ...
    'Repository root not found.');

cd(repoRoot);

%% ============================================================
% Repository paths
% =============================================================

addpath(setupDir,'-begin');
setup_paths;

atmosphereDir = fullfile(matlabDir,'atmosphere');

assert(isfolder(atmosphereDir), ...
    'Atmosphere directory not found:\n%s',atmosphereDir);

addpath(atmosphereDir,'-begin');

assert(exist('itu2145_annual50_at_location','file') == 2, ...
    'itu2145_annual50_at_location is not available.');

assert(exist('itu676_annex2_slant_gas_loss','file') == 2, ...
    'itu676_annex2_slant_gas_loss is not available.');

fprintf('\n');
fprintf('============================================================\n');
fprintf('BUILD NOMINAL LINK BUDGET - PORTABLE P.2145 + P.676-13\n');
fprintf('============================================================\n');

%% ============================================================
% Load corrected 26-Tx pool
% =============================================================

poolFile = fullfile( ...
    repoRoot, ...
    'results', ...
    'candidate_pool_26.mat');

assert(isfile(poolFile), ...
    'Corrected candidate pool not found:\n%s',poolFile);

P = load(poolFile);

Tpool = P.Tpool;

analysisTimeUTC = P.analysisTimeUTC;

userLat_deg = P.userLat_deg;
userLon_deg = P.userLon_deg;
userH_m     = P.userH_m;

Npool = height(Tpool);

assert(Npool == 26, ...
    'Expected 26 transmitters.');

Label        = string(Tpool.Label);
Architecture = string(Tpool.Architecture);
NORAD        = string(Tpool.NORAD_CAT_ID);
ObjectName   = string(Tpool.ObjectName);

Elevation_deg = Tpool.Elevation_deg;
SlantRange_m  = 1000*Tpool.SlantRange_km;

%% ============================================================
% Ranging parameters
% =============================================================

c_mps  = 299792458;
beta   = 1.023e6;
Tcoh   = 20e-3;
etaPos = 0.01;

Branging_Hz  = beta;
Branging_MHz = beta/1e6;

%% ============================================================
% Propagation parameters
% =============================================================

GasAnnualExceedance_pct = 50.0;

assert(abs(GasAnnualExceedance_pct-50.0) < 1e-12, ...
    ['The portable meteorological implementation currently uses ' ...
     'the annual 50th-percentile P.2145 data set.']);

PropagationModel = ...
    "free-space + portable ITU-R P.2145/P.676-13 annual-median gaseous attenuation";

%% ============================================================
% RF architecture parameters
% =============================================================

% HAPS
HAPS_Frequency_GHz       = 31.15;
HAPS_EIRP_PSD_dBW_MHz    = -9.1;
HAPS_ReceiverFoM_dBK     = 12.3;

% LEO / Starlink
LEO_Frequency_GHz = 20.0;

LEO_ReceiverFoM_dBK = 10.5;

LEO_ElevationControl_deg = [ ...
    10;
    25;
    90];

LEO_EIRP_Control_dBW_MHz = [ ...
    -2.7;
    17.9;
    11.3];

% MEO / O3b mPOWER
MEO_Frequency_GHz = 20.0;
EIRP_MEO_total_dBW = 49.7;
BW_MEO_MHz = 100;
MEO_ReceiverFoM_dBK = 18.4;

% GEO
GEO_Frequency_GHz = 19.95;
EIRP_GEO_total_dBW = 57.0;
BW_GEO_MHz = 81.0;
G_rx_GEO_dBi = 44.5;
Tsys_GEO_K = 250;

%% ============================================================
% Allocate arrays
% =============================================================

Frequency_GHz = nan(Npool,1);
EIRP_PSD_dBW_MHz = nan(Npool,1);
ReceiverFoM_dBK = nan(Npool,1);

ParameterSource = strings(Npool,1);
ReceiverChain   = strings(Npool,1);

%% ============================================================
% HAPS
% =============================================================

idx = Architecture == "HAPS";

Frequency_GHz(idx) = HAPS_Frequency_GHz;

EIRP_PSD_dBW_MHz(idx) = ...
    HAPS_EIRP_PSD_dBW_MHz;

ReceiverFoM_dBK(idx) = ...
    HAPS_ReceiverFoM_dBK;

ParameterSource(idx) = ...
    "ITU-R F.2439 System 6: HAPS->CPE nominal parameters";

ReceiverChain(idx) = ...
    "HAPS reference CPE";

%% ============================================================
% LEO
% =============================================================

idx = Architecture == "LEO";

Frequency_GHz(idx) = LEO_Frequency_GHz;

ReceiverFoM_dBK(idx) = ...
    LEO_ReceiverFoM_dBK;

EIRP_PSD_dBW_MHz(idx) = ...
    starlinkEIRPFromElevation( ...
        Elevation_deg(idx), ...
        LEO_ElevationControl_deg, ...
        LEO_EIRP_Control_dBW_MHz);

ParameterSource(idx) = ...
    "SpaceX Gen2 Technical Attachment Table A.4-5 + interpolation";

ReceiverChain(idx) = ...
    "Starlink Gen2 user terminal";

%% ============================================================
% MEO
% =============================================================

idx = Architecture == "MEO";

Frequency_GHz(idx) = ...
    MEO_Frequency_GHz;

EIRP_PSD_MEO_dBW_MHz = ...
    EIRP_MEO_total_dBW ...
    - 10*log10(BW_MEO_MHz);

EIRP_PSD_dBW_MHz(idx) = ...
    EIRP_PSD_MEO_dBW_MHz;

ReceiverFoM_dBK(idx) = ...
    MEO_ReceiverFoM_dBK;

ParameterSource(idx) = ...
    "O3b mPOWER + mP85L reference";

ReceiverChain(idx) = ...
    "O3b mPOWER mP85L terminal";

%% ============================================================
% GEO
% =============================================================

idx = Architecture == "GEO";

Frequency_GHz(idx) = ...
    GEO_Frequency_GHz;

EIRP_PSD_GEO_dBW_MHz = ...
    EIRP_GEO_total_dBW ...
    - 10*log10(BW_GEO_MHz);

EIRP_PSD_dBW_MHz(idx) = ...
    EIRP_PSD_GEO_dBW_MHz;

ReceiverFoM_GEO_dBK = ...
    G_rx_GEO_dBi ...
    - 10*log10(Tsys_GEO_K);

ReceiverFoM_dBK(idx) = ...
    ReceiverFoM_GEO_dBK;

ParameterSource(idx) = ...
    "ITU-R S.1328-5 representative GSO Ka system";

ReceiverChain(idx) = ...
    "Representative GSO Ka earth station";

%% ============================================================
% Check architecture parameters
% =============================================================

assert(all(isfinite(Frequency_GHz)));
assert(all(isfinite(EIRP_PSD_dBW_MHz)));
assert(all(isfinite(ReceiverFoM_dBK)));

%% ============================================================
% Portable P.2145 annual-median meteorology
% =============================================================

p2145Dir = fullfile( ...
    repoRoot, ...
    'data', ...
    'itu', ...
    'p2145', ...
    'annual');

assert(isfolder(p2145Dir), ...
    'P.2145 annual-data directory not found:\n%s',p2145Dir);

met = itu2145_annual50_at_location( ...
    p2145Dir, ...
    userLat_deg, ...
    userLon_deg, ...
    userH_m);

P2145_meteorology = met;

fprintf('\nP.2145 annual-median meteorology\n');
fprintf('P   = %.9f hPa\n',met.P_hPa);
fprintf('T   = %.9f K\n',met.T_K);
fprintf('rho = %.9f g/m^3\n',met.rho_gm3);
fprintf('V   = %.9f kg/m^2\n',met.V_kgm2);

%% ============================================================
% Atmospheric gaseous attenuation
% =============================================================

fprintf('\nComputing portable ITU-R P.676-13 gaseous attenuation...\n');

AtmosphericGasLoss_dB = computePortableGasLoss( ...
    Frequency_GHz, ...
    Elevation_deg, ...
    userH_m, ...
    met);

%% ============================================================
% Free-space path loss
% =============================================================

SlantRange_km = SlantRange_m/1e3;

FSPL_dB = ...
    92.45 ...
    + 20*log10(SlantRange_km) ...
    + 20*log10(Frequency_GHz);

%% ============================================================
% Positioning signal EIRP
% =============================================================

EIRP_positioning_dBW = ...
    EIRP_PSD_dBW_MHz ...
    + 10*log10(Branging_MHz) ...
    + 10*log10(etaPos);

%% ============================================================
% C/N0
% =============================================================

CN0_dBHz = ...
    EIRP_positioning_dBW ...
    - FSPL_dB ...
    - AtmosphericGasLoss_dB ...
    + ReceiverFoM_dBK ...
    + 228.6;

CN0_linear = ...
    10.^(CN0_dBHz/10);

%% ============================================================
% Ranging standard deviation
% =============================================================

SigmaRho_m = ...
    c_mps ./ ...
    (2*pi*beta .* ...
    sqrt(CN0_linear*Tcoh));

VarianceRho_m2 = ...
    SigmaRho_m.^2;

%% ============================================================
% Compatibility variables
% =============================================================

% Gaseous attenuation is already represented explicitly by
% AtmosphericGasLoss_dB. No additional loss is imposed here.
AdditionalLoss_dB = zeros(Npool,1);

GT_dBK = ReceiverFoM_dBK;

%% ============================================================
% Build output table
% =============================================================

PoolIndex = Tpool.PoolIndex;

Tlink = table( ...
    PoolIndex, ...
    Label, ...
    Architecture, ...
    NORAD, ...
    ObjectName, ...
    Elevation_deg, ...
    SlantRange_km, ...
    Frequency_GHz, ...
    EIRP_PSD_dBW_MHz, ...
    ReceiverFoM_dBK, ...
    AtmosphericGasLoss_dB, ...
    AdditionalLoss_dB, ...
    FSPL_dB, ...
    EIRP_positioning_dBW, ...
    CN0_dBHz, ...
    SigmaRho_m, ...
    VarianceRho_m2, ...
    ReceiverChain, ...
    ParameterSource, ...
    'VariableNames',{ ...
    'PoolIndex', ...
    'Label', ...
    'Architecture', ...
    'NORAD_CAT_ID', ...
    'ObjectName', ...
    'Elevation_deg', ...
    'SlantRange_km', ...
    'Frequency_GHz', ...
    'EIRP_PSD_dBW_MHz', ...
    'ReceiverFoM_dBK', ...
    'AtmosphericGasLoss_dB', ...
    'AdditionalLoss_dB', ...
    'FSPL_dB', ...
    'EIRP_positioning_dBW', ...
    'CN0_dBHz', ...
    'SigmaRho_m', ...
    'VarianceRho_m2', ...
    'ReceiverChain', ...
    'ParameterSource'});

%% ============================================================
% Summary
% =============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('PORTABLE LINK-BUDGET SUMMARY\n');
fprintf('============================================================\n');

architectures = ["HAPS","LEO","MEO","GEO"];

for a = 1:numel(architectures)

    arch = architectures(a);

    idx = Architecture == arch;

    fprintf('\n%s (%d)\n',arch,sum(idx));

    fprintf('Gas loss min / max   = %.6f / %.6f dB\n', ...
        min(AtmosphericGasLoss_dB(idx)), ...
        max(AtmosphericGasLoss_dB(idx)));

    fprintf('C/N0 min / max       = %.3f / %.3f dB-Hz\n', ...
        min(CN0_dBHz(idx)), ...
        max(CN0_dBHz(idx)));

    fprintf('SigmaRho min / max   = %.6f / %.6f m\n', ...
        min(SigmaRho_m(idx)), ...
        max(SigmaRho_m(idx)));

end

%% ============================================================
% Save
% =============================================================

resultsDir = fullfile(repoRoot,'results');

if ~exist(resultsDir,'dir')
    mkdir(resultsDir);
end

matFile = fullfile( ...
    resultsDir, ...
    'link_budget.mat');

csvFile = fullfile( ...
    resultsDir, ...
    'link_budget.csv');

save(matFile, ...
    'Tlink', ...
    'Label', ...
    'Architecture', ...
    'NORAD', ...
    'ObjectName', ...
    'Elevation_deg', ...
    'SlantRange_m', ...
    'Frequency_GHz', ...
    'EIRP_PSD_dBW_MHz', ...
    'ReceiverFoM_dBK', ...
    'AtmosphericGasLoss_dB', ...
    'AdditionalLoss_dB', ...
    'FSPL_dB', ...
    'EIRP_positioning_dBW', ...
    'CN0_dBHz', ...
    'CN0_linear', ...
    'SigmaRho_m', ...
    'VarianceRho_m2', ...
    'ReceiverChain', ...
    'ParameterSource', ...
    'analysisTimeUTC', ...
    'userLat_deg', ...
    'userLon_deg', ...
    'userH_m', ...
    'beta', ...
    'Tcoh', ...
    'etaPos', ...
    'Branging_Hz', ...
    'Branging_MHz', ...
    'GasAnnualExceedance_pct', ...
    'PropagationModel', ...
    'P2145_meteorology');

writetable(Tlink,csvFile);

fprintf('\nMAT file:\n%s\n',matFile);
fprintf('\nCSV file:\n%s\n',csvFile);

fprintf('\n');
fprintf('============================================================\n');
fprintf('NOMINAL PORTABLE LINK BUDGET COMPLETED\n');
fprintf('Satellite Communications Toolbox: NOT REQUIRED\n');
fprintf('============================================================\n');

%% ============================================================
% Local function - Starlink EIRP interpolation
% =============================================================

function eirp = starlinkEIRPFromElevation( ...
    elevation_deg, ...
    controlElevation_deg, ...
    controlEIRP_dBW_MHz)

    elevation_deg = elevation_deg(:);

    elevClamped = min( ...
        max(elevation_deg,min(controlElevation_deg)), ...
        max(controlElevation_deg));

    eirp = interp1( ...
        controlElevation_deg, ...
        controlEIRP_dBW_MHz, ...
        elevClamped, ...
        'linear');

    eirp = eirp(:);

end

%% ============================================================
% Local function - Portable ITU-R P.676-13 gaseous loss
% =============================================================

function Ag_dB = computePortableGasLoss( ...
    Frequency_GHz, ...
    Elevation_deg, ...
    userH_m, ...
    met)

    Frequency_GHz = Frequency_GHz(:);
    Elevation_deg = Elevation_deg(:);

    if numel(Frequency_GHz) ~= numel(Elevation_deg)
        error('Frequency and elevation arrays must have equal length.');
    end

    if any(Elevation_deg < 5) || any(Elevation_deg > 90)
        error('Portable P.676 model requires elevation between 5 and 90 deg.');
    end

    N = numel(Frequency_GHz);
    Ag_dB = nan(N,1);

    stationHeight_km = userH_m/1e3;

    for ii = 1:N

        Ag_dB(ii) = ...
            itu676_annex2_slant_gas_loss( ...
                Frequency_GHz(ii), ...
                Elevation_deg(ii), ...
                stationHeight_km, ...
                met);

    end

    if any(~isfinite(Ag_dB)) || any(Ag_dB < 0)
        error('Invalid portable ITU-R gaseous attenuation.');
    end

end
