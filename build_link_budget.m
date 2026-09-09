clear;
clc;

%% ============================================================
% BUILD LINK BUDGET - CORRECTED VALLADO CANDIDATE POOL
%
% Geometry:
%   Corrected Vallado SGP4 / WGS-72 orbital states
%
% Link-budget model:
%   Identical to the original paper implementation
%
% Atmospheric loss:
%   ITU-R P.618/P.676 using Satellite Communications Toolbox
%
% Purpose:
%   Recompute the paper results after correcting orbital
%   propagation, without changing the RF assumptions.
% =============================================================

rootDir = fileparts(mfilename('fullpath'));

if isempty(rootDir)
    rootDir = pwd;
end

cd(rootDir);
setup_paths;

fprintf('\n============================================================\n');
fprintf('BUILD LINK BUDGET - VALLADO CORRECTED POOL\n');
fprintf('============================================================\n');

%% ============================================================
% Check required MATLAB functions
% =============================================================

assert(exist('p618Config','file') ~= 0, ...
    'p618Config is not available.');

assert(exist('p618PropagationLosses','file') ~= 0, ...
    'p618PropagationLosses is not available.');

%% ============================================================
% Load corrected 26-Tx pool
% =============================================================

poolFile = fullfile( ...
    rootDir, ...
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

PropagationModel = ...
    "free-space + ITU-R P.618/P.676 median clear-sky gaseous attenuation";

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

%% Check architecture parameters

assert(all(isfinite(Frequency_GHz)));
assert(all(isfinite(EIRP_PSD_dBW_MHz)));
assert(all(isfinite(ReceiverFoM_dBK)));

%% ============================================================
% Atmospheric gaseous attenuation
% =============================================================

fprintf('\nComputing ITU-R P.618/P.676 gaseous attenuation...\n');

AtmosphericGasLoss_dB = computeP618GasLoss( ...
    Frequency_GHz, ...
    Elevation_deg, ...
    userLat_deg, ...
    userLon_deg, ...
    userH_m, ...
    GasAnnualExceedance_pct);

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

%% Compatibility variable

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
% Load old link budget for controlled comparison
% =============================================================

oldLBFile = fullfile( ...
    'C:\Repository\master-thesis-ntn-positioning', ...
    'matlab_code', ...
    'results_scenario4_link_budget', ...
    'scenario4_link_budget.mat');

if isfile(oldLBFile)

    Old = load(oldLBFile);

    Told = Old.Tlink;

    assert(height(Told) == Npool);

    dElevation_deg = ...
        Tlink.Elevation_deg ...
        - Told.Elevation_deg;

    dRange_m = ...
        1000*(Tlink.SlantRange_km ...
        - Told.SlantRange_km);

    dGas_dB = ...
        Tlink.AtmosphericGasLoss_dB ...
        - Told.AtmosphericGasLoss_dB;

    dFSPL_dB = ...
        Tlink.FSPL_dB ...
        - Told.FSPL_dB;

    dEIRP_PSD_dB = ...
        Tlink.EIRP_PSD_dBW_MHz ...
        - Told.EIRP_PSD_dBW_MHz;

    dCN0_dB = ...
        Tlink.CN0_dBHz ...
        - Told.CN0_dBHz;

    dSigmaRho_m = ...
        Tlink.SigmaRho_m ...
        - Told.SigmaRho_m;

    fprintf('\n============================================================\n');
    fprintf('OLD LINK BUDGET vs CORRECTED VALLADO LINK BUDGET\n');
    fprintf('============================================================\n');

    architectures = ["HAPS","LEO","MEO","GEO"];

    for a = 1:numel(architectures)

        arch = architectures(a);

        idx = Architecture == arch;

        fprintf('\n%s (%d)\n', ...
            arch,sum(idx));

        fprintf('max |dElevation| = %.6f deg\n', ...
            max(abs(dElevation_deg(idx))));

        fprintf('max |dRange|     = %.3f m\n', ...
            max(abs(dRange_m(idx))));

        fprintf('max |dGas|       = %.6f dB\n', ...
            max(abs(dGas_dB(idx))));

        fprintf('max |dFSPL|      = %.6f dB\n', ...
            max(abs(dFSPL_dB(idx))));

        fprintf('max |dEIRP PSD|  = %.6f dB\n', ...
            max(abs(dEIRP_PSD_dB(idx))));

        fprintf('max |dC/N0|      = %.6f dB-Hz\n', ...
            max(abs(dCN0_dB(idx))));

        fprintf('max |dSigmaRho|  = %.9f m\n', ...
            max(abs(dSigmaRho_m(idx))));

    end

else

    warning('Old link budget file was not found.');

end

%% ============================================================
% Summary corrected link budget
% =============================================================

fprintf('\n============================================================\n');
fprintf('CORRECTED LINK-BUDGET SUMMARY\n');
fprintf('============================================================\n');

architectures = ["HAPS","LEO","MEO","GEO"];

for a = 1:numel(architectures)

    arch = architectures(a);

    idx = Architecture == arch;

    fprintf('\n%s\n',arch);

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

resultsDir = fullfile(rootDir,'results');

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
    'PropagationModel');

writetable(Tlink,csvFile);

fprintf('\nMAT file:\n%s\n',matFile);

fprintf('\nCSV file:\n%s\n',csvFile);

fprintf('\n============================================================\n');
fprintf('CORRECTED LINK BUDGET COMPLETED\n');
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
% Local function - ITU-R P.618/P.676 gaseous loss
% =============================================================

function Ag_dB = computeP618GasLoss( ...
    Frequency_GHz, ...
    Elevation_deg, ...
    lat_deg, ...
    lon_deg, ...
    userH_m, ...
    GasAnnualExceedance_pct)

    Frequency_GHz = Frequency_GHz(:);
    Elevation_deg = Elevation_deg(:);

    if any(Elevation_deg < 5)

        error([ ...
            'The P.618 model used in this script ' ...
            'requires elevation >= 5 deg.']);

    end

    N = numel(Frequency_GHz);

    Ag_dB = nan(N,1);

    cfgP = p618Config;

    cfgP.Latitude = lat_deg;
    cfgP.Longitude = lon_deg;

    cfgP.GasAnnualExceedance = ...
        GasAnnualExceedance_pct;

    for ii = 1:N

        cfgP.Frequency = ...
            Frequency_GHz(ii)*1e9;

        cfgP.ElevationAngle = ...
            Elevation_deg(ii);

        pl = p618PropagationLosses( ...
            cfgP, ...
            'StationHeight', ...
            userH_m/1e3);

        Ag_dB(ii) = pl.Ag;

    end

    if any(~isfinite(Ag_dB)) || ...
       any(Ag_dB < 0)

        error('Invalid ITU-R gaseous attenuation.');

    end

end