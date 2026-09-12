%% scenario4_task1_build_temporal_link_budget.m
% TASK CODE 1 - STEP 2
% Build four link-budget / measurement-quality sets for the paper.
%
% Repository:
%   /MATLAB Drive/ntn-greedy-transmitter-selection
%
% Input geometry:
%   results/temporal_geometry_4epochs/scenario4_geometry_4epochs.mat
%
% Nominal reference:
%   results/link_budget.csv
%
% Epochs:
%   t = 0, 30, 60, 90 s
%
% For each epoch, this script recomputes for all 26 candidates:
%   - carrier frequency
%   - architecture-dependent EIRP spectral density
%   - receiver G/T
%   - ITU-R P.618/P.676 gaseous attenuation
%   - free-space path loss
%   - band-integrated EIRP
%   - link C/N0 before positioning-power allocation, q_i [dB-Hz]
%   - positioning C/N0 after eta_Pos = 0.01 [dB-Hz]
%   - gamma_i [Hz]
%   - sigma_rho_i [m]
%   - variance sigma_rho_i^2 [m^2]
%   - covariance matrices R_C and R_visible
%
% IMPORTANT:
%   This step does NOT run FTS or BTS.
%
% The RF assumptions intentionally match the current paper implementation
% used by build_link_budget.m:
%   beta = W = 1.023 MHz
%   Tcoh = 20 ms
%   etaPos = 0.01
%   gaseous attenuation = ITU-R P.618/P.676, 50% annual exceedance
%
% Outputs:
%   results/temporal_link_budget_4epochs/
%       scenario4_link_budget_4epochs.mat
%       link_budget_t000s.csv
%       link_budget_t030s.csv
%       link_budget_t060s.csv
%       link_budget_t090s.csv
%       link_budget_epoch_summary.csv
%       link_budget_architecture_summary.csv
%
% -------------------------------------------------------------------------

close all;
clear;
clc;
format long;

%% ========================================================================
% Repository and files
% =========================================================================

paperRoot = '/MATLAB Drive/ntn-greedy-transmitter-selection';

if ~isfolder(paperRoot)
    error('Paper repository not found: %s',paperRoot);
end

cd(paperRoot);
setup_paths;

geometryMat = fullfile( ...
    paperRoot, ...
    'results', ...
    'temporal_geometry_4epochs', ...
    'scenario4_geometry_4epochs.mat');

nominalLinkCsv = fullfile( ...
    paperRoot, ...
    'results', ...
    'link_budget.csv');

outDir = fullfile( ...
    paperRoot, ...
    'results', ...
    'temporal_link_budget_4epochs');

if ~isfile(geometryMat)
    error(['Temporal geometry MAT file not found:\n%s\n' ...
           'Run scenario4_task1_build_temporal_geometry.m first.'], ...
           geometryMat);
end

if ~isfile(nominalLinkCsv)
    error('Nominal link-budget CSV not found: %s',nominalLinkCsv);
end

if ~isfolder(outDir)
    mkdir(outDir);
end

%% ========================================================================
% Required ITU-R functions
% =========================================================================

assert(exist('p618Config','file') ~= 0, ...
    'p618Config is not available.');

assert(exist('p618PropagationLosses','file') ~= 0, ...
    'p618PropagationLosses is not available.');

%% ========================================================================
% Load validated geometry and nominal link-budget reference
% =========================================================================

G = load(geometryMat);
Tref = readtable(nominalLinkCsv,'TextType','string');

requiredGeometryFields = { ...
    'epochs','metadata','PoolIndex','Label','Architecture', ...
    'NORAD_CAT_ID','timeOffset_s','timeUTC'};

for ii = 1:numel(requiredGeometryFields)
    if ~isfield(G,requiredGeometryFields{ii})
        error('Geometry MAT is missing field: %s',requiredGeometryFields{ii});
    end
end

epochs = G.epochs;
metadata = G.metadata;

PoolIndex = double(G.PoolIndex(:));
Label = string(G.Label(:));
Architecture = upper(string(G.Architecture(:)));
NORAD_CAT_ID = string(G.NORAD_CAT_ID(:));

timeOffset_s = double(G.timeOffset_s(:));
timeUTC = G.timeUTC(:);

N = numel(Label);
Ntime = numel(timeOffset_s);

if N ~= 26
    error('Expected 26 candidates; found %d.',N);
end

if Ntime ~= 4 || ~isequal(timeOffset_s,[0;30;60;90])
    error('Expected epochs t = [0 30 60 90] s.');
end

if ismember("PoolIndex",string(Tref.Properties.VariableNames))
    Tref = sortrows(Tref,'PoolIndex');
end

requiredRefColumns = [ ...
    "Label","Architecture","CN0_dBHz", ...
    "AtmosphericGasLoss_dB","SigmaRho_m"];

availableRefColumns = string(Tref.Properties.VariableNames);
missingRefColumns = setdiff(requiredRefColumns,availableRefColumns);

if ~isempty(missingRefColumns)
    error('Nominal link-budget CSV is missing: %s', ...
        strjoin(missingRefColumns,", "));
end

if height(Tref) ~= N || any(string(Tref.Label(:)) ~= Label)
    error('Nominal link-budget rows do not match geometry ordering.');
end

%% ========================================================================
% Fixed ranging / propagation assumptions
% =========================================================================

c_mps = 299792458;

beta_Hz = 1.023e6;
Tcoh_s = 20e-3;
etaPos = 0.01;

Branging_Hz = beta_Hz;
Branging_MHz = Branging_Hz/1e6;

GasAnnualExceedance_pct = 50.0;

PropagationModel = ...
    "free-space + ITU-R P.618/P.676 median clear-sky gaseous attenuation";

userLat_deg = metadata.UserLatitude_deg;
userLon_deg = metadata.UserLongitude_deg;
userH_m = metadata.UserHeight_m;

%% ========================================================================
% RF architecture parameters
% =========================================================================

% HAPS
HAPS_Frequency_GHz = 31.15;
HAPS_EIRP_PSD_dBW_MHz = -9.1;
HAPS_ReceiverFoM_dBK = 12.3;

% LEO / Starlink Gen2
LEO_Frequency_GHz = 20.0;
LEO_ReceiverFoM_dBK = 10.5;

LEO_ElevationControl_deg = [10;25;90];
LEO_EIRP_Control_dBW_MHz = [-2.7;17.9;11.3];

% MEO / O3b mPOWER
MEO_Frequency_GHz = 20.0;
MEO_EIRP_total_dBW = 49.7;
MEO_BW_MHz = 100;
MEO_ReceiverFoM_dBK = 18.4;

% GEO
GEO_Frequency_GHz = 19.95;
GEO_EIRP_total_dBW = 57.0;
GEO_BW_MHz = 81.0;
GEO_G_rx_dBi = 44.5;
GEO_Tsys_K = 250;

%% ========================================================================
% Static architecture-dependent vectors
% =========================================================================

Frequency_GHz = nan(N,1);
ReceiverFoM_dBK = nan(N,1);
ReceiverChain = strings(N,1);
ParameterSource = strings(N,1);

idx = Architecture=="HAPS";
Frequency_GHz(idx) = HAPS_Frequency_GHz;
ReceiverFoM_dBK(idx) = HAPS_ReceiverFoM_dBK;
ReceiverChain(idx) = "HAPS reference CPE";
ParameterSource(idx) = ...
    "ITU-R F.2439 System 6: HAPS->CPE nominal parameters";

idx = Architecture=="LEO";
Frequency_GHz(idx) = LEO_Frequency_GHz;
ReceiverFoM_dBK(idx) = LEO_ReceiverFoM_dBK;
ReceiverChain(idx) = "Starlink Gen2 user terminal";
ParameterSource(idx) = ...
    "SpaceX Gen2 Technical Attachment Table A.4-5 + interpolation";

idx = Architecture=="MEO";
Frequency_GHz(idx) = MEO_Frequency_GHz;
ReceiverFoM_dBK(idx) = MEO_ReceiverFoM_dBK;
ReceiverChain(idx) = "O3b mPOWER mP85L terminal";
ParameterSource(idx) = "O3b mPOWER + mP85L reference";

idx = Architecture=="GEO";
Frequency_GHz(idx) = GEO_Frequency_GHz;
ReceiverFoM_dBK(idx) = ...
    GEO_G_rx_dBi-10*log10(GEO_Tsys_K);
ReceiverChain(idx) = "Representative GSO Ka earth station";
ParameterSource(idx) = ...
    "ITU-R S.1328-5 representative GSO Ka system";

if any(~isfinite(Frequency_GHz)) || any(~isfinite(ReceiverFoM_dBK))
    error('Invalid architecture-dependent RF parameters.');
end

%% ========================================================================
% Allocate output structure
% =========================================================================

linkEpochs = repmat(struct( ...
    'TimeOffset_s',[], ...
    'TimeUTC',[], ...
    'Visible',[], ...
    'VisibleIndices',[], ...
    'VisibleCount',[], ...
    'Frequency_GHz',[], ...
    'EIRP_PSD_dBW_MHz',[], ...
    'ReceiverFoM_dBK',[], ...
    'AtmosphericGasLoss_dB',[], ...
    'FSPL_dB',[], ...
    'EIRP_band_dBW',[], ...
    'EIRP_positioning_dBW',[], ...
    'q_link_dBHz',[], ...
    'CN0_positioning_dBHz',[], ...
    'Gamma_Hz',[], ...
    'SigmaRho_m',[], ...
    'VarianceRho_m2',[], ...
    'R_C',[], ...
    'R_Visible',[], ...
    'LinkTable',[]),Ntime,1);

epochSummaryRows = [];
archSummaryRows = [];

fprintf('\n============================================================\n');
fprintf('TASK CODE 1 - STEP 2: TEMPORAL LINK BUDGET\n');
fprintf('============================================================\n');
fprintf('Epochs       : %s s\n',mat2str(timeOffset_s.'));
fprintf('Candidates   : %d\n',N);
fprintf('beta = W     : %.3f MHz\n',beta_Hz/1e6);
fprintf('Tcoh         : %.3f ms\n',1e3*Tcoh_s);
fprintf('eta_Pos      : %.4f (%.2f%%)\n',etaPos,100*etaPos);
fprintf('Gas model    : ITU-R P.618/P.676, %.1f%% annual exceedance\n', ...
    GasAnnualExceedance_pct);
fprintf('============================================================\n');

%% ========================================================================
% Recompute link budget at each epoch
% =========================================================================

for k = 1:Ntime

    Elevation_deg = double(epochs(k).Elevation_deg(:));
    SlantRange_m = double(epochs(k).SlantRange_m(:));
    Visible = logical(epochs(k).Visible(:));

    if numel(Elevation_deg) ~= N || numel(SlantRange_m) ~= N
        error('Geometry dimension mismatch at t=%g s.',timeOffset_s(k));
    end

    if ~all(Visible)
        warning('%d of %d candidates are visible at t=%g s.', ...
            sum(Visible),N,timeOffset_s(k));
    end

    % --------------------------------------------------------------------
    % EIRP PSD
    % --------------------------------------------------------------------

    EIRP_PSD_dBW_MHz = nan(N,1);

    idx = Architecture=="HAPS";
    EIRP_PSD_dBW_MHz(idx) = HAPS_EIRP_PSD_dBW_MHz;

    idx = Architecture=="LEO";
    EIRP_PSD_dBW_MHz(idx) = starlinkEIRPFromElevation( ...
        Elevation_deg(idx), ...
        LEO_ElevationControl_deg, ...
        LEO_EIRP_Control_dBW_MHz);

    idx = Architecture=="MEO";
    EIRP_PSD_dBW_MHz(idx) = ...
        MEO_EIRP_total_dBW-10*log10(MEO_BW_MHz);

    idx = Architecture=="GEO";
    EIRP_PSD_dBW_MHz(idx) = ...
        GEO_EIRP_total_dBW-10*log10(GEO_BW_MHz);

    if any(~isfinite(EIRP_PSD_dBW_MHz))
        error('Invalid EIRP PSD at t=%g s.',timeOffset_s(k));
    end

    % --------------------------------------------------------------------
    % ITU-R P.618/P.676 gaseous attenuation
    % --------------------------------------------------------------------

    fprintf('\nComputing epoch t = %3d s ...\n',timeOffset_s(k));

    AtmosphericGasLoss_dB = computeP618GasLoss( ...
        Frequency_GHz, ...
        Elevation_deg, ...
        userLat_deg, ...
        userLon_deg, ...
        userH_m, ...
        GasAnnualExceedance_pct);

    % --------------------------------------------------------------------
    % Free-space path loss
    % --------------------------------------------------------------------

    SlantRange_km = SlantRange_m/1e3;

    FSPL_dB = ...
        92.45 ...
        +20*log10(SlantRange_km) ...
        +20*log10(Frequency_GHz);

    % --------------------------------------------------------------------
    % Band-integrated EIRP before / after positioning allocation
    % --------------------------------------------------------------------

    EIRP_band_dBW = ...
        EIRP_PSD_dBW_MHz ...
        +10*log10(Branging_MHz);

    EIRP_positioning_dBW = ...
        EIRP_band_dBW ...
        +10*log10(etaPos);

    % --------------------------------------------------------------------
    % q_i = link C/N0 before positioning-power allocation [dB-Hz]
    %
    % gamma_i = etaPos * 10^(q_i/10) [Hz]
    % --------------------------------------------------------------------

    q_link_dBHz = ...
        EIRP_band_dBW ...
        -FSPL_dB ...
        -AtmosphericGasLoss_dB ...
        +ReceiverFoM_dBK ...
        +228.6;

    CN0_positioning_dBHz = ...
        q_link_dBHz ...
        +10*log10(etaPos);

    Gamma_Hz = ...
        etaPos.*10.^(q_link_dBHz/10);

    % Numerical identity check for the two equivalent formulations.
    gammaFromCN0_Hz = 10.^(CN0_positioning_dBHz/10);

    if max(abs(Gamma_Hz-gammaFromCN0_Hz)./Gamma_Hz) > 1e-12
        error('gamma_i formulation mismatch at t=%g s.',timeOffset_s(k));
    end

    % --------------------------------------------------------------------
    % Pseudorange uncertainty and covariance
    % --------------------------------------------------------------------

    SigmaRho_m = ...
        c_mps ./ ...
        (2*pi*beta_Hz.*sqrt(Gamma_Hz*Tcoh_s));

    VarianceRho_m2 = SigmaRho_m.^2;

    if any(~isfinite(SigmaRho_m)) || any(SigmaRho_m <= 0)
        error('Invalid sigma_rho at t=%g s.',timeOffset_s(k));
    end

    R_C = diag(VarianceRho_m2);

    visibleIdx = find(Visible);
    R_Visible = diag(VarianceRho_m2(visibleIdx));

    % --------------------------------------------------------------------
    % Per-epoch output table
    % --------------------------------------------------------------------

    LinkTable = table( ...
        PoolIndex, ...
        Label, ...
        Architecture, ...
        NORAD_CAT_ID, ...
        Visible, ...
        Elevation_deg, ...
        SlantRange_km, ...
        Frequency_GHz, ...
        EIRP_PSD_dBW_MHz, ...
        ReceiverFoM_dBK, ...
        AtmosphericGasLoss_dB, ...
        FSPL_dB, ...
        EIRP_band_dBW, ...
        EIRP_positioning_dBW, ...
        q_link_dBHz, ...
        CN0_positioning_dBHz, ...
        Gamma_Hz, ...
        SigmaRho_m, ...
        VarianceRho_m2, ...
        ReceiverChain, ...
        ParameterSource, ...
        'VariableNames',{ ...
        'PoolIndex', ...
        'Label', ...
        'Architecture', ...
        'NORAD_CAT_ID', ...
        'Visible', ...
        'Elevation_deg', ...
        'SlantRange_km', ...
        'Frequency_GHz', ...
        'EIRP_PSD_dBW_MHz', ...
        'ReceiverFoM_dBK', ...
        'AtmosphericGasLoss_dB', ...
        'FSPL_dB', ...
        'EIRP_band_dBW', ...
        'EIRP_positioning_dBW', ...
        'q_link_dBHz', ...
        'CN0_positioning_dBHz', ...
        'Gamma_Hz', ...
        'SigmaRho_m', ...
        'VarianceRho_m2', ...
        'ReceiverChain', ...
        'ParameterSource'});

    linkEpochs(k).TimeOffset_s = timeOffset_s(k);
    linkEpochs(k).TimeUTC = timeUTC(k);
    linkEpochs(k).Visible = Visible;
    linkEpochs(k).VisibleIndices = visibleIdx;
    linkEpochs(k).VisibleCount = numel(visibleIdx);
    linkEpochs(k).Frequency_GHz = Frequency_GHz;
    linkEpochs(k).EIRP_PSD_dBW_MHz = EIRP_PSD_dBW_MHz;
    linkEpochs(k).ReceiverFoM_dBK = ReceiverFoM_dBK;
    linkEpochs(k).AtmosphericGasLoss_dB = AtmosphericGasLoss_dB;
    linkEpochs(k).FSPL_dB = FSPL_dB;
    linkEpochs(k).EIRP_band_dBW = EIRP_band_dBW;
    linkEpochs(k).EIRP_positioning_dBW = EIRP_positioning_dBW;
    linkEpochs(k).q_link_dBHz = q_link_dBHz;
    linkEpochs(k).CN0_positioning_dBHz = CN0_positioning_dBHz;
    linkEpochs(k).Gamma_Hz = Gamma_Hz;
    linkEpochs(k).SigmaRho_m = SigmaRho_m;
    linkEpochs(k).VarianceRho_m2 = VarianceRho_m2;
    linkEpochs(k).R_C = R_C;
    linkEpochs(k).R_Visible = R_Visible;
    linkEpochs(k).LinkTable = LinkTable;

    epochCsv = fullfile( ...
        outDir, ...
        sprintf('link_budget_t%03ds.csv',timeOffset_s(k)));

    writetable(LinkTable,epochCsv);

    % Overall epoch summary.
    epochRow = table( ...
        timeOffset_s(k), ...
        timeUTC(k), ...
        numel(visibleIdx), ...
        min(CN0_positioning_dBHz(Visible)), ...
        median(CN0_positioning_dBHz(Visible)), ...
        max(CN0_positioning_dBHz(Visible)), ...
        min(SigmaRho_m(Visible)), ...
        median(SigmaRho_m(Visible)), ...
        max(SigmaRho_m(Visible)), ...
        'VariableNames',{ ...
        'TimeOffset_s', ...
        'TimeUTC', ...
        'VisibleCount', ...
        'CN0_Min_dBHz', ...
        'CN0_Median_dBHz', ...
        'CN0_Max_dBHz', ...
        'SigmaRho_Min_m', ...
        'SigmaRho_Median_m', ...
        'SigmaRho_Max_m'});

    epochSummaryRows = [epochSummaryRows; epochRow]; %#ok<AGROW>

    % Architecture summary.
    architectures = ["HAPS","LEO","MEO","GEO"];

    for aa = 1:numel(architectures)

        arch = architectures(aa);
        idxArch = Architecture==arch & Visible;

        archRow = table( ...
            timeOffset_s(k), ...
            timeUTC(k), ...
            arch, ...
            sum(idxArch), ...
            min(CN0_positioning_dBHz(idxArch)), ...
            median(CN0_positioning_dBHz(idxArch)), ...
            max(CN0_positioning_dBHz(idxArch)), ...
            min(SigmaRho_m(idxArch)), ...
            median(SigmaRho_m(idxArch)), ...
            max(SigmaRho_m(idxArch)), ...
            min(AtmosphericGasLoss_dB(idxArch)), ...
            median(AtmosphericGasLoss_dB(idxArch)), ...
            max(AtmosphericGasLoss_dB(idxArch)), ...
            'VariableNames',{ ...
            'TimeOffset_s', ...
            'TimeUTC', ...
            'Architecture', ...
            'N', ...
            'CN0_Min_dBHz', ...
            'CN0_Median_dBHz', ...
            'CN0_Max_dBHz', ...
            'SigmaRho_Min_m', ...
            'SigmaRho_Median_m', ...
            'SigmaRho_Max_m', ...
            'GasLoss_Min_dB', ...
            'GasLoss_Median_dB', ...
            'GasLoss_Max_dB'});

        archSummaryRows = [archSummaryRows; archRow]; %#ok<AGROW>
    end

    fprintf('Visible            : %d / %d\n',numel(visibleIdx),N);
    fprintf('C/N0 min/median/max: %.3f / %.3f / %.3f dB-Hz\n', ...
        min(CN0_positioning_dBHz(Visible)), ...
        median(CN0_positioning_dBHz(Visible)), ...
        max(CN0_positioning_dBHz(Visible)));
    fprintf('sigma min/med/max  : %.6f / %.6f / %.6f m\n', ...
        min(SigmaRho_m(Visible)), ...
        median(SigmaRho_m(Visible)), ...
        max(SigmaRho_m(Visible)));
end

%% ========================================================================
% Validate t = 0 against the nominal link-budget CSV
% =========================================================================

refCN0_dBHz = double(Tref.CN0_dBHz(:));
refGas_dB = double(Tref.AtmosphericGasLoss_dB(:));
refSigmaRho_m = double(Tref.SigmaRho_m(:));

t0CN0_dBHz = linkEpochs(1).CN0_positioning_dBHz;
t0Gas_dB = linkEpochs(1).AtmosphericGasLoss_dB;
t0SigmaRho_m = linkEpochs(1).SigmaRho_m;

maxCN0Difference_dB = max(abs(t0CN0_dBHz-refCN0_dBHz));
maxGasDifference_dB = max(abs(t0Gas_dB-refGas_dB));
maxSigmaDifference_m = max(abs(t0SigmaRho_m-refSigmaRho_m));

fprintf('\n============================================================\n');
fprintf('t = 0 LINK-BUDGET VALIDATION\n');
fprintf('============================================================\n');
fprintf('max |Delta C/N0|     = %.6e dB\n',maxCN0Difference_dB);
fprintf('max |Delta gas loss| = %.6e dB\n',maxGasDifference_dB);
fprintf('max |Delta sigma|    = %.6e m\n',maxSigmaDifference_m);

validationPassed = ...
    maxCN0Difference_dB <= 1e-4 && ...
    maxGasDifference_dB <= 1e-4 && ...
    maxSigmaDifference_m <= 1e-6;

if validationPassed
    fprintf('Status               = VALIDATED\n');
else
    warning([ ...
        't=0 temporal link budget differs from the nominal reference ' ...
        'above the validation tolerance. Inspect before running FTS/BTS.']);
end

%% ========================================================================
% Save outputs
% =========================================================================

TepochSummary = epochSummaryRows;
TarchitectureSummary = archSummaryRows;

writetable( ...
    TepochSummary, ...
    fullfile(outDir,'link_budget_epoch_summary.csv'));

writetable( ...
    TarchitectureSummary, ...
    fullfile(outDir,'link_budget_architecture_summary.csv'));

linkMetadata = struct;
linkMetadata.Repository = paperRoot;
linkMetadata.GeometryInput = geometryMat;
linkMetadata.NominalReference = nominalLinkCsv;
linkMetadata.TimeOffset_s = timeOffset_s;
linkMetadata.beta_Hz = beta_Hz;
linkMetadata.Tcoh_s = Tcoh_s;
linkMetadata.etaPos = etaPos;
linkMetadata.Branging_Hz = Branging_Hz;
linkMetadata.GasAnnualExceedance_pct = GasAnnualExceedance_pct;
linkMetadata.PropagationModel = PropagationModel;
linkMetadata.MaxCN0DifferenceAtT0_dB = maxCN0Difference_dB;
linkMetadata.MaxGasDifferenceAtT0_dB = maxGasDifference_dB;
linkMetadata.MaxSigmaDifferenceAtT0_m = maxSigmaDifference_m;
linkMetadata.ValidationPassed = validationPassed;

matFile = fullfile( ...
    outDir, ...
    'scenario4_link_budget_4epochs.mat');

save( ...
    matFile, ...
    'linkEpochs', ...
    'linkMetadata', ...
    'TepochSummary', ...
    'TarchitectureSummary', ...
    'PoolIndex', ...
    'Label', ...
    'Architecture', ...
    'NORAD_CAT_ID', ...
    'Frequency_GHz', ...
    'ReceiverFoM_dBK', ...
    'ReceiverChain', ...
    'ParameterSource', ...
    'timeOffset_s', ...
    'timeUTC', ...
    'beta_Hz', ...
    'Tcoh_s', ...
    'etaPos', ...
    'Branging_Hz', ...
    'Branging_MHz', ...
    'GasAnnualExceedance_pct', ...
    'PropagationModel', ...
    '-v7.3');

fprintf('\n============================================================\n');
fprintf('TEMPORAL LINK BUDGET COMPLETED\n');
fprintf('============================================================\n');
disp(TepochSummary);

fprintf('MAT file:\n  %s\n',matFile);
fprintf('Validation passed: %s\n',string(validationPassed));
fprintf('============================================================\n');

%% ========================================================================
% Local function - Starlink EIRP interpolation
% =========================================================================

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

%% ========================================================================
% Local function - ITU-R P.618/P.676 gaseous loss
% =========================================================================

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
            'The P.618 model used in this script requires ' ...
            'elevation >= 5 deg.']);
    end

    N = numel(Frequency_GHz);
    Ag_dB = nan(N,1);

    cfgP = p618Config;

    cfgP.Latitude = lat_deg;
    cfgP.Longitude = lon_deg;
    cfgP.GasAnnualExceedance = GasAnnualExceedance_pct;

    for ii = 1:N

        cfgP.Frequency = Frequency_GHz(ii)*1e9;
        cfgP.ElevationAngle = Elevation_deg(ii);

        pl = p618PropagationLosses( ...
            cfgP, ...
            'StationHeight', ...
            userH_m/1e3);

        Ag_dB(ii) = pl.Ag;
    end

    if any(~isfinite(Ag_dB)) || any(Ag_dB < 0)
        error('Invalid ITU-R gaseous attenuation.');
    end
end
