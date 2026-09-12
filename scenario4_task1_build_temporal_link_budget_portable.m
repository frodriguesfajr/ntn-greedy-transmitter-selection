%% scenario4_task1_build_temporal_link_budget_portable.m
% TASK CODE 1 - STEP 2 (PORTABLE)
% Builds four temporal link-budget / measurement-quality sets without
% Satellite Communications Toolbox.
%
% Input:
%   results/temporal_geometry_4epochs/scenario4_geometry_4epochs.mat
%   results/link_budget.csv
%
% Epochs:
%   t = 0, 30, 60, 90 s
%
% Portable gaseous-loss update:
%   Lgas_i(t) = Lgas_i(0) * sin(el_i(0))/sin(el_i(t))
%
% Thus t=0 reproduces the previously validated ITU-R gaseous loss exactly,
% while later epochs use a short-term elevation-scaled approximation.
%
% Output:
%   results/temporal_link_budget_4epochs_portable/
%       scenario4_link_budget_4epochs_portable.mat
%       link_budget_t000s.csv
%       link_budget_t030s.csv
%       link_budget_t060s.csv
%       link_budget_t090s.csv
%       link_budget_epoch_summary.csv
%
% This script DOES NOT run FTS or BTS.

close all;
clear;
clc;
format long;

%% Repository
paperRoot = '/MATLAB Drive/ntn-greedy-transmitter-selection';
assert(isfolder(paperRoot),'Paper repository not found.');
cd(paperRoot);

if isfile(fullfile(paperRoot,'setup_paths.m'))
    setup_paths;
end

geometryMat = fullfile(paperRoot,'results','temporal_geometry_4epochs', ...
    'scenario4_geometry_4epochs.mat');
nominalLinkCsv = fullfile(paperRoot,'results','link_budget.csv');
outDir = fullfile(paperRoot,'results','temporal_link_budget_4epochs_portable');

assert(isfile(geometryMat),'Temporal geometry MAT not found.');
assert(isfile(nominalLinkCsv),'Nominal link_budget.csv not found.');
if ~isfolder(outDir), mkdir(outDir); end

%% Load geometry and nominal link budget
G = load(geometryMat);
Tref = readtable(nominalLinkCsv,'TextType','string');

PoolIndex = double(G.PoolIndex(:));
Label = string(G.Label(:));
Architecture = upper(string(G.Architecture(:)));
NORAD_CAT_ID = string(G.NORAD_CAT_ID(:));
timeOffset_s = double(G.timeOffset_s(:));
timeUTC = G.timeUTC(:);
epochs = G.epochs;

N = numel(Label);
Ntime = numel(timeOffset_s);
assert(N==26,'Expected 26 candidates.');
assert(isequal(timeOffset_s,[0;30;60;90]),'Expected t=[0 30 60 90] s.');

if ismember("PoolIndex",string(Tref.Properties.VariableNames))
    Tref = sortrows(Tref,'PoolIndex');
end

req = ["Label","Architecture","CN0_dBHz", ...
       "AtmosphericGasLoss_dB","SigmaRho_m"];
miss = setdiff(req,string(Tref.Properties.VariableNames));
assert(isempty(miss),'Missing nominal link-budget columns: %s',strjoin(miss,', '));
assert(height(Tref)==N,'Nominal link-budget size mismatch.');
assert(all(string(Tref.Label(:))==Label),'Nominal link-budget row ordering mismatch.');

%% Ranging assumptions
c_mps = 299792458;
beta_Hz = 1.023e6;
Tcoh_s = 20e-3;
etaPos = 0.01;
Branging_MHz = beta_Hz/1e6;

%% Architecture RF parameters used in current paper implementation
Frequency_GHz = nan(N,1);
ReceiverFoM_dBK = nan(N,1);

% HAPS
idx = Architecture=="HAPS";
Frequency_GHz(idx) = 31.15;
ReceiverFoM_dBK(idx) = 12.3;

% LEO
idx = Architecture=="LEO";
Frequency_GHz(idx) = 20.0;
ReceiverFoM_dBK(idx) = 10.5;
LEO_ElevationControl_deg = [10;25;90];
LEO_EIRP_Control_dBW_MHz = [-2.7;17.9;11.3];

% MEO
idx = Architecture=="MEO";
Frequency_GHz(idx) = 20.0;
ReceiverFoM_dBK(idx) = 18.4;
MEO_EIRP_total_dBW = 49.7;
MEO_BW_MHz = 100;

% GEO
idx = Architecture=="GEO";
Frequency_GHz(idx) = 19.95;
ReceiverFoM_dBK(idx) = 44.5 - 10*log10(250);
GEO_EIRP_total_dBW = 57.0;
GEO_BW_MHz = 81.0;

assert(all(isfinite(Frequency_GHz)) && all(isfinite(ReceiverFoM_dBK)), ...
    'Invalid RF parameter vectors.');

%% Reference values at t=0
refCN0_dBHz = double(Tref.CN0_dBHz(:));
refGas_dB = double(Tref.AtmosphericGasLoss_dB(:));
refSigmaRho_m = double(Tref.SigmaRho_m(:));
Elevation0_deg = double(epochs(1).Elevation_deg(:));

assert(all(Elevation0_deg>0),'Invalid t=0 elevation.');
assert(all(isfinite(refGas_dB)) && all(refGas_dB>=0), ...
    'Invalid stored gaseous attenuation.');

%% Allocate output structure
linkEpochs = repmat(struct( ...
    'TimeOffset_s',[], ...
    'TimeUTC',[], ...
    'Visible',[], ...
    'VisibleIndices',[], ...
    'VisibleCount',[], ...
    'EIRP_PSD_dBW_MHz',[], ...
    'AtmosphericGasLoss_dB',[], ...
    'FSPL_dB',[], ...
    'q_link_dBHz',[], ...
    'CN0_positioning_dBHz',[], ...
    'Gamma_Hz',[], ...
    'SigmaRho_m',[], ...
    'VarianceRho_m2',[], ...
    'R_C',[], ...
    'R_Visible',[], ...
    'LinkTable',[]),Ntime,1);

summaryRows = [];

fprintf('\n============================================================\n');
fprintf('TASK CODE 1 - STEP 2: PORTABLE TEMPORAL LINK BUDGET\n');
fprintf('============================================================\n');
fprintf('Epochs  : %s s\n',mat2str(timeOffset_s.'));
fprintf('beta=W  : %.3f MHz\n',beta_Hz/1e6);
fprintf('Tcoh    : %.1f ms\n',1e3*Tcoh_s);
fprintf('eta_Pos : %.2f %%\n',100*etaPos);
fprintf('Gas     : validated t0 loss + sin(el0)/sin(el(t)) scaling\n');

%% Four epochs
for k = 1:Ntime
    Elevation_deg = double(epochs(k).Elevation_deg(:));
    SlantRange_m = double(epochs(k).SlantRange_m(:));
    Visible = logical(epochs(k).Visible(:));

    assert(all(Elevation_deg>0),'Portable gas scaling requires positive elevation.');

    % EIRP spectral density by architecture
    EIRP_PSD_dBW_MHz = nan(N,1);

    idx = Architecture=="HAPS";
    EIRP_PSD_dBW_MHz(idx) = -9.1;

    idx = Architecture=="LEO";
    elevLEO = Elevation_deg(idx);
    elevLEO = min(max(elevLEO,min(LEO_ElevationControl_deg)), ...
                  max(LEO_ElevationControl_deg));
    EIRP_PSD_dBW_MHz(idx) = interp1( ...
        LEO_ElevationControl_deg,LEO_EIRP_Control_dBW_MHz, ...
        elevLEO,'linear');

    idx = Architecture=="MEO";
    EIRP_PSD_dBW_MHz(idx) = MEO_EIRP_total_dBW - 10*log10(MEO_BW_MHz);

    idx = Architecture=="GEO";
    EIRP_PSD_dBW_MHz(idx) = GEO_EIRP_total_dBW - 10*log10(GEO_BW_MHz);

    assert(all(isfinite(EIRP_PSD_dBW_MHz)),'Invalid EIRP PSD.');

    % Portable gaseous attenuation: exact at t=0 by construction
    AtmosphericGasLoss_dB = refGas_dB .* ...
        sind(Elevation0_deg) ./ sind(Elevation_deg);

    % Free-space path loss
    SlantRange_km = SlantRange_m/1e3;
    FSPL_dB = 92.45 + 20*log10(SlantRange_km) + 20*log10(Frequency_GHz);

    % EIRP integrated over adopted W = beta
    EIRP_band_dBW = EIRP_PSD_dBW_MHz + 10*log10(Branging_MHz);

    % q_i: link-budget C/N0 before positioning-power allocation [dB-Hz]
    q_link_dBHz = EIRP_band_dBW ...
        - FSPL_dB ...
        - AtmosphericGasLoss_dB ...
        + ReceiverFoM_dBK ...
        + 228.6;

    % Positioning-signal allocation
    CN0_positioning_dBHz = q_link_dBHz + 10*log10(etaPos);
    Gamma_Hz = etaPos .* 10.^(q_link_dBHz/10);

    % Pseudorange uncertainty
    SigmaRho_m = c_mps ./ ...
        (2*pi*beta_Hz .* sqrt(Gamma_Hz*Tcoh_s));
    VarianceRho_m2 = SigmaRho_m.^2;

    assert(all(isfinite(SigmaRho_m)) && all(SigmaRho_m>0), ...
        'Invalid sigma_rho.');

    R_C = diag(VarianceRho_m2);
    visibleIdx = find(Visible);
    R_Visible = diag(VarianceRho_m2(visibleIdx));

    LinkTable = table( ...
        PoolIndex,Label,Architecture,NORAD_CAT_ID,Visible, ...
        Elevation_deg,SlantRange_km,Frequency_GHz, ...
        EIRP_PSD_dBW_MHz,ReceiverFoM_dBK, ...
        AtmosphericGasLoss_dB,FSPL_dB,EIRP_band_dBW, ...
        q_link_dBHz,CN0_positioning_dBHz,Gamma_Hz, ...
        SigmaRho_m,VarianceRho_m2);

    linkEpochs(k).TimeOffset_s = timeOffset_s(k);
    linkEpochs(k).TimeUTC = timeUTC(k);
    linkEpochs(k).Visible = Visible;
    linkEpochs(k).VisibleIndices = visibleIdx;
    linkEpochs(k).VisibleCount = numel(visibleIdx);
    linkEpochs(k).EIRP_PSD_dBW_MHz = EIRP_PSD_dBW_MHz;
    linkEpochs(k).AtmosphericGasLoss_dB = AtmosphericGasLoss_dB;
    linkEpochs(k).FSPL_dB = FSPL_dB;
    linkEpochs(k).q_link_dBHz = q_link_dBHz;
    linkEpochs(k).CN0_positioning_dBHz = CN0_positioning_dBHz;
    linkEpochs(k).Gamma_Hz = Gamma_Hz;
    linkEpochs(k).SigmaRho_m = SigmaRho_m;
    linkEpochs(k).VarianceRho_m2 = VarianceRho_m2;
    linkEpochs(k).R_C = R_C;
    linkEpochs(k).R_Visible = R_Visible;
    linkEpochs(k).LinkTable = LinkTable;

    writetable(LinkTable,fullfile(outDir, ...
        sprintf('link_budget_t%03ds.csv',timeOffset_s(k))));

    row = table( ...
        timeOffset_s(k),timeUTC(k),numel(visibleIdx), ...
        min(CN0_positioning_dBHz(Visible)), ...
        median(CN0_positioning_dBHz(Visible)), ...
        max(CN0_positioning_dBHz(Visible)), ...
        min(SigmaRho_m(Visible)), ...
        median(SigmaRho_m(Visible)), ...
        max(SigmaRho_m(Visible)), ...
        'VariableNames',{'TimeOffset_s','TimeUTC','VisibleCount', ...
        'CN0_Min_dBHz','CN0_Median_dBHz','CN0_Max_dBHz', ...
        'SigmaRho_Min_m','SigmaRho_Median_m','SigmaRho_Max_m'});

    summaryRows = [summaryRows;row]; %#ok<AGROW>

    fprintf('\n--- t = %3d s ---\n',timeOffset_s(k));
    fprintf('Visible            : %d / %d\n',numel(visibleIdx),N);
    fprintf('C/N0 min/med/max   : %.3f / %.3f / %.3f dB-Hz\n', ...
        min(CN0_positioning_dBHz(Visible)), ...
        median(CN0_positioning_dBHz(Visible)), ...
        max(CN0_positioning_dBHz(Visible)));
    fprintf('sigma min/med/max  : %.6f / %.6f / %.6f m\n', ...
        min(SigmaRho_m(Visible)), ...
        median(SigmaRho_m(Visible)), ...
        max(SigmaRho_m(Visible)));
end

%% Validate t=0 against existing nominal reference
maxCN0Difference_dB = max(abs( ...
    linkEpochs(1).CN0_positioning_dBHz - refCN0_dBHz));
maxGasDifference_dB = max(abs( ...
    linkEpochs(1).AtmosphericGasLoss_dB - refGas_dB));
maxSigmaDifference_m = max(abs( ...
    linkEpochs(1).SigmaRho_m - refSigmaRho_m));

fprintf('\n============================================================\n');
fprintf('t = 0 PORTABLE LINK-BUDGET VALIDATION\n');
fprintf('============================================================\n');
fprintf('max |Delta C/N0|     = %.6e dB\n',maxCN0Difference_dB);
fprintf('max |Delta gas loss| = %.6e dB\n',maxGasDifference_dB);
fprintf('max |Delta sigma|    = %.6e m\n',maxSigmaDifference_m);

validationPassed = ...
    maxCN0Difference_dB <= 1e-4 && ...
    maxGasDifference_dB <= 1e-12 && ...
    maxSigmaDifference_m <= 1e-6;

if validationPassed
    fprintf('Status               = VALIDATED\n');
else
    warning('Portable t=0 result differs from nominal reference.');
end

%% Save
TepochSummary = summaryRows;
writetable(TepochSummary,fullfile(outDir,'link_budget_epoch_summary.csv'));

linkMetadata = struct;
linkMetadata.GeometryInput = geometryMat;
linkMetadata.NominalReference = nominalLinkCsv;
linkMetadata.TimeOffset_s = timeOffset_s;
linkMetadata.beta_Hz = beta_Hz;
linkMetadata.Tcoh_s = Tcoh_s;
linkMetadata.etaPos = etaPos;
linkMetadata.PortableGasModel = ...
    'Lgas(t)=Lgas(t0)*sin(el(t0))/sin(el(t)); t0 from validated ITU-R result';
linkMetadata.ValidationPassed = validationPassed;
linkMetadata.MaxCN0DifferenceAtT0_dB = maxCN0Difference_dB;
linkMetadata.MaxGasDifferenceAtT0_dB = maxGasDifference_dB;
linkMetadata.MaxSigmaDifferenceAtT0_m = maxSigmaDifference_m;

matFile = fullfile(outDir,'scenario4_link_budget_4epochs_portable.mat');
save(matFile,'linkEpochs','linkMetadata','TepochSummary', ...
    'PoolIndex','Label','Architecture','NORAD_CAT_ID', ...
    'Frequency_GHz','ReceiverFoM_dBK','timeOffset_s','timeUTC', ...
    'beta_Hz','Tcoh_s','etaPos','-v7.3');

fprintf('\n============================================================\n');
fprintf('PORTABLE TEMPORAL LINK BUDGET COMPLETED\n');
fprintf('============================================================\n');
disp(TepochSummary);
fprintf('MAT file:\n  %s\n',matFile);
fprintf('Validation passed: %s\n',string(validationPassed));
fprintf('============================================================\n');
