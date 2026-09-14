clear;
clc;
format long;

%% ============================================================
% EXPERIMENT 1 - SHORT-TERM TEMPORAL EVOLUTION
%
% Repository location:
%   matlab/experiments/run_experiment1_temporal.m
%
% Purpose:
%   Reproduce the temporal experiment used in the paper at
%
%       t = {0, 60, 120, 180} s
%
%   using the CURRENT validated atmospheric/link model.
%
% At each epoch:
%   1) propagate the same nominal 26 transmitter identities;
%   2) keep the four synthetic HAPS fixed in ECEF;
%   3) recompute azimuth, elevation, visibility, range;
%   4) recompute P.676-13 gaseous attenuation using the same
%      P.2145 annual-median meteorology as the nominal link budget;
%   5) recompute C/N0 and sigma_rho;
%   6) run Forward Adding (FA);
%   7) run Backward Elimination (BE).
%
% IMPORTANT:
%   - This script computes NUMERICAL RESULTS ONLY.
%   - It does NOT generate or modify the paper plots.
%   - The output contains both current names (FA/BE) and legacy
%     aliases (FTS/BTS) so the already-tuned plotting script can
%     be reused with only path/legend adjustments if desired.
%
% Required inputs:
%   results/candidate_pool_26.mat
%   results/link_budget.mat
%   results/FA_results.mat
%   results/BE_results.mat
%   archived TLE files
%
% Outputs:
%   results/experiment1_temporal/
%       experiment1_geometry.mat
%       experiment1_results.mat
%       experiment1_summary.csv
%       experiment1_epoch_t000s.csv
%       experiment1_epoch_t060s.csv
%       experiment1_epoch_t120s.csv
%       experiment1_epoch_t180s.csv
% ============================================================

%% ============================================================
% Locate repository root robustly
% =============================================================

thisFile = mfilename('fullpath');

if isempty(thisFile)
    error('Save this script before running it.');
end

repoRoot = fileparts(thisFile);

while ~isfolder(fullfile(repoRoot,'.git'))

    parentDir = fileparts(repoRoot);

    if strcmp(parentDir,repoRoot)
        error('Could not locate Git repository root from:\n%s',thisFile);
    end

    repoRoot = parentDir;
end

setupFile = fullfile( ...
    repoRoot,'matlab','setup','setup_paths.m');

assert(isfile(setupFile), ...
    'setup_paths.m not found:\n%s',setupFile);

addpath(fullfile(repoRoot,'matlab','setup'),'-begin');
setup_paths;
cd(repoRoot);

fprintf('Repository root: %s\n',repoRoot);

% Critical functions required by Experiment 1.
requiredFunctions = { ...
    'paper_epoch_eop', ...
    'propagate_tle_vallado', ...
    'itu676_annex2_slant_gas_loss'};

for kReq = 1:numel(requiredFunctions)

    if exist(requiredFunctions{kReq},'file') == 0

        error([ ...
            'Required function "%s" was not found on the MATLAB path.\n' ...
            'Repository root: %s\n' ...
            'MATLAB source root: %s'], ...
            requiredFunctions{kReq}, ...
            repoRoot, ...
            matlabRoot);
    end
end

fprintf('\n');
fprintf('============================================================\n');
fprintf('EXPERIMENT 1 - TEMPORAL EVOLUTION\n');
fprintf('============================================================\n');

%% ============================================================
% Inputs
% =============================================================

candidateFile = fullfile( ...
    repoRoot,'results','candidate_pool_26.mat');

linkBudgetFile = fullfile( ...
    repoRoot,'results','link_budget.mat');

faNominalFile = fullfile( ...
    repoRoot,'results','FA_results.mat');

beNominalFile = fullfile( ...
    repoRoot,'results','BE_results.mat');

assert(isfile(candidateFile), ...
    'Candidate pool not found:\n%s',candidateFile);

assert(isfile(linkBudgetFile), ...
    'Link budget not found:\n%s',linkBudgetFile);

assert(isfile(faNominalFile), ...
    'Nominal FA result not found:\n%s',faNominalFile);

assert(isfile(beNominalFile), ...
    'Nominal BE result not found:\n%s',beNominalFile);

S = load(candidateFile);
LB0 = load(linkBudgetFile);
FA0 = load(faNominalFile);
BE0 = load(beNominalFile);

%% ============================================================
% Candidate pool
% =============================================================

requiredCandidateVars = [ ...
    "UserPositionECEF", ...
    "PositionECEF", ...
    "Label", ...
    "Architecture", ...
    "NORAD", ...
    "ObjectName"];

for k = 1:numel(requiredCandidateVars)
    assert(isfield(S,requiredCandidateVars(k)), ...
        'Missing candidate-pool variable: %s', ...
        requiredCandidateVars(k));
end

UserPositionECEF = double(S.UserPositionECEF(:).');
PositionECEF_t0 = double(S.PositionECEF);

Label = string(S.Label(:));
Architecture = upper(string(S.Architecture(:)));
NORAD = string(S.NORAD(:));
ObjectName = string(S.ObjectName(:));

Npool = size(PositionECEF_t0,1);

assert(Npool == 26, ...
    'Expected 26 candidates; found %d.',Npool);

PoolIndex = (1:Npool).';

%% User coordinates

if isfield(S,'userLat_deg')
    userLat_deg = double(S.userLat_deg);
elseif isfield(LB0,'userLat_deg')
    userLat_deg = double(LB0.userLat_deg);
else
    userLat_deg = -22.8596582;
end

if isfield(S,'userLon_deg')
    userLon_deg = double(S.userLon_deg);
elseif isfield(LB0,'userLon_deg')
    userLon_deg = double(LB0.userLon_deg);
else
    userLon_deg = -43.2303236;
end

if isfield(S,'userH_m')
    userH_m = double(S.userH_m);
elseif isfield(LB0,'userH_m')
    userH_m = double(LB0.userH_m);
else
    userH_m = 10;
end

%% Reference epoch

if isfield(S,'analysisTimeUTC')
    analysisTimeUTC = S.analysisTimeUTC;
elseif isfield(LB0,'analysisTimeUTC')
    analysisTimeUTC = LB0.analysisTimeUTC;
else
    analysisTimeUTC = datetime( ...
        2026,8,1,12,0,0,'TimeZone','UTC');
end

if isempty(analysisTimeUTC.TimeZone)
    analysisTimeUTC.TimeZone = 'UTC';
end

%% ============================================================
% Link-budget parameters
% =============================================================

requiredLinkVars = [ ...
    "Label", ...
    "Architecture", ...
    "CN0_dBHz", ...
    "SigmaRho_m", ...
    "beta", ...
    "Tcoh", ...
    "etaPos"];

for k = 1:numel(requiredLinkVars)
    assert(isfield(LB0,requiredLinkVars(k)), ...
        'Missing nominal link-budget variable: %s', ...
        requiredLinkVars(k));
end

assert(isequal(Label,string(LB0.Label(:))), ...
    'Candidate/link-budget label order mismatch.');

assert(isequal(Architecture,upper(string(LB0.Architecture(:)))), ...
    'Candidate/link-budget architecture order mismatch.');

nominalCN0_dBHz = double(LB0.CN0_dBHz(:));
nominalSigmaRho_m = double(LB0.SigmaRho_m(:));

beta = double(LB0.beta);
Tcoh = double(LB0.Tcoh);
etaPos = double(LB0.etaPos);

c_mps = 299792458;
Branging_MHz = beta/1e6;

%% Current RF assumptions - same as nominal link-budget builder

HAPS_Frequency_GHz = 31.15;
HAPS_EIRP_PSD_dBW_MHz = -9.1;
HAPS_ReceiverFoM_dBK = 12.3;

LEO_Frequency_GHz = 20.0;
LEO_ReceiverFoM_dBK = 10.5;
LEO_ElevationControl_deg = [10;25;90];
LEO_EIRP_Control_dBW_MHz = [-2.7;17.9;11.3];

MEO_Frequency_GHz = 20.0;
EIRP_MEO_total_dBW = 49.7;
BW_MEO_MHz = 100;
MEO_ReceiverFoM_dBK = 18.4;

GEO_Frequency_GHz = 19.95;
EIRP_GEO_total_dBW = 57.0;
BW_GEO_MHz = 81.0;
G_rx_GEO_dBi = 44.5;
Tsys_GEO_K = 250;

%% ============================================================
% P.2145 meteorology
% =============================================================

if isfield(LB0,'P2145_meteorology')

    met = LB0.P2145_meteorology;

elseif isfield(LB0,'met')

    met = LB0.met;

else

    atmosphereDir = fullfile(repoRoot,'matlab','atmosphere');

    if ~isfolder(atmosphereDir)
        atmosphereDir = fullfile(repoRoot,'atmosphere');
    end

    if isfolder(atmosphereDir)
        addpath(atmosphereDir,'-begin');
    end

    p2145Dir = fullfile( ...
        repoRoot,'data','itu','p2145','annual');

    assert(isfolder(p2145Dir), ...
        'P.2145 data directory not found:\n%s',p2145Dir);

    assert(exist('itu2145_annual50_at_location','file') == 2, ...
        'itu2145_annual50_at_location.m not found.');

    met = itu2145_annual50_at_location( ...
        p2145Dir,userLat_deg,userLon_deg,userH_m);
end

assert(exist('itu676_annex2_slant_gas_loss','file') == 2, ...
    'itu676_annex2_slant_gas_loss.m not found on MATLAB path.');

%% ============================================================
% Temporal grid
% =============================================================

timeOffset_s = [0;60;120;180];

timeUTC = ...
    analysisTimeUTC + seconds(timeOffset_s);

timeUTC.Format = 'yyyy-MM-dd HH:mm:ss';

Ntime = numel(timeOffset_s);

target_m = 0.6;
targetTol = 1e-12;

fprintf('Reference UTC : %s\n',char(string(analysisTimeUTC)));
fprintf('Epochs        : %s s\n',mat2str(timeOffset_s.'));
fprintf('Candidates    : %d\n',Npool);
fprintf('Nominal FA K  : %d\n',numel(FA0.selected));
fprintf('Nominal BE K  : %d\n',numel(BE0.selected));
fprintf('beta          : %.6f MHz\n',beta/1e6);
fprintf('Tcoh          : %.3f ms\n',1e3*Tcoh);
fprintf('etaPos        : %.6f\n',etaPos);

%% ============================================================
% Elevation masks
% =============================================================

mask_deg = nan(Npool,1);

mask_deg(Architecture=="HAPS") = 15;
mask_deg(Architecture=="LEO")  = 5;
mask_deg(Architecture=="MEO")  = 5;
mask_deg(Architecture=="GEO")  = 5;

assert(all(isfinite(mask_deg)), ...
    'Unsupported architecture in candidate pool.');

%% ============================================================
% Locate TLE root
% =============================================================

tleCandidates = { ...
    fullfile(repoRoot,'tle'), ...
    fullfile(repoRoot,'data','tle')};

tleRoot = "";

for k = 1:numel(tleCandidates)

    if isfolder(tleCandidates{k})
        tleRoot = string(tleCandidates{k});
        break;
    end
end

assert(strlength(tleRoot)>0, ...
    'Could not locate TLE root under repository.');

assert(exist('paper_epoch_eop','file') ~= 0, ...
    'paper_epoch_eop.m is not on the MATLAB path.');

assert(exist('propagate_tle_vallado','file') ~= 0, ...
    'propagate_tle_vallado.m is not on the MATLAB path.');

%% ============================================================
% Propagate same 26 identities
% =============================================================

PositionECEF_m = nan(Npool,3,Ntime);

eop = paper_epoch_eop();

for i = 1:Npool

    if Architecture(i) == "HAPS"

        for tidx = 1:Ntime
            PositionECEF_m(i,:,tidx) = PositionECEF_t0(i,:);
        end

    else

        norad = strip(NORAD(i));

        assert(strlength(norad)>0, ...
            'Missing NORAD ID for %s.',Label(i));

        tleFile = fullfile( ...
            char(tleRoot), ...
            char(Architecture(i)), ...
            char(norad), ...
            [char(norad) '.tle']);

        assert(isfile(tleFile), ...
            'TLE file not found for %s:\n%s', ...
            Label(i),tleFile);

        for tidx = 1:Ntime

            state = propagate_tle_vallado( ...
                tleFile,timeUTC(tidx),eop);

            if isfield(state,'sgp4_error')
                assert(state.sgp4_error == 0, ...
                    'SGP4 error %d for %s.', ...
                    state.sgp4_error,Label(i));
            end

            PositionECEF_m(i,:,tidx) = ...
                state.rECEF_m(:).';
        end
    end
end

%% t0 geometry validation

dPos_t0_m = ...
    vecnorm(PositionECEF_m(:,:,1)-PositionECEF_t0,2,2);

maxPositionDifferenceT0_m = max(dPos_t0_m);

fprintf('Max t0 position difference: %.9f m\n', ...
    maxPositionDifferenceT0_m);

assert(maxPositionDifferenceT0_m <= 0.05, ...
    'Temporal propagation does not reproduce nominal t0 geometry.');

%% ============================================================
% Allocate results
% =============================================================

VisibleCount = zeros(Ntime,1);
NoSelection_AlphaCRB_m = nan(Ntime,1);
NoSelection_PDOP = nan(Ntime,1);

FA_K = nan(Ntime,1);
FA_AlphaCRB_m = nan(Ntime,1);
FA_PDOP = nan(Ntime,1);
FA_Evaluations = nan(Ntime,1);

BE_K = nan(Ntime,1);
BE_AlphaCRB_m = nan(Ntime,1);
BE_PDOP = nan(Ntime,1);
BE_Evaluations = nan(Ntime,1);

FA_HAPS = zeros(Ntime,1);
FA_LEO  = zeros(Ntime,1);
FA_MEO  = zeros(Ntime,1);
FA_GEO  = zeros(Ntime,1);

BE_HAPS = zeros(Ntime,1);
BE_LEO  = zeros(Ntime,1);
BE_MEO  = zeros(Ntime,1);
BE_GEO  = zeros(Ntime,1);

SelectedFA = false(Npool,Ntime);
SelectedBE = false(Npool,Ntime);

FA_Labels = strings(Ntime,1);
BE_Labels = strings(Ntime,1);

epochs = repmat(struct( ...
    'TimeOffset_s',[], ...
    'TimeUTC',[], ...
    'PositionECEF_m',[], ...
    'SlantRange_m',[], ...
    'Azimuth_deg',[], ...
    'Elevation_deg',[], ...
    'ElevationMask_deg',[], ...
    'Visible',[], ...
    'H',[], ...
    'CN0_dBHz',[], ...
    'SigmaRho_m',[], ...
    'LinkTable',[]),Ntime,1);

%% Output directory

outDir = fullfile( ...
    repoRoot,'results','experiment1_temporal');

if ~isfolder(outDir)
    mkdir(outDir);
end

%% ============================================================
% Main epoch loop
% =============================================================

for tidx = 1:Ntime

    fprintf('\n');
    fprintf('============================================================\n');
    fprintf('EPOCH t = %d s\n',timeOffset_s(tidx));
    fprintf('============================================================\n');

    positionNow_m = PositionECEF_m(:,:,tidx);

    [Hfull,slantRange_m,azimuth_deg,elevation_deg] = ...
        geometryFromECEF( ...
            positionNow_m, ...
            UserPositionECEF, ...
            userLat_deg, ...
            userLon_deg);

    visible = elevation_deg >= mask_deg;
    visibleIdx = find(visible);

    VisibleCount(tidx) = numel(visibleIdx);

    [ ...
        Frequency_GHz, ...
        EIRP_PSD_dBW_MHz, ...
        ReceiverFoM_dBK, ...
        AtmosphericGasLoss_dB, ...
        FSPL_dB, ...
        CN0_dBHz, ...
        SigmaRho_m] = ...
        computeEpochLinkBudget( ...
            Architecture, ...
            elevation_deg, ...
            slantRange_m, ...
            visible, ...
            userH_m, ...
            met, ...
            c_mps,beta,Tcoh,etaPos,Branging_MHz, ...
            HAPS_Frequency_GHz,HAPS_EIRP_PSD_dBW_MHz, ...
            HAPS_ReceiverFoM_dBK, ...
            LEO_Frequency_GHz,LEO_ReceiverFoM_dBK, ...
            LEO_ElevationControl_deg,LEO_EIRP_Control_dBW_MHz, ...
            MEO_Frequency_GHz,EIRP_MEO_total_dBW,BW_MEO_MHz, ...
            MEO_ReceiverFoM_dBK, ...
            GEO_Frequency_GHz,EIRP_GEO_total_dBW,BW_GEO_MHz, ...
            G_rx_GEO_dBi,Tsys_GEO_K);

    %% t0 link-budget validation

    if tidx == 1

        assert(all(visible), ...
            'Nominal t0 pool is expected to be fully visible.');

        maxCN0Diff = ...
            max(abs(CN0_dBHz-nominalCN0_dBHz));

        maxSigmaDiff = ...
            max(abs(SigmaRho_m-nominalSigmaRho_m));

        fprintf('t0 max |Delta C/N0|      : %.3e dB-Hz\n', ...
            maxCN0Diff);

        fprintf('t0 max |Delta sigma_rho| : %.3e m\n', ...
            maxSigmaDiff);

        assert(maxCN0Diff < 1e-10);
        assert(maxSigmaDiff < 1e-10);
    end

    %% Full visible pool

    [fullAlpha,fullPDOP,fullRank] = ...
        subsetMetrics(visibleIdx,Hfull,SigmaRho_m);

    NoSelection_AlphaCRB_m(tidx) = fullAlpha;
    NoSelection_PDOP(tidx) = fullPDOP;

    assert(fullRank == 4, ...
        'Visible pool is rank deficient at t=%d s.', ...
        timeOffset_s(tidx));

    %% Forward Adding

    [faIdx,faAlpha,faPDOP,faEval,faFeasible] = ...
        forwardAdding( ...
            visibleIdx,Hfull,SigmaRho_m, ...
            target_m,targetTol);

    %% Backward Elimination

    [beIdx,beAlpha,bePDOP,beEval,beFeasible] = ...
        backwardElimination( ...
            visibleIdx,Hfull,SigmaRho_m, ...
            target_m,targetTol);

    assert(faFeasible, ...
        'FA failed target at t=%d s.',timeOffset_s(tidx));

    assert(beFeasible, ...
        'BE failed target at t=%d s.',timeOffset_s(tidx));

    %% t0 nominal validation

    if tidx == 1

        assert(isequal(sort(faIdx(:).'),sort(double(FA0.selected(:).'))), ...
            'Experiment 1 FA does not reproduce nominal t0 FA subset.');

        assert(isequal(sort(beIdx(:).'),sort(double(BE0.selected(:).'))), ...
            'Experiment 1 BE does not reproduce nominal t0 BE subset.');

        assert(abs(faAlpha-double(FA0.finalBound)) < 1e-10, ...
            'Experiment 1 FA t0 bound mismatch.');

        assert(abs(beAlpha-double(BE0.finalBound)) < 1e-10, ...
            'Experiment 1 BE t0 bound mismatch.');
    end

    %% Store selections

    SelectedFA(faIdx,tidx) = true;
    SelectedBE(beIdx,tidx) = true;

    FA_K(tidx) = numel(faIdx);
    FA_AlphaCRB_m(tidx) = faAlpha;
    FA_PDOP(tidx) = faPDOP;
    FA_Evaluations(tidx) = faEval;

    BE_K(tidx) = numel(beIdx);
    BE_AlphaCRB_m(tidx) = beAlpha;
    BE_PDOP(tidx) = bePDOP;
    BE_Evaluations(tidx) = beEval;

    FA_HAPS(tidx) = sum(Architecture(faIdx)=="HAPS");
    FA_LEO(tidx)  = sum(Architecture(faIdx)=="LEO");
    FA_MEO(tidx)  = sum(Architecture(faIdx)=="MEO");
    FA_GEO(tidx)  = sum(Architecture(faIdx)=="GEO");

    BE_HAPS(tidx) = sum(Architecture(beIdx)=="HAPS");
    BE_LEO(tidx)  = sum(Architecture(beIdx)=="LEO");
    BE_MEO(tidx)  = sum(Architecture(beIdx)=="MEO");
    BE_GEO(tidx)  = sum(Architecture(beIdx)=="GEO");

    FA_Labels(tidx) = strjoin(Label(faIdx)," ");
    BE_Labels(tidx) = strjoin(Label(beIdx)," ");

    %% Per-epoch table

    epochTable = table( ...
        PoolIndex, ...
        Label, ...
        Architecture, ...
        NORAD, ...
        ObjectName, ...
        positionNow_m(:,1), ...
        positionNow_m(:,2), ...
        positionNow_m(:,3), ...
        slantRange_m, ...
        azimuth_deg, ...
        elevation_deg, ...
        mask_deg, ...
        visible, ...
        Frequency_GHz, ...
        EIRP_PSD_dBW_MHz, ...
        ReceiverFoM_dBK, ...
        AtmosphericGasLoss_dB, ...
        FSPL_dB, ...
        CN0_dBHz, ...
        SigmaRho_m, ...
        SelectedFA(:,tidx), ...
        SelectedBE(:,tidx), ...
        'VariableNames',{ ...
        'PoolIndex', ...
        'Label', ...
        'Architecture', ...
        'NORAD_CAT_ID', ...
        'ObjectName', ...
        'X_ECEF_m', ...
        'Y_ECEF_m', ...
        'Z_ECEF_m', ...
        'SlantRange_m', ...
        'Azimuth_deg', ...
        'Elevation_deg', ...
        'ElevationMask_deg', ...
        'Visible', ...
        'Frequency_GHz', ...
        'EIRP_PSD_dBW_MHz', ...
        'ReceiverFoM_dBK', ...
        'AtmosphericGasLoss_dB', ...
        'FSPL_dB', ...
        'CN0_dBHz', ...
        'SigmaRho_m', ...
        'SelectedFA', ...
        'SelectedBE'});

    epochCsv = fullfile( ...
        outDir, ...
        sprintf('experiment1_epoch_t%03ds.csv', ...
        timeOffset_s(tidx)));

    writetable(epochTable,epochCsv);

    %% Geometry package for plotting

    epochs(tidx).TimeOffset_s = timeOffset_s(tidx);
    epochs(tidx).TimeUTC = timeUTC(tidx);
    epochs(tidx).PositionECEF_m = positionNow_m;
    epochs(tidx).SlantRange_m = slantRange_m;
    epochs(tidx).Azimuth_deg = azimuth_deg;
    epochs(tidx).Elevation_deg = elevation_deg;
    epochs(tidx).ElevationMask_deg = mask_deg;
    epochs(tidx).Visible = visible;
    epochs(tidx).H = Hfull;
    epochs(tidx).CN0_dBHz = CN0_dBHz;
    epochs(tidx).SigmaRho_m = SigmaRho_m;
    epochs(tidx).LinkTable = epochTable;

    %% Console

    fprintf('Visible       : %d / %d\n', ...
        numel(visibleIdx),Npool);

    fprintf('No selection  : alpha_CRB=%.12f m\n', ...
        fullAlpha);

    fprintf('FA             : K=%d | alpha_CRB=%.12f m | eval=%d\n', ...
        numel(faIdx),faAlpha,faEval);

    fprintf('BE             : K=%d | alpha_CRB=%.12f m | eval=%d\n', ...
        numel(beIdx),beAlpha,beEval);

    fprintf('FA labels      : %s\n',FA_Labels(tidx));
    fprintf('BE labels      : %s\n',BE_Labels(tidx));
end

%% ============================================================
% Current summary table
% =============================================================

TimeOffset_min = timeOffset_s/60;

Tsummary = table( ...
    timeOffset_s, ...
    TimeOffset_min, ...
    timeUTC, ...
    VisibleCount, ...
    NoSelection_AlphaCRB_m, ...
    NoSelection_PDOP, ...
    FA_K, ...
    FA_AlphaCRB_m, ...
    FA_PDOP, ...
    FA_Evaluations, ...
    BE_K, ...
    BE_AlphaCRB_m, ...
    BE_PDOP, ...
    BE_Evaluations, ...
    FA_HAPS,FA_LEO,FA_MEO,FA_GEO, ...
    BE_HAPS,BE_LEO,BE_MEO,BE_GEO, ...
    FA_Labels,BE_Labels);

%% ============================================================
% Legacy aliases for the already-tuned plot interface
%
% FTS = legacy name for current FA
% BTS = legacy name for current BE
% =============================================================

SelectedFTS = SelectedFA;
SelectedBTS = SelectedBE;

Tsummary.FTS_K = Tsummary.FA_K;
Tsummary.FTS_AlphaCRB_m = Tsummary.FA_AlphaCRB_m;
Tsummary.FTS_PDOP = Tsummary.FA_PDOP;
Tsummary.FTS_Evaluations = Tsummary.FA_Evaluations;

Tsummary.BTS_K = Tsummary.BE_K;
Tsummary.BTS_AlphaCRB_m = Tsummary.BE_AlphaCRB_m;
Tsummary.BTS_PDOP = Tsummary.BE_PDOP;
Tsummary.BTS_Evaluations = Tsummary.BE_Evaluations;

%% ============================================================
% Save
% =============================================================

geometryMat = fullfile( ...
    outDir,'experiment1_geometry.mat');

resultsMat = fullfile( ...
    outDir,'experiment1_results.mat');

summaryCsv = fullfile( ...
    outDir,'experiment1_summary.csv');

metadata = struct;
metadata.ReferenceUTC = analysisTimeUTC;
metadata.TimeOffset_s = timeOffset_s;
metadata.Target_m = target_m;
metadata.CandidateCount = Npool;
metadata.Atmosphere = ...
    'ITU-R P.2145 annual median + P.676-13';
metadata.Propagation = ...
    'Vallado SGP4/WGS-72';
metadata.BaselineK_FA = numel(FA0.selected);
metadata.BaselineK_BE = numel(BE0.selected);
metadata.MaxPositionDifferenceT0_m = ...
    maxPositionDifferenceT0_m;

save( ...
    geometryMat, ...
    'epochs', ...
    'timeOffset_s', ...
    'timeUTC', ...
    'Label', ...
    'Architecture', ...
    'NORAD', ...
    'metadata', ...
    '-v7.3');

save( ...
    resultsMat, ...
    'Tsummary', ...
    'SelectedFA', ...
    'SelectedBE', ...
    'SelectedFTS', ...
    'SelectedBTS', ...
    'timeOffset_s', ...
    'timeUTC', ...
    'Label', ...
    'Architecture', ...
    'target_m', ...
    'metadata', ...
    '-v7.3');

writetable(Tsummary,summaryCsv);

fprintf('\n');
fprintf('============================================================\n');
fprintf('EXPERIMENT 1 SUMMARY\n');
fprintf('============================================================\n');

disp(Tsummary(:,{ ...
    'timeOffset_s', ...
    'VisibleCount', ...
    'NoSelection_AlphaCRB_m', ...
    'FA_K','FA_AlphaCRB_m','FA_Evaluations', ...
    'BE_K','BE_AlphaCRB_m','BE_Evaluations'}));

fprintf('\nSaved:\n');
fprintf('  %s\n',geometryMat);
fprintf('  %s\n',resultsMat);
fprintf('  %s\n',summaryCsv);

fprintf('\n');
fprintf('============================================================\n');
fprintf('EXPERIMENT 1 COMPLETED\n');
fprintf('============================================================\n');

%% ============================================================
% Local function - geometry
% =============================================================

function [H,slantRange_m,azimuth_deg,elevation_deg] = ...
    geometryFromECEF( ...
        transmitterECEF_m, ...
        userECEF_m, ...
        userLat_deg, ...
        userLon_deg)

    userECEF_m = userECEF_m(:).';

    delta_m = transmitterECEF_m-userECEF_m;

    slantRange_m = vecnorm(delta_m,2,2);

    assert(all(isfinite(slantRange_m)));
    assert(all(slantRange_m > eps));

    u = delta_m ./ slantRange_m;

    H = [-u,ones(size(transmitterECEF_m,1),1)];

    lat = deg2rad(userLat_deg);
    lon = deg2rad(userLon_deg);

    Renu = [ ...
        -sin(lon),              cos(lon),             0; ...
        -sin(lat)*cos(lon),    -sin(lat)*sin(lon),    cos(lat); ...
         cos(lat)*cos(lon),     cos(lat)*sin(lon),    sin(lat)];

    enu_m = (Renu*delta_m.').';

    azimuth_deg = ...
        mod(atan2d(enu_m(:,1),enu_m(:,2)),360);

    elevation_deg = ...
        atan2d(enu_m(:,3),hypot(enu_m(:,1),enu_m(:,2)));
end

%% ============================================================
% Local function - current temporal link budget
% =============================================================

function [ ...
    Frequency_GHz, ...
    EIRP_PSD_dBW_MHz, ...
    ReceiverFoM_dBK, ...
    AtmosphericGasLoss_dB, ...
    FSPL_dB, ...
    CN0_dBHz, ...
    SigmaRho_m] = ...
    computeEpochLinkBudget( ...
        Architecture, ...
        Elevation_deg, ...
        SlantRange_m, ...
        Visible, ...
        userH_m, ...
        met, ...
        c_mps,beta,Tcoh,etaPos,Branging_MHz, ...
        HAPS_Frequency_GHz,HAPS_EIRP_PSD_dBW_MHz, ...
        HAPS_ReceiverFoM_dBK, ...
        LEO_Frequency_GHz,LEO_ReceiverFoM_dBK, ...
        LEO_ElevationControl_deg,LEO_EIRP_Control_dBW_MHz, ...
        MEO_Frequency_GHz,EIRP_MEO_total_dBW,BW_MEO_MHz, ...
        MEO_ReceiverFoM_dBK, ...
        GEO_Frequency_GHz,EIRP_GEO_total_dBW,BW_GEO_MHz, ...
        G_rx_GEO_dBi,Tsys_GEO_K)

    N = numel(Architecture);

    Frequency_GHz = nan(N,1);
    EIRP_PSD_dBW_MHz = nan(N,1);
    ReceiverFoM_dBK = nan(N,1);

    %% HAPS

    idx = Architecture=="HAPS";

    Frequency_GHz(idx) = HAPS_Frequency_GHz;
    EIRP_PSD_dBW_MHz(idx) = HAPS_EIRP_PSD_dBW_MHz;
    ReceiverFoM_dBK(idx) = HAPS_ReceiverFoM_dBK;

    %% LEO

    idx = Architecture=="LEO";

    Frequency_GHz(idx) = LEO_Frequency_GHz;
    ReceiverFoM_dBK(idx) = LEO_ReceiverFoM_dBK;

    elev = Elevation_deg(idx);

    elevClamped = min( ...
        max(elev,min(LEO_ElevationControl_deg)), ...
        max(LEO_ElevationControl_deg));

    EIRP_PSD_dBW_MHz(idx) = ...
        interp1( ...
            LEO_ElevationControl_deg, ...
            LEO_EIRP_Control_dBW_MHz, ...
            elevClamped,'linear');

    %% MEO

    idx = Architecture=="MEO";

    Frequency_GHz(idx) = MEO_Frequency_GHz;

    EIRP_PSD_dBW_MHz(idx) = ...
        EIRP_MEO_total_dBW ...
        - 10*log10(BW_MEO_MHz);

    ReceiverFoM_dBK(idx) = MEO_ReceiverFoM_dBK;

    %% GEO

    idx = Architecture=="GEO";

    Frequency_GHz(idx) = GEO_Frequency_GHz;

    EIRP_PSD_dBW_MHz(idx) = ...
        EIRP_GEO_total_dBW ...
        - 10*log10(BW_GEO_MHz);

    ReceiverFoM_dBK(idx) = ...
        G_rx_GEO_dBi ...
        - 10*log10(Tsys_GEO_K);

    %% Allocate visible-link outputs

    AtmosphericGasLoss_dB = nan(N,1);
    FSPL_dB = nan(N,1);
    CN0_dBHz = nan(N,1);
    SigmaRho_m = nan(N,1);

    visibleIdx = find(Visible);

    stationHeight_km = userH_m/1e3;

    for jj = 1:numel(visibleIdx)

        ii = visibleIdx(jj);

        AtmosphericGasLoss_dB(ii) = ...
            itu676_annex2_slant_gas_loss( ...
                Frequency_GHz(ii), ...
                Elevation_deg(ii), ...
                stationHeight_km, ...
                met);
    end

    SlantRange_km = ...
        SlantRange_m(visibleIdx)/1e3;

    FSPL_dB(visibleIdx) = ...
        92.45 ...
        + 20*log10(SlantRange_km) ...
        + 20*log10(Frequency_GHz(visibleIdx));

    EIRP_positioning_dBW = ...
        EIRP_PSD_dBW_MHz(visibleIdx) ...
        + 10*log10(Branging_MHz) ...
        + 10*log10(etaPos);

    CN0_dBHz(visibleIdx) = ...
        EIRP_positioning_dBW ...
        - FSPL_dB(visibleIdx) ...
        - AtmosphericGasLoss_dB(visibleIdx) ...
        + ReceiverFoM_dBK(visibleIdx) ...
        + 228.6;

    CN0_linear = ...
        10.^(CN0_dBHz(visibleIdx)/10);

    SigmaRho_m(visibleIdx) = ...
        c_mps ./ ...
        (2*pi*beta .* sqrt(CN0_linear*Tcoh));
end

%% ============================================================
% Local function - Forward Adding
% =============================================================

function [selected,finalBound,finalPDOP,evaluationCount,feasible] = ...
    forwardAdding( ...
        candidateIdx,Hfull,sigmaRho_m,target_m,targetTol)

    candidateIdx = candidateIdx(:).';

    N = numel(candidateIdx);

    if N < 4
        selected = candidateIdx;
        finalBound = inf;
        finalPDOP = inf;
        evaluationCount = 0;
        feasible = false;
        return;
    end

    combLocal = nchoosek(1:N,4);

    best4 = [];
    best4Bound = inf;

    evaluationCount = 0;

    for k = 1:size(combLocal,1)

        idx = candidateIdx(combLocal(k,:));

        [b,~,rankH] = ...
            subsetMetrics(idx,Hfull,sigmaRho_m);

        evaluationCount = evaluationCount + 1;

        if rankH == 4 && b < best4Bound
            best4Bound = b;
            best4 = idx;
        end
    end

    if isempty(best4)

        selected = [];
        finalBound = inf;
        finalPDOP = inf;
        feasible = false;
        return;
    end

    selected = best4(:).';
    currentBound = best4Bound;

    while currentBound > target_m + targetTol && ...
          numel(selected) < N

        remaining = ...
            setdiff(candidateIdx,selected,'stable');

        bestCandidate = NaN;
        bestCandidateBound = inf;

        for k = 1:numel(remaining)

            trial = [selected,remaining(k)];

            [b,~,rankH] = ...
                subsetMetrics(trial,Hfull,sigmaRho_m);

            evaluationCount = evaluationCount + 1;

            if rankH == 4 && b < bestCandidateBound

                bestCandidate = remaining(k);
                bestCandidateBound = b;
            end
        end

        if isnan(bestCandidate)
            break;
        end

        selected(end+1) = bestCandidate; %#ok<AGROW>
        currentBound = bestCandidateBound;
    end

    [finalBound,finalPDOP,rankH] = ...
        subsetMetrics(selected,Hfull,sigmaRho_m);

    feasible = ...
        rankH == 4 && ...
        finalBound <= target_m + targetTol;
end

%% ============================================================
% Local function - Backward Elimination
% =============================================================

function [selected,finalBound,finalPDOP,evaluationCount,feasible] = ...
    backwardElimination( ...
        candidateIdx,Hfull,sigmaRho_m,target_m,targetTol)

    selected = candidateIdx(:).';

    evaluationCount = 0;

    [fullBound,~,fullRank] = ...
        subsetMetrics(selected,Hfull,sigmaRho_m);

    if fullRank < 4 || fullBound > target_m + targetTol

        finalBound = fullBound;
        finalPDOP = inf;
        feasible = false;
        return;
    end

    while numel(selected) > 4

        bestRemoved = NaN;
        bestTrialBound = inf;

        for k = 1:numel(selected)

            candidate = selected(k);

            trial = selected(selected ~= candidate);

            [b,~,rankH] = ...
                subsetMetrics(trial,Hfull,sigmaRho_m);

            evaluationCount = evaluationCount + 1;

            if rankH ~= 4
                continue;
            end

            if b > target_m + targetTol
                continue;
            end

            if b < bestTrialBound

                bestTrialBound = b;
                bestRemoved = candidate;

            elseif abs(b-bestTrialBound) <= 1e-14

                if isnan(bestRemoved) || candidate < bestRemoved
                    bestTrialBound = b;
                    bestRemoved = candidate;
                end
            end
        end

        if isnan(bestRemoved)
            break;
        end

        selected(selected == bestRemoved) = [];
    end

    [finalBound,finalPDOP,rankH] = ...
        subsetMetrics(selected,Hfull,sigmaRho_m);

    feasible = ...
        rankH == 4 && ...
        finalBound <= target_m + targetTol;
end

%% ============================================================
% Local function - metrics
% =============================================================

function [alphaCRB_m,PDOP,rankH] = ...
    subsetMetrics(idx,Hfull,sigmaRho_m)

    idx = idx(:);

    if isempty(idx) || ...
       any(~isfinite(sigmaRho_m(idx)))

        alphaCRB_m = inf;
        PDOP = inf;
        rankH = 0;
        return;
    end

    H = Hfull(idx,:);
    sigma = sigmaRho_m(idx);

    rankH = rank(H);

    alphaCRB_m = inf;
    PDOP = inf;

    if rankH < 4
        return;
    end

    A = H.'*H;

    if rcond(A) > 1e-14

        Cgeom = A\eye(4);

        trGeom = trace(Cgeom(1:3,1:3));

        if trGeom > 0 && isfinite(trGeom)
            PDOP = sqrt(trGeom);
        end
    end

    R = diag(sigma.^2);

    J = H.'*(R\H);

    if rcond(J) <= 1e-14
        return;
    end

    C = J\eye(4);

    trPos = trace(C(1:3,1:3));

    if trPos > 0 && isfinite(trPos)
        alphaCRB_m = sqrt(trPos);
    end
end

