clear;
clc;
format long;

%% ============================================================
% EXPERIMENT 3 - HAPS AVAILABILITY AND LEO COMPENSATION
%
% Repository location:
%   experiments/run_experiment3_availability.m
%
% Purpose:
%   Recompute the availability experiment using the CURRENT
%   validated portable physical model.
%
% PART A - HAPS availability constraint
%   HAPS limits: {4, 2, 0}
%   Base nominal pool: 4 HAPS + 8 LEO + 7 MEO + 7 GEO
%   Both Forward Adding (FA) and Backward Elimination (BE)
%   are evaluated.
%
% PART B - LEO compensation without HAPS
%   Available LEO: {8, 12, 16, 20}
%   Fixed: 7 MEO + 7 GEO
%   No HAPS
%   Both FA and BE are evaluated.
%
% Physical model:
%   - nominal 26-Tx geometry from candidate_pool_26.mat
%   - expanded orbital geometry from orbital_catalog.mat
%   - portable ITU-R P.2145 annual-median meteorology
%   - portable ITU-R P.676-13 Annex 2 gaseous attenuation
%   - same RF assumptions as build_nominal_link_budget.m
%
% IMPORTANT:
%   - This script computes NUMERICAL RESULTS ONLY.
%   - It does not generate figures.
%   - Output column aliases FTS/BTS are retained only so the
%     already-tuned Experiment-3 plotting script remains compatible.
%     Numerically:
%         FTS -> FA
%         BTS -> BE
%
% Required inputs:
%   results/candidate_pool_26.mat
%   results/link_budget.mat
%   results/FA_results.mat
%   results/BE_results.mat
%   results/orbital_catalog.mat
%
% Outputs:
%   results/experiment3_availability/
%       experiment3_haps_summary.csv
%       experiment3_leo_summary.csv
%       experiment3_haps_selected.csv
%       experiment3_leo_selected.csv
%       experiment3_availability_results.mat
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

atmosphereDir = fullfile(repoRoot,'matlab','atmosphere');

if isfolder(atmosphereDir)
    addpath(atmosphereDir,'-begin');
end

assert(exist('itu676_annex2_slant_gas_loss','file')==2, ...
    'itu676_annex2_slant_gas_loss.m not found.');

fprintf('\n');
fprintf('============================================================\n');
fprintf('EXPERIMENT 3 - HAPS AVAILABILITY AND LEO COMPENSATION\n');
fprintf('============================================================\n');

%% ============================================================
% Inputs
% =============================================================

candidateFile = fullfile( ...
    repoRoot,'results','candidate_pool_26.mat');

linkBudgetFile = fullfile( ...
    repoRoot,'results','link_budget.mat');

faFile = fullfile( ...
    repoRoot,'results','FA_results.mat');

beFile = fullfile( ...
    repoRoot,'results','BE_results.mat');

orbitalCatalogFile = fullfile( ...
    repoRoot,'results','orbital_catalog.mat');

assert(isfile(candidateFile), ...
    'Candidate pool not found:\n%s',candidateFile);

assert(isfile(linkBudgetFile), ...
    'Link budget not found:\n%s',linkBudgetFile);

assert(isfile(faFile), ...
    'FA result not found:\n%s',faFile);

assert(isfile(beFile), ...
    'BE result not found:\n%s',beFile);

assert(isfile(orbitalCatalogFile), ...
    ['Orbital catalog not found:\n%s\n' ...
     'The expanded 20-LEO experiment requires orbital_catalog.mat.'], ...
    orbitalCatalogFile);

P = load(candidateFile);
LB = load(linkBudgetFile);
FA0 = load(faFile);
BE0 = load(beFile);
OC = load(orbitalCatalogFile);

%% ============================================================
% Nominal 26-transmitter pool
% =============================================================

requiredPoolVars = [ ...
    "UserPositionECEF", ...
    "PositionECEF", ...
    "Label", ...
    "Architecture"];

for k = 1:numel(requiredPoolVars)
    assert(isfield(P,requiredPoolVars(k)), ...
        'Missing candidate-pool variable: %s', ...
        requiredPoolVars(k));
end

requiredLinkVars = [ ...
    "Label", ...
    "Architecture", ...
    "CN0_dBHz", ...
    "SigmaRho_m", ...
    "beta", ...
    "Tcoh", ...
    "etaPos"];

for k = 1:numel(requiredLinkVars)
    assert(isfield(LB,requiredLinkVars(k)), ...
        'Missing link-budget variable: %s', ...
        requiredLinkVars(k));
end

UserPositionECEF = ...
    double(P.UserPositionECEF(:).');

PositionECEF26 = ...
    double(P.PositionECEF);

Label26 = ...
    string(P.Label(:));

Architecture26 = ...
    upper(string(P.Architecture(:)));

CN0_26_dBHz = ...
    double(LB.CN0_dBHz(:));

sigma26_m = ...
    double(LB.SigmaRho_m(:));

Nnom = size(PositionECEF26,1);

assert(Nnom==26, ...
    'Expected 26 nominal candidates; found %d.',Nnom);

assert(isequal(Label26,string(LB.Label(:))));
assert(isequal(Architecture26,upper(string(LB.Architecture(:)))));

assert(all(isfinite(CN0_26_dBHz)));
assert(all(isfinite(sigma26_m)));
assert(all(sigma26_m>0));

%% User geodetic values

if isfield(LB,'userLat_deg')
    userLat_deg = double(LB.userLat_deg);
elseif isfield(P,'userLat_deg')
    userLat_deg = double(P.userLat_deg);
else
    userLat_deg = -22.8596582;
end

if isfield(LB,'userLon_deg')
    userLon_deg = double(LB.userLon_deg);
elseif isfield(P,'userLon_deg')
    userLon_deg = double(P.userLon_deg);
else
    userLon_deg = -43.2303236;
end

if isfield(LB,'userH_m')
    userH_m = double(LB.userH_m);
elseif isfield(P,'userH_m')
    userH_m = double(P.userH_m);
else
    userH_m = 10;
end

if isfield(LB,'analysisTimeUTC')
    analysisTimeUTC = LB.analysisTimeUTC;
elseif isfield(P,'analysisTimeUTC')
    analysisTimeUTC = P.analysisTimeUTC;
else
    analysisTimeUTC = datetime( ...
        2026,8,1,12,0,0,'TimeZone','UTC');
end

%% ============================================================
% Current ranging parameters and meteorology
% =============================================================

c_mps = 299792458;

beta = double(LB.beta);
Tcoh = double(LB.Tcoh);
etaPos = double(LB.etaPos);

Branging_MHz = beta/1e6;

if isfield(LB,'P2145_meteorology')

    met = LB.P2145_meteorology;

else

    p2145Dir = fullfile( ...
        repoRoot,'data','itu','p2145','annual');

    assert(isfolder(p2145Dir), ...
        'P.2145 annual-data directory not found.');

    assert(exist('itu2145_annual50_at_location','file')==2, ...
        'itu2145_annual50_at_location.m not found.');

    met = itu2145_annual50_at_location( ...
        p2145Dir, ...
        userLat_deg, ...
        userLon_deg, ...
        userH_m);
end

target_m = 0.6;
targetTol = 1e-12;

fprintf('Nominal candidates : %d\n',Nnom);
fprintf('Target             : %.3f m\n',target_m);
fprintf('beta               : %.6f MHz\n',beta/1e6);
fprintf('Tcoh               : %.3f ms\n',1e3*Tcoh);
fprintf('etaPos             : %.6f\n',etaPos);

%% ============================================================
% Nominal geometry matrix
% =============================================================

delta26 = ...
    PositionECEF26-UserPositionECEF;

range26_m = ...
    vecnorm(delta26,2,2);

u26 = ...
    delta26./range26_m;

H26 = ...
    [-u26,ones(Nnom,1)];

assert(rank(H26)==4);

%% ============================================================
% Validate nominal FA/BE baseline
% =============================================================

assert(isfield(FA0,'selected'));
assert(isfield(BE0,'selected'));

selectedFA0 = double(FA0.selected(:).');
selectedBE0 = double(BE0.selected(:).');

[fa0Bound,~,fa0Rank] = ...
    subsetMetrics(selectedFA0,H26,sigma26_m);

[be0Bound,~,be0Rank] = ...
    subsetMetrics(selectedBE0,H26,sigma26_m);

assert(fa0Rank==4 && be0Rank==4);

fprintf('Nominal FA         : K=%d | alpha_CRB=%.12f m\n', ...
    numel(selectedFA0),fa0Bound);

fprintf('Nominal BE         : K=%d | alpha_CRB=%.12f m\n', ...
    numel(selectedBE0),be0Bound);

%% ============================================================
% Output directory
% =============================================================

outDir = fullfile( ...
    repoRoot,'results','experiment3_availability');

if ~isfolder(outDir)
    mkdir(outDir);
end

%% ============================================================
% PART A - HAPS AVAILABILITY / SELECTION LIMIT
%
% HAPS_Max is interpreted exactly as the maximum number of HAPS
% allowed in the selected set.
%
% FA:
%   the best initial K=4 subset must respect the HAPS limit,
%   and every added transmitter must preserve the limit.
%
% BE:
%   if the full nominal set violates the HAPS limit, HAPS are
%   greedily removed first until the constraint is satisfied.
%   Standard BE then continues while the target remains feasible.
% =============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('PART A - HAPS AVAILABILITY\n');
fprintf('============================================================\n');

hapsCases = [4 2 0];
nH = numel(hapsCases);

HAPS_Max = hapsCases(:);

FA_K = nan(nH,1);
FA_Bound_m = nan(nH,1);
FA_PDOP = nan(nH,1);
FA_Feasible = false(nH,1);
FA_Evaluations = nan(nH,1);

FA_HAPS = zeros(nH,1);
FA_LEO = zeros(nH,1);
FA_MEO = zeros(nH,1);
FA_GEO = zeros(nH,1);

BE_K = nan(nH,1);
BE_Bound_m = nan(nH,1);
BE_PDOP = nan(nH,1);
BE_Feasible = false(nH,1);
BE_Evaluations = nan(nH,1);

BE_HAPS = zeros(nH,1);
BE_LEO = zeros(nH,1);
BE_MEO = zeros(nH,1);
BE_GEO = zeros(nH,1);

AdmissibleFullBound_m = nan(nH,1);

FA_SelectedIndices = strings(nH,1);
FA_SelectedLabels = strings(nH,1);

BE_SelectedIndices = strings(nH,1);
BE_SelectedLabels = strings(nH,1);

hapsSelectedRows = table();

FA_selected_haps = cell(nH,1);
BE_selected_haps = cell(nH,1);

for icase = 1:nH

    maxHAPS = hapsCases(icase);

    fprintf('\n');
    fprintf('------------------------------------------------------------\n');
    fprintf('HAPS maximum = %d\n',maxHAPS);
    fprintf('------------------------------------------------------------\n');

    %% Best information set admissible under the HAPS constraint

    [admissibleSet,admissibleBound] = ...
        bestAdmissibleFullSet( ...
            Architecture26, ...
            H26, ...
            sigma26_m, ...
            maxHAPS);

    AdmissibleFullBound_m(icase) = admissibleBound;

    fprintf('Best admissible full-set bound : %.12f m\n', ...
        admissibleBound);

    %% Forward Adding with HAPS cap

    [ ...
        selectedFA, ...
        boundFA, ...
        pdopFA, ...
        evalFA, ...
        feasibleFA] = ...
        forwardAddingWithHAPSCap( ...
            H26, ...
            sigma26_m, ...
            Architecture26, ...
            maxHAPS, ...
            target_m, ...
            targetTol);

    %% Backward Elimination with HAPS cap

    [ ...
        selectedBE, ...
        boundBE, ...
        pdopBE, ...
        evalBE, ...
        feasibleBE] = ...
        backwardEliminationWithHAPSCap( ...
            H26, ...
            sigma26_m, ...
            Architecture26, ...
            maxHAPS, ...
            target_m, ...
            targetTol);

    %% Nominal case must reproduce baseline

    if maxHAPS==4

        assert(isequal(sort(selectedFA),sort(selectedFA0)), ...
            'HAPS=4 FA does not reproduce nominal FA selection.');

        assert(isequal(sort(selectedBE),sort(selectedBE0)), ...
            'HAPS=4 BE does not reproduce nominal BE selection.');

        assert(abs(boundFA-fa0Bound)<1e-10);
        assert(abs(boundBE-be0Bound)<1e-10);
    end

    %% Store

    FA_selected_haps{icase} = selectedFA;
    BE_selected_haps{icase} = selectedBE;

    FA_K(icase) = numel(selectedFA);
    FA_Bound_m(icase) = boundFA;
    FA_PDOP(icase) = pdopFA;
    FA_Feasible(icase) = feasibleFA;
    FA_Evaluations(icase) = evalFA;

    FA_HAPS(icase) = sum(Architecture26(selectedFA)=="HAPS");
    FA_LEO(icase)  = sum(Architecture26(selectedFA)=="LEO");
    FA_MEO(icase)  = sum(Architecture26(selectedFA)=="MEO");
    FA_GEO(icase)  = sum(Architecture26(selectedFA)=="GEO");

    BE_K(icase) = numel(selectedBE);
    BE_Bound_m(icase) = boundBE;
    BE_PDOP(icase) = pdopBE;
    BE_Feasible(icase) = feasibleBE;
    BE_Evaluations(icase) = evalBE;

    BE_HAPS(icase) = sum(Architecture26(selectedBE)=="HAPS");
    BE_LEO(icase)  = sum(Architecture26(selectedBE)=="LEO");
    BE_MEO(icase)  = sum(Architecture26(selectedBE)=="MEO");
    BE_GEO(icase)  = sum(Architecture26(selectedBE)=="GEO");

    FA_SelectedIndices(icase) = ...
        strjoin(string(selectedFA)," ");

    FA_SelectedLabels(icase) = ...
        strjoin(Label26(selectedFA)," ");

    BE_SelectedIndices(icase) = ...
        strjoin(string(selectedBE)," ");

    BE_SelectedLabels(icase) = ...
        strjoin(Label26(selectedBE)," ");

    fprintf('FA : K=%d | alpha_CRB=%.12f m | feasible=%s | eval=%d\n', ...
        FA_K(icase),FA_Bound_m(icase), ...
        string(FA_Feasible(icase)),FA_Evaluations(icase));

    fprintf('     HAPS=%d LEO=%d MEO=%d GEO=%d\n', ...
        FA_HAPS(icase),FA_LEO(icase),FA_MEO(icase),FA_GEO(icase));

    fprintf('BE : K=%d | alpha_CRB=%.12f m | feasible=%s | eval=%d\n', ...
        BE_K(icase),BE_Bound_m(icase), ...
        string(BE_Feasible(icase)),BE_Evaluations(icase));

    fprintf('     HAPS=%d LEO=%d MEO=%d GEO=%d\n', ...
        BE_HAPS(icase),BE_LEO(icase),BE_MEO(icase),BE_GEO(icase));

    %% Long-form selected rows

    hapsSelectedRows = [ ...
        hapsSelectedRows; ...
        selectedRowsForCase( ...
            "HAPS", ...
            maxHAPS, ...
            "FA", ...
            selectedFA, ...
            Label26, ...
            Architecture26); ...
        selectedRowsForCase( ...
            "HAPS", ...
            maxHAPS, ...
            "BE", ...
            selectedBE, ...
            Label26, ...
            Architecture26)]; %#ok<AGROW>
end

%% ============================================================
% HAPS summary table
%
% Current names + legacy aliases required by the tuned plot.
% =============================================================

ThapsCurrent = table( ...
    HAPS_Max, ...
    AdmissibleFullBound_m, ...
    FA_K, ...
    FA_Bound_m, ...
    FA_PDOP, ...
    FA_Feasible, ...
    FA_Evaluations, ...
    FA_HAPS,FA_LEO,FA_MEO,FA_GEO, ...
    BE_K, ...
    BE_Bound_m, ...
    BE_PDOP, ...
    BE_Feasible, ...
    BE_Evaluations, ...
    BE_HAPS,BE_LEO,BE_MEO,BE_GEO, ...
    FA_SelectedIndices, ...
    FA_SelectedLabels, ...
    BE_SelectedIndices, ...
    BE_SelectedLabels);

% Legacy plot-compatible aliases.
FTS_K = FA_K;
FTS_Bound_m = FA_Bound_m;
FTS_PDOP = FA_PDOP;
FTS_Feasible = FA_Feasible;
FTS_Evaluations = FA_Evaluations;

FTS_HAPS = FA_HAPS;
FTS_LEO = FA_LEO;
FTS_MEO = FA_MEO;
FTS_GEO = FA_GEO;

BTS_K = BE_K;
BTS_Bound_m = BE_Bound_m;
BTS_PDOP = BE_PDOP;
BTS_Feasible = BE_Feasible;
BTS_Evaluations = BE_Evaluations;

BTS_HAPS = BE_HAPS;
BTS_LEO = BE_LEO;
BTS_MEO = BE_MEO;
BTS_GEO = BE_GEO;

Thaps = table( ...
    HAPS_Max, ...
    AdmissibleFullBound_m, ...
    FTS_K, ...
    FTS_Bound_m, ...
    FTS_PDOP, ...
    FTS_Feasible, ...
    FTS_Evaluations, ...
    FTS_HAPS,FTS_LEO,FTS_MEO,FTS_GEO, ...
    BTS_K, ...
    BTS_Bound_m, ...
    BTS_PDOP, ...
    BTS_Feasible, ...
    BTS_Evaluations, ...
    BTS_HAPS,BTS_LEO,BTS_MEO,BTS_GEO, ...
    FA_SelectedIndices, ...
    FA_SelectedLabels, ...
    BE_SelectedIndices, ...
    BE_SelectedLabels);

fprintf('\n');
fprintf('============================================================\n');
fprintf('HAPS AVAILABILITY SUMMARY\n');
fprintf('============================================================\n');

disp(Thaps(:,{ ...
    'HAPS_Max', ...
    'AdmissibleFullBound_m', ...
    'FTS_K','FTS_Bound_m','FTS_Evaluations', ...
    'FTS_HAPS','FTS_LEO','FTS_MEO','FTS_GEO', ...
    'BTS_K','BTS_Bound_m','BTS_Evaluations'}));

%% ============================================================
% PART B - NO-HAPS LEO COMPENSATION
% =============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('PART B - LEO COMPENSATION WITHOUT HAPS\n');
fprintf('============================================================\n');

assert(isfield(OC,'TOrbitalCatalog'), ...
    'TOrbitalCatalog not found in orbital_catalog.mat.');

Tcatalog = OC.TOrbitalCatalog;

if ismember('SGP4Error',Tcatalog.Properties.VariableNames)
    assert(all(Tcatalog.SGP4Error==0), ...
        'At least one orbital-catalog entry has an SGP4 error.');
end

%% Exact nested orbital ordering

leoOrder = [ ...
    "65175"
    "66864"
    "65348"
    "59340"
    "67736"
    "59960"
    "65279"
    "48482"
    "59544"
    "57470"
    "67011"
    "66267"
    "65388"
    "64047"
    "64037"
    "56537"
    "60271"
    "58736"
    "67734"
    "67950"];

meoOrder = [ ...
    "62362"
    "40081"
    "44115"
    "43231"
    "43234"
    "40080"
    "40079"];

geoOrder = [ ...
    "41589"
    "41904"
    "42692"
    "43562"
    "43175"
    "43228"
    "38087"];

TLEO20 = reorderArchitecture( ...
    Tcatalog,"LEO",leoOrder,"LEO");

TMEO7 = reorderArchitecture( ...
    Tcatalog,"MEO",meoOrder,"MEO");

TGEO7 = reorderArchitecture( ...
    Tcatalog,"GEO",geoOrder,"GEO");

assert(all(TLEO20.Elevation_deg>=5));
assert(all(TMEO7.Elevation_deg>=5));
assert(all(TGEO7.Elevation_deg>=5));

%% ============================================================
% RF parameters - same assumptions as nominal link-budget builder
% =============================================================

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

%% Portable link budget for all extended orbital candidates

TLEO20 = addPortableLinkBudget( ...
    TLEO20, ...
    userH_m,met, ...
    c_mps,beta,Tcoh,etaPos,Branging_MHz, ...
    LEO_Frequency_GHz,LEO_ReceiverFoM_dBK, ...
    LEO_ElevationControl_deg,LEO_EIRP_Control_dBW_MHz, ...
    MEO_Frequency_GHz,EIRP_MEO_total_dBW,BW_MEO_MHz, ...
    MEO_ReceiverFoM_dBK, ...
    GEO_Frequency_GHz,EIRP_GEO_total_dBW,BW_GEO_MHz, ...
    G_rx_GEO_dBi,Tsys_GEO_K);

TMEO7 = addPortableLinkBudget( ...
    TMEO7, ...
    userH_m,met, ...
    c_mps,beta,Tcoh,etaPos,Branging_MHz, ...
    LEO_Frequency_GHz,LEO_ReceiverFoM_dBK, ...
    LEO_ElevationControl_deg,LEO_EIRP_Control_dBW_MHz, ...
    MEO_Frequency_GHz,EIRP_MEO_total_dBW,BW_MEO_MHz, ...
    MEO_ReceiverFoM_dBK, ...
    GEO_Frequency_GHz,EIRP_GEO_total_dBW,BW_GEO_MHz, ...
    G_rx_GEO_dBi,Tsys_GEO_K);

TGEO7 = addPortableLinkBudget( ...
    TGEO7, ...
    userH_m,met, ...
    c_mps,beta,Tcoh,etaPos,Branging_MHz, ...
    LEO_Frequency_GHz,LEO_ReceiverFoM_dBK, ...
    LEO_ElevationControl_deg,LEO_EIRP_Control_dBW_MHz, ...
    MEO_Frequency_GHz,EIRP_MEO_total_dBW,BW_MEO_MHz, ...
    MEO_ReceiverFoM_dBK, ...
    GEO_Frequency_GHz,EIRP_GEO_total_dBW,BW_GEO_MHz, ...
    G_rx_GEO_dBi,Tsys_GEO_K);

%% ============================================================
% Consistency check:
% common 8 LEO + 7 MEO + 7 GEO must reproduce nominal links
% =============================================================

Tcommon = [ ...
    TLEO20(1:8,:); ...
    TMEO7; ...
    TGEO7];

nominalOrbital = ...
    Architecture26~="HAPS";

nominalNORAD = getNominalNORAD(P,LB);

nominalNORADOrbital = ...
    nominalNORAD(nominalOrbital);

assert(isequal( ...
    string(Tcommon.NORAD_CAT_ID), ...
    nominalNORADOrbital), ...
    'Extended orbital ordering differs from nominal orbital ordering.');

sigmaNomOrbital = ...
    sigma26_m(nominalOrbital);

cn0NomOrbital = ...
    CN0_26_dBHz(nominalOrbital);

sigmaOverlapDiff = ...
    Tcommon.SigmaRho_m-sigmaNomOrbital;

cn0OverlapDiff = ...
    Tcommon.CN0_dBHz-cn0NomOrbital;

fprintf('\nNominal-overlap validation\n');
fprintf('Max |Delta sigma_rho| : %.3e m\n', ...
    max(abs(sigmaOverlapDiff)));
fprintf('Max |Delta C/N0|      : %.3e dB-Hz\n', ...
    max(abs(cn0OverlapDiff)));

assert(max(abs(sigmaOverlapDiff))<1e-10, ...
    'Extended-pool sigma_rho does not reproduce nominal links.');

assert(max(abs(cn0OverlapDiff))<1e-10, ...
    'Extended-pool C/N0 does not reproduce nominal links.');

%% ============================================================
% LEO cases
% =============================================================

leoCases = [8 12 16 20];
nL = numel(leoCases);

LEO_Available = leoCases(:);
CandidateCount_LEO = zeros(nL,1);

FullPoolBound_m = nan(nL,1);
FullPoolPDOP = nan(nL,1);
FullPoolFeasible = false(nL,1);

LEO_FA_K = nan(nL,1);
LEO_FA_Bound_m = nan(nL,1);
LEO_FA_PDOP = nan(nL,1);
LEO_FA_Feasible = false(nL,1);
LEO_FA_Evaluations = nan(nL,1);

LEO_FA_LEO = zeros(nL,1);
LEO_FA_MEO = zeros(nL,1);
LEO_FA_GEO = zeros(nL,1);

LEO_BE_K = nan(nL,1);
LEO_BE_Bound_m = nan(nL,1);
LEO_BE_PDOP = nan(nL,1);
LEO_BE_Feasible = false(nL,1);
LEO_BE_Evaluations = nan(nL,1);

LEO_BE_LEO = zeros(nL,1);
LEO_BE_MEO = zeros(nL,1);
LEO_BE_GEO = zeros(nL,1);

LEO_FA_SelectedIndices = strings(nL,1);
LEO_FA_SelectedLabels = strings(nL,1);

LEO_BE_SelectedIndices = strings(nL,1);
LEO_BE_SelectedLabels = strings(nL,1);

leoSelectedRows = table();

FA_selected_leo = cell(nL,1);
BE_selected_leo = cell(nL,1);

for icase = 1:nL

    nLEO = leoCases(icase);

    fprintf('\n');
    fprintf('------------------------------------------------------------\n');
    fprintf('NO HAPS | %d LEO + 7 MEO + 7 GEO\n',nLEO);
    fprintf('------------------------------------------------------------\n');

    Tcase = [ ...
        TLEO20(1:nLEO,:); ...
        TMEO7; ...
        TGEO7];

    Ncase = height(Tcase);

    CandidateCount_LEO(icase) = Ncase;

    Architecture = ...
        upper(string(Tcase.Architecture));

    Label = ...
        string(Tcase.Label);

    PositionECEF = [ ...
        Tcase.X_ECEF_m, ...
        Tcase.Y_ECEF_m, ...
        Tcase.Z_ECEF_m];

    sigmaRho_m = ...
        double(Tcase.SigmaRho_m(:));

    delta = ...
        PositionECEF-UserPositionECEF;

    ranges = ...
        vecnorm(delta,2,2);

    u = ...
        delta./ranges;

    H = ...
        [-u,ones(Ncase,1)];

    allIdx = 1:Ncase;

    [fullBound,fullPDOP,fullRank] = ...
        subsetMetrics( ...
            allIdx,H,sigmaRho_m);

    FullPoolBound_m(icase) = fullBound;
    FullPoolPDOP(icase) = fullPDOP;
    FullPoolFeasible(icase) = ...
        fullRank==4 && ...
        fullBound<=target_m+targetTol;

    fprintf('Full pool : N=%d | alpha_CRB=%.12f m | feasible=%s\n', ...
        Ncase,fullBound,string(FullPoolFeasible(icase)));

    %% If the full pool is infeasible, no subset can satisfy the target.
    % Return the full pool for both methods so the plot shows the
    % physically infeasible bound rather than an artificial algorithm error.

    if ~FullPoolFeasible(icase)

        selectedFA = allIdx;
        selectedBE = allIdx;

        boundFA = fullBound;
        boundBE = fullBound;

        pdopFA = fullPDOP;
        pdopBE = fullPDOP;

        feasibleFA = false;
        feasibleBE = false;

        evalFA = nchoosek(Ncase,4);
        evalBE = 0;

    else

        [ ...
            selectedFA, ...
            boundFA, ...
            pdopFA, ...
            evalFA, ...
            feasibleFA] = ...
            forwardAddingStandard( ...
                H,sigmaRho_m, ...
                target_m,targetTol);

        [ ...
            selectedBE, ...
            boundBE, ...
            pdopBE, ...
            evalBE, ...
            feasibleBE] = ...
            backwardEliminationStandard( ...
                H,sigmaRho_m, ...
                target_m,targetTol);
    end

    %% Store

    FA_selected_leo{icase} = selectedFA;
    BE_selected_leo{icase} = selectedBE;

    LEO_FA_K(icase) = numel(selectedFA);
    LEO_FA_Bound_m(icase) = boundFA;
    LEO_FA_PDOP(icase) = pdopFA;
    LEO_FA_Feasible(icase) = feasibleFA;
    LEO_FA_Evaluations(icase) = evalFA;

    LEO_FA_LEO(icase) = sum(Architecture(selectedFA)=="LEO");
    LEO_FA_MEO(icase) = sum(Architecture(selectedFA)=="MEO");
    LEO_FA_GEO(icase) = sum(Architecture(selectedFA)=="GEO");

    LEO_BE_K(icase) = numel(selectedBE);
    LEO_BE_Bound_m(icase) = boundBE;
    LEO_BE_PDOP(icase) = pdopBE;
    LEO_BE_Feasible(icase) = feasibleBE;
    LEO_BE_Evaluations(icase) = evalBE;

    LEO_BE_LEO(icase) = sum(Architecture(selectedBE)=="LEO");
    LEO_BE_MEO(icase) = sum(Architecture(selectedBE)=="MEO");
    LEO_BE_GEO(icase) = sum(Architecture(selectedBE)=="GEO");

    LEO_FA_SelectedIndices(icase) = ...
        strjoin(string(selectedFA)," ");

    LEO_FA_SelectedLabels(icase) = ...
        strjoin(Label(selectedFA)," ");

    LEO_BE_SelectedIndices(icase) = ...
        strjoin(string(selectedBE)," ");

    LEO_BE_SelectedLabels(icase) = ...
        strjoin(Label(selectedBE)," ");

    fprintf('FA : K=%d | alpha_CRB=%.12f m | feasible=%s | eval=%d\n', ...
        LEO_FA_K(icase),LEO_FA_Bound_m(icase), ...
        string(LEO_FA_Feasible(icase)),LEO_FA_Evaluations(icase));

    fprintf('BE : K=%d | alpha_CRB=%.12f m | feasible=%s | eval=%d\n', ...
        LEO_BE_K(icase),LEO_BE_Bound_m(icase), ...
        string(LEO_BE_Feasible(icase)),LEO_BE_Evaluations(icase));

    %% Long-form selected rows

    leoSelectedRows = [ ...
        leoSelectedRows; ...
        selectedRowsForCase( ...
            "LEO", ...
            nLEO, ...
            "FA", ...
            selectedFA, ...
            Label, ...
            Architecture); ...
        selectedRowsForCase( ...
            "LEO", ...
            nLEO, ...
            "BE", ...
            selectedBE, ...
            Label, ...
            Architecture)]; %#ok<AGROW>
end

%% ============================================================
% Cross-check no-HAPS 8-LEO case with Part A
% =============================================================

idxH0 = find(HAPS_Max==0,1);

assert(~isempty(idxH0));

deltaNoHAPS_FA = ...
    LEO_FA_Bound_m(1)-FA_Bound_m(idxH0);

deltaNoHAPS_BE = ...
    LEO_BE_Bound_m(1)-BE_Bound_m(idxH0);

fprintf('\nNo-HAPS cross-check\n');
fprintf('FA difference : %.3e m\n',deltaNoHAPS_FA);
fprintf('BE difference : %.3e m\n',deltaNoHAPS_BE);

assert(abs(deltaNoHAPS_FA)<1e-10);
assert(abs(deltaNoHAPS_BE)<1e-10);

%% ============================================================
% LEO summary table
% Current names + legacy aliases for tuned plot.
% =============================================================

LEO_FTS_K = LEO_FA_K;
LEO_FTS_Bound_m = LEO_FA_Bound_m;
LEO_FTS_PDOP = LEO_FA_PDOP;
LEO_FTS_Feasible = LEO_FA_Feasible;
LEO_FTS_Evaluations = LEO_FA_Evaluations;
LEO_FTS_LEO = LEO_FA_LEO;
LEO_FTS_MEO = LEO_FA_MEO;
LEO_FTS_GEO = LEO_FA_GEO;

LEO_BTS_K = LEO_BE_K;
LEO_BTS_Bound_m = LEO_BE_Bound_m;
LEO_BTS_PDOP = LEO_BE_PDOP;
LEO_BTS_Feasible = LEO_BE_Feasible;
LEO_BTS_Evaluations = LEO_BE_Evaluations;
LEO_BTS_LEO = LEO_BE_LEO;
LEO_BTS_MEO = LEO_BE_MEO;
LEO_BTS_GEO = LEO_BE_GEO;

Tleo = table( ...
    LEO_Available, ...
    CandidateCount_LEO, ...
    FullPoolBound_m, ...
    FullPoolPDOP, ...
    FullPoolFeasible, ...
    LEO_FTS_K, ...
    LEO_FTS_Bound_m, ...
    LEO_FTS_PDOP, ...
    LEO_FTS_Feasible, ...
    LEO_FTS_Evaluations, ...
    LEO_FTS_LEO,LEO_FTS_MEO,LEO_FTS_GEO, ...
    LEO_BTS_K, ...
    LEO_BTS_Bound_m, ...
    LEO_BTS_PDOP, ...
    LEO_BTS_Feasible, ...
    LEO_BTS_Evaluations, ...
    LEO_BTS_LEO,LEO_BTS_MEO,LEO_BTS_GEO, ...
    LEO_FA_SelectedIndices, ...
    LEO_FA_SelectedLabels, ...
    LEO_BE_SelectedIndices, ...
    LEO_BE_SelectedLabels);

fprintf('\n');
fprintf('============================================================\n');
fprintf('LEO COMPENSATION SUMMARY\n');
fprintf('============================================================\n');

disp(Tleo(:,{ ...
    'LEO_Available', ...
    'CandidateCount_LEO', ...
    'FullPoolBound_m', ...
    'LEO_FTS_K','LEO_FTS_Bound_m','LEO_FTS_Evaluations', ...
    'LEO_BTS_K','LEO_BTS_Bound_m','LEO_BTS_Evaluations'}));

%% ============================================================
% Save outputs
% =============================================================

hapsCSV = fullfile( ...
    outDir,'experiment3_haps_summary.csv');

leoCSV = fullfile( ...
    outDir,'experiment3_leo_summary.csv');

hapsSelectedCSV = fullfile( ...
    outDir,'experiment3_haps_selected.csv');

leoSelectedCSV = fullfile( ...
    outDir,'experiment3_leo_selected.csv');

matFile = fullfile( ...
    outDir,'experiment3_availability_results.mat');

writetable(Thaps,hapsCSV);
writetable(Tleo,leoCSV);
writetable(hapsSelectedRows,hapsSelectedCSV);
writetable(leoSelectedRows,leoSelectedCSV);

metadata = struct;
metadata.AnalysisTimeUTC = analysisTimeUTC;
metadata.Target_m = target_m;
metadata.AtmosphereModel = ...
    'Portable ITU-R P.2145 annual median + P.676-13 Annex 2';
metadata.Beta_Hz = beta;
metadata.Tcoh_s = Tcoh;
metadata.EtaPos = etaPos;
metadata.HAPSConstraint = ...
    'Maximum HAPS allowed in selected subset';
metadata.LEOCases = leoCases;
metadata.HAPSCases = hapsCases;

save( ...
    matFile, ...
    'Thaps', ...
    'ThapsCurrent', ...
    'Tleo', ...
    'hapsSelectedRows', ...
    'leoSelectedRows', ...
    'FA_selected_haps', ...
    'BE_selected_haps', ...
    'FA_selected_leo', ...
    'BE_selected_leo', ...
    'TLEO20', ...
    'TMEO7', ...
    'TGEO7', ...
    'target_m', ...
    'met', ...
    'metadata', ...
    '-v7.3');

fprintf('\nSaved:\n');
fprintf('  %s\n',hapsCSV);
fprintf('  %s\n',leoCSV);
fprintf('  %s\n',hapsSelectedCSV);
fprintf('  %s\n',leoSelectedCSV);
fprintf('  %s\n',matFile);

fprintf('\n');
fprintf('============================================================\n');
fprintf('EXPERIMENT 3 COMPLETED\n');
fprintf('============================================================\n');

%% ============================================================
% LOCAL FUNCTIONS
% =============================================================

function [alphaCRB_m,PDOP,rankH] = ...
    subsetMetrics(idx,Hfull,sigmaRho_m)

    idx = idx(:);

    if isempty(idx)
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

    if rankH<4
        return;
    end

    A = H.'*H;

    if rcond(A)>1e-14

        Cgeom = A\eye(4);

        trGeom = trace(Cgeom(1:3,1:3));

        if trGeom>0 && isfinite(trGeom)
            PDOP = sqrt(trGeom);
        end
    end

    R = diag(sigma.^2);

    J = H.'*(R\H);

    if rcond(J)<=1e-14
        return;
    end

    C = J\eye(4);

    trPos = trace(C(1:3,1:3));

    if trPos>0 && isfinite(trPos)
        alphaCRB_m = sqrt(trPos);
    end
end

function [bestSet,bestBound] = ...
    bestAdmissibleFullSet( ...
        Architecture,Hfull,sigmaRho_m,maxHAPS)

    hapsIdx = find(Architecture=="HAPS").';
    orbitalIdx = find(Architecture~="HAPS").';

    bestBound = inf;
    bestSet = [];

    nUse = min(maxHAPS,numel(hapsIdx));

    if nUse==0

        trial = orbitalIdx;

        [b,~,rk] = ...
            subsetMetrics(trial,Hfull,sigmaRho_m);

        if rk==4
            bestBound = b;
            bestSet = trial;
        end

        return;
    end

    comb = nchoosek(hapsIdx,nUse);

    for k = 1:size(comb,1)

        trial = [orbitalIdx,comb(k,:)];

        [b,~,rk] = ...
            subsetMetrics(trial,Hfull,sigmaRho_m);

        if rk==4 && b<bestBound
            bestBound = b;
            bestSet = trial;
        end
    end
end

function [selected,finalBound,finalPDOP,evaluationCount,feasible] = ...
    forwardAddingWithHAPSCap( ...
        Hfull,sigmaRho_m,Architecture,maxHAPS,target_m,targetTol)

    N = size(Hfull,1);

    comb4 = nchoosek(1:N,4);

    best4 = [];
    best4Bound = inf;

    evaluationCount = 0;

    for k = 1:size(comb4,1)

        idx = comb4(k,:);

        if sum(Architecture(idx)=="HAPS")>maxHAPS
            continue;
        end

        [b,~,rk] = ...
            subsetMetrics(idx,Hfull,sigmaRho_m);

        evaluationCount = evaluationCount+1;

        if rk==4 && b<best4Bound
            best4Bound = b;
            best4 = idx;
        end
    end

    assert(~isempty(best4), ...
        'No valid K=4 subset under HAPS constraint.');

    selected = best4(:).';
    currentBound = best4Bound;

    while currentBound>target_m+targetTol && ...
          numel(selected)<N

        remaining = ...
            setdiff(1:N,selected,'stable');

        bestCandidate = NaN;
        bestCandidateBound = inf;

        for k = 1:numel(remaining)

            candidate = remaining(k);
            trial = [selected,candidate];

            if sum(Architecture(trial)=="HAPS")>maxHAPS
                continue;
            end

            [b,~,rk] = ...
                subsetMetrics(trial,Hfull,sigmaRho_m);

            evaluationCount = evaluationCount+1;

            if rk==4 && b<bestCandidateBound

                bestCandidate = candidate;
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
        rankH==4 && ...
        finalBound<=target_m+targetTol;
end

function [selected,finalBound,finalPDOP,evaluationCount,feasible] = ...
    backwardEliminationWithHAPSCap( ...
        Hfull,sigmaRho_m,Architecture,maxHAPS,target_m,targetTol)

    N = size(Hfull,1);

    selected = 1:N;
    evaluationCount = 0;

    %% First enforce HAPS cap greedily.
    while sum(Architecture(selected)=="HAPS")>maxHAPS

        removableHAPS = ...
            selected(Architecture(selected)=="HAPS");

        bestRemove = NaN;
        bestBound = inf;

        for k = 1:numel(removableHAPS)

            candidate = removableHAPS(k);

            trial = ...
                selected(selected~=candidate);

            [b,~,rk] = ...
                subsetMetrics(trial,Hfull,sigmaRho_m);

            evaluationCount = evaluationCount+1;

            if rk==4 && b<bestBound
                bestBound = b;
                bestRemove = candidate;
            end
        end

        assert(~isnan(bestRemove), ...
            'Could not enforce HAPS constraint in BE.');

        selected(selected==bestRemove) = [];
    end

    [currentBound,currentPDOP,currentRank] = ...
        subsetMetrics(selected,Hfull,sigmaRho_m);

    %% If constraint itself makes the complete admissible set infeasible,
    % no further removal can restore information.
    if currentRank<4 || currentBound>target_m+targetTol

        finalBound = currentBound;
        finalPDOP = currentPDOP;
        feasible = false;
        return;
    end

    %% Standard BE once HAPS constraint is satisfied.
    while numel(selected)>4

        bestRemove = NaN;
        bestBound = inf;

        for k = 1:numel(selected)

            candidate = selected(k);

            trial = ...
                selected(selected~=candidate);

            if sum(Architecture(trial)=="HAPS")>maxHAPS
                continue;
            end

            [b,~,rk] = ...
                subsetMetrics(trial,Hfull,sigmaRho_m);

            evaluationCount = evaluationCount+1;

            if rk~=4
                continue;
            end

            if b>target_m+targetTol
                continue;
            end

            if b<bestBound

                bestBound = b;
                bestRemove = candidate;

            elseif abs(b-bestBound)<=1e-14

                if isnan(bestRemove) || candidate<bestRemove
                    bestBound = b;
                    bestRemove = candidate;
                end
            end
        end

        if isnan(bestRemove)
            break;
        end

        selected(selected==bestRemove) = [];
    end

    [finalBound,finalPDOP,rankH] = ...
        subsetMetrics(selected,Hfull,sigmaRho_m);

    feasible = ...
        rankH==4 && ...
        finalBound<=target_m+targetTol;
end

function [selected,finalBound,finalPDOP,evaluationCount,feasible] = ...
    forwardAddingStandard( ...
        Hfull,sigmaRho_m,target_m,targetTol)

    N = size(Hfull,1);

    comb4 = nchoosek(1:N,4);

    best4 = [];
    best4Bound = inf;

    evaluationCount = size(comb4,1);

    for k = 1:size(comb4,1)

        idx = comb4(k,:);

        [b,~,rk] = ...
            subsetMetrics(idx,Hfull,sigmaRho_m);

        if rk==4 && b<best4Bound
            best4Bound = b;
            best4 = idx;
        end
    end

    assert(~isempty(best4));

    selected = best4(:).';
    currentBound = best4Bound;

    while currentBound>target_m+targetTol && ...
          numel(selected)<N

        remaining = ...
            setdiff(1:N,selected,'stable');

        bestCandidate = NaN;
        bestCandidateBound = inf;

        for k = 1:numel(remaining)

            candidate = remaining(k);
            trial = [selected,candidate];

            [b,~,rk] = ...
                subsetMetrics(trial,Hfull,sigmaRho_m);

            evaluationCount = evaluationCount+1;

            if rk==4 && b<bestCandidateBound
                bestCandidate = candidate;
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
        rankH==4 && ...
        finalBound<=target_m+targetTol;
end

function [selected,finalBound,finalPDOP,evaluationCount,feasible] = ...
    backwardEliminationStandard( ...
        Hfull,sigmaRho_m,target_m,targetTol)

    N = size(Hfull,1);

    selected = 1:N;
    evaluationCount = 0;

    [fullBound,fullPDOP,fullRank] = ...
        subsetMetrics(selected,Hfull,sigmaRho_m);

    if fullRank<4 || fullBound>target_m+targetTol

        finalBound = fullBound;
        finalPDOP = fullPDOP;
        feasible = false;
        return;
    end

    while numel(selected)>4

        bestRemove = NaN;
        bestBound = inf;

        for k = 1:numel(selected)

            candidate = selected(k);

            trial = ...
                selected(selected~=candidate);

            [b,~,rk] = ...
                subsetMetrics(trial,Hfull,sigmaRho_m);

            evaluationCount = evaluationCount+1;

            if rk~=4
                continue;
            end

            if b>target_m+targetTol
                continue;
            end

            if b<bestBound

                bestBound = b;
                bestRemove = candidate;

            elseif abs(b-bestBound)<=1e-14

                if isnan(bestRemove) || candidate<bestRemove
                    bestBound = b;
                    bestRemove = candidate;
                end
            end
        end

        if isnan(bestRemove)
            break;
        end

        selected(selected==bestRemove) = [];
    end

    [finalBound,finalPDOP,rankH] = ...
        subsetMetrics(selected,Hfull,sigmaRho_m);

    feasible = ...
        rankH==4 && ...
        finalBound<=target_m+targetTol;
end

function Tout = ...
    reorderArchitecture( ...
        T,arch,noradOrder,labelPrefix)

    Tout = T([],:);

    for k = 1:numel(noradOrder)

        idx = ...
            upper(string(T.Architecture))==upper(arch) & ...
            string(T.NORAD_CAT_ID)==noradOrder(k);

        if sum(idx)~=1
            error( ...
                'Could not uniquely locate %s NORAD %s.', ...
                arch,noradOrder(k));
        end

        Tout = [Tout;T(idx,:)]; %#ok<AGROW>
    end

    Tout.Label = ...
        labelPrefix + string((1:height(Tout)).');
end

function T = ...
    addPortableLinkBudget( ...
        T, ...
        userH_m,met, ...
        c_mps,beta,Tcoh,etaPos,Branging_MHz, ...
        LEO_Frequency_GHz,LEO_ReceiverFoM_dBK, ...
        LEO_ElevationControl_deg,LEO_EIRP_Control_dBW_MHz, ...
        MEO_Frequency_GHz,EIRP_MEO_total_dBW,BW_MEO_MHz, ...
        MEO_ReceiverFoM_dBK, ...
        GEO_Frequency_GHz,EIRP_GEO_total_dBW,BW_GEO_MHz, ...
        G_rx_GEO_dBi,Tsys_GEO_K)

    N = height(T);

    Architecture = ...
        upper(string(T.Architecture));

    Frequency_GHz = nan(N,1);
    EIRP_PSD_dBW_MHz = nan(N,1);
    ReceiverFoM_dBK = nan(N,1);

    %% LEO

    idx = Architecture=="LEO";

    Frequency_GHz(idx) = ...
        LEO_Frequency_GHz;

    ReceiverFoM_dBK(idx) = ...
        LEO_ReceiverFoM_dBK;

    elev = ...
        T.Elevation_deg(idx);

    elevClamped = min( ...
        max(elev,min(LEO_ElevationControl_deg)), ...
        max(LEO_ElevationControl_deg));

    EIRP_PSD_dBW_MHz(idx) = ...
        interp1( ...
            LEO_ElevationControl_deg, ...
            LEO_EIRP_Control_dBW_MHz, ...
            elevClamped, ...
            'linear');

    %% MEO

    idx = Architecture=="MEO";

    Frequency_GHz(idx) = ...
        MEO_Frequency_GHz;

    EIRP_PSD_dBW_MHz(idx) = ...
        EIRP_MEO_total_dBW ...
        - 10*log10(BW_MEO_MHz);

    ReceiverFoM_dBK(idx) = ...
        MEO_ReceiverFoM_dBK;

    %% GEO

    idx = Architecture=="GEO";

    Frequency_GHz(idx) = ...
        GEO_Frequency_GHz;

    EIRP_PSD_dBW_MHz(idx) = ...
        EIRP_GEO_total_dBW ...
        - 10*log10(BW_GEO_MHz);

    ReceiverFoM_dBK(idx) = ...
        G_rx_GEO_dBi ...
        - 10*log10(Tsys_GEO_K);

    assert(all(isfinite(Frequency_GHz)));
    assert(all(isfinite(EIRP_PSD_dBW_MHz)));
    assert(all(isfinite(ReceiverFoM_dBK)));

    %% Portable gaseous attenuation

    AtmosphericGasLoss_dB = ...
        nan(N,1);

    stationHeight_km = ...
        userH_m/1e3;

    for ii = 1:N

        AtmosphericGasLoss_dB(ii) = ...
            itu676_annex2_slant_gas_loss( ...
                Frequency_GHz(ii), ...
                T.Elevation_deg(ii), ...
                stationHeight_km, ...
                met);
    end

    assert(all(isfinite(AtmosphericGasLoss_dB)));
    assert(all(AtmosphericGasLoss_dB>=0));

    %% FSPL

    FSPL_dB = ...
        92.45 ...
        + 20*log10(T.SlantRange_km) ...
        + 20*log10(Frequency_GHz);

    %% Positioning EIRP

    EIRP_positioning_dBW = ...
        EIRP_PSD_dBW_MHz ...
        + 10*log10(Branging_MHz) ...
        + 10*log10(etaPos);

    %% C/N0

    CN0_dBHz = ...
        EIRP_positioning_dBW ...
        - FSPL_dB ...
        - AtmosphericGasLoss_dB ...
        + ReceiverFoM_dBK ...
        + 228.6;

    CN0_linear = ...
        10.^(CN0_dBHz/10);

    %% Pseudorange uncertainty

    SigmaRho_m = ...
        c_mps ./ ...
        (2*pi*beta .* ...
        sqrt(CN0_linear*Tcoh));

    %% Store

    T.Frequency_GHz = Frequency_GHz;
    T.EIRP_PSD_dBW_MHz = EIRP_PSD_dBW_MHz;
    T.ReceiverFoM_dBK = ReceiverFoM_dBK;
    T.AtmosphericGasLoss_dB = AtmosphericGasLoss_dB;
    T.FSPL_dB = FSPL_dB;
    T.EIRP_positioning_dBW = EIRP_positioning_dBW;
    T.CN0_dBHz = CN0_dBHz;
    T.SigmaRho_m = SigmaRho_m;
end

function norad = ...
    getNominalNORAD(P,LB)

    if isfield(LB,'NORAD')

        norad = string(LB.NORAD(:));

    elseif isfield(P,'NORAD')

        norad = string(P.NORAD(:));

    elseif isfield(P,'Tpool') && ...
           istable(P.Tpool) && ...
           ismember('NORAD_CAT_ID',P.Tpool.Properties.VariableNames)

        norad = string(P.Tpool.NORAD_CAT_ID);

    else
        error('Could not locate nominal NORAD identifiers.');
    end
end

function T = ...
    selectedRowsForCase( ...
        experimentPart, ...
        caseValue, ...
        method, ...
        selected, ...
        Label, ...
        Architecture)

    n = numel(selected);

    T = table( ...
        repmat(string(experimentPart),n,1), ...
        repmat(caseValue,n,1), ...
        repmat(string(method),n,1), ...
        (1:n).', ...
        selected(:), ...
        Label(selected(:)), ...
        Architecture(selected(:)), ...
        'VariableNames',{ ...
        'ExperimentPart', ...
        'CaseValue', ...
        'Method', ...
        'SelectionOrder', ...
        'PoolIndex', ...
        'Label', ...
        'Architecture'});
end
