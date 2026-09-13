clear;
clc;
format long;

%% ============================================================
% ETA_POS SENSITIVITY - NOMINAL TRANSMITTER SELECTION
%
% Purpose:
%   Quantify how the common positioning-power allocation etaPos
%   affects the positional bound and the minimum subset size
%   returned by Forward Adding (FA) and Backward Elimination (BE).
%
% Important:
%   etaPos is a COMMON multiplicative power factor for every
%   candidate. Therefore, relative measurement-quality ordering is
%   preserved and all nominal pseudorange standard deviations scale as
%
%       sigma_rho(eta) = sigma_rho(eta0)*sqrt(eta0/eta).
%
%   This script does NOT modify results/link_budget.mat.
%
% Default sensitivity grid:
%   0.25%, 0.5%, 1%, 2%, 4%
%   corresponding to approximately -6, -3, 0, +3, +6 dB
%   relative to the nominal etaPos = 1%.
%
% Suggested location:
%   matlab/experiments/run_etaPos_sensitivity.m
% =============================================================

%% Repository paths

thisFile = mfilename('fullpath');

if isempty(thisFile)
    error('This script must be executed from a saved .m file.');
end

experimentsDir = fileparts(thisFile);
matlabDir = fileparts(experimentsDir);
repoRoot = fileparts(matlabDir);

addpath(repoRoot,'-begin');
setup_paths;

fprintf('\n');
fprintf('============================================================\n');
fprintf('ETA_POS SENSITIVITY\n');
fprintf('============================================================\n');

%% Input files

candidateFile = fullfile( ...
    repoRoot,'results','candidate_pool_26.mat');

linkBudgetFile = fullfile( ...
    repoRoot,'results','link_budget.mat');

assert(isfile(candidateFile), ...
    'Candidate pool not found:\n%s',candidateFile);

assert(isfile(linkBudgetFile), ...
    'Link budget not found:\n%s',linkBudgetFile);

S  = load(candidateFile);
LB = load(linkBudgetFile);

%% Required variables

requiredCandidateVars = [ ...
    "UserPositionECEF", ...
    "PositionECEF", ...
    "Label", ...
    "Architecture"];

for k = 1:numel(requiredCandidateVars)
    assert(isfield(S,requiredCandidateVars(k)), ...
        'Missing candidate-pool variable: %s', ...
        requiredCandidateVars(k));
end

requiredLinkVars = [ ...
    "SigmaRho_m", ...
    "etaPos"];

for k = 1:numel(requiredLinkVars)
    assert(isfield(LB,requiredLinkVars(k)), ...
        'Missing link-budget variable: %s', ...
        requiredLinkVars(k));
end

UserPositionECEF = S.UserPositionECEF(:).';
PositionECEF = S.PositionECEF;

Label = string(S.Label(:));
Architecture = string(S.Architecture(:));

sigmaNominal_m = LB.SigmaRho_m(:);
etaNominal = double(LB.etaPos);

Npool = size(PositionECEF,1);

assert(Npool == numel(sigmaNominal_m));
assert(Npool == numel(Label));
assert(Npool == numel(Architecture));

assert(etaNominal > 0 && etaNominal <= 1);
assert(all(isfinite(sigmaNominal_m)));
assert(all(sigmaNominal_m > 0));

%% Geometry matrix

r = PositionECEF - UserPositionECEF;
d = sqrt(sum(r.^2,2));

assert(all(isfinite(d)) && all(d > eps), ...
    'Invalid transmitter-user geometry.');

u = r ./ d;

% Same convention used by run_FA.m and run_BE.m
Hfull = [-u,ones(Npool,1)];

assert(rank(Hfull) == 4, ...
    'Full candidate geometry does not have rank(H)=4.');

%% Selection target

target_m = 0.6;
targetTol = 1e-12;

%% Sensitivity values
%
% Multipliers relative to nominal etaPos:
% 0.25, 0.5, 1, 2, 4.
%
% If etaNominal = 0.01, this gives:
% 0.0025, 0.005, 0.01, 0.02, 0.04.

etaMultiplier = [0.25;0.5;1;2;4];
etaGrid = etaNominal*etaMultiplier;

if any(etaGrid > 1)
    error('Sensitivity grid produced etaPos > 1.');
end

%% Baseline full-pool metric

[fullBoundNominal,fullPDOP,fullRankH] = ...
    subsetMetrics( ...
        1:Npool, ...
        Hfull, ...
        sigmaNominal_m);

assert(fullRankH == 4);

etaFullPoolThreshold = ...
    etaNominal*(fullBoundNominal/target_m)^2;

fprintf('Nominal etaPos          : %.8f (%.4f %%)\n', ...
    etaNominal,100*etaNominal);

fprintf('Nominal full-pool bound : %.12f m\n', ...
    fullBoundNominal);

fprintf('Target                  : %.6f m\n', ...
    target_m);

fprintf('etaPos required for full-pool target (scaling law): %.8f (%.4f %%)\n', ...
    etaFullPoolThreshold,100*etaFullPoolThreshold);

%% Allocate outputs

nCases = numel(etaGrid);

Allocation_pct = 100*etaGrid;
OffsetRelativeNominal_dB = 10*log10(etaGrid/etaNominal);
SigmaScale = sqrt(etaNominal./etaGrid);

FullBound_m = nan(nCases,1);
PredictedFullBound_m = nan(nCases,1);
ScalingResidual_m = nan(nCases,1);

FA_Feasible = false(nCases,1);
FA_K = nan(nCases,1);
FA_Bound_m = nan(nCases,1);
FA_PDOP = nan(nCases,1);
FA_HAPS = nan(nCases,1);
FA_LEO = nan(nCases,1);
FA_MEO = nan(nCases,1);
FA_GEO = nan(nCases,1);
FA_Evaluations = nan(nCases,1);
FA_Selected = strings(nCases,1);

BE_Feasible = false(nCases,1);
BE_K = nan(nCases,1);
BE_Bound_m = nan(nCases,1);
BE_PDOP = nan(nCases,1);
BE_HAPS = nan(nCases,1);
BE_LEO = nan(nCases,1);
BE_MEO = nan(nCases,1);
BE_GEO = nan(nCases,1);
BE_Evaluations = nan(nCases,1);
BE_Selected = strings(nCases,1);

%% Evaluate cases

for caseIdx = 1:nCases

    etaTest = etaGrid(caseIdx);

    sigmaTest_m = ...
        sigmaNominal_m*sqrt(etaNominal/etaTest);

    [fullBound,~,rankH] = ...
        subsetMetrics( ...
            1:Npool, ...
            Hfull, ...
            sigmaTest_m);

    assert(rankH == 4);

    predictedFullBound = ...
        fullBoundNominal*sqrt(etaNominal/etaTest);

    FullBound_m(caseIdx) = fullBound;
    PredictedFullBound_m(caseIdx) = predictedFullBound;
    ScalingResidual_m(caseIdx) = ...
        fullBound-predictedFullBound;

    %% Forward Adding

    fa = runFA( ...
        Hfull, ...
        sigmaTest_m, ...
        target_m, ...
        Label, ...
        Architecture);

    FA_Feasible(caseIdx) = fa.feasible;
    FA_K(caseIdx) = numel(fa.selected);
    FA_Bound_m(caseIdx) = fa.bound_m;
    FA_PDOP(caseIdx) = fa.PDOP;
    FA_Evaluations(caseIdx) = fa.evaluations;

    FA_HAPS(caseIdx) = sum(Architecture(fa.selected)=="HAPS");
    FA_LEO(caseIdx)  = sum(Architecture(fa.selected)=="LEO");
    FA_MEO(caseIdx)  = sum(Architecture(fa.selected)=="MEO");
    FA_GEO(caseIdx)  = sum(Architecture(fa.selected)=="GEO");

    FA_Selected(caseIdx) = ...
        strjoin(Label(fa.selected)," ");

    %% Backward Elimination

    be = runBE( ...
        Hfull, ...
        sigmaTest_m, ...
        target_m, ...
        targetTol, ...
        Label, ...
        Architecture);

    BE_Feasible(caseIdx) = be.feasible;
    BE_K(caseIdx) = numel(be.selected);
    BE_Bound_m(caseIdx) = be.bound_m;
    BE_PDOP(caseIdx) = be.PDOP;
    BE_Evaluations(caseIdx) = be.evaluations;

    BE_HAPS(caseIdx) = sum(Architecture(be.selected)=="HAPS");
    BE_LEO(caseIdx)  = sum(Architecture(be.selected)=="LEO");
    BE_MEO(caseIdx)  = sum(Architecture(be.selected)=="MEO");
    BE_GEO(caseIdx)  = sum(Architecture(be.selected)=="GEO");

    BE_Selected(caseIdx) = ...
        strjoin(Label(be.selected)," ");

end

%% Validate common-scaling law

maxScalingResidual = max(abs(ScalingResidual_m));

fprintf('\nMaximum full-pool scaling residual: %.3e m\n', ...
    maxScalingResidual);

assert(maxScalingResidual < 1e-9, ...
    'Common etaPos scaling law validation failed.');

%% Summary table

Teta = table( ...
    etaGrid, ...
    Allocation_pct, ...
    OffsetRelativeNominal_dB, ...
    SigmaScale, ...
    FullBound_m, ...
    PredictedFullBound_m, ...
    ScalingResidual_m, ...
    FA_Feasible, ...
    FA_K, ...
    FA_Bound_m, ...
    FA_PDOP, ...
    FA_HAPS, ...
    FA_LEO, ...
    FA_MEO, ...
    FA_GEO, ...
    FA_Evaluations, ...
    FA_Selected, ...
    BE_Feasible, ...
    BE_K, ...
    BE_Bound_m, ...
    BE_PDOP, ...
    BE_HAPS, ...
    BE_LEO, ...
    BE_MEO, ...
    BE_GEO, ...
    BE_Evaluations, ...
    BE_Selected, ...
    'VariableNames',{ ...
    'EtaPos', ...
    'Allocation_pct', ...
    'OffsetRelativeNominal_dB', ...
    'SigmaScale', ...
    'FullBound_m', ...
    'PredictedFullBound_m', ...
    'ScalingResidual_m', ...
    'FA_Feasible', ...
    'FA_K', ...
    'FA_Bound_m', ...
    'FA_PDOP', ...
    'FA_HAPS', ...
    'FA_LEO', ...
    'FA_MEO', ...
    'FA_GEO', ...
    'FA_Evaluations', ...
    'FA_Selected', ...
    'BE_Feasible', ...
    'BE_K', ...
    'BE_Bound_m', ...
    'BE_PDOP', ...
    'BE_HAPS', ...
    'BE_LEO', ...
    'BE_MEO', ...
    'BE_GEO', ...
    'BE_Evaluations', ...
    'BE_Selected'});

fprintf('\n');
fprintf('============================================================\n');
fprintf('ETA_POS SENSITIVITY SUMMARY\n');
fprintf('============================================================\n');

disp(Teta(:,{ ...
    'EtaPos', ...
    'Allocation_pct', ...
    'OffsetRelativeNominal_dB', ...
    'FullBound_m', ...
    'FA_Feasible', ...
    'FA_K', ...
    'FA_Bound_m', ...
    'BE_Feasible', ...
    'BE_K', ...
    'BE_Bound_m'}));

%% Save

resultsDir = fullfile( ...
    repoRoot, ...
    'results', ...
    'etaPos_sensitivity');

if ~isfolder(resultsDir)
    mkdir(resultsDir);
end

csvFile = fullfile( ...
    resultsDir, ...
    'etaPos_sensitivity.csv');

matFile = fullfile( ...
    resultsDir, ...
    'etaPos_sensitivity.mat');

writetable(Teta,csvFile);

save(matFile, ...
    'Teta', ...
    'etaNominal', ...
    'etaGrid', ...
    'etaMultiplier', ...
    'etaFullPoolThreshold', ...
    'target_m', ...
    'fullBoundNominal', ...
    'fullPDOP', ...
    'Hfull', ...
    'sigmaNominal_m', ...
    'Label', ...
    'Architecture');

fprintf('\nSaved:\n%s\n%s\n', ...
    csvFile,matFile);

fprintf('\n');
fprintf('============================================================\n');
fprintf('ETA_POS SENSITIVITY COMPLETED\n');
fprintf('============================================================\n');

%% ============================================================
% Local function: Forward Adding
% =============================================================

function out = runFA( ...
    Hfull, ...
    sigmaRho_m, ...
    target_m, ...
    Label, ...
    Architecture)

    Npool = size(Hfull,1);

    comb4 = nchoosek(1:Npool,4);

    best4Bound = inf;
    best4PDOP = inf;
    best4 = [];

    for k = 1:size(comb4,1)

        idx = comb4(k,:);

        [b,pdop,rankH] = ...
            subsetMetrics( ...
                idx, ...
                Hfull, ...
                sigmaRho_m);

        if rankH == 4 && b < best4Bound
            best4Bound = b;
            best4PDOP = pdop;
            best4 = idx;
        end

    end

    if isempty(best4)
        error('FA: no full-rank K=4 subset found.');
    end

    selected = best4(:).';
    currentBound = best4Bound;
    currentPDOP = best4PDOP;

    evaluations = size(comb4,1);

    while currentBound > target_m && ...
          numel(selected) < Npool

        remaining = setdiff( ...
            1:Npool, ...
            selected, ...
            'stable');

        bestCandidate = NaN;
        bestCandidateBound = inf;
        bestCandidatePDOP = inf;

        for k = 1:numel(remaining)

            candidate = remaining(k);
            idxTrial = [selected,candidate];

            [b,pdop,rankH] = ...
                subsetMetrics( ...
                    idxTrial, ...
                    Hfull, ...
                    sigmaRho_m);

            evaluations = evaluations + 1;

            if rankH == 4 && ...
               b < bestCandidateBound

                bestCandidateBound = b;
                bestCandidatePDOP = pdop;
                bestCandidate = candidate;

            end
        end

        if isnan(bestCandidate)
            break;
        end

        selected(end+1) = bestCandidate; %#ok<AGROW>
        currentBound = bestCandidateBound;
        currentPDOP = bestCandidatePDOP;

    end

    [finalBound,finalPDOP,finalRankH] = ...
        subsetMetrics( ...
            selected, ...
            Hfull, ...
            sigmaRho_m);

    out.selected = selected;
    out.bound_m = finalBound;
    out.PDOP = finalPDOP;
    out.rankH = finalRankH;
    out.feasible = ...
        finalRankH == 4 && ...
        finalBound <= target_m;
    out.evaluations = evaluations;

    %#ok<NASGU>
    Label = Label;
    Architecture = Architecture;

end

%% ============================================================
% Local function: Backward Elimination
% =============================================================

function out = runBE( ...
    Hfull, ...
    sigmaRho_m, ...
    target_m, ...
    targetTol, ...
    Label, ...
    Architecture)

    Npool = size(Hfull,1);

    selected = 1:Npool;

    [fullBound,fullPDOP,fullRankH] = ...
        subsetMetrics( ...
            selected, ...
            Hfull, ...
            sigmaRho_m);

    if fullRankH < 4
        error('BE: full candidate set does not have rank(H)=4.');
    end

    evaluations = 0;

    % If even the complete pool misses the target, BE is infeasible.
    if fullBound > target_m + targetTol

        out.selected = selected;
        out.bound_m = fullBound;
        out.PDOP = fullPDOP;
        out.rankH = fullRankH;
        out.feasible = false;
        out.evaluations = evaluations;
        return;

    end

    currentBound = fullBound;
    currentPDOP = fullPDOP;

    while numel(selected) > 4

        removable = selected;

        bestRemoved = NaN;
        bestTrialBound = inf;
        bestTrialPDOP = inf;

        for k = 1:numel(removable)

            candidateToRemove = removable(k);

            idxTrial = ...
                selected(selected ~= candidateToRemove);

            [trialBound,trialPDOP,trialRankH] = ...
                subsetMetrics( ...
                    idxTrial, ...
                    Hfull, ...
                    sigmaRho_m);

            evaluations = evaluations + 1;

            if trialRankH ~= 4
                continue;
            end

            if trialBound > target_m + targetTol
                continue;
            end

            if trialBound < bestTrialBound

                bestTrialBound = trialBound;
                bestTrialPDOP = trialPDOP;
                bestRemoved = candidateToRemove;

            elseif abs(trialBound-bestTrialBound) <= 1e-14

                if isnan(bestRemoved) || ...
                   candidateToRemove < bestRemoved

                    bestTrialBound = trialBound;
                    bestTrialPDOP = trialPDOP;
                    bestRemoved = candidateToRemove;

                end
            end
        end

        if isnan(bestRemoved)
            break;
        end

        selected(selected == bestRemoved) = [];

        currentBound = bestTrialBound;
        currentPDOP = bestTrialPDOP;

    end

    [finalBound,finalPDOP,finalRankH] = ...
        subsetMetrics( ...
            selected, ...
            Hfull, ...
            sigmaRho_m);

    out.selected = selected;
    out.bound_m = finalBound;
    out.PDOP = finalPDOP;
    out.rankH = finalRankH;
    out.feasible = ...
        finalRankH == 4 && ...
        finalBound <= target_m + targetTol;
    out.evaluations = evaluations;

    %#ok<NASGU>
    currentBound = currentBound;
    currentPDOP = currentPDOP;
    Label = Label;
    Architecture = Architecture;

end

%% ============================================================
% Local function: subset metrics
% =============================================================

function [bound_m,PDOP,rankH] = ...
    subsetMetrics( ...
        idx, ...
        Hfull, ...
        sigmaRho_m)

    idx = idx(:);

    H = Hfull(idx,:);
    sigma = sigmaRho_m(idx);

    rankH = rank(H);

    bound_m = inf;
    PDOP = inf;

    if rankH < 4
        return;
    end

    %% Geometric PDOP

    A = H.'*H;

    if rcond(A) > 1e-14

        Cgeom = A\eye(4);

        trGeom = trace(Cgeom(1:3,1:3));

        if trGeom > 0 && isfinite(trGeom)
            PDOP = sqrt(trGeom);
        end

    end

    %% Weighted positional bound

    R = diag(sigma.^2);

    J = H.'*(R\H);

    if rcond(J) <= 1e-14
        return;
    end

    C = J\eye(4);

    trPos = trace(C(1:3,1:3));

    if trPos > 0 && isfinite(trPos)
        bound_m = sqrt(trPos);
    end

end
