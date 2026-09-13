clear;
clc;
format long;

%% ============================================================
% EXACT SUBSET VALIDATION
%
% Purpose:
%   Validate the nominal greedy result by exhaustive enumeration.
%
% The script evaluates:
%   1) every subset of size K_greedy - 1;
%   2) every subset of size K_greedy.
%
% If no subset of size K_greedy - 1 satisfies the positional-bound
% target and at least one subset of size K_greedy does, then the
% greedy solution is globally optimal in NUMBER OF TRANSMITTERS.
%
% The K_greedy enumeration also gives:
%   - exact number/fraction of feasible subsets;
%   - globally best bound at that K;
%   - number of subsets strictly better than the FA/BE solution;
%   - top subsets by positional bound.
%
% Implementation:
%   - combinations are generated in blocks (no giant nchoosek matrix);
%   - each transmitter contributes one 4x4 Fisher-information matrix;
%   - page-wise 4x4 solves use pagemldivide (MATLAB R2022a+).
%
% Suggested location:
%   matlab/validate/run_exact_subset_validation.m
%
% Inputs:
%   results/candidate_pool_26.mat
%   results/link_budget.mat
%   results/FA_results.mat
%   results/BE_results.mat
%
% Outputs:
%   results/exact_subset_validation/
%       exact_subset_summary.csv
%       exact_top_subsets_K*.csv
%       exact_subset_validation.mat
% =============================================================

%% Repository paths

thisFile = mfilename('fullpath');

if isempty(thisFile)
    error('This script must be executed from a saved .m file.');
end

validateDir = fileparts(thisFile);
matlabDir   = fileparts(validateDir);
repoRoot    = fileparts(matlabDir);

addpath(repoRoot,'-begin');
setup_paths;

fprintf('\n');
fprintf('============================================================\n');
fprintf('EXACT SUBSET VALIDATION\n');
fprintf('============================================================\n');

%% Inputs

candidateFile = fullfile( ...
    repoRoot,'results','candidate_pool_26.mat');

linkBudgetFile = fullfile( ...
    repoRoot,'results','link_budget.mat');

faFile = fullfile( ...
    repoRoot,'results','FA_results.mat');

beFile = fullfile( ...
    repoRoot,'results','BE_results.mat');

assert(isfile(candidateFile), ...
    'Candidate pool not found:\n%s',candidateFile);

assert(isfile(linkBudgetFile), ...
    'Link budget not found:\n%s',linkBudgetFile);

assert(isfile(faFile), ...
    'FA result not found:\n%s',faFile);

assert(isfile(beFile), ...
    'BE result not found:\n%s',beFile);

S  = load(candidateFile);
LB = load(linkBudgetFile);
FA = load(faFile);
BE = load(beFile);

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
    "SigmaRho_m"];

for k = 1:numel(requiredLinkVars)
    assert(isfield(LB,requiredLinkVars(k)), ...
        'Missing link-budget variable: %s', ...
        requiredLinkVars(k));
end

assert(isfield(FA,'selected'), ...
    'FA_results.mat does not contain selected.');

assert(isfield(FA,'finalBound'), ...
    'FA_results.mat does not contain finalBound.');

assert(isfield(FA,'target_m'), ...
    'FA_results.mat does not contain target_m.');

assert(isfield(BE,'selected'), ...
    'BE_results.mat does not contain selected.');

assert(isfield(BE,'finalBound'), ...
    'BE_results.mat does not contain finalBound.');

%% Nominal data

UserPositionECEF = double(S.UserPositionECEF(:).');
PositionECEF = double(S.PositionECEF);

Label = string(S.Label(:));
Architecture = string(S.Architecture(:));

sigmaRho_m = double(LB.SigmaRho_m(:));

Npool = size(PositionECEF,1);

assert(numel(Label) == Npool);
assert(numel(Architecture) == Npool);
assert(numel(sigmaRho_m) == Npool);
assert(all(isfinite(sigmaRho_m)));
assert(all(sigmaRho_m > 0));

%% Geometry matrix

delta = PositionECEF - UserPositionECEF;
range_m = sqrt(sum(delta.^2,2));

assert(all(isfinite(range_m)));
assert(all(range_m > eps));

u = delta ./ range_m;

% Same convention used by run_FA.m and run_BE.m.
Hfull = [-u,ones(Npool,1)];

assert(rank(Hfull) == 4, ...
    'Full candidate geometry does not have rank(H)=4.');

%% Greedy result consistency

selectedFA = sort(double(FA.selected(:).'));
selectedBE = sort(double(BE.selected(:).'));

assert(numel(selectedFA) == numel(selectedBE), ...
    'FA and BE returned different subset sizes.');

assert(isequal(selectedFA,selectedBE), ...
    ['FA and BE returned different subsets. ' ...
     'Exact validation can still be done, but this script expects ' ...
     'the nominal FA/BE coincidence.']);

Kgreedy = numel(selectedFA);
Klower  = Kgreedy - 1;

assert(Klower >= 4, ...
    'Kgreedy-1 must be at least 4.');

target_m = double(FA.target_m);
targetTol = 1e-12;
comparisonTol = 1e-12;

faBound = double(FA.finalBound);
beBound = double(BE.finalBound);

assert(abs(faBound-beBound) < 1e-10, ...
    'FA and BE bounds differ unexpectedly.');

%% Direct validation of the FA/BE selected set

[greedyBoundDirect,greedyPDOP,greedyRank] = ...
    subsetMetricsDirect( ...
        selectedFA, ...
        Hfull, ...
        sigmaRho_m);

assert(greedyRank == 4);

assert(abs(greedyBoundDirect-faBound) < 1e-10, ...
    ['Direct recomputation of the selected-set bound does not ' ...
     'match FA_results.mat.']);

fprintf('Candidate pool          : %d\n',Npool);
fprintf('Greedy subset size      : %d\n',Kgreedy);
fprintf('Lower size to test      : %d\n',Klower);
fprintf('Target                  : %.12f m\n',target_m);
fprintf('FA/BE selected bound    : %.12f m\n',greedyBoundDirect);
fprintf('FA/BE selected PDOP     : %.12f\n',greedyPDOP);
fprintf('FA/BE selected indices  : ');
fprintf('%d ',selectedFA);
fprintf('\n');

fprintf('FA/BE selected labels   : ');
fprintf('%s ',Label(selectedFA));
fprintf('\n');

%% Fisher-information contributions
%
% For row h_i and variance sigma_i^2,
%
%   J_i = h_i^T h_i / sigma_i^2.
%
% Any subset FIM is the sum of its transmitter contributions.

JcontribFlat = zeros(Npool,16);

for i = 1:Npool

    h = Hfull(i,:);
    Ji = (h.'*h)/(sigmaRho_m(i)^2);

    JcontribFlat(i,:) = reshape(Ji,1,16);

end

%% Page-wise solver availability

assert(~isempty(which('pagemldivide')), ...
    ['pagemldivide is required for the block-vectorized exact ' ...
     'enumeration (MATLAB R2022a or newer).']);

%% Enumeration settings

batchSize = 100000;
topN = 20;

fprintf('\nEnumeration batch size  : %d\n',batchSize);
fprintf('Top subsets retained    : %d per K\n',topN);

%% Exact K = Kgreedy - 1

fprintf('\n');
fprintf('============================================================\n');
fprintf('EXACT ENUMERATION: K=%d\n',Klower);
fprintf('============================================================\n');

lowerResult = enumerateExactK( ...
    Npool, ...
    Klower, ...
    JcontribFlat, ...
    target_m, ...
    targetTol, ...
    NaN, ...
    comparisonTol, ...
    batchSize, ...
    topN);

printExactResult(lowerResult,Label);

%% Exact K = Kgreedy

fprintf('\n');
fprintf('============================================================\n');
fprintf('EXACT ENUMERATION: K=%d\n',Kgreedy);
fprintf('============================================================\n');

greedyKResult = enumerateExactK( ...
    Npool, ...
    Kgreedy, ...
    JcontribFlat, ...
    target_m, ...
    targetTol, ...
    greedyBoundDirect, ...
    comparisonTol, ...
    batchSize, ...
    topN);

printExactResult(greedyKResult,Label);

%% Global-cardinality conclusion

globalMinimumProved = ...
    lowerResult.NFeasible == 0 && ...
    greedyKResult.NFeasible > 0;

if globalMinimumProved

    globalMinimumK = Kgreedy;

    fprintf('\n');
    fprintf('============================================================\n');
    fprintf('GLOBAL CARDINALITY RESULT\n');
    fprintf('============================================================\n');
    fprintf('No K=%d subset satisfies the %.3f m target.\n', ...
        Klower,target_m);
    fprintf('At least one K=%d subset satisfies the target.\n', ...
        Kgreedy);
    fprintf('Therefore the global minimum number of transmitters is K=%d.\n', ...
        globalMinimumK);
    fprintf(['FA and BE are globally optimal in subset size for this ' ...
             'nominal instance.\n']);

else

    globalMinimumK = NaN;

    fprintf('\n');
    fprintf('============================================================\n');
    fprintf('GLOBAL CARDINALITY RESULT\n');
    fprintf('============================================================\n');

    if lowerResult.NFeasible > 0

        fprintf(['At least one K=%d subset is feasible. Therefore the ' ...
                 'greedy K=%d result is NOT proven globally minimal.\n'], ...
                 Klower,Kgreedy);

    elseif greedyKResult.NFeasible == 0

        fprintf(['No K=%d subset was found feasible, which is inconsistent ' ...
                 'with the nominal greedy result and requires investigation.\n'], ...
                 Kgreedy);

    end
end

%% Greedy quality within the fixed K

bestKBound = greedyKResult.BestBound_m;

greedyIsBestBoundAtK = ...
    abs(greedyBoundDirect-bestKBound) <= comparisonTol;

fprintf('\n');
fprintf('============================================================\n');
fprintf('FIXED-K QUALITY OF FA/BE\n');
fprintf('============================================================\n');

fprintf('Globally best K=%d bound       : %.12f m\n', ...
    Kgreedy,bestKBound);

fprintf('FA/BE K=%d bound               : %.12f m\n', ...
    Kgreedy,greedyBoundDirect);

fprintf('Subsets strictly better than FA/BE : %d\n', ...
    greedyKResult.NStrictlyBetterThanReference);

fprintf('Subsets <= FA/BE bound (tol.)      : %d\n', ...
    greedyKResult.NAtOrBelowReference);

fprintf('FA/BE globally best bound at K?    : %s\n', ...
    string(greedyIsBestBoundAtK));

fprintf('Exact feasible fraction at K=%d    : %.12g %%\n', ...
    Kgreedy,100*greedyKResult.FeasibleFraction);

%% Summary table

K = [Klower;Kgreedy];

TotalSubsets = [ ...
    lowerResult.TotalSubsets; ...
    greedyKResult.TotalSubsets];

NumericallyValid = [ ...
    lowerResult.NNumericallyValid; ...
    greedyKResult.NNumericallyValid];

FeasibleSubsets = [ ...
    lowerResult.NFeasible; ...
    greedyKResult.NFeasible];

FeasibleFraction = [ ...
    lowerResult.FeasibleFraction; ...
    greedyKResult.FeasibleFraction];

FeasiblePercent = 100*FeasibleFraction;

BestBound_m = [ ...
    lowerResult.BestBound_m; ...
    greedyKResult.BestBound_m];

BestIndices = [ ...
    combinationToString(lowerResult.BestCombination); ...
    combinationToString(greedyKResult.BestCombination)];

BestLabels = [ ...
    labelsToString(lowerResult.BestCombination,Label); ...
    labelsToString(greedyKResult.BestCombination,Label)];

ReferenceBound_m = [ ...
    NaN; ...
    greedyBoundDirect];

StrictlyBetterThanReference = [ ...
    NaN; ...
    greedyKResult.NStrictlyBetterThanReference];

AtOrBelowReference = [ ...
    NaN; ...
    greedyKResult.NAtOrBelowReference];

Tsummary = table( ...
    K, ...
    TotalSubsets, ...
    NumericallyValid, ...
    FeasibleSubsets, ...
    FeasibleFraction, ...
    FeasiblePercent, ...
    BestBound_m, ...
    BestIndices, ...
    BestLabels, ...
    ReferenceBound_m, ...
    StrictlyBetterThanReference, ...
    AtOrBelowReference);

fprintf('\n');
disp(Tsummary);

%% Top-subset tables

TtopLower = makeTopTable( ...
    lowerResult, ...
    Label, ...
    Architecture);

TtopGreedyK = makeTopTable( ...
    greedyKResult, ...
    Label, ...
    Architecture);

%% Save

resultsDir = fullfile( ...
    repoRoot, ...
    'results', ...
    'exact_subset_validation');

if ~isfolder(resultsDir)
    mkdir(resultsDir);
end

summaryCSV = fullfile( ...
    resultsDir, ...
    'exact_subset_summary.csv');

topLowerCSV = fullfile( ...
    resultsDir, ...
    sprintf('exact_top_subsets_K%d.csv',Klower));

topGreedyCSV = fullfile( ...
    resultsDir, ...
    sprintf('exact_top_subsets_K%d.csv',Kgreedy));

matFile = fullfile( ...
    resultsDir, ...
    'exact_subset_validation.mat');

writetable(Tsummary,summaryCSV);
writetable(TtopLower,topLowerCSV);
writetable(TtopGreedyK,topGreedyCSV);

save(matFile, ...
    'Tsummary', ...
    'TtopLower', ...
    'TtopGreedyK', ...
    'lowerResult', ...
    'greedyKResult', ...
    'globalMinimumProved', ...
    'globalMinimumK', ...
    'greedyIsBestBoundAtK', ...
    'selectedFA', ...
    'selectedBE', ...
    'greedyBoundDirect', ...
    'greedyPDOP', ...
    'target_m', ...
    'targetTol', ...
    'comparisonTol', ...
    'batchSize', ...
    'topN', ...
    'Label', ...
    'Architecture');

fprintf('\nSaved:\n');
fprintf('%s\n',summaryCSV);
fprintf('%s\n',topLowerCSV);
fprintf('%s\n',topGreedyCSV);
fprintf('%s\n',matFile);

fprintf('\n');
fprintf('============================================================\n');
fprintf('EXACT SUBSET VALIDATION COMPLETED\n');
fprintf('============================================================\n');

%% ============================================================
% Local function: exact enumeration for one K
% =============================================================

function result = enumerateExactK( ...
    N, ...
    K, ...
    JcontribFlat, ...
    target_m, ...
    targetTol, ...
    referenceBound_m, ...
    comparisonTol, ...
    batchSize, ...
    topN)

    totalSubsets = nchoosek(N,K);

    fprintf('Total subsets          : %.0f\n',totalSubsets);

    currentCombination = 1:K;
    done = false;

    processed = 0;
    nNumericallyValid = 0;
    nFeasible = 0;

    nStrictlyBetterThanReference = 0;
    nAtOrBelowReference = 0;

    topBounds = zeros(0,1);
    topCombinations = zeros(0,K,'uint8');

    nextProgress = min(500000,totalSubsets);

    while ~done

        [batch,currentCombination,done] = ...
            nextCombinationBatch( ...
                N, ...
                K, ...
                currentCombination, ...
                batchSize);

        nBatch = size(batch,1);

        bounds = ...
            boundsForCombinationBatch( ...
                batch, ...
                JcontribFlat);

        processed = processed + nBatch;

        finiteMask = isfinite(bounds);

        nNumericallyValid = ...
            nNumericallyValid + sum(finiteMask);

        nFeasible = ...
            nFeasible + sum(bounds <= target_m + targetTol);

        if isfinite(referenceBound_m)

            nStrictlyBetterThanReference = ...
                nStrictlyBetterThanReference + ...
                sum(bounds < referenceBound_m - comparisonTol);

            nAtOrBelowReference = ...
                nAtOrBelowReference + ...
                sum(bounds <= referenceBound_m + comparisonTol);

        end

        %% Maintain exact top-N subsets

        finiteIdx = find(finiteMask);

        if ~isempty(finiteIdx)

            finiteBounds = bounds(finiteIdx);

            nTake = min(topN,numel(finiteBounds));

            [batchTopBounds,batchOrder] = ...
                mink(finiteBounds,nTake);

            batchTopCombinations = ...
                batch(finiteIdx(batchOrder),:);

            candidateBounds = [ ...
                topBounds; ...
                batchTopBounds(:)];

            candidateCombinations = [ ...
                topCombinations; ...
                batchTopCombinations];

            nKeep = min(topN,numel(candidateBounds));

            [topBounds,keepOrder] = ...
                mink(candidateBounds,nKeep);

            topCombinations = ...
                candidateCombinations(keepOrder,:);

        end

        %% Progress

        if processed >= nextProgress || done

            fprintf( ...
                '  %9.0f / %9.0f (%.2f %%) | feasible=%d\n', ...
                processed, ...
                totalSubsets, ...
                100*processed/totalSubsets, ...
                nFeasible);

            while nextProgress <= processed
                nextProgress = nextProgress + 500000;
            end

        end

    end

    assert(processed == totalSubsets, ...
        ['Combination generator count mismatch: processed %.0f, ' ...
         'expected %.0f.'], ...
         processed,totalSubsets);

    assert(~isempty(topBounds), ...
        'No numerically valid subsets were found.');

    result = struct();

    result.K = K;
    result.TotalSubsets = totalSubsets;
    result.NNumericallyValid = nNumericallyValid;
    result.NFeasible = nFeasible;
    result.FeasibleFraction = nFeasible/totalSubsets;

    result.BestBound_m = topBounds(1);
    result.BestCombination = double(topCombinations(1,:));

    result.TopBounds_m = topBounds(:);
    result.TopCombinations = double(topCombinations);

    result.ReferenceBound_m = referenceBound_m;
    result.NStrictlyBetterThanReference = ...
        nStrictlyBetterThanReference;
    result.NAtOrBelowReference = ...
        nAtOrBelowReference;

end

%% ============================================================
% Local function: generate next lexicographic batch
% =============================================================

function [batch,nextCombination,done] = ...
    nextCombinationBatch( ...
        N, ...
        K, ...
        currentCombination, ...
        batchSize)

    batch = zeros(batchSize,K,'uint8');

    combination = double(currentCombination(:).');

    nWritten = 0;
    done = false;

    while nWritten < batchSize

        nWritten = nWritten + 1;

        batch(nWritten,:) = uint8(combination);

        %% Advance lexicographically

        p = K;

        while p >= 1 && ...
              combination(p) == N-K+p

            p = p - 1;

        end

        if p == 0

            done = true;
            nextCombination = [];
            break;

        end

        combination(p) = combination(p) + 1;

        for q = p+1:K
            combination(q) = combination(q-1) + 1;
        end

    end

    batch = batch(1:nWritten,:);

    if ~done
        nextCombination = combination;
    end

end

%% ============================================================
% Local function: vectorized positional bounds for a batch
% =============================================================

function bounds = ...
    boundsForCombinationBatch( ...
        combinations, ...
        JcontribFlat)

    nBatch = size(combinations,1);
    K = size(combinations,2);

    Jflat = zeros(nBatch,16);

    for j = 1:K

        idx = double(combinations(:,j));

        Jflat = Jflat + ...
            JcontribFlat(idx,:);

    end

    % Each row contains one 4x4 FIM in MATLAB column-major order.
    Jpages = reshape(Jflat.',4,4,nBatch);

    % Page-wise inverse through linear solves.
    [Cpages,rcondJ] = ...
        pagemldivide( ...
            Jpages, ...
            eye(4));

    trPos = reshape( ...
        Cpages(1,1,:) + ...
        Cpages(2,2,:) + ...
        Cpages(3,3,:), ...
        [],1);

    rcondJ = reshape(rcondJ,[],1);

    valid = ...
        rcondJ > 1e-14 & ...
        isfinite(trPos) & ...
        trPos > 0;

    bounds = inf(nBatch,1);

    bounds(valid) = ...
        sqrt(trPos(valid));

end

%% ============================================================
% Local function: direct metrics, identical convention to FA/BE
% =============================================================

function [bound_m,PDOP,rankH] = ...
    subsetMetricsDirect( ...
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

%% ============================================================
% Local function: print exact result
% =============================================================

function printExactResult(result,Label)

    fprintf('Numerically valid      : %d\n', ...
        result.NNumericallyValid);

    fprintf('Feasible subsets       : %d\n', ...
        result.NFeasible);

    fprintf('Feasible fraction      : %.12g %%\n', ...
        100*result.FeasibleFraction);

    fprintf('Best bound             : %.12f m\n', ...
        result.BestBound_m);

    fprintf('Best indices           : ');
    fprintf('%d ',result.BestCombination);
    fprintf('\n');

    fprintf('Best labels            : ');
    fprintf('%s ',Label(result.BestCombination));
    fprintf('\n');

    if isfinite(result.ReferenceBound_m)

        fprintf('Reference bound        : %.12f m\n', ...
            result.ReferenceBound_m);

        fprintf('Strictly better than reference : %d\n', ...
            result.NStrictlyBetterThanReference);

        fprintf('At/below reference (tol.)      : %d\n', ...
            result.NAtOrBelowReference);

    end

end

%% ============================================================
% Local function: top-subset table
% =============================================================

function T = makeTopTable( ...
    result, ...
    Label, ...
    Architecture)

    n = numel(result.TopBounds_m);

    Rank = (1:n).';
    K = repmat(result.K,n,1);
    PositionalBound_m = result.TopBounds_m(:);

    Indices = strings(n,1);
    Labels = strings(n,1);

    HAPS = zeros(n,1);
    LEO  = zeros(n,1);
    MEO  = zeros(n,1);
    GEO  = zeros(n,1);

    for i = 1:n

        idx = result.TopCombinations(i,:);

        Indices(i) = combinationToString(idx);
        Labels(i) = labelsToString(idx,Label);

        HAPS(i) = sum(Architecture(idx)=="HAPS");
        LEO(i)  = sum(Architecture(idx)=="LEO");
        MEO(i)  = sum(Architecture(idx)=="MEO");
        GEO(i)  = sum(Architecture(idx)=="GEO");

    end

    T = table( ...
        Rank, ...
        K, ...
        PositionalBound_m, ...
        HAPS, ...
        LEO, ...
        MEO, ...
        GEO, ...
        Indices, ...
        Labels);

end

%% ============================================================
% Local string helpers
% =============================================================

function s = combinationToString(idx)

    idx = double(idx(:).');

    if isempty(idx)
        s = "";
        return;
    end

    pieces = string(idx);

    s = strjoin(pieces," ");

end

function s = labelsToString(idx,Label)

    idx = double(idx(:).');

    if isempty(idx)
        s = "";
        return;
    end

    s = strjoin(Label(idx)," ");

end
