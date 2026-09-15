clear;
clc;

format long;

%% ============================================================
% EXHAUSTIVE EVALUATION OF ALL K = 11 SUBSETS FROM 26 Tx
%
% Purpose
%   Evaluate every nchoosek(26,11) = 7,726,160 subset once and save
%   the complete result for later figures, tables, and verification.
%
% Independence
%   This script does not call the Experiment-2 script or any external
%   project function. It requires only these consolidated MAT files:
%       results/candidate_pool_26.mat
%       results/link_budget.mat
%       results/FA_results.mat
%       results/BE_results.mat
%
% Output
%   results/experiment2_selection_robustness/
%       experiment2_exhaustive_K11_all_subsets.mat
%       experiment2_exhaustive_K11_summary.csv
%       experiment2_exhaustive_K11_checkpoint.mat
%
% The large MAT file stores, for every subset:
%   CombinationIndices : 7,726,160 x 11 uint8
%   AlphaCRB_m         : 7,726,160 x 1 double
%   HAPS/LEO/MEO/GEO   : architecture counts as uint8
%
% Safe restart
%   Results are written in blocks. If execution stops, run this same
%   script again. It resumes from the last completed block.
% =============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('EXHAUSTIVE K=11 SUBSET EVALUATION\n');
fprintf('============================================================\n');

%% ============================================================
% User settings
% =============================================================

Krequested = 11;
target_m = 0.6;
targetTol = 1e-12;

% Number of combinations committed per checkpoint.
% A smaller value loses less work after an interruption; a larger value
% reduces disk-write overhead. 25,000 is a conservative default.
blockSize = 25000;

% Must match the numerical safeguard in subsetMetrics from Experiment 2.
rcondThreshold = 1e-14;

%% ============================================================
% Locate the project root without requiring a .git directory
% =============================================================

thisFile = mfilename('fullpath');

if isempty(thisFile)
    error('Save this script before running it.');
end

searchDir = fileparts(thisFile);
repoRoot = '';

while true

    candidateTest = fullfile( ...
        searchDir,'results','candidate_pool_26.mat');

    linkBudgetTest = fullfile( ...
        searchDir,'results','link_budget.mat');

    if isfile(candidateTest) && isfile(linkBudgetTest)
        repoRoot = searchDir;
        break;
    end

    parentDir = fileparts(searchDir);

    if strcmp(parentDir,searchDir)
        break;
    end

    searchDir = parentDir;
end

if isempty(repoRoot)
    error([ ...
        'Could not find results/candidate_pool_26.mat and ' ...
        'results/link_budget.mat in this script''s directory ' ...
        'or any parent directory.']);
end

candidateFile = fullfile( ...
    repoRoot,'results','candidate_pool_26.mat');

linkBudgetFile = fullfile( ...
    repoRoot,'results','link_budget.mat');

faFile = fullfile( ...
    repoRoot,'results','FA_results.mat');

beFile = fullfile( ...
    repoRoot,'results','BE_results.mat');

assert(isfile(faFile), ...
    'FA result not found:\n%s',faFile);

assert(isfile(beFile), ...
    'BE result not found:\n%s',beFile);

outDir = fullfile( ...
    repoRoot, ...
    'results', ...
    'experiment2_selection_robustness');

if ~isfolder(outDir)
    mkdir(outDir);
end

outputMat = fullfile( ...
    outDir, ...
    'experiment2_exhaustive_K11_all_subsets.mat');

summaryCSV = fullfile( ...
    outDir, ...
    'experiment2_exhaustive_K11_summary.csv');

checkpointFile = fullfile( ...
    outDir, ...
    'experiment2_exhaustive_K11_checkpoint.mat');

checkpointTemp = fullfile( ...
    outDir, ...
    'experiment2_exhaustive_K11_checkpoint_tmp.mat');

%% ============================================================
% Load and validate nominal data
% =============================================================

S = load(candidateFile);
LB = load(linkBudgetFile);
FA = load(faFile);
BE = load(beFile);

requiredCandidateVars = [ ...
    "UserPositionECEF", ...
    "PositionECEF", ...
    "Label", ...
    "Architecture"];

for q = 1:numel(requiredCandidateVars)
    assert(isfield(S,requiredCandidateVars(q)), ...
        'Missing candidate-pool variable: %s', ...
        requiredCandidateVars(q));
end

requiredLinkVars = [ ...
    "Label", ...
    "Architecture", ...
    "SigmaRho_m"];

for q = 1:numel(requiredLinkVars)
    assert(isfield(LB,requiredLinkVars(q)), ...
        'Missing link-budget variable: %s', ...
        requiredLinkVars(q));
end

assert(isfield(FA,'selected'), ...
    'FA_results.mat does not contain selected.');

assert(isfield(BE,'selected'), ...
    'BE_results.mat does not contain selected.');

UserPositionECEF = double(S.UserPositionECEF(:));
PositionECEF = double(S.PositionECEF);
Label = string(S.Label(:));
Architecture = upper(string(S.Architecture(:)));
sigmaRho_m = double(LB.SigmaRho_m(:));

Npool = size(PositionECEF,1);
K = Krequested;

assert(Npool==26, ...
    'Expected 26 nominal candidates; found %d.',Npool);

assert(K>=4 && K<=Npool, ...
    'K must satisfy 4 <= K <= Npool.');

assert(isequal(Label,string(LB.Label(:))), ...
    'Candidate/link-budget label order mismatch.');

assert(isequal( ...
    Architecture, ...
    upper(string(LB.Architecture(:)))), ...
    'Candidate/link-budget architecture order mismatch.');

assert(all(isfinite(sigmaRho_m))); 
assert(all(sigmaRho_m>0));

selectedFA = sort(double(FA.selected(:).'));
selectedBE = sort(double(BE.selected(:).'));

assert(numel(selectedFA)==K, ...
    'The nominal FA solution has K=%d, not K=%d.', ...
    numel(selectedFA),K);

assert(numel(selectedBE)==K, ...
    'The nominal BE solution has K=%d, not K=%d.', ...
    numel(selectedBE),K);

%% ============================================================
% Geometry and additive information contributions
% =============================================================

diffTrue = PositionECEF - UserPositionECEF.';
trueRange_m = sqrt(sum(diffTrue.^2,2));

assert(all(trueRange_m>eps));

u = diffTrue ./ trueRange_m;
Hfull = [-u,ones(Npool,1)];

assert(rank(Hfull)==4, ...
    'Nominal candidate geometry is rank deficient.');

% With independent pseudorange errors,
% J(S) = sum_{i in S} (h_i^T h_i)/sigma_i^2.
% Precomputing each transmitter contribution avoids rebuilding R and
% repeatedly solving R\H for all 7.7 million subsets.
Jcontribution = zeros(4,4,Npool);

for i = 1:Npool
    hi = Hfull(i,:);
    Jcontribution(:,:,i) = ...
        (hi.'*hi)/(sigmaRho_m(i)^2);
end

architectureCode = zeros(Npool,1,'uint8');
architectureCode(Architecture=="HAPS") = uint8(1);
architectureCode(Architecture=="LEO")  = uint8(2);
architectureCode(Architecture=="MEO")  = uint8(3);
architectureCode(Architecture=="GEO")  = uint8(4);

assert(all(architectureCode>0), ...
    'Unknown architecture label found.');

totalCombinations = nchoosek(Npool,K);

assert(totalCombinations==7726160, ...
    'Unexpected number of K=11 subsets: %d.',totalCombinations);

fprintf('Project root       : %s\n',repoRoot);
fprintf('Candidate count    : %d\n',Npool);
fprintf('Subset size K      : %d\n',K);
fprintf('Total subsets      : %s\n', ...
    strrep(sprintf('%d',totalCombinations), ...
    '7726160','7,726,160'));
fprintf('Checkpoint block   : %d subsets\n',blockSize);

%% ============================================================
% Create a new result or resume a prior execution
% =============================================================

if isfile(checkpointFile)

    checkpoint = load(checkpointFile);

    requiredCheckpointVars = { ...
        'CompletedCount', ...
        'NextCombination', ...
        'RunCompleted', ...
        'NpoolSaved', ...
        'KSaved', ...
        'totalCombinationsSaved', ...
        'HfullSaved', ...
        'sigmaRhoSaved'};

    for q = 1:numel(requiredCheckpointVars)
        assert(isfield(checkpoint,requiredCheckpointVars{q}), ...
            'Checkpoint is missing %s.', ...
            requiredCheckpointVars{q});
    end

    assert(isfile(outputMat), ...
        'Checkpoint exists, but the large result MAT file is missing.');

    assert(checkpoint.NpoolSaved==Npool);
    assert(checkpoint.KSaved==K);
    assert(checkpoint.totalCombinationsSaved==totalCombinations);
    assert(isequaln(checkpoint.HfullSaved,Hfull), ...
        'Geometry differs from the checkpointed run.');
    assert(isequaln(checkpoint.sigmaRhoSaved,sigmaRho_m), ...
        'Measurement uncertainties differ from the checkpointed run.');

    CompletedCount = double(checkpoint.CompletedCount);
    NextCombination = double(checkpoint.NextCombination);
    RunCompleted = logical(checkpoint.RunCompleted);

    fprintf('\nExisting checkpoint found.\n');
    fprintf('Completed subsets  : %d / %d (%.3f %%)\n', ...
        CompletedCount,totalCombinations, ...
        100*CompletedCount/totalCombinations);

else

    if isfile(outputMat)
        error([ ...
            'The output MAT file exists without a checkpoint:\n%s\n' ...
            'Rename it or remove it manually before starting a new run.'], ...
            outputMat);
    end

    CompletedCount = 0;
    NextCombination = 1:K;
    RunCompleted = false;

    creationUTC = datetime('now','TimeZone','UTC');
    description = [ ...
        "Exact exhaustive evaluation of every K=11 subset " ...
        "from the nominal 26-transmitter candidate pool."];

    save( ...
        outputMat, ...
        'Npool', ...
        'K', ...
        'totalCombinations', ...
        'target_m', ...
        'targetTol', ...
        'rcondThreshold', ...
        'blockSize', ...
        'Hfull', ...
        'sigmaRho_m', ...
        'Label', ...
        'Architecture', ...
        'selectedFA', ...
        'selectedBE', ...
        'creationUTC', ...
        'description', ...
        '-v7.3');

    resultFile = matfile(outputMat,'Writable',true);

    % Preallocate directly on disk. Indices fit in uint8 because Npool=26.
    resultFile.CombinationIndices(totalCombinations,K) = uint8(0);
    resultFile.AlphaCRB_m(totalCombinations,1) = NaN;
    resultFile.HAPS(totalCombinations,1) = uint8(0);
    resultFile.LEO(totalCombinations,1) = uint8(0);
    resultFile.MEO(totalCombinations,1) = uint8(0);
    resultFile.GEO(totalCombinations,1) = uint8(0);

    clear resultFile;

    writeCheckpoint( ...
        checkpointTemp, ...
        checkpointFile, ...
        CompletedCount, ...
        NextCombination, ...
        RunCompleted, ...
        Npool, ...
        K, ...
        totalCombinations, ...
        Hfull, ...
        sigmaRho_m);

    fprintf('\nNew output initialized:\n%s\n',outputMat);
end

if RunCompleted
    fprintf('\nThe exhaustive run is already complete.\n');
    fprintf('Summary: %s\n',summaryCSV);
    fprintf('Results: %s\n',outputMat);
    return;
end

%% ============================================================
% Exhaustive evaluation in checkpointed blocks
% =============================================================

resultFile = matfile(outputMat,'Writable',true);
runTimer = tic;
lastReportTimer = tic;
completedAtRunStart = CompletedCount;

while CompletedCount < totalCombinations

    rowsThisBlock = min( ...
        blockSize, ...
        totalCombinations-CompletedCount);

    [combinationBlock,nextAfterBlock] = ...
        generateCombinationBlock( ...
            NextCombination, ...
            rowsThisBlock, ...
            Npool, ...
            K);

    actualRows = size(combinationBlock,1);

    if actualRows==0
        error('Combination generator stopped before the expected total.');
    end

    alphaBlock = inf(actualRows,1);

    for q = 1:actualRows

        idx = combinationBlock(q,:);

        J = sum(Jcontribution(:,:,idx),3);

        if rcond(J)>rcondThreshold

            C = J\eye(4);
            trPos = trace(C(1:3,1:3));

            if trPos>0 && isfinite(trPos)
                alphaBlock(q) = sqrt(trPos);
            end
        end
    end

    archBlock = architectureCode(combinationBlock);

    hapsBlock = uint8(sum(archBlock==1,2));
    leoBlock  = uint8(sum(archBlock==2,2));
    meoBlock  = uint8(sum(archBlock==3,2));
    geoBlock  = uint8(sum(archBlock==4,2));

    firstRow = CompletedCount+1;
    lastRow = CompletedCount+actualRows;

    % Commit the numerical block before advancing the checkpoint.
    resultFile.CombinationIndices(firstRow:lastRow,:) = ...
        uint8(combinationBlock);

    resultFile.AlphaCRB_m(firstRow:lastRow,1) = alphaBlock;
    resultFile.HAPS(firstRow:lastRow,1) = hapsBlock;
    resultFile.LEO(firstRow:lastRow,1) = leoBlock;
    resultFile.MEO(firstRow:lastRow,1) = meoBlock;
    resultFile.GEO(firstRow:lastRow,1) = geoBlock;

    CompletedCount = lastRow;
    NextCombination = nextAfterBlock;

    writeCheckpoint( ...
        checkpointTemp, ...
        checkpointFile, ...
        CompletedCount, ...
        NextCombination, ...
        RunCompleted, ...
        Npool, ...
        K, ...
        totalCombinations, ...
        Hfull, ...
        sigmaRho_m);

    if toc(lastReportTimer)>=10 || ...
       CompletedCount==totalCombinations

        elapsed_s = toc(runTimer);
        processedThisRun = ...
            CompletedCount-completedAtRunStart;
        rate_per_s = processedThisRun/max(elapsed_s,eps);
        remaining_s = ...
            (totalCombinations-CompletedCount)/max(rate_per_s,eps);

        fprintf([ ...
            'Completed %d / %d (%.3f %%) | ' ...
            '%.1f subsets/s | ETA this run: %s\n'], ...
            CompletedCount, ...
            totalCombinations, ...
            100*CompletedCount/totalCombinations, ...
            rate_per_s, ...
            formatDuration(remaining_s));

        lastReportTimer = tic;
    end
end

clear resultFile;

%% ============================================================
% Exact population summary
% =============================================================

fprintf('\nComputing exact population statistics...\n');

resultFile = matfile(outputMat);
AlphaCRB_m = resultFile.AlphaCRB_m;

finiteMask = isfinite(AlphaCRB_m);
finiteBounds = AlphaCRB_m(finiteMask);

assert(numel(AlphaCRB_m)==totalCombinations);
assert(~isempty(finiteBounds));

NFinite = sum(finiteMask);
NInvalid = totalCombinations-NFinite;
NTargetReached = sum( ...
    finiteMask & AlphaCRB_m<=target_m+targetTol);

TargetReached_pct = ...
    100*NTargetReached/totalCombinations;

MeanBound_m = mean(finiteBounds);
MedianBound_m = median(finiteBounds);
P01Bound_m = empiricalPercentile(finiteBounds,1);
P05Bound_m = empiricalPercentile(finiteBounds,5);
P25Bound_m = empiricalPercentile(finiteBounds,25);
P75Bound_m = empiricalPercentile(finiteBounds,75);
P95Bound_m = empiricalPercentile(finiteBounds,95);
P99Bound_m = empiricalPercentile(finiteBounds,99);
[MinimumBound_m,bestRow] = min(AlphaCRB_m);
MaximumFiniteBound_m = max(finiteBounds);

BestIndices = double( ...
    resultFile.CombinationIndices(bestRow,:));

faRow = lexicographicCombinationRank(selectedFA,Npool,K);
beRow = lexicographicCombinationRank(selectedBE,Npool,K);

FAAlphaCRB_m = AlphaCRB_m(faRow);
BEAlphaCRB_m = AlphaCRB_m(beRow);

NStrictlyBetterThanFA = sum( ...
    finiteMask & AlphaCRB_m<FAAlphaCRB_m-targetTol);

NAtOrBelowFA = sum( ...
    finiteMask & AlphaCRB_m<=FAAlphaCRB_m+targetTol);

NStrictlyBetterThanBE = sum( ...
    finiteMask & AlphaCRB_m<BEAlphaCRB_m-targetTol);

NAtOrBelowBE = sum( ...
    finiteMask & AlphaCRB_m<=BEAlphaCRB_m+targetTol);

Tsummary = table( ...
    Npool, ...
    K, ...
    totalCombinations, ...
    target_m, ...
    NFinite, ...
    NInvalid, ...
    NTargetReached, ...
    TargetReached_pct, ...
    MeanBound_m, ...
    MedianBound_m, ...
    P01Bound_m, ...
    P05Bound_m, ...
    P25Bound_m, ...
    P75Bound_m, ...
    P95Bound_m, ...
    P99Bound_m, ...
    MinimumBound_m, ...
    MaximumFiniteBound_m, ...
    FAAlphaCRB_m, ...
    BEAlphaCRB_m, ...
    NStrictlyBetterThanFA, ...
    NAtOrBelowFA, ...
    NStrictlyBetterThanBE, ...
    NAtOrBelowBE);

writetable(Tsummary,summaryCSV);

completionUTC = datetime('now','TimeZone','UTC');

save( ...
    outputMat, ...
    'Tsummary', ...
    'BestIndices', ...
    'bestRow', ...
    'faRow', ...
    'beRow', ...
    'completionUTC', ...
    '-append');

RunCompleted = true;
NextCombination = zeros(1,0);

writeCheckpoint( ...
    checkpointTemp, ...
    checkpointFile, ...
    CompletedCount, ...
    NextCombination, ...
    RunCompleted, ...
    Npool, ...
    K, ...
    totalCombinations, ...
    Hfull, ...
    sigmaRho_m);

disp(Tsummary);

fprintf('Best subset indices: %s\n', ...
    strjoin(string(BestIndices),', '));

fprintf('\nSaved summary:\n%s\n',summaryCSV);
fprintf('\nSaved complete population:\n%s\n',outputMat);
fprintf('\nCheckpoint retained as completion record:\n%s\n', ...
    checkpointFile);

fprintf('\n');
fprintf('============================================================\n');
fprintf('EXHAUSTIVE K=11 EVALUATION COMPLETED\n');
fprintf('============================================================\n');

%% ============================================================
% Local functions
% =============================================================

function [block,nextCombination] = ...
    generateCombinationBlock( ...
        firstCombination, ...
        requestedRows, ...
        N, ...
        K)

    if isempty(firstCombination)
        block = zeros(0,K);
        nextCombination = zeros(1,0);
        return;
    end

    block = zeros(requestedRows,K,'uint8');
    current = double(firstCombination(:).');
    actualRows = 0;

    for row = 1:requestedRows

        actualRows = actualRows+1;
        block(actualRows,:) = uint8(current);

        [current,hasNext] = nextCombinationLexicographic( ...
            current,N,K);

        if ~hasNext
            current = zeros(1,0);
            break;
        end
    end

    block = block(1:actualRows,:);
    nextCombination = current;
end

function [combination,hasNext] = ...
    nextCombinationLexicographic(combination,N,K)

    position = K;

    while position>=1 && ...
          combination(position)==N-K+position
        position = position-1;
    end

    if position==0
        hasNext = false;
        return;
    end

    combination(position) = combination(position)+1;

    for q = position+1:K
        combination(q) = combination(q-1)+1;
    end

    hasNext = true;
end

function row = ...
    lexicographicCombinationRank(combination,N,K)

    combination = sort(double(combination(:).'));

    assert(numel(combination)==K);
    assert(all(diff(combination)>0));
    assert(combination(1)>=1 && combination(end)<=N);

    zeroBasedRank = 0;
    previous = 0;

    for position = 1:K

        for candidate = previous+1:combination(position)-1
            zeroBasedRank = zeroBasedRank + ...
                nchoosek(N-candidate,K-position);
        end

        previous = combination(position);
    end

    row = zeroBasedRank+1;
end

function p = empiricalPercentile(x,percent)

    x = sort(x(isfinite(x)));
    x = x(:);

    if isempty(x)
        p = NaN;
        return;
    end

    if numel(x)==1
        p = x;
        return;
    end

    position = 1+(numel(x)-1)*(percent/100);
    lowerIndex = floor(position);
    upperIndex = ceil(position);

    if lowerIndex==upperIndex
        p = x(lowerIndex);
    else
        weight = position-lowerIndex;
        p = ...
            (1-weight)*x(lowerIndex) ...
            + weight*x(upperIndex);
    end
end

function writeCheckpoint( ...
    temporaryFile, ...
    checkpointFile, ...
    CompletedCount, ...
    NextCombination, ...
    RunCompleted, ...
    NpoolSaved, ...
    KSaved, ...
    totalCombinationsSaved, ...
    HfullSaved, ...
    sigmaRhoSaved)

    checkpointUTC = datetime('now','TimeZone','UTC');

    save( ...
        temporaryFile, ...
        'CompletedCount', ...
        'NextCombination', ...
        'RunCompleted', ...
        'NpoolSaved', ...
        'KSaved', ...
        'totalCombinationsSaved', ...
        'HfullSaved', ...
        'sigmaRhoSaved', ...
        'checkpointUTC');

    [moveSucceeded,moveMessage] = ...
        movefile(temporaryFile,checkpointFile,'f');

    if ~moveSucceeded
        error('Could not update checkpoint: %s',moveMessage);
    end
end

function text = formatDuration(seconds)

    if ~isfinite(seconds)
        text = 'unknown';
        return;
    end

    seconds = max(0,round(seconds));
    hours = floor(seconds/3600);
    minutes = floor(mod(seconds,3600)/60);
    remainingSeconds = mod(seconds,60);

    text = sprintf( ...
        '%02d:%02d:%02d', ...
        hours,minutes,remainingSeconds);
end
