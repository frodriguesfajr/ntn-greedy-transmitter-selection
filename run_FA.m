clear;
clc;
format long;

%% ============================================================
% FORWARD ADDING (FA)
%
% Objective:
%   Apply Forward Adding (FA) to the nominal 26-transmitter candidate
%   pool and determine a small transmitter subset that satisfies the
%   prescribed positional-bound target.
%
% The measurement-quality chain is
%
%   FSPL_i -> C/N0_i -> sigma_rho_i -> R -> positional bound.
%
% Ranging-signal assumptions:
%   beta   = 1.023 MHz
%   Tcoh   = 20 ms
%   etaPos = 0.01, fraction of link power allocated to the positioning
%            signal (1%, equivalent to -20 dB), common to all architectures.
%
% Selection procedure:
%   1) exhaustively evaluate all full-rank four-transmitter subsets;
%   2) initialize FA with the subset having the smallest positional bound;
%   3) at each iteration, add the remaining candidate that minimizes the
%      positional bound;
%   4) stop when the positional bound reaches 0.6 m or no candidate remains.
%
% Geometry is read from the candidate-pool file, while nominal
% measurement quality is read from the link-budget file.
%
% Outputs:
%   results/FA_results.mat
%   results/FA_selected.csv
%   results/FA_trace.csv
%% ============================================================

rootDir = fileparts(mfilename('fullpath'));

if isempty(rootDir)
    rootDir = pwd;
end

cd(rootDir);
setup_paths;

fprintf('\n============================================================\n');
fprintf('FORWARD ADDING\n');
fprintf('============================================================\n');

%% ============================================================
% Input files
%% ============================================================

candidateMat = fullfile( ...
    rootDir, ...
    'results', ...
    'candidate_pool_26.mat');

linkBudgetMat = fullfile( ...
    rootDir, ...
    'results', ...
    'link_budget.mat');

assert(isfile(candidateMat), ...
    'Candidate-pool MAT file not found:\n%s',candidateMat);

assert(isfile(linkBudgetMat), ...
    'Link-budget MAT file not found:\n%s',linkBudgetMat);

%% ============================================================
% Candidate pool
%% ============================================================

S = load(candidateMat);

requiredCandidateVars = [ ...
    "UserPositionECEF", ...
    "PositionECEF", ...
    "Label", ...
    "Architecture", ...
    "NORAD", ...
    "ObjectName", ...
    "Elevation_deg", ...
    "SlantRange_m"];

for k = 1:numel(requiredCandidateVars)

    if ~isfield(S,requiredCandidateVars(k))

        error( ...
            'Missing candidate-pool variable: %s', ...
            requiredCandidateVars(k));

    end
end

UserPositionECEF = S.UserPositionECEF;

PositionECEF = S.PositionECEF;

Label = string(S.Label(:));

Architecture = string(S.Architecture(:));

NORAD = string(S.NORAD(:));

ObjectName = string(S.ObjectName(:));

Elevation_deg = S.Elevation_deg(:);

SlantRange_m = S.SlantRange_m(:);

Npool = size(PositionECEF,1);

assert(Npool == 26, ...
    'Expected 26 candidates; obtained %d.',Npool);

%% ============================================================
% Link budget
%% ============================================================

LB = load(linkBudgetMat);

requiredLinkVars = [ ...
    "Label", ...
    "Architecture", ...
    "CN0_dBHz", ...
    "SigmaRho_m", ...
    "VarianceRho_m2", ...
    "beta", ...
    "Tcoh", ...
    "etaPos"];

for k = 1:numel(requiredLinkVars)

    if ~isfield(LB,requiredLinkVars(k))

        error( ...
            'Missing link-budget variable: %s', ...
            requiredLinkVars(k));

    end
end

LabelLB = string(LB.Label(:));

ArchitectureLB = string(LB.Architecture(:));

assert(isequal(Label,LabelLB), ...
    'Link-budget transmitter order differs from candidate pool.');

assert(isequal(Architecture,ArchitectureLB), ...
    'Link-budget architecture order differs from candidate pool.');

CN0_dBHz = LB.CN0_dBHz(:);

sigmaRho_m = LB.SigmaRho_m(:);

VarianceRho_m2 = LB.VarianceRho_m2(:);

beta = LB.beta;
Tcoh = LB.Tcoh;
etaPos = LB.etaPos;

assert(all(isfinite(CN0_dBHz)));
assert(all(isfinite(sigmaRho_m)));
assert(all(sigmaRho_m > 0));

%% Optional link-budget fields

FSPL_dB = optionalVector(LB,'FSPL_dB',Npool);

Frequency_GHz = ...
    optionalVector(LB,'Frequency_GHz',Npool);

EIRP_PSD_dBW_MHz = ...
    optionalVector(LB,'EIRP_PSD_dBW_MHz',Npool);

EIRP_positioning_dBW = ...
    optionalVector(LB,'EIRP_positioning_dBW',Npool);

AdditionalLoss_dB = ...
    optionalVector(LB,'AdditionalLoss_dB',Npool);

if isfield(LB,'ReceiverChain')
    ReceiverChain = string(LB.ReceiverChain(:));
else
    ReceiverChain = repmat("",Npool,1);
end

%% ============================================================
% Geometry matrix H
%% ============================================================

r = PositionECEF - UserPositionECEF(:).';

d = sqrt(sum(r.^2,2));

if any(d <= eps)
    error('Invalid transmitter-user distance.');
end

u = r ./ d;

% Geometry-matrix convention used by the selection algorithm
Hfull = [ ...
    -u, ...
    ones(Npool,1)];

%% ============================================================
% Quality table
%% ============================================================

PoolIndex = (1:Npool).';

Tquality = table( ...
    PoolIndex, ...
    Label, ...
    Architecture, ...
    NORAD, ...
    ObjectName, ...
    Elevation_deg, ...
    SlantRange_m/1e3, ...
    Frequency_GHz, ...
    EIRP_PSD_dBW_MHz, ...
    AdditionalLoss_dB, ...
    FSPL_dB, ...
    EIRP_positioning_dBW, ...
    CN0_dBHz, ...
    sigmaRho_m, ...
    VarianceRho_m2, ...
    ReceiverChain, ...
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
    'AdditionalLoss_dB', ...
    'FSPL_dB', ...
    'EIRP_positioning_dBW', ...
    'CN0_dBHz', ...
    'SigmaRho_m', ...
    'VarianceRho_m2', ...
    'ReceiverChain'});

%% ============================================================
% Full-pool metrics
%% ============================================================

allIdx = (1:Npool).';

[fullBound,fullPDOP,fullRankH] = ...
    subsetMetrics( ...
    allIdx, ...
    Hfull, ...
    sigmaRho_m);

fprintf('\nFull 26-Tx pool\n');
fprintf('Positional bound : %.12f m\n',fullBound);
fprintf('PDOP             : %.12f\n',fullPDOP);
fprintf('rank(H)          : %d\n',fullRankH);

%% ============================================================
% Selection target
%% ============================================================

target_m = 0.6;

fprintf('\nSelection target: positional bound <= %.3f m\n', ...
    target_m);

%% ============================================================
% Exhaustive best-K=4 initialization
%% ============================================================

fprintf('\n============================================================\n');
fprintf('EXHAUSTIVE SEARCH FOR BEST K=4 SUBSET\n');
fprintf('============================================================\n');

comb4 = nchoosek(1:Npool,4);

nComb4 = size(comb4,1);

best4Bound = inf;
best4PDOP  = inf;
best4 = [];

for k = 1:nComb4

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
    error('No full-rank K=4 subset was found.');
end

fprintf('Subsets evaluated : %d\n',nComb4);
fprintf('Best K=4 bound    : %.12f m\n',best4Bound);
fprintf('Best K=4 PDOP     : %.12f\n',best4PDOP);

fprintf('Best K=4 indices  : ');
fprintf('%d ',best4);
fprintf('\n');

fprintf('Best K=4 labels   : ');
for k = 1:numel(best4)
    fprintf('%s ',char(Label(best4(k))));
end
fprintf('\n');

fprintf('Composition        : HAPS=%d | LEO=%d | MEO=%d | GEO=%d\n', ...
    sum(Architecture(best4)=="HAPS"), ...
    sum(Architecture(best4)=="LEO"), ...
    sum(Architecture(best4)=="MEO"), ...
    sum(Architecture(best4)=="GEO"));

%% ============================================================
% Forward Adding
%% ============================================================

selected = best4(:).';

currentBound = best4Bound;

traceK = numel(selected);
traceAddedIndex = NaN;
traceAddedLabel = "";
traceAddedArchitecture = "";
traceBound = currentBound;

% Evaluation counter:
% includes exhaustive K=4 search and every FA trial.
evaluationCount = nComb4;

fprintf('\n============================================================\n');
fprintf('FORWARD ADDING\n');
fprintf('============================================================\n');

fprintf( ...
    'K=%d | positional bound = %.12f m | initialization\n', ...
    numel(selected), ...
    currentBound);

while currentBound > target_m && ...
      numel(selected) < Npool

    remaining = setdiff( ...
        1:Npool, ...
        selected, ...
        'stable');

    bestCandidate = NaN;
    bestCandidateBound = inf;

    for k = 1:numel(remaining)

        candidate = remaining(k);

        idxTrial = [ ...
            selected, ...
            candidate];

        [b,~,rankH] = ...
            subsetMetrics( ...
            idxTrial, ...
            Hfull, ...
            sigmaRho_m);

        evaluationCount = ...
            evaluationCount + 1;

        if rankH == 4 && ...
           b < bestCandidateBound

            bestCandidateBound = b;
            bestCandidate = candidate;

        end
    end

    if isnan(bestCandidate)

        warning( ...
            'No valid additional candidate was found.');

        break;

    end

    selected(end+1) = ...
        bestCandidate; %#ok<SAGROW>

    currentBound = ...
        bestCandidateBound;

    traceK(end+1,1) = ...
        numel(selected); %#ok<SAGROW>

    traceAddedIndex(end+1,1) = ...
        bestCandidate; %#ok<SAGROW>

    traceAddedLabel(end+1,1) = ...
        Label(bestCandidate); %#ok<SAGROW>

    traceAddedArchitecture(end+1,1) = ...
        Architecture(bestCandidate); %#ok<SAGROW>

    traceBound(end+1,1) = ...
        currentBound; %#ok<SAGROW>

    fprintf( ...
        'K=%d | add %-6s (%s) | positional bound = %.12f m\n', ...
        numel(selected), ...
        char(Label(bestCandidate)), ...
        char(Architecture(bestCandidate)), ...
        currentBound);

end

%% ============================================================
% Final result
%% ============================================================

[finalBound,finalPDOP,finalRankH] = ...
    subsetMetrics( ...
    selected, ...
    Hfull, ...
    sigmaRho_m);

targetReached = ...
    finalBound <= target_m;

Tselected = ...
    Tquality(selected,:);

Tselected.SelectionOrder = ...
    (1:height(Tselected)).';

Tselected = movevars( ...
    Tselected, ...
    'SelectionOrder', ...
    'Before', ...
    'PoolIndex');

fprintf('\n============================================================\n');
fprintf('FA RESULT\n');
fprintf('============================================================\n');

fprintf('Final K             : %d\n',numel(selected));

fprintf('Positional bound    : %.12f m\n', ...
    finalBound);

fprintf('PDOP                : %.12f\n', ...
    finalPDOP);

fprintf('rank(H)             : %d\n', ...
    finalRankH);

fprintf('Target <= %.3f m    : %s\n', ...
    target_m,string(targetReached));

fprintf('Evaluations         : %d\n', ...
    evaluationCount);

fprintf('Composition         : HAPS=%d | LEO=%d | MEO=%d | GEO=%d\n', ...
    sum(Architecture(selected)=="HAPS"), ...
    sum(Architecture(selected)=="LEO"), ...
    sum(Architecture(selected)=="MEO"), ...
    sum(Architecture(selected)=="GEO"));

fprintf('\nSelected indices:\n');

disp(selected);

fprintf('\nSelected transmitters:\n');

disp(Tselected(:,{ ...
    'SelectionOrder', ...
    'PoolIndex', ...
    'Label', ...
    'Architecture', ...
    'NORAD_CAT_ID', ...
    'ObjectName', ...
    'SigmaRho_m'}));

%% ============================================================
% Trace
%% ============================================================

Ttrace = table( ...
    traceK(:), ...
    traceAddedIndex(:), ...
    traceAddedLabel(:), ...
    traceAddedArchitecture(:), ...
    traceBound(:), ...
    'VariableNames',{ ...
    'K', ...
    'AddedPoolIndex', ...
    'AddedLabel', ...
    'AddedArchitecture', ...
    'PositionalBound_m'});
%% ============================================================
% Save results
%% ============================================================

resultsDir = fullfile( ...
    rootDir, ...
    'results');

if ~exist(resultsDir,'dir')
    mkdir(resultsDir);
end

matFile = fullfile( ...
    resultsDir, ...
    'FA_results.mat');

selectedCSV = fullfile( ...
    resultsDir, ...
    'FA_selected.csv');

traceCSV = fullfile( ...
    resultsDir, ...
    'FA_trace.csv');

save(matFile, ...
    'candidateMat', ...
    'linkBudgetMat', ...
    'target_m', ...
    'Hfull', ...
    'sigmaRho_m', ...
    'VarianceRho_m2', ...
    'CN0_dBHz', ...
    'best4', ...
    'best4Bound', ...
    'selected', ...
    'finalBound', ...
    'finalPDOP', ...
    'finalRankH', ...
    'targetReached', ...
    'evaluationCount', ...
    'Tselected', ...
    'Ttrace', ...
    'fullBound', ...
    'fullPDOP', ...
    'fullRankH', ...
    'beta', ...
    'Tcoh', ...
    'etaPos');

writetable( ...
    Tselected, ...
    selectedCSV);

writetable( ...
    Ttrace, ...
    traceCSV);

fprintf('\nSaved:\n%s\n',matFile);
fprintf('%s\n',selectedCSV);
fprintf('%s\n',traceCSV);

fprintf('\n============================================================\n');
fprintf('FORWARD ADDING COMPLETED\n');
fprintf('============================================================\n');


%% ============================================================
% Local weighted positional metric
%% ============================================================

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

    %% Pure geometric PDOP

    A = H.'*H;

    if rcond(A) > 1e-14

        Cgeom = A\eye(4);

        PDOP = sqrt( ...
            trace(Cgeom(1:3,1:3)));

    end

    %% Weighted positional bound

    R = diag(sigma.^2);

    J = H.'*(R\H);

    if rcond(J) <= 1e-14
        return;
    end

    C = J\eye(4);

    trPos = trace( ...
        C(1:3,1:3));

    if trPos > 0 && ...
       isfinite(trPos)

        bound_m = sqrt(trPos);

    end

end


%% ============================================================
% Optional vector helper
%% ============================================================

function v = optionalVector(S,fieldName,N)

    if isfield(S,fieldName)

        v = S.(fieldName);

        v = v(:);

        if numel(v) ~= N

            error( ...
                'Field %s has %d elements; expected %d.', ...
                fieldName,numel(v),N);

        end

    else

        v = nan(N,1);

    end

end