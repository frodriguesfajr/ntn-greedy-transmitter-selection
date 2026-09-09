clear;
clc;
format long;

%% ============================================================
% BACKWARD ELIMINATION - GEOMETRY
%
% Same BE algorithm used in the paper baseline.
%
%
%   MATLAB satelliteScenario SGP4
%                  ->
%   Vallado SGP4 / WGS-72
%
% Target:
%   positional bound <= 0.6 m
%% ============================================================

rootDir = fileparts(mfilename('fullpath'));

if isempty(rootDir)
    rootDir = pwd;
end

cd(rootDir);
setup_paths;

fprintf('\n============================================================\n');
fprintf('BACKWARD ELIMINATION - GEOMETRY\n');
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

% Force row form for compatibility with 26x3 PositionECEF
UserPositionECEF = S.UserPositionECEF(:).';

PositionECEF = S.PositionECEF;

Label = string(S.Label(:));
Architecture = string(S.Architecture(:));
NORAD = string(S.NORAD(:));
ObjectName = string(S.ObjectName(:));

Elevation_deg = S.Elevation_deg(:);
SlantRange_m  = S.SlantRange_m(:);

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
    "VarianceRho_m2"];

for k = 1:numel(requiredLinkVars)

    if ~isfield(LB,requiredLinkVars(k))
        error( ...
            'Missing link-budget variable: %s', ...
            requiredLinkVars(k));
    end

end

LabelLB = string(LB.Label(:));
ArchitectureLB = string(LB.Architecture(:));

if ~isequal(LabelLB,Label)
    error('Label order differs between candidate pool and link budget.');
end

if ~isequal(ArchitectureLB,Architecture)
    error('Architecture order differs between candidate pool and link budget.');
end

CN0_dBHz = LB.CN0_dBHz(:);
sigmaRho_m = LB.SigmaRho_m(:);

VarianceRho_m2 = sigmaRho_m.^2;

if any(~isfinite(sigmaRho_m)) || ...
   any(sigmaRho_m <= 0)

    error('SigmaRho_m contains invalid values.');

end

%% ============================================================
% Parameters
%% ============================================================

target_m = 0.6;

% Same numerical tolerance as original BE
targetTol = 1e-12;

%% ============================================================
% Full geometry matrix
%% ============================================================

r = PositionECEF - UserPositionECEF;

d = sqrt(sum(r.^2,2));

if any(~isfinite(d)) || any(d <= eps)
    error('Invalid transmitter-user distance.');
end

u = r ./ d;

% Exact convention used by the original algorithm
Hfull = [-u,ones(Npool,1)];

%% ============================================================
% Measurement-quality table
%% ============================================================

PoolIndex = (1:Npool).';

if isfield(LB,'Tlink') && ...
   istable(LB.Tlink) && ...
   height(LB.Tlink) == Npool

    Tquality = LB.Tlink;

    if ~ismember( ...
            'PoolIndex', ...
            Tquality.Properties.VariableNames)

        Tquality.PoolIndex = PoolIndex;

        Tquality = movevars( ...
            Tquality, ...
            'PoolIndex', ...
            'Before',1);

    end

else

    Tquality = table( ...
        PoolIndex, ...
        Label, ...
        Architecture, ...
        NORAD, ...
        ObjectName, ...
        Elevation_deg, ...
        SlantRange_m/1e3, ...
        CN0_dBHz, ...
        sigmaRho_m, ...
        VarianceRho_m2, ...
        'VariableNames',{ ...
        'PoolIndex', ...
        'Label', ...
        'Architecture', ...
        'NORAD_CAT_ID', ...
        'ObjectName', ...
        'Elevation_deg', ...
        'SlantRange_km', ...
        'CN0_dBHz', ...
        'SigmaRho_m', ...
        'VarianceRho_m2'});

end

%% ============================================================
% Initial full-set diagnostics
%% ============================================================

allIdx = 1:Npool;

[fullBound,fullPDOP,fullRankH] = ...
    subsetMetrics( ...
    allIdx, ...
    Hfull, ...
    sigmaRho_m);

fprintf('\n===== Full initial set =====\n');

fprintf('Initial K            : %d\n',Npool);
fprintf('Positional bound     : %.12f m\n',fullBound);
fprintf('PDOP                 : %.12f\n',fullPDOP);
fprintf('rank(H)              : %d\n',fullRankH);
fprintf('Initial target met   : %s\n', ...
    string(fullBound <= target_m + targetTol));

if fullRankH < 4
    error('The full candidate set does not have rank(H)=4.');
end

if fullBound > target_m + targetTol

    error([ ...
        'The full candidate set does not satisfy the %.3f m target. ' ...
        'Backward Elimination cannot be initialized.'], ...
        target_m);

end

%% ============================================================
% Backward Elimination
%% ============================================================

selected = allIdx;

currentBound = fullBound;
currentPDOP  = fullPDOP;

traceK = Npool;

traceRemovedPoolIndex = NaN;
traceRemovedLabel = "";
traceRemovedArchitecture = "";

traceBound = fullBound;
tracePDOP  = fullPDOP;

traceCandidatesEvaluated = 0;

totalSubsetEvaluations = 0;
nAcceptedRemovals = 0;

fprintf('\n============================================================\n');
fprintf('BACKWARD ELIMINATION\n');
fprintf('============================================================\n');

fprintf( ...
    'K=%d | bound = %.12f m | full set\n', ...
    numel(selected), ...
    currentBound);

while numel(selected) > 4

    removable = selected;

    bestRemoved = NaN;
    bestTrialBound = inf;
    bestTrialPDOP = inf;

    nEvaluatedThisIteration = 0;

    for k = 1:numel(removable)

        candidateToRemove = removable(k);

        idxTrial = ...
            selected(selected ~= candidateToRemove);

        [trialBound,trialPDOP,trialRankH] = ...
            subsetMetrics( ...
            idxTrial, ...
            Hfull, ...
            sigmaRho_m);

        totalSubsetEvaluations = ...
            totalSubsetEvaluations + 1;

        nEvaluatedThisIteration = ...
            nEvaluatedThisIteration + 1;

        %% Removal must preserve rank and target

        if trialRankH ~= 4
            continue;
        end

        if trialBound > target_m + targetTol
            continue;
        end

        %% Select least detrimental removal

        if trialBound < bestTrialBound

            bestTrialBound = trialBound;
            bestTrialPDOP = trialPDOP;
            bestRemoved = candidateToRemove;

        elseif abs(trialBound-bestTrialBound) <= 1e-14

            % Deterministic tie break:
            % remove smallest PoolIndex.

            if isnan(bestRemoved) || ...
               candidateToRemove < bestRemoved

                bestTrialBound = trialBound;
                bestTrialPDOP = trialPDOP;
                bestRemoved = candidateToRemove;

            end
        end
    end

    %% Stopping criterion

    if isnan(bestRemoved)

        fprintf([ ...
            'Stop: no additional removal preserves ' ...
            'rank(H)=4 and bound <= %.3f m.\n'], ...
            target_m);

        break;

    end

    %% Accept removal

    selected(selected == bestRemoved) = [];

    currentBound = bestTrialBound;
    currentPDOP  = bestTrialPDOP;

    nAcceptedRemovals = ...
        nAcceptedRemovals + 1;

    traceK(end+1,1) = ...
        numel(selected); %#ok<SAGROW>

    traceRemovedPoolIndex(end+1,1) = ...
        bestRemoved; %#ok<SAGROW>

    traceRemovedLabel(end+1,1) = ...
        Label(bestRemoved); %#ok<SAGROW>

    traceRemovedArchitecture(end+1,1) = ...
        Architecture(bestRemoved); %#ok<SAGROW>

    traceBound(end+1,1) = ...
        currentBound; %#ok<SAGROW>

    tracePDOP(end+1,1) = ...
        currentPDOP; %#ok<SAGROW>

    traceCandidatesEvaluated(end+1,1) = ...
        nEvaluatedThisIteration; %#ok<SAGROW>

    fprintf([ ...
        'K=%d | remove %-6s (%-4s) | ' ...
        'bound = %.12f m | PDOP = %.12f\n'], ...
        numel(selected), ...
        char(Label(bestRemoved)), ...
        char(Architecture(bestRemoved)), ...
        currentBound, ...
        currentPDOP);

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
    finalBound <= target_m + targetTol;

nHAPS = sum(Architecture(selected)=="HAPS");
nLEO  = sum(Architecture(selected)=="LEO");
nMEO  = sum(Architecture(selected)=="MEO");
nGEO  = sum(Architecture(selected)=="GEO");

Tselected = ...
    Tquality(selected,:);

Tselected.FinalSetOrder = ...
    (1:height(Tselected)).';

Tselected = movevars( ...
    Tselected, ...
    'FinalSetOrder', ...
    'Before',1);

fprintf('\n============================================================\n');
fprintf('BE RESULT\n');
fprintf('============================================================\n');

fprintf('Initial K               : %d\n',Npool);
fprintf('Final K                 : %d\n',numel(selected));
fprintf('Accepted removals       : %d\n',nAcceptedRemovals);
fprintf('Subsets evaluated       : %d\n',totalSubsetEvaluations);

fprintf('Final positional bound  : %.12f m\n', ...
    finalBound);

fprintf('Final PDOP              : %.12f\n', ...
    finalPDOP);

fprintf('Final rank(H)           : %d\n', ...
    finalRankH);

fprintf('Target <= %.3f m        : %s\n', ...
    target_m, ...
    string(targetReached));

fprintf( ...
    'Final composition        : HAPS=%d | LEO=%d | MEO=%d | GEO=%d\n', ...
    nHAPS,nLEO,nMEO,nGEO);

fprintf('\nSelected indices:\n');
disp(selected);

fprintf('\nSelected transmitters:\n');

displayVars = { ...
    'PoolIndex', ...
    'Label', ...
    'Architecture', ...
    'NORAD_CAT_ID', ...
    'ObjectName'};

if ismember( ...
        'SigmaRho_m', ...
        Tselected.Properties.VariableNames)

    displayVars{end+1} = ...
        'SigmaRho_m';

end

disp(Tselected(:,displayVars));

%% ============================================================
% Trace table
%% ============================================================

TtraceBE = table( ...
    traceK(:), ...
    traceRemovedPoolIndex(:), ...
    traceRemovedLabel(:), ...
    traceRemovedArchitecture(:), ...
    traceBound(:), ...
    tracePDOP(:), ...
    traceCandidatesEvaluated(:), ...
    'VariableNames',{ ...
    'K', ...
    'RemovedPoolIndex', ...
    'RemovedLabel', ...
    'RemovedArchitecture', ...
    'PositionalBound_m', ...
    'PDOP', ...
    'CandidatesEvaluated'});
%% ============================================================

resultsDir = fullfile(rootDir,'results');

if ~exist(resultsDir,'dir')
    mkdir(resultsDir);
end

matFile = fullfile( ...
    resultsDir, ...
    'BE_results.mat');

selectedCSV = fullfile( ...
    resultsDir, ...
    'BE_selected.csv');

traceCSV = fullfile( ...
    resultsDir, ...
    'BE_trace.csv');

save(matFile, ...
    'candidateMat', ...
    'linkBudgetMat', ...
    'target_m', ...
    'targetTol', ...
    'Hfull', ...
    'sigmaRho_m', ...
    'VarianceRho_m2', ...
    'CN0_dBHz', ...
    'selected', ...
    'fullBound', ...
    'fullPDOP', ...
    'fullRankH', ...
    'finalBound', ...
    'finalPDOP', ...
    'finalRankH', ...
    'targetReached', ...
    'totalSubsetEvaluations', ...
    'nAcceptedRemovals', ...
    'Tselected', ...
    'TtraceBE');

writetable( ...
    Tselected, ...
    selectedCSV);

writetable( ...
    TtraceBE, ...
    traceCSV);

fprintf('\nSaved:\n%s\n',matFile);
fprintf('%s\n',selectedCSV);
fprintf('%s\n',traceCSV);

fprintf('\n============================================================\n');
fprintf('BACKWARD ELIMINATION COMPLETED\n');
fprintf('============================================================\n');


%% ============================================================
% Local function
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

    %% Geometric PDOP

    A = H.'*H;

    if rcond(A) > 1e-14

        Cgeom = A\eye(4);

        trGeom = ...
            trace(Cgeom(1:3,1:3));

        if trGeom > 0 && ...
           isfinite(trGeom)

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

    trPos = ...
        trace(C(1:3,1:3));

    if trPos > 0 && ...
       isfinite(trPos)

        bound_m = sqrt(trPos);

    end

end