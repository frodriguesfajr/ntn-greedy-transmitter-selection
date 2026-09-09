clear;
clc;
format long;

%% ============================================================
% HAPS AVAILABILITY RESTRICTIONS
%
% Cases:
%   maxHAPS = 4
%   maxHAPS = 2
%   maxHAPS = 0
%
% Same algorithm and target used in the paper.
% Only change:
%   orbital propagation -> Vallado SGP4 / WGS-72
%% ============================================================

rootDir = fileparts(mfilename('fullpath'));

if isempty(rootDir)
    rootDir = pwd;
end

cd(rootDir);
setup_paths;

fprintf('\n============================================================\n');
fprintf('HAPS RESTRICTIONS - GEOMETRY\n');
fprintf('============================================================\n');

%% ============================================================
% Inputs
%% ============================================================

candidateMat = fullfile( ...
    rootDir,'results','candidate_pool_26.mat');

linkBudgetMat = fullfile( ...
    rootDir,'results','link_budget.mat');

assert(isfile(candidateMat), ...
    'Candidate pool not found:\n%s',candidateMat);

assert(isfile(linkBudgetMat), ...
    'Link budget not found:\n%s',linkBudgetMat);

%% ============================================================
% Candidate pool
%% ============================================================

S = load(candidateMat);

UserPositionECEF = S.UserPositionECEF(:).';
PositionECEF     = S.PositionECEF;

Label        = string(S.Label(:));
Architecture = string(S.Architecture(:));
NORAD        = string(S.NORAD(:));
ObjectName   = string(S.ObjectName(:));

Elevation_deg = S.Elevation_deg(:);
SlantRange_m  = S.SlantRange_m(:);

Npool = size(PositionECEF,1);

assert(Npool == 26);

%% ============================================================
% Link budget
%% ============================================================

LB = load(linkBudgetMat);

LabelLB        = string(LB.Label(:));
ArchitectureLB = string(LB.Architecture(:));

assert(isequal(Label,LabelLB), ...
    'Candidate/link-budget labels do not match.');

assert(isequal(Architecture,ArchitectureLB), ...
    'Candidate/link-budget architectures do not match.');

CN0_dBHz       = LB.CN0_dBHz(:);
sigmaRho_m     = LB.SigmaRho_m(:);
VarianceRho_m2 = LB.VarianceRho_m2(:);

assert(all(isfinite(sigmaRho_m)));
assert(all(sigmaRho_m > 0));

%% Use prepared corrected link-budget table when available

if isfield(LB,'Tlink') && ...
   istable(LB.Tlink) && ...
   height(LB.Tlink) == Npool

    Tquality = LB.Tlink;

else

    PoolIndex = (1:Npool).';

    Tquality = table( ...
        PoolIndex,Label,Architecture,NORAD,ObjectName, ...
        Elevation_deg,SlantRange_m/1e3,CN0_dBHz, ...
        sigmaRho_m,VarianceRho_m2, ...
        'VariableNames',{ ...
        'PoolIndex','Label','Architecture','NORAD_CAT_ID', ...
        'ObjectName','Elevation_deg','SlantRange_km', ...
        'CN0_dBHz','SigmaRho_m','VarianceRho_m2'});

end

%% ============================================================
% Geometry
%% ============================================================

r = PositionECEF - UserPositionECEF;

d = sqrt(sum(r.^2,2));

if any(~isfinite(d)) || any(d <= eps)
    error('Invalid transmitter-user distance.');
end

u = r ./ d;

Hfull = [-u,ones(Npool,1)];

%% ============================================================
% Cases
%% ============================================================

target_m = 0.6;

maxHAPS_cases = [4 2 0];

caseNames = [ ...
    "Reference_max4"
    "Limited_max2"
    "No_HAPS"];

nCases = numel(maxHAPS_cases);

%% Outputs

MaxHAPS = maxHAPS_cases(:);

CandidatePoolCount = ...
    repmat(Npool,nCases,1);

MaxFeasibleSetSize = nan(nCases,1);

BestMaxFeasibleBound_m = nan(nCases,1);
BestMaxFeasiblePDOP    = nan(nCases,1);

Best4Bound_m = nan(nCases,1);

TargetReached = false(nCases,1);

Kfinal      = nan(nCases,1);
FinalBound_m = nan(nCases,1);
FinalPDOP    = nan(nCases,1);

CountHAPS = zeros(nCases,1);
CountLEO  = zeros(nCases,1);
CountMEO  = zeros(nCases,1);
CountGEO  = zeros(nCases,1);

InitialEvaluations = zeros(nCases,1);
FAEvaluations      = zeros(nCases,1);

selectedPerCase   = cell(nCases,1);
bestMaxSetPerCase = cell(nCases,1);
tracePerCase      = cell(nCases,1);

%% ============================================================
% Case loop
%% ============================================================

for caseIdx = 1:nCases

    maxHAPS = maxHAPS_cases(caseIdx);
    caseName = caseNames(caseIdx);

    fprintf('\n\n############################################################\n');
    fprintf('CASE %s | maxHAPS = %d\n',caseName,maxHAPS);
    fprintf('############################################################\n');

    %% --------------------------------------------------------
    % Largest feasible set under HAPS restriction
    %% --------------------------------------------------------

    idxNonHAPS = find(Architecture ~= "HAPS");
    idxHAPS    = find(Architecture == "HAPS");

    nHAPSuse = min(maxHAPS,numel(idxHAPS));

    MaxFeasibleSetSize(caseIdx) = ...
        numel(idxNonHAPS) + nHAPSuse;

    if nHAPSuse == 0

        bestMaxSet = idxNonHAPS(:).';

        [bestMaxBound,bestMaxPDOP,bestMaxRank] = ...
            subsetMetrics(bestMaxSet,Hfull,sigmaRho_m);

    elseif nHAPSuse == numel(idxHAPS)

        bestMaxSet = ...
            [idxNonHAPS(:).',idxHAPS(:).'];

        [bestMaxBound,bestMaxPDOP,bestMaxRank] = ...
            subsetMetrics(bestMaxSet,Hfull,sigmaRho_m);

    else

        combHAPS = nchoosek(idxHAPS,nHAPSuse);

        bestMaxBound = inf;
        bestMaxPDOP  = inf;
        bestMaxRank  = 0;
        bestMaxSet   = [];

        for kk = 1:size(combHAPS,1)

            idxTrial = ...
                [idxNonHAPS(:).',combHAPS(kk,:)];

            [b,p,rk] = ...
                subsetMetrics(idxTrial,Hfull,sigmaRho_m);

            if rk == 4 && b < bestMaxBound

                bestMaxBound = b;
                bestMaxPDOP  = p;
                bestMaxRank  = rk;
                bestMaxSet   = idxTrial;

            end
        end
    end

    BestMaxFeasibleBound_m(caseIdx) = bestMaxBound;
    BestMaxFeasiblePDOP(caseIdx)    = bestMaxPDOP;

    bestMaxSetPerCase{caseIdx} = bestMaxSet;

    fprintf('\nLargest feasible set\n');
    fprintf('Maximum feasible N : %d\n', ...
        MaxFeasibleSetSize(caseIdx));

    fprintf('HAPS                : %d\n', ...
        sum(Architecture(bestMaxSet)=="HAPS"));

    fprintf('Bound               : %.12f m\n', ...
        bestMaxBound);

    fprintf('PDOP                : %.12f\n', ...
        bestMaxPDOP);

    fprintf('rank(H)             : %d\n', ...
        bestMaxRank);

    fprintf('Target achievable   : %s\n', ...
        string(bestMaxBound <= target_m));

    %% --------------------------------------------------------
    % Allowed candidate set
    %% --------------------------------------------------------

    if maxHAPS == 0
        allowed = idxNonHAPS(:).';
    else
        allowed = 1:Npool;
    end

    %% --------------------------------------------------------
    % Best initial K=4 subset
    %% --------------------------------------------------------

    comb4 = nchoosek(allowed,4);

    best4Bound = inf;
    best4 = [];

    nEvaluatedFeasible = 0;

    for k = 1:size(comb4,1)

        idx = comb4(k,:);

        if sum(Architecture(idx)=="HAPS") > maxHAPS
            continue;
        end

        nEvaluatedFeasible = ...
            nEvaluatedFeasible + 1;

        [b,~,rankH] = ...
            subsetMetrics(idx,Hfull,sigmaRho_m);

        if rankH == 4 && b < best4Bound

            best4Bound = b;
            best4 = idx;

        end
    end

    InitialEvaluations(caseIdx) = ...
        nEvaluatedFeasible;

    if isempty(best4)

        warning( ...
            'No valid initial K=4 subset for maxHAPS=%d.', ...
            maxHAPS);

        continue;
    end

    Best4Bound_m(caseIdx) = best4Bound;

    fprintf('\nBest initial K=4\n');
    fprintf('Feasible subsets evaluated : %d\n', ...
        nEvaluatedFeasible);

    fprintf('Bound                      : %.12f m\n', ...
        best4Bound);

    fprintf('Indices                    : ');
    fprintf('%d ',best4);
    fprintf('\n');

    fprintf('Labels                     : ');
    fprintf('%s ',Label(best4));
    fprintf('\n');

    %% --------------------------------------------------------
    % Forward Adding
    %% --------------------------------------------------------

    selected = best4(:).';
    currentBound = best4Bound;

    traceK = numel(selected);
    traceAddedIndex = NaN;
    traceAddedLabel = "";
    traceAddedArchitecture = "";
    traceHAPS = sum(Architecture(selected)=="HAPS");
    traceBound = currentBound;

    nFAEvaluations = 0;

    fprintf('\nForward Adding\n');

    fprintf( ...
        'K=%d | HAPS=%d | bound=%.12f m | initialization\n', ...
        numel(selected),traceHAPS,currentBound);

    while currentBound > target_m && ...
          numel(selected) < numel(allowed)

        remaining = ...
            setdiff(allowed,selected,'stable');

        currentHAPS = ...
            sum(Architecture(selected)=="HAPS");

        if currentHAPS >= maxHAPS

            remaining = ...
                remaining(Architecture(remaining) ~= "HAPS");

        end

        if isempty(remaining)
            break;
        end

        bestCandidate = NaN;
        bestCandidateBound = inf;

        for k = 1:numel(remaining)

            candidate = remaining(k);

            idxTrial = [selected,candidate];

            if sum(Architecture(idxTrial)=="HAPS") > maxHAPS
                continue;
            end

            [b,~,rankH] = ...
                subsetMetrics(idxTrial,Hfull,sigmaRho_m);

            nFAEvaluations = ...
                nFAEvaluations + 1;

            if rankH == 4 && ...
               b < bestCandidateBound

                bestCandidateBound = b;
                bestCandidate = candidate;

            end
        end

        if isnan(bestCandidate)
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

        traceHAPS(end+1,1) = ...
            sum(Architecture(selected)=="HAPS"); %#ok<SAGROW>

        traceBound(end+1,1) = ...
            currentBound; %#ok<SAGROW>

        fprintf( ...
            'K=%d | HAPS=%d | add %-6s (%s) | bound=%.12f m\n', ...
            numel(selected), ...
            sum(Architecture(selected)=="HAPS"), ...
            Label(bestCandidate), ...
            Architecture(bestCandidate), ...
            currentBound);

    end

    FAEvaluations(caseIdx) = nFAEvaluations;

    %% --------------------------------------------------------
    % Final result
    %% --------------------------------------------------------

    [finalBound,finalPDOP,finalRank] = ...
        subsetMetrics(selected,Hfull,sigmaRho_m);

    reached = finalBound <= target_m;

    TargetReached(caseIdx) = reached;

    Kfinal(caseIdx)       = numel(selected);
    FinalBound_m(caseIdx) = finalBound;
    FinalPDOP(caseIdx)    = finalPDOP;

    CountHAPS(caseIdx) = ...
        sum(Architecture(selected)=="HAPS");

    CountLEO(caseIdx) = ...
        sum(Architecture(selected)=="LEO");

    CountMEO(caseIdx) = ...
        sum(Architecture(selected)=="MEO");

    CountGEO(caseIdx) = ...
        sum(Architecture(selected)=="GEO");

    selectedPerCase{caseIdx} = selected;

    Ttrace = table( ...
        traceK(:), ...
        traceAddedIndex(:), ...
        traceAddedLabel(:), ...
        traceAddedArchitecture(:), ...
        traceHAPS(:), ...
        traceBound(:), ...
        'VariableNames',{ ...
        'K', ...
        'AddedPoolIndex', ...
        'AddedLabel', ...
        'AddedArchitecture', ...
        'HAPSCount', ...
        'PositionalBound_m'});

    tracePerCase{caseIdx} = Ttrace;

    fprintf('\nRESULT maxHAPS=%d\n',maxHAPS);

    fprintf('Final K       : %d\n', ...
        numel(selected));

    fprintf('Final bound   : %.12f m\n', ...
        finalBound);

    fprintf('PDOP          : %.12f\n', ...
        finalPDOP);

    fprintf('rank(H)       : %d\n', ...
        finalRank);

    fprintf('Target met    : %s\n', ...
        string(reached));

    fprintf( ...
        'Composition   : HAPS=%d | LEO=%d | MEO=%d | GEO=%d\n', ...
        CountHAPS(caseIdx), ...
        CountLEO(caseIdx), ...
        CountMEO(caseIdx), ...
        CountGEO(caseIdx));

    fprintf('Selected      : ');
    fprintf('%d ',selected);
    fprintf('\n');

end

%% ============================================================
% Corrected comparison table
%% ============================================================

Tcomparison = table( ...
    caseNames, ...
    MaxHAPS, ...
    CandidatePoolCount, ...
    MaxFeasibleSetSize, ...
    BestMaxFeasibleBound_m, ...
    BestMaxFeasiblePDOP, ...
    Best4Bound_m, ...
    TargetReached, ...
    Kfinal, ...
    FinalBound_m, ...
    FinalPDOP, ...
    CountHAPS, ...
    CountLEO, ...
    CountMEO, ...
    CountGEO, ...
    InitialEvaluations, ...
    FAEvaluations, ...
    'VariableNames',{ ...
    'Case', ...
    'MaxHAPS', ...
    'CandidatePoolCount', ...
    'MaxFeasibleSetSize', ...
    'BestMaxFeasibleBound_m', ...
    'BestMaxFeasiblePDOP', ...
    'Best4Bound_m', ...
    'TargetReached', ...
    'Kfinal', ...
    'FinalBound_m', ...
    'FinalPDOP', ...
    'HAPS', ...
    'LEO', ...
    'MEO', ...
    'GEO', ...
    'InitialEvaluations', ...
    'FAEvaluations'});

fprintf('\n============================================================\n');
fprintf('HAPS-RESTRICTION SUMMARY\n');
fprintf('============================================================\n');

disp(Tcomparison);

%% ============================================================
% Compare with old paper results
%% ============================================================
% Save
%% ============================================================

resultsDir = fullfile(rootDir,'results');

matFile = fullfile( ...
    resultsDir, ...
    'HAPS_restrictions_results.mat');

csvFile = fullfile( ...
    resultsDir, ...
    'HAPS_restrictions_summary.csv');

save(matFile, ...
    'Tcomparison', ...
    'selectedPerCase', ...
    'bestMaxSetPerCase', ...
    'tracePerCase', ...
    'target_m', ...
    'maxHAPS_cases', ...
    'caseNames', ...
    'Hfull', ...
    'sigmaRho_m');

writetable(Tcomparison,csvFile);

fprintf('\nSaved:\n%s\n%s\n', ...
    matFile,csvFile);

fprintf('\n============================================================\n');
fprintf('HAPS RESTRICTIONS COMPLETED\n');
fprintf('============================================================\n');


%% ============================================================
% Local subset metric
%% ============================================================

function [bound_m,PDOP,rankH] = ...
    subsetMetrics(idx,Hfull,sigmaRho_m)

    idx = idx(:);

    H = Hfull(idx,:);
    sigma = sigmaRho_m(idx);

    rankH = rank(H);

    bound_m = inf;
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
        bound_m = sqrt(trPos);
    end

end