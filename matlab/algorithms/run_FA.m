clear;
clc;
format long;

%% ============================================================
% FORWARD ADDING (FA)
%
% Canonical location:
%   matlab/algorithms/run_FA.m
%
% Inputs:
%   results/candidate_pool_26.mat
%   results/link_budget.mat
%
% Outputs:
%   results/FA_results.mat
%   results/FA_selected.csv
%   results/FA_trace.csv
% ============================================================

repoRoot = locateRepoRoot(mfilename('fullpath'));
addpath(genpath(fullfile(repoRoot,'matlab')),'-begin');

setupFile = fullfile(repoRoot,'matlab','setup','setup_paths.m');
if isfile(setupFile)
    run(setupFile);
end

candidateMat = fullfile(repoRoot,'results','candidate_pool_26.mat');
linkBudgetMat = fullfile(repoRoot,'results','link_budget.mat');

assert(isfile(candidateMat), ...
    'Candidate pool not found:\n%s',candidateMat);

assert(isfile(linkBudgetMat), ...
    'Link budget not found:\n%s',linkBudgetMat);

S = load(candidateMat);
LB = load(linkBudgetMat);

requiredCandidateVars = [ ...
    "UserPositionECEF","PositionECEF","Label","Architecture", ...
    "NORAD","ObjectName","Elevation_deg","SlantRange_m"];

for k = 1:numel(requiredCandidateVars)
    assert(isfield(S,requiredCandidateVars(k)), ...
        'Missing candidate-pool variable: %s',requiredCandidateVars(k));
end

requiredLinkVars = [ ...
    "Label","Architecture","CN0_dBHz","SigmaRho_m", ...
    "VarianceRho_m2","beta","Tcoh","etaPos"];

for k = 1:numel(requiredLinkVars)
    assert(isfield(LB,requiredLinkVars(k)), ...
        'Missing link-budget variable: %s',requiredLinkVars(k));
end

UserPositionECEF = double(S.UserPositionECEF(:).');
PositionECEF = double(S.PositionECEF);

Label = string(S.Label(:));
Architecture = upper(string(S.Architecture(:)));
NORAD = string(S.NORAD(:));
ObjectName = string(S.ObjectName(:));
Elevation_deg = double(S.Elevation_deg(:));
SlantRange_m = double(S.SlantRange_m(:));

LabelLB = string(LB.Label(:));
ArchitectureLB = upper(string(LB.Architecture(:)));

assert(isequal(Label,LabelLB), ...
    'Candidate/link-budget label order mismatch.');

assert(isequal(Architecture,ArchitectureLB), ...
    'Candidate/link-budget architecture order mismatch.');

CN0_dBHz = double(LB.CN0_dBHz(:));
sigmaRho_m = double(LB.SigmaRho_m(:));
VarianceRho_m2 = double(LB.VarianceRho_m2(:));

beta = double(LB.beta);
Tcoh = double(LB.Tcoh);
etaPos = double(LB.etaPos);

Npool = size(PositionECEF,1);

assert(Npool==26,'Expected 26 candidates.');
assert(all(isfinite(sigmaRho_m)) && all(sigmaRho_m>0));

delta = PositionECEF-UserPositionECEF;
range_m = vecnorm(delta,2,2);
u = delta./range_m;
Hfull = [-u,ones(Npool,1)];

assert(rank(Hfull)==4,'Full candidate set must have rank(H)=4.');

target_m = 0.6;
targetTol = 1e-12;

allIdx = 1:Npool;

[fullBound,fullPDOP,fullRankH] = ...
    subsetMetrics(allIdx,Hfull,sigmaRho_m);

fprintf('\n============================================================\n');
fprintf('FORWARD ADDING (FA)\n');
fprintf('============================================================\n');
fprintf('Candidates        : %d\n',Npool);
fprintf('Full bound        : %.12f m\n',fullBound);
fprintf('Full PDOP         : %.12f\n',fullPDOP);
fprintf('rank(H)           : %d\n',fullRankH);
fprintf('Target            : %.3f m\n',target_m);

%% Best initial 4-transmitter subset

comb4 = nchoosek(1:Npool,4);

best4 = [];
best4Bound = inf;
best4PDOP = inf;

evaluationCount = 0;

for k = 1:size(comb4,1)

    idx = comb4(k,:);

    [b,pdop,rankH] = ...
        subsetMetrics(idx,Hfull,sigmaRho_m);

    evaluationCount = evaluationCount+1;

    if rankH==4 && b<best4Bound
        best4 = idx;
        best4Bound = b;
        best4PDOP = pdop;
    end
end

assert(~isempty(best4),'No valid rank-4 initial subset was found.');

selected = best4(:).';
currentBound = best4Bound;

traceK = numel(selected);
traceAddedIndex = NaN;
traceAddedLabel = "";
traceAddedArchitecture = "";
traceBound = currentBound;

while currentBound>target_m+targetTol && numel(selected)<Npool

    remaining = setdiff(allIdx,selected,'stable');

    bestCandidate = NaN;
    bestCandidateBound = inf;

    for k = 1:numel(remaining)

        candidate = remaining(k);
        trial = [selected,candidate];

        [b,~,rankH] = ...
            subsetMetrics(trial,Hfull,sigmaRho_m);

        evaluationCount = evaluationCount+1;

        if rankH==4 && b<bestCandidateBound
            bestCandidate = candidate;
            bestCandidateBound = b;
        end
    end

    if isnan(bestCandidate)
        break;
    end

    selected(end+1) = bestCandidate; %#ok<SAGROW>
    currentBound = bestCandidateBound;

    traceK(end+1,1) = numel(selected); %#ok<SAGROW>
    traceAddedIndex(end+1,1) = bestCandidate; %#ok<SAGROW>
    traceAddedLabel(end+1,1) = Label(bestCandidate); %#ok<SAGROW>
    traceAddedArchitecture(end+1,1) = Architecture(bestCandidate); %#ok<SAGROW>
    traceBound(end+1,1) = currentBound; %#ok<SAGROW>
end

[finalBound,finalPDOP,finalRankH] = ...
    subsetMetrics(selected,Hfull,sigmaRho_m);

targetReached = finalBound<=target_m+targetTol;

PoolIndex = (1:Npool).';

Tquality = table( ...
    PoolIndex,Label,Architecture,NORAD,ObjectName, ...
    Elevation_deg,SlantRange_m/1e3,CN0_dBHz,sigmaRho_m,VarianceRho_m2, ...
    'VariableNames',{ ...
    'PoolIndex','Label','Architecture','NORAD_CAT_ID','ObjectName', ...
    'Elevation_deg','SlantRange_km','CN0_dBHz','SigmaRho_m','VarianceRho_m2'});

Tselected = Tquality(selected,:);
Tselected.SelectionOrder = (1:height(Tselected)).';
Tselected = movevars(Tselected,'SelectionOrder','Before','PoolIndex');

Ttrace = table( ...
    traceK(:),traceAddedIndex(:),traceAddedLabel(:), ...
    traceAddedArchitecture(:),traceBound(:), ...
    'VariableNames',{ ...
    'K','AddedPoolIndex','AddedLabel','AddedArchitecture','PositionalBound_m'});

fprintf('\nFA RESULT\n');
fprintf('Final K           : %d\n',numel(selected));
fprintf('Positional bound  : %.12f m\n',finalBound);
fprintf('PDOP              : %.12f\n',finalPDOP);
fprintf('rank(H)           : %d\n',finalRankH);
fprintf('Target reached    : %s\n',string(targetReached));
fprintf('Evaluations       : %d\n',evaluationCount);
fprintf('Composition       : HAPS=%d | LEO=%d | MEO=%d | GEO=%d\n', ...
    sum(Architecture(selected)=="HAPS"), ...
    sum(Architecture(selected)=="LEO"), ...
    sum(Architecture(selected)=="MEO"), ...
    sum(Architecture(selected)=="GEO"));
fprintf('Labels            : %s\n',strjoin(Label(selected)," "));

resultsDir = fullfile(repoRoot,'results');

if ~isfolder(resultsDir)
    mkdir(resultsDir);
end

matFile = fullfile(resultsDir,'FA_results.mat');
selectedCSV = fullfile(resultsDir,'FA_selected.csv');
traceCSV = fullfile(resultsDir,'FA_trace.csv');

save(matFile, ...
    'candidateMat','linkBudgetMat','target_m','Hfull', ...
    'sigmaRho_m','VarianceRho_m2','CN0_dBHz', ...
    'best4','best4Bound','best4PDOP', ...
    'selected','finalBound','finalPDOP','finalRankH', ...
    'targetReached','evaluationCount','Tselected','Ttrace', ...
    'fullBound','fullPDOP','fullRankH','beta','Tcoh','etaPos');

writetable(Tselected,selectedCSV);
writetable(Ttrace,traceCSV);

fprintf('\nSaved:\n  %s\n',matFile);
fprintf('  %s\n',selectedCSV);
fprintf('  %s\n',traceCSV);
fprintf('============================================================\n');

%% ============================================================
% Local functions
% =============================================================

function repoRoot = locateRepoRoot(thisFile)

    repoRoot = fileparts(thisFile);

    while ~isfolder(fullfile(repoRoot,'.git'))

        parentDir = fileparts(repoRoot);

        if strcmp(parentDir,repoRoot)
            error('Could not locate Git repository root from:\n%s',thisFile);
        end

        repoRoot = parentDir;
    end
end

function [bound_m,PDOP,rankH] = subsetMetrics(idx,Hfull,sigmaRho_m)

    idx = idx(:);
    H = Hfull(idx,:);
    sigma = sigmaRho_m(idx);

    rankH = rank(H);
    bound_m = inf;
    PDOP = inf;

    if rankH<4
        return;
    end

    A = H.'*H;

    if rcond(A)>1e-14
        Cgeom = A\eye(4);
        PDOP = sqrt(trace(Cgeom(1:3,1:3)));
    end

    R = diag(sigma.^2);
    J = H.'*(R\H);

    if rcond(J)<=1e-14
        return;
    end

    C = J\eye(4);
    bound_m = sqrt(trace(C(1:3,1:3)));
end
