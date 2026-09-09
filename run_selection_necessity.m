clear;
clc;
format long;

%% ============================================================
% SELECTION NECESSITY - FIXED-SIZE SUBSET COMPARISON
%
% Objective:
%   Show why transmitter selection is necessary when only a limited number
%   of measurements can be used.
%
% The nominal full pool contains 26 transmitters. The FA and BE algorithms
% both return a 12-transmitter subset for the 0.6 m positional-bound target.
% This experiment therefore fixes the subset size at K=12 and compares:
%
%   1) the subset selected by FA;
%   2) the subset selected by BE;
%   3) the K transmitters with the largest C/N0;
%   4) the K transmitters with the largest elevation;
%   5) uniformly random K-transmitter subsets.
%
% The full 26-transmitter pool is included only as a reference. Under the
% nominal independent-measurement model, removing valid measurements cannot
% improve the Fisher-information bound relative to the complete pool.
%
% The key quantity is
%
%   P( B_CRB(S_random) <= epsilon ),
%
% which measures how often an arbitrary subset with the SAME number of
% transmitters satisfies the positioning requirement. A low probability
% indicates that reducing the number of measurements is not sufficient by
% itself: the identities of the selected transmitters matter.
%
% Random experiment:
%   N_random = 50000
%   rng(2,'twister')
%
% Outputs:
%   results/selection_necessity/
%       selection_necessity_methods.csv
%       selection_necessity_random_summary.csv
%       selection_necessity_random_subsets.csv
%       selection_necessity_results.mat
%
% The script computes numerical results only. The corresponding paper
% figure should be generated separately by generate_paper_figures.m.
%% ============================================================

rootDir = fileparts(mfilename('fullpath'));

if isempty(rootDir)
    rootDir = pwd;
end

cd(rootDir);
setup_paths;

fprintf('\n============================================================\n');
fprintf('SELECTION NECESSITY - FIXED-SIZE SUBSET COMPARISON\n');
fprintf('============================================================\n');

%% ============================================================
% Input files
%% ============================================================

candidateFile = fullfile( ...
    rootDir,'results','candidate_pool_26.mat');

linkBudgetFile = fullfile( ...
    rootDir,'results','link_budget.mat');

faFile = fullfile( ...
    rootDir,'results','FA_results.mat');

beFile = fullfile( ...
    rootDir,'results','BE_results.mat');

assert(isfile(candidateFile), ...
    'Candidate pool not found:\n%s\nRun build_candidate_pool.m first.', ...
    candidateFile);

assert(isfile(linkBudgetFile), ...
    'Link budget not found:\n%s\nRun build_link_budget.m first.', ...
    linkBudgetFile);

assert(isfile(faFile), ...
    'FA result not found:\n%s\nRun run_FA.m first.', ...
    faFile);

assert(isfile(beFile), ...
    'BE result not found:\n%s\nRun run_BE.m first.', ...
    beFile);

%% ============================================================
% Load nominal geometry and measurement quality
%% ============================================================

S  = load(candidateFile);
LB = load(linkBudgetFile);
FA = load(faFile);
BE = load(beFile);

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
    "Label", ...
    "Architecture", ...
    "CN0_dBHz", ...
    "SigmaRho_m"];

for k = 1:numel(requiredLinkVars)

    assert(isfield(LB,requiredLinkVars(k)), ...
        'Missing link-budget variable: %s', ...
        requiredLinkVars(k));

end

assert(isfield(FA,'selected'), ...
    'FA_results.mat does not contain selected.');

assert(isfield(BE,'selected'), ...
    'BE_results.mat does not contain selected.');

UserPositionECEF = S.UserPositionECEF(:).';
PositionECEF = double(S.PositionECEF);

Label = string(S.Label(:));
Architecture = string(S.Architecture(:));

LabelLB = string(LB.Label(:));
ArchitectureLB = string(LB.Architecture(:));

CN0_dBHz = double(LB.CN0_dBHz(:));
sigmaRho_m = double(LB.SigmaRho_m(:));

Npool = size(PositionECEF,1);

assert(Npool == 26, ...
    'Expected 26 nominal candidates; found %d.',Npool);

assert(numel(Label) == Npool);
assert(numel(Architecture) == Npool);
assert(numel(CN0_dBHz) == Npool);
assert(numel(sigmaRho_m) == Npool);

assert(isequal(Label,LabelLB), ...
    'Candidate-pool and link-budget label order differs.');

assert(isequal(Architecture,ArchitectureLB), ...
    'Candidate-pool and link-budget architecture order differs.');

assert(all(isfinite(CN0_dBHz)));
assert(all(isfinite(sigmaRho_m)));
assert(all(sigmaRho_m > 0));

%% ============================================================
% Geometry matrix
%% ============================================================

delta = PositionECEF - UserPositionECEF;

range_m = vecnorm(delta,2,2);

assert(all(isfinite(range_m)));
assert(all(range_m > 0));

u = delta ./ range_m;

Hfull = [-u,ones(Npool,1)];

assert(rank(Hfull) == 4, ...
    'The nominal candidate pool does not have rank(H)=4.');

%% ============================================================
% Experiment parameters
%% ============================================================

target_m = 0.6;

selectedFA = FA.selected(:).';
selectedBE = BE.selected(:).';

K = numel(selectedFA);

assert(K >= 4 && K <= Npool);

assert(numel(selectedBE) == K, ...
    ['FA and BE must have the same final subset size for this ' ...
     'fixed-K comparison. Obtained FA=%d and BE=%d.'], ...
    K,numel(selectedBE));

N_random = 50000;

rng(2,'twister');

fprintf('Nominal candidates : %d\n',Npool);
fprintf('Fixed subset size  : %d\n',K);
fprintf('Random subsets     : %d\n',N_random);
fprintf('Target             : %.3f m\n',target_m);

%% ============================================================
% Deterministic comparison methods
%% ============================================================

% Highest-C/N0 heuristic.
[~,orderCN0] = sort(CN0_dBHz,'descend');
selectedTopCN0 = sort(orderCN0(1:K)).';

% Highest-elevation heuristic.
if isfield(S,'Elevation_deg')
    elevation_deg = double(S.Elevation_deg(:));
elseif isfield(S,'Tpool') && ...
        istable(S.Tpool) && ...
        ismember('Elevation_deg',S.Tpool.Properties.VariableNames)
    elevation_deg = double(S.Tpool.Elevation_deg(:));
else
    error('Candidate pool does not contain Elevation_deg.');
end

[~,orderElevation] = sort(elevation_deg,'descend');
selectedTopElevation = sort(orderElevation(1:K)).';

fullSet = (1:Npool).';

[fullBound,fullPDOP,fullRank] = ...
    subsetMetrics(fullSet,Hfull,sigmaRho_m);

[faBound,faPDOP,faRank] = ...
    subsetMetrics(selectedFA,Hfull,sigmaRho_m);

[beBound,bePDOP,beRank] = ...
    subsetMetrics(selectedBE,Hfull,sigmaRho_m);

[topCN0Bound,topCN0PDOP,topCN0Rank] = ...
    subsetMetrics(selectedTopCN0,Hfull,sigmaRho_m);

[topElevationBound,topElevationPDOP,topElevationRank] = ...
    subsetMetrics(selectedTopElevation,Hfull,sigmaRho_m);

%% ============================================================
% Random fixed-K subsets
%% ============================================================

RandomIndices = zeros(N_random,K);

RandomBound_m = inf(N_random,1);
RandomPDOP = inf(N_random,1);
RandomRankH = zeros(N_random,1);
RandomTargetMet = false(N_random,1);

RandomHAPS = zeros(N_random,1);
RandomLEO  = zeros(N_random,1);
RandomMEO  = zeros(N_random,1);
RandomGEO  = zeros(N_random,1);

fprintf('\nSampling uniformly random K=%d subsets...\n',K);

for trial = 1:N_random

    idx = sort(randperm(Npool,K));

    RandomIndices(trial,:) = idx;

    [b,p,rk] = ...
        subsetMetrics(idx,Hfull,sigmaRho_m);

    RandomBound_m(trial) = b;
    RandomPDOP(trial) = p;
    RandomRankH(trial) = rk;

    RandomTargetMet(trial) = ...
        rk == 4 && b <= target_m;

    RandomHAPS(trial) = ...
        sum(Architecture(idx)=="HAPS");

    RandomLEO(trial) = ...
        sum(Architecture(idx)=="LEO");

    RandomMEO(trial) = ...
        sum(Architecture(idx)=="MEO");

    RandomGEO(trial) = ...
        sum(Architecture(idx)=="GEO");

    if mod(trial,10000)==0 || trial==N_random

        fprintf('  %6d / %6d subsets evaluated\n', ...
            trial,N_random);

    end
end

%% ============================================================
% Random-subset statistics
%% ============================================================

finiteRandom = ...
    isfinite(RandomBound_m) & ...
    RandomRankH==4;

finiteBounds = ...
    RandomBound_m(finiteRandom);

if isempty(finiteBounds)
    error('No full-rank random subset was obtained.');
end

FullRankProb = mean(finiteRandom);

TargetMetProb = mean(RandomTargetMet);

RandomMeanBound_m = mean(finiteBounds);
RandomMedianBound_m = median(finiteBounds);

RandomP05Bound_m = empiricalPercentile(finiteBounds,5);
RandomP25Bound_m = empiricalPercentile(finiteBounds,25);
RandomP75Bound_m = empiricalPercentile(finiteBounds,75);
RandomP95Bound_m = empiricalPercentile(finiteBounds,95);

RandomMinBound_m = min(finiteBounds);
RandomMaxBound_m = max(finiteBounds);

ProbRandomAtOrBelowFA = ...
    mean(finiteRandom & RandomBound_m <= faBound);

ProbRandomAtOrBelowBE = ...
    mean(finiteRandom & RandomBound_m <= beBound);

FAImprovementVsRandomMedian_pct = ...
    100*(RandomMedianBound_m-faBound)/RandomMedianBound_m;

BEImprovementVsRandomMedian_pct = ...
    100*(RandomMedianBound_m-beBound)/RandomMedianBound_m;

%% ============================================================
% Method comparison table
%% ============================================================

Method = [ ...
    "Full pool"
    "FA selected"
    "BE selected"
    "Top C/N0"
    "Top elevation"];

Kmethod = [ ...
    Npool
    K
    K
    K
    K];

PositionalBound_m = [ ...
    fullBound
    faBound
    beBound
    topCN0Bound
    topElevationBound];

PDOP = [ ...
    fullPDOP
    faPDOP
    bePDOP
    topCN0PDOP
    topElevationPDOP];

RankH = [ ...
    fullRank
    faRank
    beRank
    topCN0Rank
    topElevationRank];

TargetReached = ...
    RankH==4 & ...
    PositionalBound_m<=target_m;

MethodSets = { ...
    fullSet
    selectedFA(:)
    selectedBE(:)
    selectedTopCN0(:)
    selectedTopElevation(:)};

HAPS = zeros(numel(Method),1);
LEO  = zeros(numel(Method),1);
MEO  = zeros(numel(Method),1);
GEO  = zeros(numel(Method),1);

SelectedIndices = strings(numel(Method),1);

for m = 1:numel(Method)

    idx = MethodSets{m};

    HAPS(m) = sum(Architecture(idx)=="HAPS");
    LEO(m)  = sum(Architecture(idx)=="LEO");
    MEO(m)  = sum(Architecture(idx)=="MEO");
    GEO(m)  = sum(Architecture(idx)=="GEO");

    SelectedIndices(m) = ...
        strjoin(string(idx(:).'),";");

end

Tmethods = table( ...
    Method, ...
    Kmethod, ...
    PositionalBound_m, ...
    PDOP, ...
    RankH, ...
    TargetReached, ...
    HAPS, ...
    LEO, ...
    MEO, ...
    GEO, ...
    SelectedIndices, ...
    'VariableNames',{ ...
    'Method', ...
    'K', ...
    'PositionalBound_m', ...
    'PDOP', ...
    'RankH', ...
    'TargetReached', ...
    'HAPS', ...
    'LEO', ...
    'MEO', ...
    'GEO', ...
    'SelectedIndices'});

%% ============================================================
% Random-subset summary
%% ============================================================

TrandomSummary = table( ...
    N_random, ...
    K, ...
    target_m, ...
    FullRankProb, ...
    TargetMetProb, ...
    RandomMeanBound_m, ...
    RandomMedianBound_m, ...
    RandomP05Bound_m, ...
    RandomP25Bound_m, ...
    RandomP75Bound_m, ...
    RandomP95Bound_m, ...
    RandomMinBound_m, ...
    RandomMaxBound_m, ...
    ProbRandomAtOrBelowFA, ...
    ProbRandomAtOrBelowBE, ...
    FAImprovementVsRandomMedian_pct, ...
    BEImprovementVsRandomMedian_pct);

%% ============================================================
% Random-subset trial table
%% ============================================================

TrialID = (1:N_random).';

Trandom = table( ...
    TrialID, ...
    RandomBound_m, ...
    RandomPDOP, ...
    RandomRankH, ...
    RandomTargetMet, ...
    RandomHAPS, ...
    RandomLEO, ...
    RandomMEO, ...
    RandomGEO, ...
    'VariableNames',{ ...
    'TrialID', ...
    'PositionalBound_m', ...
    'PDOP', ...
    'RankH', ...
    'TargetReached', ...
    'HAPS', ...
    'LEO', ...
    'MEO', ...
    'GEO'});

%% ============================================================
% Console summary
%% ============================================================

fprintf('\n============================================================\n');
fprintf('DETERMINISTIC METHODS\n');
fprintf('============================================================\n');

disp(Tmethods(:,{ ...
    'Method', ...
    'K', ...
    'PositionalBound_m', ...
    'PDOP', ...
    'TargetReached', ...
    'HAPS', ...
    'LEO', ...
    'MEO', ...
    'GEO'}));

fprintf('\n============================================================\n');
fprintf('RANDOM K=%d SUBSETS\n',K);
fprintf('============================================================\n');

fprintf('Full-rank probability             : %.4f %%\n', ...
    100*FullRankProb);

fprintf('P(random bound <= %.3f m)         : %.4f %%\n', ...
    target_m,100*TargetMetProb);

fprintf('Random median bound               : %.12f m\n', ...
    RandomMedianBound_m);

fprintf('Random 5th / 95th percentile      : %.12f / %.12f m\n', ...
    RandomP05Bound_m,RandomP95Bound_m);

fprintf('Best random sampled bound         : %.12f m\n', ...
    RandomMinBound_m);

fprintf('Worst finite random sampled bound : %.12f m\n', ...
    RandomMaxBound_m);

fprintf('FA selected bound                 : %.12f m\n', ...
    faBound);

fprintf('BE selected bound                 : %.12f m\n', ...
    beBound);

fprintf('P(random bound <= FA bound)       : %.4f %%\n', ...
    100*ProbRandomAtOrBelowFA);

fprintf('P(random bound <= BE bound)       : %.4f %%\n', ...
    100*ProbRandomAtOrBelowBE);

fprintf('FA reduction vs random median     : %.4f %%\n', ...
    FAImprovementVsRandomMedian_pct);

fprintf('BE reduction vs random median     : %.4f %%\n', ...
    BEImprovementVsRandomMedian_pct);

%% ============================================================
% Save
%% ============================================================

outDir = fullfile( ...
    rootDir, ...
    'results', ...
    'selection_necessity');

if ~exist(outDir,'dir')
    mkdir(outDir);
end

methodsFile = fullfile( ...
    outDir, ...
    'selection_necessity_methods.csv');

randomSummaryFile = fullfile( ...
    outDir, ...
    'selection_necessity_random_summary.csv');

randomTrialsFile = fullfile( ...
    outDir, ...
    'selection_necessity_random_subsets.csv');

matFile = fullfile( ...
    outDir, ...
    'selection_necessity_results.mat');

writetable(Tmethods,methodsFile);
writetable(TrandomSummary,randomSummaryFile);
writetable(Trandom,randomTrialsFile);

save(matFile, ...
    'Npool', ...
    'K', ...
    'N_random', ...
    'target_m', ...
    'Hfull', ...
    'sigmaRho_m', ...
    'CN0_dBHz', ...
    'elevation_deg', ...
    'Label', ...
    'Architecture', ...
    'selectedFA', ...
    'selectedBE', ...
    'selectedTopCN0', ...
    'selectedTopElevation', ...
    'RandomIndices', ...
    'RandomBound_m', ...
    'RandomPDOP', ...
    'RandomRankH', ...
    'RandomTargetMet', ...
    'Tmethods', ...
    'TrandomSummary', ...
    'Trandom');

fprintf('\nSaved:\n');
fprintf('  %s\n',methodsFile);
fprintf('  %s\n',randomSummaryFile);
fprintf('  %s\n',randomTrialsFile);
fprintf('  %s\n',matFile);

fprintf('\n============================================================\n');
fprintf('SELECTION-NECESSITY EXPERIMENT COMPLETED\n');
fprintf('============================================================\n');


%% ============================================================
% Local weighted positional metric
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


%% ============================================================
% Empirical percentile without Statistics Toolbox
%% ============================================================

function p = empiricalPercentile(x,percent)

    x = x(isfinite(x));
    x = sort(x(:));

    if isempty(x)

        p = NaN;
        return;
    end

    if numel(x)==1

        p = x;
        return;
    end

    pos = 1 + ...
        (numel(x)-1)*(percent/100);

    lo = floor(pos);
    hi = ceil(pos);

    if lo == hi

        p = x(lo);

    else

        w = pos-lo;

        p = ...
            (1-w)*x(lo) ...
            + w*x(hi);

    end
end
