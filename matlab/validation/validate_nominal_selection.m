clear;
clc;

%% ============================================================
%  NOMINAL FTS/BTS VALIDATION
%  ============================================================

%% Repository paths

thisFile  = mfilename('fullpath');
setupDir  = fileparts(thisFile);
matlabDir = fileparts(setupDir);
repoRoot  = fileparts(matlabDir);

addpath(fullfile(repoRoot,'matlab','core'));
addpath(fullfile(repoRoot,'matlab','algorithms'));

%% Positioning target used in the paper
epsilon_m = 0.6;

%% ============================================================
%  LOAD NOMINAL CANDIDATE POOL
%  ============================================================

poolFile = fullfile( ...
    repoRoot, ...
    'data','generated', ...
    'nominal_candidate_pool.mat');

P = load(poolFile);

candidatePool = P.candidatePool;
config        = P.config;

%% ============================================================
%  LOAD NOMINAL LINK BUDGET
%  ============================================================

budgetFile = fullfile( ...
    repoRoot, ...
    'data','generated', ...
    'nominal_link_budget.mat');

B = load(budgetFile);

R = B.R;

%% ============================================================
%  FULL-POOL GEOMETRY
%  ============================================================

H = build_geometry_matrix( ...
    candidatePool, ...
    config.receiver.rECEF_m);

fullMetrics = compute_position_metrics(H,R);

%% ============================================================
%  FTS
%  ============================================================

fprintf('\nRunning FTS...\n');

tic;
fts = forward_transmitter_selection( ...
    H,R,epsilon_m);
ftsTime_s = toc;

%% ============================================================
%  BTS
%  ============================================================

fprintf('Running BTS...\n');

tic;
bts = backward_transmitter_selection( ...
    H,R,epsilon_m);
btsTime_s = toc;

%% ============================================================
%  SELECTED SUBSETS
%  ============================================================

ftsPool = candidatePool(fts.selectedIdx,:);
btsPool = candidatePool(bts.selectedIdx,:);

%% Architecture counts

archNames = ["HAPS","LEO","MEO","GEO"];

ftsCounts = zeros(size(archNames));
btsCounts = zeros(size(archNames));

for a = 1:numel(archNames)

    ftsCounts(a) = sum( ...
        ftsPool.Architecture == archNames(a));

    btsCounts(a) = sum( ...
        btsPool.Architecture == archNames(a));
end

%% ============================================================
%  REPORT
%  ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('NOMINAL TRANSMITTER-SELECTION VALIDATION\n');
fprintf('============================================================\n');

fprintf('\nNO SELECTION\n');
fprintf('------------------------------------------------------------\n');
fprintf('K              : %d\n',height(candidatePool));
fprintf('rank(H)        : %d\n',fullMetrics.rankH);
fprintf('Position bound : %.6f m\n',fullMetrics.alphaCRB_m);

fprintf('\nFTS\n');
fprintf('------------------------------------------------------------\n');
fprintf('Feasible       : %d\n',fts.feasible);
fprintf('K              : %d\n',fts.K);
fprintf('Position bound : %.6f m\n',fts.positionBound_m);
fprintf('Evaluations    : %d\n',fts.evaluations);
fprintf('Runtime        : %.3f s\n',ftsTime_s);

fprintf('Composition    : ');
fprintf('HAPS=%d  LEO=%d  MEO=%d  GEO=%d\n', ...
    ftsCounts(1), ...
    ftsCounts(2), ...
    ftsCounts(3), ...
    ftsCounts(4));

fprintf('\nSelected FTS transmitters:\n');

disp(ftsPool(:,{ ...
    'Label', ...
    'Architecture', ...
    'NORAD', ...
    'Elevation_deg'}));

fprintf('\nBTS\n');
fprintf('------------------------------------------------------------\n');
fprintf('Feasible       : %d\n',bts.feasible);
fprintf('K              : %d\n',bts.K);
fprintf('Position bound : %.6f m\n',bts.positionBound_m);
fprintf('Evaluations    : %d\n',bts.evaluations);
fprintf('Runtime        : %.3f s\n',btsTime_s);

fprintf('Composition    : ');
fprintf('HAPS=%d  LEO=%d  MEO=%d  GEO=%d\n', ...
    btsCounts(1), ...
    btsCounts(2), ...
    btsCounts(3), ...
    btsCounts(4));

fprintf('\nSelected BTS transmitters:\n');

disp(btsPool(:,{ ...
    'Label', ...
    'Architecture', ...
    'NORAD', ...
    'Elevation_deg'}));

%% ============================================================
%  COMPARE FTS AND BTS
%  ============================================================

sameSubset = isequal( ...
    sort(fts.selectedIdx), ...
    sort(bts.selectedIdx));

fprintf('\n------------------------------------------------------------\n');
fprintf('Same FTS/BTS subset : %d\n',sameSubset);

fprintf('Evaluation ratio FTS/BTS : %.2f\n', ...
    fts.evaluations/bts.evaluations);

fprintf('============================================================\n');

%% ============================================================
%  SAVE VALIDATION RESULT
%  ============================================================

outputDir = fullfile( ...
    repoRoot, ...
    'data','generated');

outputFile = fullfile( ...
    outputDir, ...
    'nominal_selection.mat');

save(outputFile, ...
    'fts', ...
    'bts', ...
    'ftsPool', ...
    'btsPool', ...
    'fullMetrics', ...
    'epsilon_m');

fprintf('\nSaved:\n%s\n',outputFile);
fprintf('============================================================\n');