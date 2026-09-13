clear;
clc;

%% ============================================================
%  VALIDATE NOMINAL FULL-POOL POSITIONING PERFORMANCE
%  ============================================================

%% Repository paths
thisFile  = mfilename('fullpath');
setupDir  = fileparts(thisFile);
matlabDir = fileparts(setupDir);
repoRoot  = fileparts(matlabDir);

addpath(fullfile(repoRoot,'matlab','core'));

%% Candidate pool
poolFile = fullfile( ...
    repoRoot, ...
    'data','generated', ...
    'nominal_candidate_pool.mat');

assert(isfile(poolFile), ...
    'Nominal candidate-pool file not found.');

P = load(poolFile);

candidatePool = P.candidatePool;
config        = P.config;

%% Link budget
budgetFile = fullfile( ...
    repoRoot, ...
    'data','generated', ...
    'nominal_link_budget.mat');

assert(isfile(budgetFile), ...
    'Nominal link-budget file not found.');

B = load(budgetFile);

R = B.R;

%% Visible nominal pool
idx = candidatePool.Visible;

pool = candidatePool(idx,:);
Rvis = R(idx,idx);

%% Geometry matrix
H = build_geometry_matrix( ...
    pool, ...
    config.receiver.rECEF_m);

%% Metrics
metrics = compute_position_metrics(H,Rvis);

%% ============================================================
%  REPORT
%  ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('NOMINAL FULL-POOL POSITIONING VALIDATION\n');
fprintf('============================================================\n');

fprintf('Transmitters    : %d\n',height(pool));
fprintf('rank(H)         : %d\n',metrics.rankH);
fprintf('PDOP            : %.6f\n',metrics.PDOP);
fprintf('Position bound  : %.6f m\n',metrics.alphaCRB_m);

%% Previous paper-reference result
referenceBound_m = 0.495580716;

difference_m = ...
    metrics.alphaCRB_m - referenceBound_m;

fprintf('\n');
fprintf('Previous reference : %.6f m\n',referenceBound_m);
fprintf('Difference         : %+.6f m\n',difference_m);
fprintf('Relative difference: %+.3f %%\n', ...
    100*difference_m/referenceBound_m);

fprintf('============================================================\n');

if metrics.rankH == 4
    fprintf('FULL-RANK GEOMETRY: OK\n');
else
    fprintf('FULL-RANK GEOMETRY: FAILED\n');
end

fprintf('============================================================\n');