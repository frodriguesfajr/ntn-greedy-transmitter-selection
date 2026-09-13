clear;
clc;

%% ============================================================
%  VALIDATE TRANSMITTER SELECTION AGAINST THE ORIGINAL
%  MATLAB-TOOLBOX GASEOUS ATTENUATION VALUES
%
%  This script is for validation only.
%  The main reproducible pipeline remains based on P.2145/P.676.
%  ============================================================

%% Repository paths

thisFile  = mfilename('fullpath');
validationDir = fileparts(thisFile);
matlabDir = fileparts(validationDir);
repoRoot  = fileparts(matlabDir);

addpath(fullfile(repoRoot,'matlab','core'));
addpath(fullfile(repoRoot,'matlab','algorithms'));

%% Load candidate pool

P = load(fullfile( ...
    repoRoot,'data','generated', ...
    'nominal_candidate_pool.mat'));

candidatePool = P.candidatePool;
config        = P.config;

%% Load current budget
%
% We reuse all current RF quantities except gaseous attenuation.

B = load(fullfile( ...
    repoRoot,'data','generated', ...
    'nominal_link_budget.mat'));

linkBudget = B.linkBudget;
linkConfig = B.linkConfig;

%% ============================================================
%  ORIGINAL MATLAB TOOLBOX GAS LOSSES
%  ============================================================

gasLossReference_dB = [ ...
    0.464872051720006
    0.915619201385095
    0.915619201385093
    1.083269801941140
    0.579644683798001
    0.746529809423433
    0.823309107604283
    1.011147561592230
    1.014773052147900
    1.223282587959560
    1.307988823149580
    1.451125653190190
    0.736578957477030
    0.740967662299461
    0.879755276729633
    1.095074353889510
    2.077132816727250
    2.077451994238410
    4.049861773076530
    5.131094883426020
    0.924362770672161
    0.784988441756104
    0.682794565160070
    0.630299567781078
    0.651756984573601
    0.692143534373514 ];

assert(numel(gasLossReference_dB)==height(candidatePool));

%% ============================================================
%  REBUILD ORIGINAL MEASUREMENT COVARIANCE
%  ============================================================

etaPos  = linkConfig.etaPos;
beta_Hz = linkConfig.beta_Hz;
Tcoh_s  = linkConfig.Tcoh_s;

c_mps = 299792458;

CN0_Link_ref_dBHz = ...
    linkBudget.EIRP_Ranging_dBW ...
    - linkBudget.FSPL_dB ...
    - gasLossReference_dB ...
    + linkBudget.ReceiverGT_dBK ...
    + 228.6;

gamma_ref_Hz = ...
    etaPos .* 10.^(CN0_Link_ref_dBHz/10);

sigma_ref_m = ...
    c_mps ./ ...
    (2*pi*beta_Hz .* sqrt(gamma_ref_Hz*Tcoh_s));

Rref = diag(sigma_ref_m.^2);

%% Geometry

H = build_geometry_matrix( ...
    candidatePool, ...
    config.receiver.rECEF_m);

%% Full-pool reference

fullRef = compute_position_metrics(H,Rref);

%% Selection

epsilon_m = 0.6;

ftsRef = forward_transmitter_selection( ...
    H,Rref,epsilon_m);

btsRef = backward_transmitter_selection( ...
    H,Rref,epsilon_m);

%% Selected pools

ftsPool = candidatePool(ftsRef.selectedIdx,:);
btsPool = candidatePool(btsRef.selectedIdx,:);

arch = ["HAPS","LEO","MEO","GEO"];

ftsCounts = zeros(1,4);
btsCounts = zeros(1,4);

for k = 1:4
    ftsCounts(k) = sum(ftsPool.Architecture==arch(k));
    btsCounts(k) = sum(btsPool.Architecture==arch(k));
end

%% ============================================================
%  REPORT
%  ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('REFERENCE-SELECTION VALIDATION\n');
fprintf('============================================================\n');

fprintf('\nORIGINAL MATLAB GAS MODEL\n');
fprintf('------------------------------------------------------------\n');
fprintf('Full-pool bound : %.6f m\n', ...
    fullRef.alphaCRB_m);

fprintf('\nFTS\n');
fprintf('K               : %d\n',ftsRef.K);
fprintf('Bound           : %.6f m\n',ftsRef.positionBound_m);
fprintf('Evaluations     : %d\n',ftsRef.evaluations);
fprintf('Composition     : HAPS=%d LEO=%d MEO=%d GEO=%d\n', ...
    ftsCounts);

fprintf('\nBTS\n');
fprintf('K               : %d\n',btsRef.K);
fprintf('Bound           : %.6f m\n',btsRef.positionBound_m);
fprintf('Evaluations     : %d\n',btsRef.evaluations);
fprintf('Composition     : HAPS=%d LEO=%d MEO=%d GEO=%d\n', ...
    btsCounts);

fprintf('\nSame subset     : %d\n', ...
    isequal(sort(ftsRef.selectedIdx), ...
            sort(btsRef.selectedIdx)));

fprintf('\nSelected FTS:\n');
disp(ftsPool(:,{'Label','Architecture','NORAD'}));

fprintf('\n============================================================\n');