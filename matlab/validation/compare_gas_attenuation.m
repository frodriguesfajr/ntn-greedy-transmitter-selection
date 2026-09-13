clear;
clc;

%% Paths

thisFile = mfilename('fullpath');
atmDir   = fileparts(thisFile);
matlabDir = fileparts(atmDir);
repoRoot = fileparts(matlabDir);

addpath(atmDir);

%% Candidate pool

poolFile = fullfile( ...
    repoRoot,'data','generated', ...
    'nominal_candidate_pool.mat');

S = load(poolFile);

candidatePool = S.candidatePool;

%% Old MATLAB-toolbox result

linkFile = fullfile( ...
    repoRoot,'data','generated', ...
    'nominal_link_budget.mat');

L = load(linkFile);

oldLoss_dB = ...
    L.linkBudget.AtmosphericGasLoss_dB;

N = height(candidatePool);

newLoss_dB = nan(N,1);

%% Compute independently

for k = 1:N

    fGHz = candidatePool.Frequency_GHz(k);
    el   = candidatePool.Elevation_deg(k);

    if candidatePool.Architecture(k) == "HAPS"
        topHeight_km = 20;
    else
        topHeight_km = 100;
    end

    newLoss_dB(k) = itu676_slant_path( ...
        fGHz, ...
        el, ...
        0.010, ...
        topHeight_km);
end

%% Comparison

delta_dB = newLoss_dB-oldLoss_dB;

relative_pct = ...
    100*delta_dB./oldLoss_dB;

comparison = table( ...
    candidatePool.Label, ...
    candidatePool.Architecture, ...
    candidatePool.Elevation_deg, ...
    candidatePool.Frequency_GHz, ...
    oldLoss_dB, ...
    newLoss_dB, ...
    delta_dB, ...
    relative_pct, ...
    'VariableNames',{ ...
    'Label', ...
    'Architecture', ...
    'Elevation_deg', ...
    'Frequency_GHz', ...
    'MATLAB_dB', ...
    'ITU676_P835_dB', ...
    'Difference_dB', ...
    'Difference_pct'});

disp(comparison);

fprintf('\n');
fprintf('==============================================\n');
fprintf('GASEOUS ATTENUATION COMPARISON\n');
fprintf('==============================================\n');

fprintf('Mean absolute difference : %.6f dB\n', ...
    mean(abs(delta_dB)));

fprintf('Maximum absolute difference: %.6f dB\n', ...
    max(abs(delta_dB)));

fprintf('Mean relative difference : %.3f %%\n', ...
    mean(abs(relative_pct)));

fprintf('==============================================\n');