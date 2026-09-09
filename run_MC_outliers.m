%% Monte Carlo with Gaussian-impulsive contamination and robust BE (v3)
%
% Objective:
%   Evaluate the robustness of Backward Elimination (BE) under a
%   DISTRIBUTION-CONTAMINATION MODEL, following the formulation
%   adopted for the paper.
%
% MIXTURE / CONTAMINATION MODEL
% --------------------------------
% For pseudorange i, the nominal error is:
%
%       n_i ~ N(0,sigma_i^2)
%
% In each Monte Carlo realization, exactly B measurements are randomly
% selected to belong to the contaminated component. Define:
%
%       z_i = 1, if measurement i is contaminated
%       z_i = 0, otherwise,
%
% and the total error is:
%
%       e_i = n_i + z_i * o_i
%
% with impulsive component:
%
%       o_i = s_i * alpha_i * sigma_i
%
%       alpha_i ~ U(20,100)
%       s_i in {-1,+1}.
%
% Therefore, the observed pseudorange is:
%
%       rho_i = rho_i_true + n_i + z_i*o_i.
%
% In terms of the marginal distribution, this experiment can be
% interpreted as a MIXTURE/CONTAMINATION MODEL:
%
%       f_e(e) = (1-epsilon) f_G(e) + epsilon f_out(e),
%
% com:
%
%       epsilon = B/N,
%
%       f_G      = nominal Gaussian component,
%       f_out    = Gaussian component contaminated by a symmetric
%                  impulse with magnitude U(20 sigma_i,100 sigma_i).
%
% Therefore, the second component is NOT Cauchy and is NOT a second
% Gaussian distribution. It is the impulsive component specified by the
% MATLAB model adopted for this paper.
%
% IMPORTANT - model fidelity:
%   - The impulse is ADDED to the Gaussian pseudorange realization.
%   - The R matrix remains NOMINAL.
%   - sigma_i in R is not increased when contamination occurs.
%   - The algorithm knows B as a controlled simulation parameter,
%     but does NOT know the contaminated indices.
%   - Each sigma_i is obtained from the Scenario 4 clear-sky link budget.
%
% -------------------------------------------------------------------------
% MODEL
% -------------------------------------------------------------------------
%
% Deterministic BE in Scenario 4d depends only on H and R; therefore,
% it cannot directly observe an impulsive pseudorange realization.
%
% In this version, robustness to outliers is introduced through a
% NORMALIZED-RESIDUAL elimination stage:
%
%   1) Start with the full candidate pool.
%   2) Estimate [x y z b]^T by WLS using nominal R.
%   3) Compute the pseudorange residuals.
%   4) Normalize the residuals using the residual covariance:
%
%          Q_v = R - H (H'R^{-1}H)^{-1} H'
%
%          z_i = |v_i| / sqrt(Q_v(i,i))
%
%   5) Remove the measurement with the largest z_i while preserving rank(H)=4.
%   6) Repeat B times.
%
% In this version, B is known because it is a CONTROLLED PARAMETER
% of the Monte Carlo simulation. The algorithm knows the value of B but
% does NOT know which measurements were contaminated.
%
% After rejecting B measurements, the geometry/measurement-weighted BE
% from Scenario 4d is applied to the remaining set, using the same target:
%
%       B_pos <= 0.6 m
%
% This explicitly separates:
%
%   Phase A: outlier detection/removal using residuals;
%   Phase B: subset-size reduction using the positional bound.
%
% Future work may replace the known value B by a statistical threshold
% applied to the normalized residuals.
%
% -------------------------------------------------------------------------
% METRICS
% -------------------------------------------------------------------------
%
% For each value of B:
%   - RMSE of WLS without rejection (full pool);
%   - RMSE after robust BE;
%   - 95th and 99th percentiles of the position error;
%   - operational probability P(e_pos <= target_m);
%   - mean outlier-identification rate;
%   - probability of identifying exactly the B contaminated indices;
%   - RMSE conditioned on exact detection;
%   - RMSE conditioned on inexact detection;
%   - mean final subset size;
%   - mean positional bound of the final subset;
%   - probability that the final subset satisfies the nominal BOUND.
%
% IMPORTANT:
%   BoundFeasibilityProb measures theoretical-bound feasibility after
%   selection. It should not be confused with P(e_pos <= target_m), which is
%   an empirical metric obtained directly from the Monte Carlo realizations.
%
%
%
%
% Outputs:
%   results_scenario4_BE_outliers_MC/
%       scenario4e_MC_summary_v3.csv
%       scenario4e_MC_trials_v3.csv
%       scenario4e_MC_top10_errors_v3.csv
%       scenario4e_MC_results_v3.mat

close all;
clear;
clc;
format long;

%% ##################### Input files #######################################

scriptDir = fileparts(mfilename('fullpath'));

if strlength(string(scriptDir)) == 0
    scriptDir = pwd;
end

candidateMat = fullfile( ...
    scriptDir, ...
    "results", ...
    "candidate_pool_26.mat");

linkBudgetMat = fullfile( ...
    scriptDir, ...
    "results", ...
    "link_budget.mat");

if ~isfile(candidateMat)
    error("Candidate pool not found: %s",candidateMat);
end

if ~isfile(linkBudgetMat)
    error("Link budget not found: %s",linkBudgetMat);
end

outDir = fullfile(scriptDir,"results","MC_outliers");

if ~exist(outDir,'dir')
    mkdir(outDir);
end

%% ##################### Load data #########################################

S = load(candidateMat);
LB = load(linkBudgetMat);

requiredPoolVars = [ ...
    "UserPositionECEF", ...
    "PositionECEF", ...
    "Label", ...
    "Architecture" ...
];

for k = 1:numel(requiredPoolVars)
    if ~isfield(S,requiredPoolVars(k))
        error("Missing variable in candidate pool: %s",requiredPoolVars(k));
    end
end

requiredLBVars = [ ...
    "Label", ...
    "Architecture", ...
    "SigmaRho_m", ...
    "CN0_dBHz" ...
];

for k = 1:numel(requiredLBVars)
    if ~isfield(LB,requiredLBVars(k))
        error("Missing variable in link budget: %s",requiredLBVars(k));
    end
end

UserPositionECEF = S.UserPositionECEF(:);
PositionECEF     = S.PositionECEF;
Label            = string(S.Label(:));
Architecture     = string(S.Architecture(:));

sigmaRho_m = LB.SigmaRho_m(:);
CN0_dBHz   = LB.CN0_dBHz(:);

LabelLB        = string(LB.Label(:));
ArchitectureLB = string(LB.Architecture(:));

Npool = size(PositionECEF,1);

if numel(sigmaRho_m) ~= Npool
    error("Number of sigma values does not match the candidate pool.");
end

if any(LabelLB ~= Label)
    error("Label order differs between candidate pool and link budget.");
end

if any(ArchitectureLB ~= Architecture)
    error("Architecture order differs between candidate pool and link budget.");
end

if any(~isfinite(sigmaRho_m)) || any(sigmaRho_m <= 0)
    error("SigmaRho_m contains invalid values.");
end

%% ##################### Monte Carlo parameters ############################
%
% Start with a reduced number of realizations for runtime validation.
% Use N_MC = 10000 for the final paper results.

N_MC = 10000;
% Runtime test setting
% N_MC = 100;


% Even B values preserve exact positive/negative balancing
% in the contamination model. The function below also supports odd B values.
B_values = [0 2 4 6];

mu_gauss = 0;

outlierLowerSigma = 20;
outlierUpperSigma = 100;

% Positional-bound target, consistent with Scenarios 4a-4d.
target_m = 0.6;
targetTol = 1e-12;

% True state.
trueBias_m = 0;
thetaTrue = [UserPositionECEF; trueBias_m];

% Deterministic WLS initialization offset.
%
% The iterative WLS estimator is intentionally initialized away from the
% true user position to avoid an artificially favorable initialization.
% This offset is used only as the numerical starting point of the estimator;
% it is not a pseudorange error, a receiver uncertainty, or a link-budget
% parameter.
initialOffset_m = [100; -100; 50];
theta0 = [UserPositionECEF + initialOffset_m; 0];

% Maximum number of WLS iterations.
maxWLSIterations = 30;

% WLS convergence criterion.
wlsTolerance_m = 1e-6;

% Reproducibility.
rng(2,'twister');

%% ##################### Geometry and true pseudorange ######################

diffTrue = PositionECEF - UserPositionECEF.';
trueRange_m = sqrt(sum(diffTrue.^2,2));

rhoTrue_m = trueRange_m + trueBias_m;

u = diffTrue ./ trueRange_m;

Hfull = [-u,ones(Npool,1)];

if rank(Hfull) < 4
    error("The full candidate pool has rank(H)<4.");
end

[fullBound_m,fullPDOP,fullRankH] = ...
    subsetMetrics((1:Npool).',Hfull,sigmaRho_m);

fprintf('\n');
fprintf('============================================================\n');
fprintf('ROBUST BE UNDER IMPULSIVE OUTLIERS - GEOMETRY\n');
fprintf('============================================================\n');
fprintf('Number of candidates  : %d\n',Npool);
if N_MC == 10000
    fprintf('Monte Carlo trials/case: %d\n',N_MC);
else
    fprintf('Monte Carlo trials/case: %d (Validation case)\n',N_MC);
end
fprintf('B values              : ');
fprintf('%d ',B_values);
fprintf('\n');
fprintf('Outlier amplitude     : U(%.0f,%.0f) x sigma_i\n', ...
    outlierLowerSigma,outlierUpperSigma);
fprintf('Positional-bound target: %.3f m\n',target_m);
fprintf('Full-pool bound       : %.6f m\n',fullBound_m);
fprintf('\n');

%% ##################### Result preallocation ##############################

nB = numel(B_values);
nTrials = nB*N_MC;

TrialID              = zeros(nTrials,1);
BColumn              = zeros(nTrials,1);
ErrorFullWLS_m       = NaN(nTrials,1);
ErrorRobustBE_m      = NaN(nTrials,1);
KfinalBE             = NaN(nTrials,1);
FinalBoundBE_m       = NaN(nTrials,1);
FinalPDOPBE          = NaN(nTrials,1);
BoundFeasible         = false(nTrials,1);
CorrectOutliers      = NaN(nTrials,1);
DetectionRate        = NaN(nTrials,1);
ExactDetection       = false(nTrials,1);
WLSFullConverged     = false(nTrials,1);
WLSRobustConverged   = false(nTrials,1);

% Per-realization operational metrics.
PositionErrorTarget_m = target_m;

ContaminationModel = ...
    "Gaussian-impulsive contamination mixture: N(0,sigma_i^2) + sparse symmetric U(20,100)*sigma_i impulses";

FullErrorWithinTarget = false(nTrials,1);
BEErrorWithinTarget   = false(nTrials,1);

trialRow = 0;

%% ##################### Monte Carlo simulation ############################

for ib = 1:nB

    B = B_values(ib);

    if B > Npool-4
        error("B=%d leaves fewer than four measurements.",B);
    end

    fprintf('------------------------------------------------------------\n');
    fprintf('Running B = %d contaminated measurements ...\n',B);

    for mc = 1:N_MC

        trialRow = trialRow + 1;

        %% --------- Nominal Gaussian noise --------------------------------
        noise_gauss = mu_gauss + sigmaRho_m.*randn(Npool,1);

        %% --------- Contaminated mixture component ------------------------
        %
        % The impulsive-outlier generation used in this section is based on
        % MATLAB code provided by Prof. Apolinario and adapted here to the
        % heterogeneous pseudorange uncertainties sigma_i obtained from the
        % Scenario 4 link budget. The adaptation preserves the original idea
        % of randomly selecting B contaminated measurements and adding
        % symmetric impulses with amplitudes between 20 and 100 times the
        % corresponding nominal sigma_i.
        %
        % Implementation of the contaminated component:
        %   - randomly select B indices;
        %   - use balanced positive/negative signs;
        %   - use amplitudes between 20 and 100 sigma.
        %
        % For heterogeneous sigma_i, the multiplier is applied to the nominal
        % sigma of the corresponding pseudorange.

        % Mixture-component indicator:
        %   contaminationIndicator(i)=1 -> contaminated component
        %   contaminationIndicator(i)=0 -> nominal Gaussian component
        contaminationIndicator = false(Npool,1);

        impulsiveComponent = zeros(Npool,1);
        trueOutlierIdx = zeros(0,1);

        if B > 0

            trueOutlierIdx = randperm(Npool,B).';
            contaminationIndicator(trueOutlierIdx) = true;

            multipliers = balancedOutlierMultipliers( ...
                B,outlierLowerSigma,outlierUpperSigma);

            % Shuffle signs/magnitudes among the selected indices.
            % No additional random call is introduced here,
            % preserving the intended pseudorandom sequence.
            multipliers = multipliers(randperm(B));

            impulsiveComponent(trueOutlierIdx) = ...
                multipliers(:).*sigmaRho_m(trueOutlierIdx);
        end

        % Mixture/contamination model:
        %   error = Gaussian component + sparse impulsive component.
        measurementError = noise_gauss + impulsiveComponent;

        rhoMeasured_m = rhoTrue_m + measurementError;

        %% --------- Baseline: WLS on the full pool -------------------------
        [thetaFull,~,~,convFull] = iterativeWLS( ...
            rhoMeasured_m, ...
            PositionECEF, ...
            sigmaRho_m, ...
            theta0, ...
            maxWLSIterations, ...
            wlsTolerance_m);

        errFull = norm(thetaFull(1:3)-UserPositionECEF);

        %% --------- Phase A: normalized-residual elimination ---------------
        active = (1:Npool).';
        removedByResidualBE = zeros(0,1);

        for j = 1:B

            rhoActive = rhoMeasured_m(active);
            posActive = PositionECEF(active,:);
            sigActive = sigmaRho_m(active);

            [thetaActive,residualActive,Hactive,~] = iterativeWLS( ...
                rhoActive, ...
                posActive, ...
                sigActive, ...
                theta0, ...
                maxWLSIterations, ...
                wlsTolerance_m);

            normalizedResidual = computeNormalizedResidual( ...
                residualActive,Hactive,sigActive);

            % Sort normalized residuals from largest to smallest.
            [~,order] = sort(normalizedResidual,'descend');

            removedThisStep = NaN;

            % Remove the largest residual whose removal preserves rank(H)=4.
            for q = 1:numel(order)

                localCandidate = order(q);
                globalCandidate = active(localCandidate);

                trialActive = active(active ~= globalCandidate);

                if numel(trialActive) < 4
                    continue;
                end

                if rank(Hfull(trialActive,:)) < 4
                    continue;
                end

                removedThisStep = globalCandidate;
                break;
            end

            if isnan(removedThisStep)
                % This should not occur for the B values used here, but the safeguard
                % prevents a removal that would destroy full rank.
                break;
            end

            active(active == removedThisStep) = [];
            removedByResidualBE(end+1,1) = removedThisStep; %#ok<SAGROW>
        end

        %% --------- Detection metrics --------------------------------------
        if B == 0

            nCorrect = 0;
            detectionRate = 1;
            exactDetection = true;

        else

            nCorrect = numel(intersect( ...
                trueOutlierIdx,removedByResidualBE));

            detectionRate = nCorrect/B;

            exactDetection = ...
                numel(removedByResidualBE)==B && ...
                isequal(sort(removedByResidualBE),sort(trueOutlierIdx));
        end

        %% --------- Phase B: BE using the positional bound -----------------
        %
        % After rejecting the B suspected indices, apply
        % the Scenario 4d BE logic to the remaining pool.

        [selectedBE,boundBE,pdopBE,rankBE,targetBE] = ...
            backwardSelectByBound( ...
                active,Hfull,sigmaRho_m,target_m,targetTol);

        %% --------- Final WLS on the robust subset -------------------------
        [thetaRobust,~,~,convRobust] = iterativeWLS( ...
            rhoMeasured_m(selectedBE), ...
            PositionECEF(selectedBE,:), ...
            sigmaRho_m(selectedBE), ...
            theta0, ...
            maxWLSIterations, ...
            wlsTolerance_m);

        errRobust = norm(thetaRobust(1:3)-UserPositionECEF);

        %% --------- Store results -----------------------------------------
        TrialID(trialRow)            = mc;
        BColumn(trialRow)            = B;
        ErrorFullWLS_m(trialRow)     = errFull;
        ErrorRobustBE_m(trialRow)    = errRobust;
        KfinalBE(trialRow)           = numel(selectedBE);
        FinalBoundBE_m(trialRow)     = boundBE;
        FinalPDOPBE(trialRow)        = pdopBE;
        BoundFeasible(trialRow)       = targetBE && rankBE==4;
        CorrectOutliers(trialRow)    = nCorrect;
        DetectionRate(trialRow)      = detectionRate;
        ExactDetection(trialRow)     = exactDetection;
        WLSFullConverged(trialRow)   = convFull;
        WLSRobustConverged(trialRow) = convRobust;

        FullErrorWithinTarget(trialRow) = ...
            isfinite(errFull) && errFull <= PositionErrorTarget_m;

        BEErrorWithinTarget(trialRow) = ...
            isfinite(errRobust) && errRobust <= PositionErrorTarget_m;
    end

    fprintf('Completed B = %d.\n',B);
end

%% ##################### Per-realization table ############################

Ttrials = table( ...
    TrialID, ...
    BColumn, ...
    ErrorFullWLS_m, ...
    ErrorRobustBE_m, ...
    FullErrorWithinTarget, ...
    BEErrorWithinTarget, ...
    KfinalBE, ...
    FinalBoundBE_m, ...
    FinalPDOPBE, ...
    BoundFeasible, ...
    CorrectOutliers, ...
    DetectionRate, ...
    ExactDetection, ...
    WLSFullConverged, ...
    WLSRobustConverged, ...
    'VariableNames',{ ...
    'TrialID', ...
    'B', ...
    'ErrorFullWLS_m', ...
    'ErrorRobustBE_m', ...
    'FullErrorWithinTarget', ...
    'BEErrorWithinTarget', ...
    'KfinalBE', ...
    'FinalBoundBE_m', ...
    'FinalPDOPBE', ...
    'BoundFeasible', ...
    'CorrectOutliers', ...
    'DetectionRate', ...
    'ExactDetection', ...
    'WLSFullConverged', ...
    'WLSRobustConverged'});

%% ##################### Summary by B ######################################

Bsummary = B_values(:);

RMSE_FullWLS_m            = NaN(nB,1);
RMSE_RobustBE_m           = NaN(nB,1);
P95_FullWLS_m             = NaN(nB,1);
P95_RobustBE_m            = NaN(nB,1);
P99_FullWLS_m             = NaN(nB,1);
P99_RobustBE_m            = NaN(nB,1);
Median_RobustBE_m         = NaN(nB,1);
Prob_Full_ErrorLETarget       = NaN(nB,1);
Prob_BE_ErrorLETarget         = NaN(nB,1);
RMSE_BE_ExactDetection_m  = NaN(nB,1);
RMSE_BE_InexactDetection_m= NaN(nB,1);
MeanK_BE                  = NaN(nB,1);
MedianK_BE                = NaN(nB,1);
MeanFinalBound_BE_m       = NaN(nB,1);
BoundFeasibilityProb      = NaN(nB,1);
MeanDetectionRate         = NaN(nB,1);
ExactDetectionProb        = NaN(nB,1);
WLSFullConvergenceProb    = NaN(nB,1);
WLSBEConvergenceProb      = NaN(nB,1);
MaxFullWLS_Error_m        = NaN(nB,1);
MaxRobustBE_Error_m       = NaN(nB,1);

for ib = 1:nB

    B = B_values(ib);
    mask = BColumn==B;

    eFull = ErrorFullWLS_m(mask);
    eBE   = ErrorRobustBE_m(mask);

    RMSE_FullWLS_m(ib) = sqrt(mean(eFull.^2,'omitnan'));
    RMSE_RobustBE_m(ib) = sqrt(mean(eBE.^2,'omitnan'));

    P95_FullWLS_m(ib) = empiricalPercentile(eFull,95);
    P95_RobustBE_m(ib) = empiricalPercentile(eBE,95);

    P99_FullWLS_m(ib) = empiricalPercentile(eFull,99);
    P99_RobustBE_m(ib) = empiricalPercentile(eBE,99);

    Median_RobustBE_m(ib) = median(eBE,'omitnan');

    Prob_Full_ErrorLETarget(ib) = ...
        mean(FullErrorWithinTarget(mask));

    Prob_BE_ErrorLETarget(ib) = ...
        mean(BEErrorWithinTarget(mask));

    exactMaskLocal = ExactDetection(mask);
    inexactMaskLocal = ~ExactDetection(mask);

    if any(exactMaskLocal)
        eExact = eBE(exactMaskLocal);
        RMSE_BE_ExactDetection_m(ib) = ...
            sqrt(mean(eExact.^2,'omitnan'));
    end

    if any(inexactMaskLocal)
        eInexact = eBE(inexactMaskLocal);
        RMSE_BE_InexactDetection_m(ib) = ...
            sqrt(mean(eInexact.^2,'omitnan'));
    end

    MeanK_BE(ib) = mean(KfinalBE(mask),'omitnan');
    MedianK_BE(ib) = median(KfinalBE(mask),'omitnan');

    MeanFinalBound_BE_m(ib) = ...
        mean(FinalBoundBE_m(mask),'omitnan');

    BoundFeasibilityProb(ib) = ...
        mean(BoundFeasible(mask));

    MeanDetectionRate(ib) = ...
        mean(DetectionRate(mask),'omitnan');

    ExactDetectionProb(ib) = ...
        mean(ExactDetection(mask));

    WLSFullConvergenceProb(ib) = ...
        mean(WLSFullConverged(mask));

    WLSBEConvergenceProb(ib) = ...
        mean(WLSRobustConverged(mask));

    MaxFullWLS_Error_m(ib) = max(eFull,[],'omitnan');
    MaxRobustBE_Error_m(ib) = max(eBE,[],'omitnan');
end

Tsummary = table( ...
    Bsummary, ...
    RMSE_FullWLS_m, ...
    RMSE_RobustBE_m, ...
    P95_FullWLS_m, ...
    P95_RobustBE_m, ...
    P99_FullWLS_m, ...
    P99_RobustBE_m, ...
    Median_RobustBE_m, ...
    Prob_Full_ErrorLETarget, ...
    Prob_BE_ErrorLETarget, ...
    RMSE_BE_ExactDetection_m, ...
    RMSE_BE_InexactDetection_m, ...
    MeanK_BE, ...
    MedianK_BE, ...
    MeanFinalBound_BE_m, ...
    BoundFeasibilityProb, ...
    MeanDetectionRate, ...
    ExactDetectionProb, ...
    WLSFullConvergenceProb, ...
    WLSBEConvergenceProb, ...
    MaxFullWLS_Error_m, ...
    MaxRobustBE_Error_m, ...
    'VariableNames',{ ...
    'B', ...
    'RMSE_FullWLS_m', ...
    'RMSE_RobustBE_m', ...
    'P95_FullWLS_m', ...
    'P95_RobustBE_m', ...
    'P99_FullWLS_m', ...
    'P99_RobustBE_m', ...
    'Median_RobustBE_m', ...
    'Prob_Full_ErrorLETarget', ...
    'Prob_BE_ErrorLETarget', ...
    'RMSE_BE_ExactDetection_m', ...
    'RMSE_BE_InexactDetection_m', ...
    'MeanK_BE', ...
    'MedianK_BE', ...
    'MeanFinalBound_BE_m', ...
    'BoundFeasibilityProb', ...
    'MeanDetectionRate', ...
    'ExactDetectionProb', ...
    'WLSFullConvergenceProb', ...
    'WLSBEConvergenceProb', ...
    'MaxFullWLS_Error_m', ...
    'MaxRobustBE_Error_m'});

fprintf('\n');
fprintf('============================================================\n');
fprintf('MONTE CARLO SUMMARY - SCENARIO 4e\n');
fprintf('============================================================\n');

disp(' ');
disp('===== Positioning accuracy =====');

TaccuracyDisplay = table( ...
    Bsummary, ...
    round(RMSE_FullWLS_m,3), ...
    round(RMSE_RobustBE_m,3), ...
    round(P95_RobustBE_m,3), ...
    round(P99_RobustBE_m,3), ...
    'VariableNames',{ ...
    'B', ...
    'WLS_RMSE_m', ...
    'RobustBE_RMSE_m', ...
    'RobustBE_P95_m', ...
    'RobustBE_P99_m'});

disp(TaccuracyDisplay);

disp(' ');
disp('===== Outlier detection and final subset =====');

MeanDetectionDisplay = 100*MeanDetectionRate;
ExactDetectionDisplay = 100*ExactDetectionProb;
BoundFeasibleDisplay = 100*BoundFeasibilityProb;

% No outliers are present for B=0, so detection metrics are not applicable.
MeanDetectionDisplay(Bsummary==0) = NaN;
ExactDetectionDisplay(Bsummary==0) = NaN;

TdetectionDisplay = table( ...
    Bsummary, ...
    round(MeanDetectionDisplay,2), ...
    round(ExactDetectionDisplay,2), ...
    round(MeanK_BE,2), ...
    round(BoundFeasibleDisplay,2), ...
    'VariableNames',{ ...
    'B', ...
    'MeanDetection_pct', ...
    'ExactDetection_pct', ...
    'MeanFinalK', ...
    'BoundFeasible_pct'});

disp(TdetectionDisplay);

fprintf(['ExactDetection_pct: trials in which all B injected outliers were ' ...
         'correctly identified.\n']);
fprintf(['BoundFeasible_pct : trials in which the final subset satisfied the ' ...
         'nominal %.3f m positional-bound target.\n'],target_m);

disp(' ');
disp('===== Robust-BE RMSE conditioned on detection outcome =====');

TconditionalDisplay = table( ...
    Bsummary, ...
    round(RMSE_BE_ExactDetection_m,3), ...
    round(RMSE_BE_InexactDetection_m,3), ...
    'VariableNames',{ ...
    'B', ...
    'RMSE_IfExact_m', ...
    'RMSE_IfInexact_m'});

disp(TconditionalDisplay);

fprintf(['NaN in RMSE_IfInexact_m means that no inexact-detection trial ' ...
         'occurred for that value of B.\n']);

%% ##################### Ten largest robust-BE errors by B #################
%
% This table helps diagnose why RMSE may increase even when
% the median and P95 remain relatively low.

TopRows = table();

for ib = 1:nB

    B = B_values(ib);
    idxB = find(BColumn==B);

    e = ErrorRobustBE_m(idxB);

    [~,ord] = sort(e,'descend');

    nTop = min(10,numel(ord));
    use = idxB(ord(1:nTop));

    Ttmp = table( ...
        repmat(B,nTop,1), ...
        (1:nTop).', ...
        TrialID(use), ...
        ErrorRobustBE_m(use), ...
        ErrorFullWLS_m(use), ...
        ExactDetection(use), ...
        DetectionRate(use), ...
        KfinalBE(use), ...
        FinalBoundBE_m(use), ...
        BoundFeasible(use), ...
        'VariableNames',{ ...
        'B', ...
        'TailRank', ...
        'TrialID', ...
        'ErrorRobustBE_m', ...
        'ErrorFullWLS_m', ...
        'ExactDetection', ...
        'DetectionRate', ...
        'KfinalBE', ...
        'FinalBoundBE_m', ...
        'BoundFeasible'});

    TopRows = [TopRows;Ttmp]; %#ok<AGROW>
end


%% ##################### Figures ############################################
% Figures are intentionally generated only by the final paper-output script.
% Scenario 4e computes and saves the numerical Monte Carlo results only.

%% ##################### Save outputs #####################################

writetable(Ttrials, ...
    fullfile(outDir,'scenario4e_MC_trials_v3.csv'));

writetable(Tsummary, ...
    fullfile(outDir,'scenario4e_MC_summary_v3.csv'));

writetable(TopRows, ...
    fullfile(outDir,'scenario4e_MC_top10_errors_v3.csv'));

save(fullfile(outDir,'scenario4e_MC_results_v3.mat'), ...
    'N_MC','B_values', ...
    'outlierLowerSigma','outlierUpperSigma', ...
    'target_m', ...
    'thetaTrue','theta0', ...
    'rhoTrue_m', ...
    'Hfull','sigmaRho_m','CN0_dBHz', ...
    'Label','Architecture', ...
    'PositionErrorTarget_m','ContaminationModel', ...
    'Ttrials','Tsummary','TopRows');

fprintf('\nResults saved to:\n  %s\n',outDir);

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function multipliers = balancedOutlierMultipliers( ...
    B,lowerSigma,upperSigma)
% Generate positive/negative multipliers according to the adopted
% contamination model:
%
%   positive = U(lowerSigma,upperSigma)
%   negative = -positive
%
% For even B, the signs are exactly balanced.
% For odd B, one additional outlier with random sign is generated.

    if B == 0
        multipliers = zeros(0,1);
        return;
    end

    nPairs = floor(B/2);

    positive = lowerSigma + ...
        (upperSigma-lowerSigma).*rand(nPairs,1);

    negative = -positive;

    multipliers = [positive;negative];

    if mod(B,2)==1

        extraMagnitude = lowerSigma + ...
            (upperSigma-lowerSigma).*rand;

        if rand < 0.5
            extraMagnitude = -extraMagnitude;
        end

        multipliers = [multipliers;extraMagnitude]; %#ok<AGROW>
    end
end

function [theta,residual,H,converged] = iterativeWLS( ...
    rho_m,txECEF,sigma_m,theta0,maxIter,tol)
% Iterative WLS estimation of:
%
%   theta = [x y z b]^T
%
% using the model:
%
%   rho_i = ||p_i-p|| + b + erro_i
%
% R remains diag(sigma_i^2), even when the realization contains an outlier.

    theta = theta0(:);
    converged = false;

    rho_m = rho_m(:);
    sigma_m = sigma_m(:);

    for iter = 1:maxIter

        p = theta(1:3).';
        b = theta(4);

        diff = txECEF - p;
        ranges = sqrt(sum(diff.^2,2));

        if any(ranges <= eps)
            error("Invalid distance during WLS.");
        end

        u = diff./ranges;

        predicted = ranges + b;
        residual = rho_m - predicted;

        H = [-u,ones(size(txECEF,1),1)];

        % Numerically simple weighted form:
        % minimize ||R^{-1/2}(rho-h(theta))||^2.
        Hw = H./sigma_m;
        rw = residual./sigma_m;

        Nmat = Hw.'*Hw;

        if rcond(Nmat) <= 1e-14
            break;
        end

        delta = Nmat\(Hw.'*rw);

        theta = theta + delta;

        if norm(delta(1:3)) <= tol && abs(delta(4)) <= tol
            converged = true;
            break;
        end
    end

    % Recompute residuals/H at the final state.
    p = theta(1:3).';
    b = theta(4);

    diff = txECEF - p;
    ranges = sqrt(sum(diff.^2,2));
    u = diff./ranges;

    predicted = ranges + b;
    residual = rho_m - predicted;
    H = [-u,ones(size(txECEF,1),1)];
end

function z = computeNormalizedResidual(residual,H,sigma_m)
% Normalized residuals using the residual covariance:
%
%   Qv = R - H (H'R^{-1}H)^(-1) H'
%
%   z_i = |v_i| / sqrt(Qv_ii)

    residual = residual(:);
    sigma_m = sigma_m(:);

    R = diag(sigma_m.^2);

    J = H.'*(R\H);

    if rcond(J) <= 1e-14
        z = abs(residual)./sigma_m;
        return;
    end

    C = J\eye(size(J));

    Qv = R - H*C*H.';

    q = real(diag(Qv));

    % Protect against small negative values caused by roundoff.
    q(q < eps) = eps;

    z = abs(residual)./sqrt(q);
end

function [selected,bound_m,pdop,rankH,targetReached] = ...
    backwardSelectByBound(initialSet,Hfull,sigmaRho_m,target_m,targetTol)
% Deterministic Scenario 4d BE applied to an initial set that has already
% passed through the outlier-rejection phase.

    selected = initialSet(:);

    [bound_m,pdop,rankH] = ...
        subsetMetrics(selected,Hfull,sigmaRho_m);

    targetReached = ...
        rankH==4 && bound_m <= target_m + targetTol;

    % If the remaining set can no longer satisfy the target,
    % do not remove additional measurements.
    if ~targetReached
        return;
    end

    while numel(selected) > 4

        bestRemoved = NaN;
        bestBound = inf;
        bestPDOP = inf;

        for k = 1:numel(selected)

            candidate = selected(k);
            trial = selected(selected ~= candidate);

            [trialBound,trialPDOP,trialRank] = ...
                subsetMetrics(trial,Hfull,sigmaRho_m);

            if trialRank ~= 4
                continue;
            end

            if trialBound > target_m + targetTol
                continue;
            end

            if trialBound < bestBound

                bestRemoved = candidate;
                bestBound = trialBound;
                bestPDOP = trialPDOP;

            elseif abs(trialBound-bestBound) <= 1e-14

                if isnan(bestRemoved) || candidate < bestRemoved
                    bestRemoved = candidate;
                    bestBound = trialBound;
                    bestPDOP = trialPDOP;
                end
            end
        end

        if isnan(bestRemoved)
            break;
        end

        selected(selected == bestRemoved) = [];

        bound_m = bestBound;
        pdop = bestPDOP;
    end

    [bound_m,pdop,rankH] = ...
        subsetMetrics(selected,Hfull,sigmaRho_m);

    targetReached = ...
        rankH==4 && bound_m <= target_m + targetTol;
end

function [bound_m,PDOP,rankH] = subsetMetrics( ...
    idx,Hfull,sigmaRho_m)
% Same metric as in Scenario 4d:
%
%   J = H'R^{-1}H
%   B_pos = sqrt(trace(J^{-1}(1:3,1:3)))

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

function p = empiricalPercentile(x,percent)
% Simple percentile calculation without Statistics Toolbox dependency.

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

    pos = 1 + (numel(x)-1)*(percent/100);

    lo = floor(pos);
    hi = ceil(pos);

    if lo==hi
        p = x(lo);
    else
        w = pos-lo;
        p = (1-w)*x(lo) + w*x(hi);
    end
end
