%% scenario4_task1_run_fts_bts_4epochs.m
% TASK CODE 1 - STEP 3
% Compare No Selection, FTS, and BTS at four temporal epochs.
%
% Repository:
%   /MATLAB Drive/ntn-greedy-transmitter-selection
%
% Inputs:
%   results/temporal_geometry_4epochs/
%       scenario4_geometry_4epochs.mat
%
%   results/temporal_link_budget_4epochs_portable/
%       scenario4_link_budget_4epochs_portable.mat
%
% Optional nominal-validation inputs:
%   results/BE_selected.csv
%   results/BE_results.mat
%
% Epochs:
%   t = 0, 30, 60, 90 s
%
% At each epoch:
%   1) use the visible candidate set;
%   2) compute the No Selection B_CRB;
%   3) run Forward Transmitter Selection (FTS);
%   4) run Backward Transmitter Selection (BTS);
%   5) record K, B_CRB, PDOP, composition, selected labels,
%      feasibility, and candidate-subset assessments;
%   6) compare whether FTS and BTS return the same subset.
%
% Figures:
%   - temporal_bcrb_fts_bts.pdf/png
%   - skyplot_t000s.pdf/png
%   - skyplot_t030s.pdf/png
%   - skyplot_t060s.pdf/png
%   - skyplot_t090s.pdf/png
%
% Main outputs:
%   results/temporal_fts_bts_4epochs/
%       temporal_fts_bts_summary.csv
%       temporal_fts_bts_selected_sets.csv
%       scenario4_temporal_fts_bts_4epochs.mat
%
% Nomenclature:
%   No Selection = all visible candidates
%   FTS          = Forward Transmitter Selection
%   BTS          = Backward Transmitter Selection
%   B_CRB        = scalar position-domain Cramer-Rao bound
%
% -------------------------------------------------------------------------

close all;
clear;
clc;
format long;

%% ========================================================================
% Repository and input files
% =========================================================================

paperRoot = '/MATLAB Drive/ntn-greedy-transmitter-selection';

if ~isfolder(paperRoot)
    error('Paper repository not found: %s',paperRoot);
end

cd(paperRoot);

geometryMat = fullfile( ...
    paperRoot, ...
    'results', ...
    'temporal_geometry_4epochs', ...
    'scenario4_geometry_4epochs.mat');

linkMat = fullfile( ...
    paperRoot, ...
    'results', ...
    'temporal_link_budget_4epochs_portable', ...
    'scenario4_link_budget_4epochs_portable.mat');

beSelectedCsv = fullfile( ...
    paperRoot,'results','BE_selected.csv');

beResultsMat = fullfile( ...
    paperRoot,'results','BE_results.mat');

outDir = fullfile( ...
    paperRoot,'results','temporal_fts_bts_4epochs');

if ~isfile(geometryMat)
    error('Geometry MAT not found: %s',geometryMat);
end

if ~isfile(linkMat)
    error(['Portable temporal link-budget MAT not found:\n%s\n' ...
           'Run scenario4_task1_build_temporal_link_budget_portable.m first.'], ...
           linkMat);
end

if ~isfolder(outDir)
    mkdir(outDir);
end

%% ========================================================================
% Load validated temporal geometry and link budget
% =========================================================================

G = load(geometryMat);
L = load(linkMat);

requiredGeometryFields = { ...
    'epochs','PoolIndex','Label','Architecture', ...
    'NORAD_CAT_ID','timeOffset_s','timeUTC'};

for ii = 1:numel(requiredGeometryFields)
    if ~isfield(G,requiredGeometryFields{ii})
        error('Geometry MAT missing field: %s',requiredGeometryFields{ii});
    end
end

if ~isfield(L,'linkEpochs') || ~isfield(L,'linkMetadata')
    error('Portable link-budget MAT is incomplete.');
end

if isfield(L.linkMetadata,'ValidationPassed') && ...
        ~logical(L.linkMetadata.ValidationPassed)
    error('Portable temporal link budget was not validated at t=0.');
end

epochs = G.epochs;
linkEpochs = L.linkEpochs;

PoolIndex = double(G.PoolIndex(:));
Label = string(G.Label(:));
Architecture = upper(string(G.Architecture(:)));
NORAD_CAT_ID = string(G.NORAD_CAT_ID(:));

timeOffset_s = double(G.timeOffset_s(:));
timeUTC = G.timeUTC(:);

N = numel(Label);
Ntime = numel(timeOffset_s);

if N ~= 26
    error('Expected 26 candidate identities; found %d.',N);
end

if Ntime ~= 4 || ~isequal(timeOffset_s,[0;30;60;90])
    error('Expected epochs t = [0 30 60 90] s.');
end

if numel(linkEpochs) ~= Ntime
    error('Geometry and link-budget epoch counts do not match.');
end

for k = 1:Ntime
    if linkEpochs(k).TimeOffset_s ~= timeOffset_s(k)
        error('Geometry/link-budget time mismatch at epoch %d.',k);
    end
end

%% ========================================================================
% Selection configuration
% =========================================================================

target_m = 0.6;
targetTol = 1e-12;

fprintf('\n============================================================\n');
fprintf('TASK CODE 1 - STEP 3: FTS vs BTS AT FOUR EPOCHS\n');
fprintf('============================================================\n');
fprintf('Epochs        : %s s\n',mat2str(timeOffset_s.'));
fprintf('Candidate IDs : %d\n',N);
fprintf('Target        : B_CRB <= %.3f m\n',target_m);
fprintf('============================================================\n');

%% ========================================================================
% Allocate results
% =========================================================================

VisibleCount = zeros(Ntime,1);

NoSelection_BCRB_m = nan(Ntime,1);
NoSelection_PDOP = nan(Ntime,1);
NoSelectionRank = zeros(Ntime,1);
NoSelectionFeasible = false(Ntime,1);

FTS_K = zeros(Ntime,1);
FTS_BCRB_m = nan(Ntime,1);
FTS_PDOP = nan(Ntime,1);
FTS_Feasible = false(Ntime,1);
FTS_Evaluations = zeros(Ntime,1);
FTS_Selected = strings(Ntime,1);

BTS_K = zeros(Ntime,1);
BTS_BCRB_m = nan(Ntime,1);
BTS_PDOP = nan(Ntime,1);
BTS_Feasible = false(Ntime,1);
BTS_Evaluations = zeros(Ntime,1);
BTS_Selected = strings(Ntime,1);

SameK = false(Ntime,1);
SameSubset = false(Ntime,1);

FTS_HAPS = zeros(Ntime,1);
FTS_LEO  = zeros(Ntime,1);
FTS_MEO  = zeros(Ntime,1);
FTS_GEO  = zeros(Ntime,1);

BTS_HAPS = zeros(Ntime,1);
BTS_LEO  = zeros(Ntime,1);
BTS_MEO  = zeros(Ntime,1);
BTS_GEO  = zeros(Ntime,1);

SelectedFTS = false(N,Ntime);
SelectedBTS = false(N,Ntime);

%% ========================================================================
% Run No Selection, FTS, and BTS independently at each epoch
% =========================================================================

for k = 1:Ntime

    Hfull = double(epochs(k).H);
    visible = logical(epochs(k).Visible(:));
    sigmaRho_m = double(linkEpochs(k).SigmaRho_m(:));

    if size(Hfull,1) ~= N || size(Hfull,2) ~= 4
        error('Invalid H dimensions at t=%g s.',timeOffset_s(k));
    end

    if numel(sigmaRho_m) ~= N
        error('Invalid sigma_rho dimensions at t=%g s.',timeOffset_s(k));
    end

    if any(~isfinite(sigmaRho_m(visible))) || ...
       any(sigmaRho_m(visible) <= 0)
        error('Invalid visible sigma_rho at t=%g s.',timeOffset_s(k));
    end

    visibleIdx = find(visible).';
    VisibleCount(k) = numel(visibleIdx);

    % --------------------------------------------------------------------
    % No Selection
    % --------------------------------------------------------------------

    [NoSelection_BCRB_m(k),NoSelection_PDOP(k),NoSelectionRank(k)] = ...
        subsetMetrics(visibleIdx,Hfull,sigmaRho_m);

    NoSelectionFeasible(k) = ...
        NoSelectionRank(k) == 4 && ...
        NoSelection_BCRB_m(k) <= target_m + targetTol;

    % --------------------------------------------------------------------
    % Forward Transmitter Selection
    % --------------------------------------------------------------------

    [ftsIdx, ...
     FTS_BCRB_m(k), ...
     FTS_PDOP(k), ...
     FTS_Feasible(k), ...
     FTS_Evaluations(k)] = ...
        forwardTransmitterSelection( ...
            visibleIdx,Hfull,sigmaRho_m,target_m,targetTol);

    FTS_K(k) = numel(ftsIdx);
    SelectedFTS(ftsIdx,k) = true;
    FTS_Selected(k) = strjoin(Label(ftsIdx),';');

    [FTS_HAPS(k),FTS_LEO(k),FTS_MEO(k),FTS_GEO(k)] = ...
        architectureCounts(Architecture,ftsIdx);

    % --------------------------------------------------------------------
    % Backward Transmitter Selection
    % --------------------------------------------------------------------

    [btsIdx, ...
     BTS_BCRB_m(k), ...
     BTS_PDOP(k), ...
     BTS_Feasible(k), ...
     BTS_Evaluations(k)] = ...
        backwardTransmitterSelection( ...
            visibleIdx,Hfull,sigmaRho_m,target_m,targetTol);

    BTS_K(k) = numel(btsIdx);
    SelectedBTS(btsIdx,k) = true;
    BTS_Selected(k) = strjoin(Label(btsIdx),';');

    [BTS_HAPS(k),BTS_LEO(k),BTS_MEO(k),BTS_GEO(k)] = ...
        architectureCounts(Architecture,btsIdx);

    % --------------------------------------------------------------------
    % Direct comparison
    % --------------------------------------------------------------------

    SameK(k) = FTS_K(k) == BTS_K(k);
    SameSubset(k) = isequal(sort(ftsIdx(:)),sort(btsIdx(:)));

    fprintf('\n--- t = %3d s ---\n',timeOffset_s(k));
    fprintf('Visible       : %d / %d\n',VisibleCount(k),N);
    fprintf('No Selection  : K=%2d, B_CRB=%.6f m, PDOP=%.6f, feasible=%d\n', ...
        VisibleCount(k), ...
        NoSelection_BCRB_m(k), ...
        NoSelection_PDOP(k), ...
        NoSelectionFeasible(k));

    fprintf('FTS           : K=%2d, B_CRB=%.6f m, PDOP=%.6f, evals=%d, feasible=%d\n', ...
        FTS_K(k), ...
        FTS_BCRB_m(k), ...
        FTS_PDOP(k), ...
        FTS_Evaluations(k), ...
        FTS_Feasible(k));

    fprintf('              : H/L/M/G = %d/%d/%d/%d\n', ...
        FTS_HAPS(k),FTS_LEO(k),FTS_MEO(k),FTS_GEO(k));

    fprintf('BTS           : K=%2d, B_CRB=%.6f m, PDOP=%.6f, evals=%d, feasible=%d\n', ...
        BTS_K(k), ...
        BTS_BCRB_m(k), ...
        BTS_PDOP(k), ...
        BTS_Evaluations(k), ...
        BTS_Feasible(k));

    fprintf('              : H/L/M/G = %d/%d/%d/%d\n', ...
        BTS_HAPS(k),BTS_LEO(k),BTS_MEO(k),BTS_GEO(k));

    fprintf('Same K        : %d\n',SameK(k));
    fprintf('Same subset   : %d\n',SameSubset(k));

    fprintf('FTS set       : %s\n',FTS_Selected(k));
    fprintf('BTS set       : %s\n',BTS_Selected(k));
end

%% ========================================================================
% Validate nominal t = 0 against existing BTS result, when available
% =========================================================================

fprintf('\n============================================================\n');
fprintf('t = 0 NOMINAL VALIDATION\n');
fprintf('============================================================\n');

if isfile(beSelectedCsv)

    Tbe = readtable(beSelectedCsv,'TextType','string');

    if ~ismember("Label",string(Tbe.Properties.VariableNames))
        warning('BE_selected.csv has no Label column; subset validation skipped.');
    else
        nominalBTSLabels = string(Tbe.Label(:));
        nominalBTSIdx = find(ismember(Label,nominalBTSLabels));

        if isequal(sort(find(SelectedBTS(:,1))),sort(nominalBTSIdx(:)))
            fprintf('BTS subset            : MATCHES nominal BE_selected.csv\n');
        else
            error('t=0 BTS subset does not match nominal BE_selected.csv.');
        end
    end
else
    fprintf('BE_selected.csv       : not found; subset validation skipped.\n');
end

if isfile(beResultsMat)

    BE = load(beResultsMat);

    expectedBound = NaN;

    if isfield(BE,'finalBound')
        expectedBound = double(BE.finalBound);
    elseif isfield(BE,'finalBound_m')
        expectedBound = double(BE.finalBound_m);
    end

    if isfinite(expectedBound)

        deltaBound = abs(BTS_BCRB_m(1)-expectedBound);

        fprintf('Nominal BTS B_CRB     : %.12f m\n',expectedBound);
        fprintf('Temporal t=0 BTS      : %.12f m\n',BTS_BCRB_m(1));
        fprintf('|Delta|               : %.3e m\n',deltaBound);

        if deltaBound > 1e-9
            error('t=0 BTS bound does not reproduce nominal BE result.');
        end
    end
else
    fprintf('BE_results.mat        : not found; bound validation skipped.\n');
end

% These two counts are the current paper convention for N=26, K=12.
if FTS_K(1) == 12 && FTS_Evaluations(1) ~= 15098
    warning('t=0 FTS evaluation count is %d, expected 15098.', ...
        FTS_Evaluations(1));
end

if BTS_K(1) == 12 && BTS_Evaluations(1) ~= 285
    warning('t=0 BTS evaluation count is %d, expected 285.', ...
        BTS_Evaluations(1));
end

fprintf('FTS/BTS same subset   : %d\n',SameSubset(1));
fprintf('FTS evaluations       : %d\n',FTS_Evaluations(1));
fprintf('BTS evaluations       : %d\n',BTS_Evaluations(1));

%% ========================================================================
% Save numerical outputs
% =========================================================================

Tsummary = table( ...
    timeOffset_s, ...
    timeUTC, ...
    VisibleCount, ...
    NoSelection_BCRB_m, ...
    NoSelection_PDOP, ...
    NoSelectionFeasible, ...
    FTS_K, ...
    FTS_BCRB_m, ...
    FTS_PDOP, ...
    FTS_Feasible, ...
    FTS_Evaluations, ...
    FTS_HAPS, ...
    FTS_LEO, ...
    FTS_MEO, ...
    FTS_GEO, ...
    BTS_K, ...
    BTS_BCRB_m, ...
    BTS_PDOP, ...
    BTS_Feasible, ...
    BTS_Evaluations, ...
    BTS_HAPS, ...
    BTS_LEO, ...
    BTS_MEO, ...
    BTS_GEO, ...
    SameK, ...
    SameSubset, ...
    'VariableNames',{ ...
    'TimeOffset_s', ...
    'TimeUTC', ...
    'VisibleCount', ...
    'NoSelection_BCRB_m', ...
    'NoSelection_PDOP', ...
    'NoSelection_Feasible', ...
    'FTS_K', ...
    'FTS_BCRB_m', ...
    'FTS_PDOP', ...
    'FTS_Feasible', ...
    'FTS_Evaluations', ...
    'FTS_HAPS', ...
    'FTS_LEO', ...
    'FTS_MEO', ...
    'FTS_GEO', ...
    'BTS_K', ...
    'BTS_BCRB_m', ...
    'BTS_PDOP', ...
    'BTS_Feasible', ...
    'BTS_Evaluations', ...
    'BTS_HAPS', ...
    'BTS_LEO', ...
    'BTS_MEO', ...
    'BTS_GEO', ...
    'SameK', ...
    'SameSubset'});

Tsets = table( ...
    timeOffset_s, ...
    timeUTC, ...
    FTS_Selected, ...
    BTS_Selected, ...
    'VariableNames',{ ...
    'TimeOffset_s', ...
    'TimeUTC', ...
    'FTS_Selected', ...
    'BTS_Selected'});

writetable( ...
    Tsummary, ...
    fullfile(outDir,'temporal_fts_bts_summary.csv'));

writetable( ...
    Tsets, ...
    fullfile(outDir,'temporal_fts_bts_selected_sets.csv'));

matFile = fullfile( ...
    outDir, ...
    'scenario4_temporal_fts_bts_4epochs.mat');

save( ...
    matFile, ...
    'Tsummary', ...
    'Tsets', ...
    'SelectedFTS', ...
    'SelectedBTS', ...
    'timeOffset_s', ...
    'timeUTC', ...
    'target_m', ...
    'PoolIndex', ...
    'Label', ...
    'Architecture', ...
    'NORAD_CAT_ID', ...
    '-v7.3');

%% ========================================================================
% Figure 1 - B_CRB versus time
% =========================================================================

makeTemporalBoundFigure( ...
    timeOffset_s, ...
    NoSelection_BCRB_m, ...
    FTS_BCRB_m, ...
    BTS_BCRB_m, ...
    target_m, ...
    outDir);

%% ========================================================================
% Figures 2-5 - skyplots at 0, 30, 60, 90 s
% =========================================================================

for k = 1:Ntime

    makeSelectionSkyplot( ...
        epochs(k).Azimuth_deg(:), ...
        epochs(k).Elevation_deg(:), ...
        Architecture, ...
        epochs(k).Visible(:), ...
        SelectedFTS(:,k), ...
        SelectedBTS(:,k), ...
        timeOffset_s(k), ...
        FTS_K(k), ...
        FTS_BCRB_m(k), ...
        BTS_K(k), ...
        BTS_BCRB_m(k), ...
        outDir);
end

%% ========================================================================
% Final console summary
% =========================================================================

fprintf('\n============================================================\n');
fprintf('FTS/BTS FOUR-EPOCH COMPARISON COMPLETED\n');
fprintf('============================================================\n');

disp(Tsummary);

fprintf('All epochs FTS feasible : %s\n',string(all(FTS_Feasible)));
fprintf('All epochs BTS feasible : %s\n',string(all(BTS_Feasible)));
fprintf('Same K at all epochs    : %s\n',string(all(SameK)));
fprintf('Same subset all epochs  : %s\n',string(all(SameSubset)));

fprintf('\nResults folder:\n  %s\n',outDir);
fprintf('MAT file:\n  %s\n',matFile);
fprintf('============================================================\n');

%% ========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function [selected,finalBound_m,finalPDOP,targetReached,evaluationCount] = ...
    forwardTransmitterSelection( ...
        initialIdx,Hfull,sigmaRho_m,target_m,targetTol)

    initialIdx = initialIdx(:).';

    selected = [];
    finalBound_m = inf;
    finalPDOP = inf;
    targetReached = false;
    evaluationCount = 0;

    if numel(initialIdx) < 4
        return;
    end

    % --------------------------------------------------------------------
    % Exhaustive initialization over all full-rank 4-transmitter subsets.
    % --------------------------------------------------------------------

    combinations4 = nchoosek(initialIdx,4);

    bestBound_m = inf;
    bestSet = [];

    for row = 1:size(combinations4,1)

        trialIdx = combinations4(row,:);
        evaluationCount = evaluationCount + 1;

        [trialBound_m,~,trialRank] = ...
            subsetMetrics(trialIdx,Hfull,sigmaRho_m);

        if trialRank ~= 4
            continue;
        end

        if trialBound_m < bestBound_m-1e-14

            bestBound_m = trialBound_m;
            bestSet = trialIdx;

        elseif abs(trialBound_m-bestBound_m) <= 1e-14

            if isempty(bestSet) || ...
                    lexicographicallySmaller(trialIdx,bestSet)
                bestSet = trialIdx;
            end
        end
    end

    if isempty(bestSet)
        return;
    end

    selected = bestSet(:).';
    currentBound_m = bestBound_m;

    % --------------------------------------------------------------------
    % Greedy additions until the target is met.
    % --------------------------------------------------------------------

    while currentBound_m > target_m+targetTol && ...
            numel(selected) < numel(initialIdx)

        remaining = initialIdx(~ismember(initialIdx,selected));

        bestAdded = NaN;
        bestTrialBound_m = inf;

        for jj = 1:numel(remaining)

            candidateToAdd = remaining(jj);
            trialIdx = [selected candidateToAdd];

            evaluationCount = evaluationCount + 1;

            [trialBound_m,~,trialRank] = ...
                subsetMetrics(trialIdx,Hfull,sigmaRho_m);

            if trialRank ~= 4
                continue;
            end

            if trialBound_m < bestTrialBound_m-1e-14

                bestTrialBound_m = trialBound_m;
                bestAdded = candidateToAdd;

            elseif abs(trialBound_m-bestTrialBound_m) <= 1e-14

                if isnan(bestAdded) || candidateToAdd < bestAdded
                    bestAdded = candidateToAdd;
                end
            end
        end

        if isnan(bestAdded)
            break;
        end

        selected(end+1) = bestAdded; %#ok<AGROW>
        currentBound_m = bestTrialBound_m;
    end

    [finalBound_m,finalPDOP,finalRank] = ...
        subsetMetrics(selected,Hfull,sigmaRho_m);

    targetReached = ...
        finalRank == 4 && ...
        finalBound_m <= target_m+targetTol;
end

function [selected,finalBound_m,finalPDOP,targetReached,evaluationCount] = ...
    backwardTransmitterSelection( ...
        initialIdx,Hfull,sigmaRho_m,target_m,targetTol)

    selected = initialIdx(:).';
    evaluationCount = 0;

    [currentBound_m,currentPDOP,currentRank] = ...
        subsetMetrics(selected,Hfull,sigmaRho_m);

    % Initial No Selection feasibility check is intentionally NOT counted,
    % matching the current paper convention for BTS assessments.
    if currentRank ~= 4 || currentBound_m > target_m+targetTol

        finalBound_m = currentBound_m;
        finalPDOP = currentPDOP;
        targetReached = false;
        return;
    end

    while numel(selected) > 4

        bestRemoved = NaN;
        bestTrialBound_m = inf;

        for jj = 1:numel(selected)

            candidateToRemove = selected(jj);
            trialIdx = selected(selected ~= candidateToRemove);

            evaluationCount = evaluationCount + 1;

            [trialBound_m,~,trialRank] = ...
                subsetMetrics(trialIdx,Hfull,sigmaRho_m);

            if trialRank ~= 4 || ...
                    trialBound_m > target_m+targetTol
                continue;
            end

            if trialBound_m < bestTrialBound_m-1e-14

                bestTrialBound_m = trialBound_m;
                bestRemoved = candidateToRemove;

            elseif abs(trialBound_m-bestTrialBound_m) <= 1e-14

                if isnan(bestRemoved) || ...
                        candidateToRemove < bestRemoved
                    bestRemoved = candidateToRemove;
                end
            end
        end

        if isnan(bestRemoved)
            break;
        end

        selected(selected == bestRemoved) = [];
        currentBound_m = bestTrialBound_m;
    end

    [finalBound_m,finalPDOP,finalRank] = ...
        subsetMetrics(selected,Hfull,sigmaRho_m);

    targetReached = ...
        finalRank == 4 && ...
        finalBound_m <= target_m+targetTol;
end

function [bound_m,PDOP,rankH] = ...
    subsetMetrics(idx,Hfull,sigmaRho_m)

    idx = idx(:);

    H = Hfull(idx,:);
    sigma = sigmaRho_m(idx);

    rankH = rank(H);
    bound_m = inf;
    PDOP = inf;

    if rankH < 4 || ...
       any(~isfinite(sigma)) || ...
       any(sigma <= 0)
        return;
    end

    % PDOP
    normalMatrix = H.'*H;

    if rcond(normalMatrix) > 1e-14

        covarianceGeometry = normalMatrix\eye(4);
        traceGeometry = trace(covarianceGeometry(1:3,1:3));

        if traceGeometry > 0 && isfinite(traceGeometry)
            PDOP = sqrt(traceGeometry);
        end
    end

    % Fisher information matrix:
    % J = H' R^{-1} H, with diagonal R.
    weights = 1./(sigma.^2);
    J = H.'*(H.*weights);

    if rcond(J) <= 1e-14
        return;
    end

    covarianceBound = J\eye(4);
    tracePosition = trace(covarianceBound(1:3,1:3));

    if tracePosition > 0 && isfinite(tracePosition)
        bound_m = sqrt(tracePosition);
    end
end

function tf = lexicographicallySmaller(a,b)

    a = sort(a(:).');
    b = sort(b(:).');

    firstDifferent = find(a~=b,1,'first');

    if isempty(firstDifferent)
        tf = false;
    else
        tf = a(firstDifferent) < b(firstDifferent);
    end
end

function [nHAPS,nLEO,nMEO,nGEO] = ...
    architectureCounts(architecture,idx)

    arch = architecture(idx);

    nHAPS = sum(arch=="HAPS");
    nLEO  = sum(arch=="LEO");
    nMEO  = sum(arch=="MEO");
    nGEO  = sum(arch=="GEO");
end

function makeTemporalBoundFigure( ...
    timeOffset_s, ...
    noSelectionBound_m, ...
    ftsBound_m, ...
    btsBound_m, ...
    target_m, ...
    outDir)

    fig = figure( ...
        'Name','Temporal FTS-BTS B_CRB', ...
        'NumberTitle','off', ...
        'Color','w', ...
        'Units','inches', ...
        'Position',[1 1 5.4 3.3]);

    ax = axes('Parent',fig);
    hold(ax,'on');
    grid(ax,'on');
    box(ax,'on');

    plot( ...
        ax,timeOffset_s,noSelectionBound_m, ...
        '-o', ...
        'LineWidth',1.4, ...
        'MarkerSize',5, ...
        'DisplayName','No Selection');

    plot( ...
        ax,timeOffset_s,ftsBound_m, ...
        '-s', ...
        'LineWidth',1.4, ...
        'MarkerSize',5, ...
        'DisplayName','FTS');

    plot( ...
        ax,timeOffset_s,btsBound_m, ...
        '-d', ...
        'LineWidth',1.4, ...
        'MarkerSize',5, ...
        'DisplayName','BTS');

    yline( ...
        ax,target_m,'--', ...
        'LineWidth',1.1, ...
        'DisplayName','Target');

    xlabel(ax,'Time offset [s]');
    ylabel(ax,'B_{CRB} [m]');

    xlim(ax,[min(timeOffset_s) max(timeOffset_s)]);
    xticks(ax,timeOffset_s);

    set( ...
        ax, ...
        'FontName','Times New Roman', ...
        'FontSize',10, ...
        'LineWidth',0.7);

    legend( ...
        ax, ...
        'Location','best', ...
        'Box','off', ...
        'FontName','Times New Roman', ...
        'FontSize',9);

    exportgraphics( ...
        fig, ...
        fullfile(outDir,'temporal_bcrb_fts_bts.pdf'), ...
        'ContentType','vector');

    exportgraphics( ...
        fig, ...
        fullfile(outDir,'temporal_bcrb_fts_bts.png'), ...
        'Resolution',300);

    savefig( ...
        fig, ...
        fullfile(outDir,'temporal_bcrb_fts_bts.fig'));
end

function makeSelectionSkyplot( ...
    azimuth_deg, ...
    elevation_deg, ...
    architecture, ...
    visible, ...
    selectedFTS, ...
    selectedBTS, ...
    timeOffset_s, ...
    ftsK, ...
    ftsBound_m, ...
    btsK, ...
    btsBound_m, ...
    outDir)

    architecture = string(architecture(:));
    visible = logical(visible(:));
    selectedFTS = logical(selectedFTS(:));
    selectedBTS = logical(selectedBTS(:));

    archList = ["HAPS","LEO","MEO","GEO"];

    % Architecture colors only. Selection status is encoded by marker shape
    % and fill. No "x" markers are used.
    archColors = [ ...
        0.0000 0.4470 0.7410
        0.8500 0.3250 0.0980
        0.9290 0.6940 0.1250
        0.4940 0.1840 0.5560];

    theta = deg2rad(mod(azimuth_deg,360));
    rho = 90-elevation_deg;

    selectedBoth = visible & selectedFTS & selectedBTS;
    selectedFTSOnly = visible & selectedFTS & ~selectedBTS;
    selectedBTSOnly = visible & ~selectedFTS & selectedBTS;
    selectedNeither = visible & ~selectedFTS & ~selectedBTS;

    fig = figure( ...
        'Name',sprintf('Skyplot t=%d s',timeOffset_s), ...
        'NumberTitle','off', ...
        'Color','w', ...
        'Units','inches', ...
        'Position',[1 1 4.3 4.0]);

    pax = polaraxes('Parent',fig);
    hold(pax,'on');

    % Plot each architecture while preserving selection status.
    for aa = 1:numel(archList)

        idxArch = architecture==archList(aa);

        idx = idxArch & selectedNeither;
        if any(idx)
            polarscatter( ...
                pax,theta(idx),rho(idx),32, ...
                archColors(aa,:), ...
                'o', ...
                'LineWidth',1.0, ...
                'HandleVisibility','off');
        end

        idx = idxArch & selectedBoth;
        if any(idx)
            polarscatter( ...
                pax,theta(idx),rho(idx),42, ...
                archColors(aa,:), ...
                'o', ...
                'filled', ...
                'MarkerEdgeColor','k', ...
                'LineWidth',0.7, ...
                'HandleVisibility','off');
        end

        idx = idxArch & selectedFTSOnly;
        if any(idx)
            polarscatter( ...
                pax,theta(idx),rho(idx),48, ...
                archColors(aa,:), ...
                '^', ...
                'filled', ...
                'MarkerEdgeColor','k', ...
                'LineWidth',0.7, ...
                'HandleVisibility','off');
        end

        idx = idxArch & selectedBTSOnly;
        if any(idx)
            polarscatter( ...
                pax,theta(idx),rho(idx),48, ...
                archColors(aa,:), ...
                's', ...
                'filled', ...
                'MarkerEdgeColor','k', ...
                'LineWidth',0.7, ...
                'HandleVisibility','off');
        end
    end

    set( ...
        pax, ...
        'ThetaZeroLocation','top', ...
        'ThetaDir','clockwise', ...
        'ThetaTick',[0 90 180 270], ...
        'ThetaTickLabel',{'N','E','S','W'}, ...
        'RLim',[0 90], ...
        'RTick',[30 60 90], ...
        'RTickLabel',{'60^\circ','30^\circ','0^\circ'}, ...
        'RAxisLocation',245, ...
        'FontName','Times New Roman', ...
        'FontSize',9, ...
        'GridAlpha',0.22);

    title( ...
        pax, ...
        sprintf(['t = %d s   |   FTS: K=%d, B_{CRB}=%.3f m   |   ' ...
                 'BTS: K=%d, B_{CRB}=%.3f m'], ...
        timeOffset_s,ftsK,ftsBound_m,btsK,btsBound_m), ...
        'FontName','Times New Roman', ...
        'FontSize',9, ...
        'FontWeight','normal');

    % --------------------------------------------------------------------
    % Compact legend: architectures + selection status
    % --------------------------------------------------------------------

    archHandles = gobjects(4,1);

    for aa = 1:4
        archHandles(aa) = polarplot( ...
            pax,nan,nan,'o', ...
            'LineStyle','none', ...
            'MarkerEdgeColor',archColors(aa,:), ...
            'MarkerFaceColor',archColors(aa,:), ...
            'MarkerSize',5);
    end

    hNot = polarplot( ...
        pax,nan,nan,'ko', ...
        'LineStyle','none', ...
        'MarkerFaceColor','none', ...
        'MarkerSize',5, ...
        'LineWidth',0.9);

    hBoth = polarplot( ...
        pax,nan,nan,'ko', ...
        'LineStyle','none', ...
        'MarkerFaceColor','k', ...
        'MarkerSize',5);

    hFTS = polarplot( ...
        pax,nan,nan,'k^', ...
        'LineStyle','none', ...
        'MarkerFaceColor','k', ...
        'MarkerSize',5);

    hBTS = polarplot( ...
        pax,nan,nan,'ks', ...
        'LineStyle','none', ...
        'MarkerFaceColor','k', ...
        'MarkerSize',5);

    lg = legend( ...
        pax, ...
        [archHandles(:);hNot;hBoth;hFTS;hBTS], ...
        {'HAPS','LEO','MEO','GEO', ...
         'Not selected','FTS & BTS','FTS only','BTS only'}, ...
        'NumColumns',2, ...
        'Location','southoutside', ...
        'Box','off', ...
        'FontName','Times New Roman', ...
        'FontSize',8);

    lg.ItemTokenSize = [12 8];

    baseName = sprintf('skyplot_t%03ds',timeOffset_s);

    exportgraphics( ...
        fig, ...
        fullfile(outDir,[baseName '.pdf']), ...
        'ContentType','vector');

    exportgraphics( ...
        fig, ...
        fullfile(outDir,[baseName '.png']), ...
        'Resolution',300);

    savefig( ...
        fig, ...
        fullfile(outDir,[baseName '.fig']));
end
