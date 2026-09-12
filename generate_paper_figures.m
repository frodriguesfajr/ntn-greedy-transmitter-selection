%% Paper figures - nominal geometry, availability, robustness, and short-term sensitivity
%
% Objective:
%   Generate the numerical short-term sensitivity result and the separate
%   figures used in the manuscript.
%
% The script produces six independent paper figures:
%   1) nominal selection skyplot;
%   2) 3-minute temporal sensitivity;
%   3) HAPS availability;
%   4) LEO compensation;
%   5) Gaussian-impulsive robustness;
%   6) fixed-size transmitter-selection necessity.
%
% The figures are intentionally exported as separate files so that their
% placement, side-by-side arrangement, and subfigure labels can be
% controlled directly in LaTeX.
%
% Inputs:
%   results/candidate_pool_26.csv
%   results/link_budget.csv
%   results/BE_selected.csv
%   results/BE_results.mat
%   results/HAPS_restrictions_summary.csv
%   results/LEO_availability_summary.csv
%   results/MC_outliers/MC_outliers_summary.csv
%   results/selection_necessity/selection_necessity_methods.csv
%   results/selection_necessity/selection_necessity_random_subsets.csv
%   tle/<architecture>/<NORAD>/<NORAD>.tle
%
% Orbital positions in the 180-s temporal experiment are propagated from
% the archived TLEs using the included Vallado SGP4 implementation with
% WGS-72 constants.
%
% The 180-s interval is intentionally short so that the experiment isolates
% geometry and link-quality evolution while retaining the same nominal
% transmitter identities.
%
% Numerical outputs:
%   results/paper_figures/temporal_sensitivity.csv
%   results/paper_figures/temporal_summary.csv
%
% Figure outputs:
%   results/paper_figures/selection_visibility.pdf
%   results/paper_figures/time_sensitivity.pdf
%   results/paper_figures/haps_availability.pdf
%   results/paper_figures/leo_compensation.pdf
%   results/paper_figures/robustness.pdf
%   results/paper_figures/selection_necessity.pdf
%
close all;
clear;
clc;
format long;

%% ##################### Paths and configuration #########################

rootDir = fileparts(mfilename('fullpath'));

if strlength(string(rootDir)) == 0
    rootDir = pwd;
end

cd(rootDir);
setup_paths;

cfg = paperFigureConfig(rootDir);

candidateCsv = fullfile(rootDir,"results","candidate_pool_26.csv");

linkBudgetCsv = fullfile(rootDir,"results","link_budget.csv");

beSelectedCsv = fullfile(rootDir,"results","BE_selected.csv");

beResultsMat = fullfile(rootDir,"results","BE_results.mat");

hapsAvailabilityCsv = fullfile( ...
    rootDir,"results","HAPS_restrictions_summary.csv");

leoCompensationCsv = fullfile( ...
    rootDir,"results","LEO_availability_summary.csv");

outlierResultsCsv = fullfile( ...
    rootDir,"results","MC_outliers","MC_outliers_summary.csv");

selectionMethodsCsv = fullfile( ...
    rootDir,"results","selection_necessity", ...
    "selection_necessity_methods.csv");

selectionRandomCsv = fullfile( ...
    rootDir,"results","selection_necessity", ...
    "selection_necessity_random_subsets.csv");

outDir = fullfile(rootDir,"results","paper_figures");

inputFiles = string({ ...
    candidateCsv,linkBudgetCsv,beSelectedCsv,beResultsMat, ...
    hapsAvailabilityCsv,leoCompensationCsv,outlierResultsCsv, ...
    selectionMethodsCsv,selectionRandomCsv});

for k = 1:numel(inputFiles)
    if ~isfile(inputFiles(k))
        error("Required input file not found: %s",inputFiles(k));
    end
end

if ~exist(outDir,'dir')
    mkdir(outDir);
end

%% ##################### Load input results ###############################

Tpool = readtable(candidateCsv,'TextType','string');
Tlink = readtable(linkBudgetCsv,'TextType','string');
Tbe = readtable(beSelectedCsv,'TextType','string');
BE = load(beResultsMat);
Thaps = readtable(hapsAvailabilityCsv,'TextType','string');
Tleo = readtable(leoCompensationCsv,'TextType','string');
Toutlier = readtable(outlierResultsCsv,'TextType','string');
TselectionMethods = readtable( ...
    selectionMethodsCsv, ...
    'Delimiter',',', ...
    'TextType','string');

TselectionRandom = readtable( ...
    selectionRandomCsv, ...
    'Delimiter',',', ...
    'TextType','string');

% Normalize candidate-pool NORAD column if needed.
if ~ismember("NORAD_CAT_ID",string(Tpool.Properties.VariableNames))
    if ismember("NORAD",string(Tpool.Properties.VariableNames))
        Tpool.NORAD_CAT_ID = string(Tpool.NORAD);
    else
        error("Candidate pool has no NORAD identifier column.");
    end
end

requiredPool = [ ...
    "PoolIndex","Label","Architecture","NORAD_CAT_ID", ...
    "X_ECEF_m","Y_ECEF_m","Z_ECEF_m", ...
    "Azimuth_deg","Elevation_deg"];

requiredLink = [ ...
    "Label","Architecture","CN0_dBHz", ...
    "AtmosphericGasLoss_dB","Elevation_deg","SigmaRho_m"];

requiredBE = "Label";

requiredHaps = [ ...
    "MaxHAPS","Kfinal","FinalBound_m","TargetReached", ...
    "HAPS","LEO","MEO","GEO"];

requiredLeo = [ ...
    "LEOAvailable","AllCandidateBound_m","FinalBound_m", ...
    "TargetReached"];

requiredOutlier = [ ...
    "B","RMSE_FullWLS_m","RMSE_RobustBE_m", ...
    "ExactDetectionProb","RMSE_BE_ExactDetection_m"];

requiredSelectionMethods = [ ...
    "Method","K","PositionalBound_m","TargetReached"];

requiredSelectionRandom = [ ...
    "TrialID","PositionalBound_m","RankH","TargetReached"];

assertTableVariables(Tpool,requiredPool,"candidate-pool table");
assertTableVariables(Tlink,requiredLink,"link-budget table");
assertTableVariables(Tbe,requiredBE,"BE selected-set table");
assertTableVariables(Thaps,requiredHaps,"HAPS-availability table");
assertTableVariables(Tleo,requiredLeo,"LEO-compensation table");
assertTableVariables(Toutlier,requiredOutlier,"outlier-results table");
assertTableVariables( ...
    TselectionMethods,requiredSelectionMethods, ...
    "selection-necessity methods table");
assertTableVariables( ...
    TselectionRandom,requiredSelectionRandom, ...
    "selection-necessity random-subset table");

Tpool = sortrows(Tpool,"PoolIndex");

if ismember("PoolIndex",string(Tlink.Properties.VariableNames))
    Tlink = sortrows(Tlink,"PoolIndex");
end

Label = string(Tpool.Label(:));
Architecture = string(Tpool.Architecture(:));
Npool = height(Tpool);

if Npool ~= 26
    error("Expected 26 nominal candidates; found %d.",Npool);
end

if height(Tlink) ~= Npool || any(string(Tlink.Label(:)) ~= Label)
    error("Link-budget rows do not match the candidate-pool ordering.");
end

selectedT0 = find(ismember(Label,string(Tbe.Label(:))));

if numel(selectedT0) ~= 12
    error("Expected 12 transmitters in the nominal BE subset.");
end

if isfield(BE,'finalBound')
    expectedBoundT0_m = double(BE.finalBound);
elseif isfield(BE,'finalBound_m')
    expectedBoundT0_m = double(BE.finalBound_m);
else
    error("BE_results.mat does not contain the expected final-bound variable.");
end

target_m = double(cfg.selectionTarget_m);
targetTol = 1e-12;

%% ##################### Temporal window #################################
%
% Ten-second samples over 180 s. This is a short-term sensitivity study,
% not a constellation-arrival/departure experiment.

timeOffset_s = (0:10:180).';
timeUTC = cfg.analysisTimeUTC + seconds(timeOffset_s);
Ntime = numel(timeOffset_s);

PositionECEF = nan(Npool,3,Ntime);

eop = paper_epoch_eop();

fprintf('\n============================================================\n');
fprintf('PAPER FIGURES - SHORT-TERM TEMPORAL PROPAGATION\n');
fprintf('============================================================\n');

for i = 1:Npool

    nominalXYZ = [ ...
        double(Tpool.X_ECEF_m(i)), ...
        double(Tpool.Y_ECEF_m(i)), ...
        double(Tpool.Z_ECEF_m(i))];

    if Architecture(i) == "HAPS"

        positionHistory = repmat(nominalXYZ,Ntime,1);

    else

        norad = strip(string(Tpool.NORAD_CAT_ID(i)));

        if ismissing(norad) || strlength(norad) == 0
            error("Missing NORAD ID for %s.",Label(i));
        end

        tleFile = fullfile( ...
            cfg.paths.tleRoot,Architecture(i),norad,norad+".tle");

        if ~isfile(tleFile)
            error("Archived TLE missing for %s: %s",Label(i),tleFile);
        end

        positionHistory = nan(Ntime,3);

        fprintf('Propagating %-6s (%s) with Vallado SGP4...\n', ...
            Label(i),Architecture(i));

        for tidx = 1:Ntime

            state = propagate_tle_vallado(tleFile,timeUTC(tidx),eop);

            if isfield(state,'sgp4_error') && state.sgp4_error ~= 0
                error("SGP4 error %d for %s at %s.", ...
                    state.sgp4_error,Label(i),string(timeUTC(tidx)));
            end

            positionHistory(tidx,:) = state.rECEF_m(:).';
        end
    end

    PositionECEF(i,:,:) = permute(positionHistory,[3 2 1]);
end

%% ##################### Receiver and link model #########################

UserPositionECEF = geodeticToECEFWGS84( ...
    cfg.user.lat_deg,cfg.user.lon_deg,cfg.user.h_m);

positionAtT0Csv = [ ...
    double(Tpool.X_ECEF_m), ...
    double(Tpool.Y_ECEF_m), ...
    double(Tpool.Z_ECEF_m)];

positionAtT0E = PositionECEF(:,:,1);

positionAgreement_m = ...
    max(vecnorm(positionAtT0E-positionAtT0Csv,2,2));

if positionAgreement_m > 0.05
    error([ ...
        "Candidate-pool and Vallado positions differ by up to %.6f m " ...
        "at t0."],positionAgreement_m);
end

archList = ["HAPS";"LEO";"MEO";"GEO"];
zenithGasLoss_dB = nan(numel(archList),1);

for aidx = 1:numel(archList)

    idx = Architecture == archList(aidx);

    if ~any(idx)
        error("No %s candidates found in the nominal pool.",archList(aidx));
    end

    zenithSamples = ...
        double(Tlink.AtmosphericGasLoss_dB(idx)).* ...
        sind(double(Tlink.Elevation_deg(idx)));

    zenithGasLoss_dB(aidx) = median(zenithSamples);

    reconstructedGasLoss = ...
        zenithGasLoss_dB(aidx)./ ...
        sind(double(Tlink.Elevation_deg(idx)));

    gasMismatch_dB = max(abs(reconstructedGasLoss - ...
        double(Tlink.AtmosphericGasLoss_dB(idx))));

    if gasMismatch_dB > 1e-6
        warning(['Stored %s gas losses differ from the temporal secant ', ...
            'approximation by up to %.6f dB.'], ...
            archList(aidx),gasMismatch_dB);
    end
end

% Per-link numerical calibration makes t0 C/N0 exactly equal to the
% stored link-budget values. Dynamic changes still come from range,
% elevation, gaseous loss, and LEO elevation-dependent EIRP.
[~,~,~,rawCN0T0_dBHz] = epochGeometryAndQuality( ...
    positionAtT0E,Architecture,UserPositionECEF,cfg, ...
    archList,zenithGasLoss_dB,zeros(Npool,1));

cn0Calibration_dB = ...
    double(Tlink.CN0_dBHz(:))-rawCN0T0_dBHz;

% Recompute the nominal BE-subset bound directly from stored inputs.
[H0,sigmaRho0_m,visible0,~] = epochGeometryAndQuality( ...
    positionAtT0E,Architecture,UserPositionECEF,cfg, ...
    archList,zenithGasLoss_dB,cn0Calibration_dB);

if ~all(visible0(selectedT0))
    error("Nominal BE subset is not fully visible at t0.");
end

[expectedBoundT0_m,~,rankT0] = ...
    subsetMetrics(selectedT0,H0,sigmaRho0_m);

if rankT0 ~= 4
    error("Nominal BE subset is rank deficient at t0.");
end

%% ##################### Temporal BE analysis ############################

FullVisibleCount = zeros(Ntime,1);
FullBound_m = nan(Ntime,1);

AdaptiveFeasible = false(Ntime,1);
AdaptiveK = nan(Ntime,1);
AdaptiveBound_m = nan(Ntime,1);
AdaptiveLabels = strings(Ntime,1);
AdaptiveSetChanged = false(Ntime,1);

FixedSubsetVisible = false(Ntime,1);
FixedSubsetBound_m = nan(Ntime,1);
FixedSubsetTargetMet = false(Ntime,1);

previousAdaptive = [];

for tidx = 1:Ntime

    positionNow = PositionECEF(:,:,tidx);

    [Hfull,sigmaRho_m,visible,~] = epochGeometryAndQuality( ...
        positionNow,Architecture,UserPositionECEF,cfg, ...
        archList,zenithGasLoss_dB,cn0Calibration_dB);

    visibleIdx = find(visible);

    FullVisibleCount(tidx) = numel(visibleIdx);

    [FullBound_m(tidx),~,fullRank] = ...
        subsetMetrics(visibleIdx,Hfull,sigmaRho_m);

    if fullRank == 4 && FullBound_m(tidx) <= target_m + targetTol

        [adaptiveIdx,AdaptiveBound_m(tidx),AdaptiveFeasible(tidx)] = ...
            backwardElimination( ...
            visibleIdx,Hfull,sigmaRho_m,target_m,targetTol);

        AdaptiveK(tidx) = numel(adaptiveIdx);
        AdaptiveLabels(tidx) = strjoin(Label(adaptiveIdx),";");

        if tidx > 1
            AdaptiveSetChanged(tidx) = ...
                ~isequal(adaptiveIdx(:),previousAdaptive(:));
        end

        previousAdaptive = adaptiveIdx;
    end

    FixedSubsetVisible(tidx) = all(visible(selectedT0));

    if FixedSubsetVisible(tidx)

        FixedSubsetBound_m(tidx) = ...
            subsetMetrics(selectedT0,Hfull,sigmaRho_m);

        FixedSubsetTargetMet(tidx) = ...
            FixedSubsetBound_m(tidx) <= target_m + targetTol;
    end
end

if any(FullVisibleCount ~= Npool)
    error([ ...
        "The 180-s short-term window contains a visibility transition. " ...
        "This figure is intended to isolate geometry/link evolution."]);
end

if abs(AdaptiveBound_m(1)-expectedBoundT0_m) > 1e-9
    error([ ...
        "The t0 BE bound (%.12f m) differs from the nominal value " ...
        "(%.12f m)."], ...
        AdaptiveBound_m(1),expectedBoundT0_m);
end

adaptiveT0 = split(AdaptiveLabels(1),";");

if ~isequal(adaptiveT0(:),Label(selectedT0))
    error("Temporal BE does not reproduce the nominal t0 subset.");
end

firstFixedViolationIdx = find( ...
    FixedSubsetVisible & ~FixedSubsetTargetMet,1,'first');

if isempty(firstFixedViolationIdx)
    firstFixedViolation_s = NaN;
else
    firstFixedViolation_s = timeOffset_s(firstFixedViolationIdx);
end

firstAdaptiveChangeIdx = find(AdaptiveSetChanged,1,'first');

if isempty(firstAdaptiveChangeIdx)
    firstAdaptiveChange_s = NaN;
else
    firstAdaptiveChange_s = timeOffset_s(firstAdaptiveChangeIdx);
end

%% ##################### Save numerical results ##########################

Ttemporal = table( ...
    timeOffset_s,timeUTC,FullVisibleCount,FullBound_m, ...
    AdaptiveFeasible,AdaptiveK,AdaptiveBound_m,AdaptiveSetChanged, ...
    FixedSubsetVisible,FixedSubsetBound_m,FixedSubsetTargetMet, ...
    AdaptiveLabels, ...
    'VariableNames',{ ...
    'TimeOffset_s','TimeUTC','FullVisibleCount','FullBound_m', ...
    'AdaptiveFeasible','AdaptiveK','AdaptiveBound_m','AdaptiveSetChanged', ...
    'FixedSubsetVisible','FixedSubsetBound_m','FixedSubsetTargetMet', ...
    'AdaptiveLabels'});

Tsummary = table( ...
    Npool,numel(selectedT0),target_m, ...
    FullBound_m(1),AdaptiveBound_m(1), ...
    firstAdaptiveChange_s,firstFixedViolation_s, ...
    min(AdaptiveK),max(AdaptiveK),sum(AdaptiveSetChanged), ...
    positionAgreement_m,max(abs(cn0Calibration_dB)), ...
    'VariableNames',{ ...
    'NominalCandidates','NominalSelected','Target_m', ...
    'FullBoundAtT0_m','SelectedBoundAtT0_m', ...
    'FirstSampledAdaptiveSetChange_s', ...
    'FirstSampledFixedTargetViolation_s', ...
    'MinAdaptiveK','MaxAdaptiveK','NumberOfAdaptiveSetChanges', ...
    'MaxPositionDifferenceAtT0_m','MaxCN0Calibration_dB'});

writetable(Ttemporal, ...
    fullfile(outDir,"temporal_sensitivity.csv"));

writetable(Tsummary, ...
    fullfile(outDir,"temporal_summary.csv"));


%% ##################### Paper figures ###################################
%
% Each panel is exported as an independent PDF/PNG/FIG.
% The intended final width is approximately half of one ICASSP column.
% LaTeX should control side-by-side placement and (a)/(b)/(c) labels.
%
% Final figure files:
%   selection_visibility.pdf
%   time_sensitivity.pdf
%   haps_availability.pdf
%   leo_compensation.pdf
%   robustness.pdf

halfColumnWidth_in = 1.68;

fontName = 'Times New Roman';
fontSize = 7;
annotationFontSize = 6.5;
legendFontSize = 5.7;
lineWidth = 0.9;
markerSize = 3;
infeasibleColor = [0.6350 0.0780 0.1840];

archColors = [ ...
    0.0000 0.4470 0.7410
    0.8500 0.3250 0.0980
    0.9290 0.6940 0.1250
    0.4940 0.1840 0.5560];

%% ----------------------------------------------------------------------
% Figure 1: nominal selection visibility
%% ----------------------------------------------------------------------

fig = paperFigure( ...
    'Nominal selection visibility', ...
    halfColumnWidth_in,1.68);

skyPos = [0.10 0.25 0.80 0.69];

pax = polaraxes( ...
    'Parent',fig, ...
    'Units','normalized', ...
    'Position',skyPos);

hold(pax,'on');

theta = deg2rad(mod(Tpool.Azimuth_deg,360));
rho = 90-Tpool.Elevation_deg;

for aidx = 1:numel(archList)

    idxArch = Architecture == archList(aidx);

    idxSelected = ...
        idxArch & ismember((1:Npool).',selectedT0);

    idxNotSelected = ...
        idxArch & ~idxSelected;

    if any(idxNotSelected)

        polarscatter( ...
            pax, ...
            theta(idxNotSelected), ...
            rho(idxNotSelected), ...
            10, ...
            archColors(aidx,:), ...
            'o', ...
            'LineWidth',0.6, ...
            'HandleVisibility','off');
    end

    if any(idxSelected)

        polarscatter( ...
            pax, ...
            theta(idxSelected), ...
            rho(idxSelected), ...
            17, ...
            archColors(aidx,:), ...
            'o', ...
            'filled', ...
            'MarkerEdgeColor','k', ...
            'LineWidth',0.5, ...
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
    'RTickLabel',{'','',''}, ...
    'RAxisLocation',255, ...
    'FontName',fontName, ...
    'FontSize',fontSize, ...
    'GridAlpha',0.22);

% Radial coordinate is zenith angle; labels below show elevation.
text( ...
    pax,deg2rad(255),30,'60^\circ', ...
    'FontName',fontName, ...
    'FontSize',fontSize, ...
    'HorizontalAlignment','right', ...
    'VerticalAlignment','top');

text( ...
    pax,deg2rad(255),60,'30^\circ', ...
    'FontName',fontName, ...
    'FontSize',fontSize, ...
    'HorizontalAlignment','right', ...
    'VerticalAlignment','top');

hArch = gobjects(numel(archList),1);

for aidx = 1:numel(archList)

    hArch(aidx) = polarplot( ...
        pax,nan,nan,'o', ...
        'LineStyle','none', ...
        'Color',archColors(aidx,:), ...
        'MarkerFaceColor',archColors(aidx,:), ...
        'MarkerSize',3, ...
        'LineWidth',0.5);
end

hSelected = polarplot( ...
    pax,nan,nan,'ko', ...
    'LineStyle','none', ...
    'MarkerFaceColor','k', ...
    'MarkerSize',3, ...
    'LineWidth',0.5);

hNotSelected = polarplot( ...
    pax,nan,nan,'ko', ...
    'LineStyle','none', ...
    'MarkerFaceColor','none', ...
    'MarkerSize',3, ...
    'LineWidth',0.6);

skyLegend = legend( ...
    pax, ...
    [hArch;hSelected;hNotSelected], ...
    {'HAPS','LEO','MEO','GEO','Selected','Unselected'}, ...
    'NumColumns',3, ...
    'FontName',fontName, ...
    'FontSize',5.5, ...
    'Box','off', ...
    'AutoUpdate','off');

skyLegend.ItemTokenSize = [7 6];

set( ...
    skyLegend, ...
    'Units','normalized', ...
    'Position',[0.02 0.015 0.96 0.16]);

set(pax,'Position',skyPos);

exportFigurePair( ...
    fig,outDir,"selection_visibility");


%% ----------------------------------------------------------------------
% Figure 2: short-term temporal sensitivity
%% ----------------------------------------------------------------------

fig = paperFigure( ...
    'Short-term temporal sensitivity', ...
    halfColumnWidth_in,1.52);

timePos = [0.23 0.20 0.73 0.75];

ax = paperAxes( ...
    fig,timePos,fontName,fontSize);

hAdaptive = plot( ...
    ax, ...
    timeOffset_s/60, ...
    AdaptiveBound_m, ...
    '-o', ...
    'Color',archColors(1,:), ...
    'LineWidth',lineWidth, ...
    'MarkerSize',markerSize);

hold(ax,'on');

hFixed = plot( ...
    ax, ...
    timeOffset_s/60, ...
    FixedSubsetBound_m, ...
    '-s', ...
    'Color',archColors(2,:), ...
    'LineWidth',lineWidth, ...
    'MarkerSize',markerSize);

yline( ...
    ax,target_m,'--', ...
    'Color',[0.4 0.4 0.4], ...
    'LineWidth',0.7, ...
    'HandleVisibility','off');

text( ...
    ax, ...
    timeOffset_s(end)/60-0.02, ...
    target_m+0.004, ...
    'Target', ...
    'FontName',fontName, ...
    'FontSize',annotationFontSize, ...
    'HorizontalAlignment','right', ...
    'VerticalAlignment','bottom');

if isfinite(firstFixedViolation_s)

    xViolation_min = firstFixedViolation_s/60;
    yViolation_m = FixedSubsetBound_m(firstFixedViolationIdx);

    plot( ...
        ax, ...
        xViolation_min, ...
        yViolation_m, ...
        'o', ...
        'MarkerFaceColor',infeasibleColor, ...
        'MarkerEdgeColor','k', ...
        'MarkerSize',markerSize+1, ...
        'HandleVisibility','off');

    text( ...
        ax, ...
        xViolation_min+0.16, ...
        yViolation_m+0.001, ...
        sprintf('%d s',firstFixedViolation_s), ...
        'FontName',fontName, ...
        'FontSize',annotationFontSize, ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','bottom');
end

xlim(ax,[-0.07 timeOffset_s(end)/60+0.07]);

ylim( ...
    ax, ...
    [min([AdaptiveBound_m;FixedSubsetBound_m])-0.004 ...
     max([AdaptiveBound_m;FixedSubsetBound_m])+0.010]);

xticks(ax,0:1:timeOffset_s(end)/60);
yticks(ax,[0.60 0.65 0.70]);

xlabel(ax,'Time [min]');
ylabel(ax,'Position bound [m]');

timeLegend = legend( ...
    ax, ...
    [hAdaptive,hFixed], ...
    {'BE updated','Fixed set'}, ...
    'Location','northwest', ...
    'FontName',fontName, ...
    'FontSize',legendFontSize, ...
    'Box','off', ...
    'AutoUpdate','off');

timeLegend.ItemTokenSize = [8 6];

set(ax,'Position',timePos);

exportFigurePair( ...
    fig,outDir,"time_sensitivity");


%% ----------------------------------------------------------------------
% Figure 3: HAPS availability
%% ----------------------------------------------------------------------

fig = paperFigure( ...
    'HAPS availability', ...
    halfColumnWidth_in,1.58);

hapsPos = [0.22 0.20 0.60 0.74];

axHaps = paperAxes( ...
    fig,hapsPos,fontName,fontSize);

xHaps = 1:height(Thaps);

composition = [ ...
    Thaps.HAPS, ...
    Thaps.LEO, ...
    Thaps.MEO, ...
    Thaps.GEO];

feasibleHaps = logical(Thaps.TargetReached);

yyaxis(axHaps,'left');

hBars = bar( ...
    axHaps, ...
    xHaps, ...
    composition, ...
    'stacked', ...
    'BarWidth',0.65);

for aidx = 1:numel(hBars)

    hBars(aidx).FaceColor = archColors(aidx,:);
    hBars(aidx).LineWidth = 0.5;
end

ylabel( ...
    axHaps, ...
    'Transmitters used', ...
    'FontName',fontName, ...
    'FontSize',fontSize);

ylim(axHaps,[0 max(24,max(Thaps.Kfinal)+2)]);
yticks(axHaps,[0 12 24]);

yyaxis(axHaps,'right');

plot( ...
    axHaps, ...
    xHaps, ...
    Thaps.FinalBound_m, ...
    '-ko', ...
    'LineWidth',lineWidth, ...
    'MarkerSize',markerSize, ...
    'MarkerFaceColor','k');

hold(axHaps,'on');

plot( ...
    axHaps, ...
    xHaps(~feasibleHaps), ...
    Thaps.FinalBound_m(~feasibleHaps), ...
    'x', ...
    'Color',infeasibleColor, ...
    'MarkerSize',markerSize+2, ...
    'LineWidth',1.1, ...
    'HandleVisibility','off');

yline( ...
    axHaps,target_m,'--', ...
    'Color',[0.4 0.4 0.4], ...
    'LineWidth',0.7, ...
    'HandleVisibility','off');

ylabel( ...
    axHaps, ...
    'Bound [m]', ...
    'FontName',fontName, ...
    'FontSize',fontSize);

ylim(axHaps,[0.56 0.63]);
yticks(axHaps,[0.56 0.60 0.62]);

for k = 1:height(Thaps)

    text( ...
        axHaps, ...
        xHaps(k), ...
        Thaps.FinalBound_m(k)+0.0028, ...
        sprintf('%.3f',Thaps.FinalBound_m(k)), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', ...
        'FontName',fontName, ...
        'FontSize',6);
end

axHaps.YAxis(1).Color = 'k';
axHaps.YAxis(2).Color = 'k';

xlim(axHaps,[0.4 height(Thaps)+0.6]);
xticks(axHaps,xHaps);

hapsLabels = cellstr(string(Thaps.MaxHAPS));

for k = 1:height(Thaps)

    if ~feasibleHaps(k)
        hapsLabels{k} = [hapsLabels{k} '*'];
    end
end

xticklabels(axHaps,hapsLabels);

xlabel(axHaps,'HAPS limit');

hapsLegend = legend( ...
    axHaps, ...
    hBars(:), ...
    {'HAPS','LEO','MEO','GEO'}, ...
    'NumColumns',2, ...
    'FontName',fontName, ...
    'FontSize',5, ...
    'Box','off', ...
    'AutoUpdate','off');

hapsLegend.ItemTokenSize = [4.5 5];

set( ...
    hapsLegend, ...
    'Units','normalized', ...
    'Position',[0.24 0.73 0.42 0.15]);

set(axHaps,'Position',hapsPos);

exportFigurePair( ...
    fig,outDir,"haps_availability");


%% ----------------------------------------------------------------------
% Figure 4: LEO compensation
%% ----------------------------------------------------------------------

fig = paperFigure( ...
    'LEO compensation', ...
    halfColumnWidth_in,1.52);

leoPos = [0.23 0.20 0.72 0.75];

axLeo = paperAxes( ...
    fig,leoPos,fontName,fontSize);

hFull = plot( ...
    axLeo, ...
    Tleo.LEOAvailable, ...
    Tleo.AllCandidateBound_m, ...
    '-o', ...
    'Color',archColors(1,:), ...
    'LineWidth',lineWidth, ...
    'MarkerSize',markerSize);

hold(axLeo,'on');

feasibleLeo = logical(Tleo.TargetReached);
selectedBound = Tleo.FinalBound_m;

hSelectedLeo = plot( ...
    axLeo, ...
    Tleo.LEOAvailable, ...
    selectedBound, ...
    '--s', ...
    'Color',archColors(2,:), ...
    'LineWidth',lineWidth, ...
    'MarkerSize',markerSize);

hInfeasible = plot( ...
    axLeo, ...
    Tleo.LEOAvailable(~feasibleLeo), ...
    Tleo.AllCandidateBound_m(~feasibleLeo), ...
    'x', ...
    'Color',infeasibleColor, ...
    'LineWidth',1.1, ...
    'MarkerSize',markerSize+2);

yline( ...
    axLeo,target_m,'--', ...
    'Color',[0.4 0.4 0.4], ...
    'LineWidth',0.7, ...
    'HandleVisibility','off');

xlabel(axLeo,'Available LEO');
ylabel(axLeo,'Bound [m]');

xticks(axLeo,Tleo.LEOAvailable);

xlim( ...
    axLeo, ...
    [min(Tleo.LEOAvailable)-0.7 ...
     max(Tleo.LEOAvailable)+0.7]);

ylim(axLeo,[0.56 0.63]);
yticks(axLeo,[0.56 0.60 0.62]);

leoLegend = legend( ...
    axLeo, ...
    [hFull,hSelectedLeo,hInfeasible], ...
    {'Full pool','FS result','Infeasible'}, ...
    'NumColumns',1, ...
    'FontName',fontName, ...
    'FontSize',legendFontSize, ...
    'Box','off', ...
    'AutoUpdate','off', ...
    'Location','best');

leoLegend.ItemTokenSize = [7 5];

set(axLeo,'Position',leoPos);

exportFigurePair( ...
    fig,outDir,"leo_compensation");


%% ----------------------------------------------------------------------
% Figure 5: robustness to Gaussian-impulsive errors
%% ----------------------------------------------------------------------

fig = paperFigure( ...
    'Robustness to impulsive errors', ...
    halfColumnWidth_in,1.52);

outlierPos = [0.24 0.20 0.72 0.75];

axOutlier = paperAxes( ...
    fig,outlierPos,fontName,fontSize);

axOutlier.YScale = 'log';

hWLS = semilogy( ...
    axOutlier, ...
    Toutlier.B, ...
    Toutlier.RMSE_FullWLS_m, ...
    '-o', ...
    'Color',archColors(1,:), ...
    'LineWidth',lineWidth, ...
    'MarkerSize',markerSize);

hold(axOutlier,'on');

hRejected = semilogy( ...
    axOutlier, ...
    Toutlier.B, ...
    Toutlier.RMSE_RejectedWLS_m, ...
    '-d', ...
    'Color',archColors(3,:), ...
    'LineWidth',lineWidth, ...
    'MarkerSize',markerSize);

hRobust = semilogy( ...
    axOutlier, ...
    Toutlier.B, ...
    Toutlier.RMSE_RobustBE_m, ...
    '-s', ...
    'Color',archColors(2,:), ...
    'LineWidth',lineWidth, ...
    'MarkerSize',markerSize);

yline( ...
    axOutlier,target_m,'--', ...
    'Color',[0.4 0.4 0.4], ...
    'LineWidth',0.7, ...
    'HandleVisibility','off');

ylim( ...
    axOutlier, ...
    [0.42 max(20,1.2*max(Toutlier.RMSE_FullWLS_m))]);

xlim( ...
    axOutlier, ...
    [min(Toutlier.B)-0.3 ...
     max(Toutlier.B)+0.5]);

xticks(axOutlier,Toutlier.B);

axOutlier.YMinorGrid = 'off';

xlabel(axOutlier,'Outliers, B');
ylabel(axOutlier,'Position RMSE [m]');

outlierLegend = legend( ...
    axOutlier, ...
    [hWLS,hRejected,hRobust], ...
    {'Full WLS','Rejected WLS','Residual-aided BE'}, ...
    'NumColumns',1, ...
    'FontName',fontName, ...
    'FontSize',legendFontSize, ...
    'Box','off', ...
    'AutoUpdate','off', ...
    'Location','northwest');

outlierLegend.ItemTokenSize = [7 5];

set(axOutlier,'Position',outlierPos);

exportFigurePair( ...
    fig,outDir,"robustness");


%% ----------------------------------------------------------------------
% Figure 6: necessity of transmitter selection at fixed K
%% ----------------------------------------------------------------------
%
% The empirical CDF compares 50,000 uniformly sampled K=12 subsets with
% the FA/BE solution, the positioning target, and a simple highest-C/N0
% heuristic. The full 26-transmitter pool is intentionally omitted because
% this figure addresses the fixed-size selection problem.

methodNames = string(TselectionMethods.Method);

idxFA = methodNames == "FA selected";
idxBE = methodNames == "BE selected";
idxTopCN0 = methodNames == "Top C/N0";

if sum(idxFA) ~= 1 || ...
   sum(idxBE) ~= 1 || ...
   sum(idxTopCN0) ~= 1

    error([ ...
        'Selection-necessity methods table must contain exactly one ' ...
        'FA selected, BE selected, and Top C/N0 row.']);
end

%% Fixed subset size

selectionK = double(TselectionMethods.K(idxFA));

if double(TselectionMethods.K(idxBE)) ~= selectionK || ...
   double(TselectionMethods.K(idxTopCN0)) ~= selectionK

    error('Selection-necessity methods must use the same fixed K.');
end

%% Deterministic bounds

faSelectionBound_m = ...
    double(TselectionMethods.PositionalBound_m(idxFA));

beSelectionBound_m = ...
    double(TselectionMethods.PositionalBound_m(idxBE));

if abs(faSelectionBound_m-beSelectionBound_m) > 1e-10

    error('FA and BE bounds differ in the fixed-K comparison.');
end

selectedSelectionBound_m = ...
    0.5*(faSelectionBound_m+beSelectionBound_m);

topCN0SelectionBound_m = ...
    double(TselectionMethods.PositionalBound_m(idxTopCN0));

%% Random-subset distribution

randomSelectionBound_m = ...
    double(TselectionRandom.PositionalBound_m(:));

randomSelectionRank = ...
    double(TselectionRandom.RankH(:));

validRandom = ...
    randomSelectionRank == 4 & ...
    isfinite(randomSelectionBound_m);

randomSelectionBound_m = ...
    sort(randomSelectionBound_m(validRandom));

if isempty(randomSelectionBound_m)

    error('No finite full-rank random subsets are available.');
end

randomCDF = ...
    (1:numel(randomSelectionBound_m)).' ...
    / numel(randomSelectionBound_m);

randomMedianBound_m = ...
    median(randomSelectionBound_m);

randomTargetProbability = ...
    mean(randomSelectionBound_m <= target_m);

%% Plot

fig = paperFigure( ...
    'Fixed-size transmitter-selection necessity', ...
    halfColumnWidth_in,1.58);

selectionPos = [0.24 0.20 0.72 0.75];

axSelection = paperAxes( ...
    fig,selectionPos,fontName,fontSize);

hRandomCDF = plot( ...
    axSelection, ...
    randomSelectionBound_m, ...
    randomCDF, ...
    '-', ...
    'Color',archColors(1,:), ...
    'LineWidth',1.0);

hold(axSelection,'on');

hSelectedBound = xline( ...
    axSelection, ...
    selectedSelectionBound_m, ...
    '-', ...
    'Color',archColors(2,:), ...
    'LineWidth',1.0);

hTargetBound = xline( ...
    axSelection, ...
    target_m, ...
    '--', ...
    'Color',[0.35 0.35 0.35], ...
    'LineWidth',0.9);

hTopCN0Bound = xline( ...
    axSelection, ...
    topCN0SelectionBound_m, ...
    ':', ...
    'Color',archColors(4,:), ...
    'LineWidth',1.0);

xlim(axSelection,[0.55 1.30]);
ylim(axSelection,[0 1]);

xlabel(axSelection,'Position bound [m]');
ylabel(axSelection,'Empirical CDF');

selectionLegend = legend( ...
    axSelection, ...
    [hRandomCDF,hSelectedBound,hTargetBound,hTopCN0Bound], ...
    {sprintf('Random subsets, K=%d',selectionK), ...
     'FS/BE','Target','Highest-C/N0 subset'}, ...
    'Location','northwest', ...
    'FontName',fontName, ...
    'FontSize',5.2, ...
    'Box','off', ...
    'AutoUpdate','off');

selectionLegend.ItemTokenSize = [7 5];

set(axSelection,'Position',selectionPos);

exportFigurePair( ...
    fig,outDir,"selection_necessity");

fprintf('\nSelection-necessity figure\n');
fprintf('K                         : %d\n',selectionK);
fprintf('FA/BE bound               : %.12f m\n', ...
    selectedSelectionBound_m);
fprintf('Top-C/N0 bound            : %.12f m\n', ...
    topCN0SelectionBound_m);
fprintf('Random median bound       : %.12f m\n', ...
    randomMedianBound_m);
fprintf('P(random bound <= %.3f m): %.4f %%\n', ...
    target_m,100*randomTargetProbability);

fprintf('\nSeparate paper figures saved in:\n  %s\n',outDir);
fprintf('LaTeX should control side-by-side placement and panel labels.\n');


%% ##################### Console summary #################################

fprintf('\n');
fprintf('============================================================\n');
fprintf('PAPER FIGURES - GEOMETRY AND TEMPORAL SENSITIVITY\n');
fprintf('============================================================\n');
fprintf('Nominal candidates                 : %d\n',Npool);
fprintf('Nominal selected subset            : %d\n',numel(selectedT0));
fprintf('Selected bound at t0               : %.12f m\n', ...
    AdaptiveBound_m(1));
fprintf('First sampled adaptive-set change  : %.0f s\n', ...
    firstAdaptiveChange_s);
fprintf('First sampled fixed-set violation  : %.0f s\n', ...
    firstFixedViolation_s);
fprintf('Adaptive K range                   : %.0f to %.0f\n', ...
    min(AdaptiveK),max(AdaptiveK));
fprintf('All nominal candidates visible     : %s\n', ...
    string(all(FullVisibleCount == Npool)));
fprintf('Selection-necessity fixed K         : %d\n',selectionK);
fprintf('Random subsets meeting target       : %.4f %%\n', ...
    100*randomTargetProbability);
fprintf('FA/BE fixed-K bound                 : %.12f m\n', ...
    selectedSelectionBound_m);
fprintf('Random fixed-K median bound         : %.12f m\n', ...
    randomMedianBound_m);
fprintf('Results saved to:\n  %s\n',outDir);

%% ##################### Local functions #################################



function cfg = paperFigureConfig(rootDir)
% Fixed configuration used by the paper-figure calculations.

    cfg.analysisTimeUTC = datetime( ...
        2026,8,1,12,0,0,'TimeZone','UTC');

    cfg.user.lat_deg = -22.8596582;
    cfg.user.lon_deg = -43.2303236;
    cfg.user.h_m = 10;

    cfg.mask.HAPS_deg = 15;
    cfg.mask.LEO_deg = 5;
    cfg.mask.MEO_deg = 5;
    cfg.mask.GEO_deg = 5;

    cfg.selectionTarget_m = 0.6;

    cfg.ranging.c_mps = 299792458;
    cfg.ranging.beta_Hz = 1.023e6;
    cfg.ranging.Tcoh_s = 20e-3;
    cfg.ranging.etaPos = 0.01;

    cfg.rf.HAPS.Frequency_GHz = 31.15;
    cfg.rf.HAPS.EIRP_PSD_dBW_MHz = -9.1;
    cfg.rf.HAPS.ReceiverFoM_dBK = 12.3;

    cfg.rf.LEO.Frequency_GHz = 20.0;
    cfg.rf.LEO.ReceiverFoM_dBK = 10.5;
    cfg.rf.LEO.ElevationControl_deg = [10;25;90];
    cfg.rf.LEO.EIRP_Control_dBW_MHz = [-2.7;17.9;11.3];

    cfg.rf.MEO.Frequency_GHz = 20.0;
    cfg.rf.MEO.EIRP_total_dBW = 49.7;
    cfg.rf.MEO.BW_MHz = 100;
    cfg.rf.MEO.ReceiverFoM_dBK = 18.4;

    cfg.rf.GEO.Frequency_GHz = 19.95;
    cfg.rf.GEO.EIRP_total_dBW = 57.0;
    cfg.rf.GEO.BW_MHz = 81;
    cfg.rf.GEO.G_rx_dBi = 44.5;
    cfg.rf.GEO.Tsys_K = 250;

    cfg.paths.tleRoot = fullfile(rootDir,'tle');
end

function assertTableVariables(T,required,description)

    available = string(T.Properties.VariableNames);
    missing = setdiff(required,available);

    if ~isempty(missing)
        error("Missing variable(s) in %s: %s", ...
            description,strjoin(missing,", "));
    end
end


function positionECEF = geodeticToECEFWGS84(lat_deg,lon_deg,h_m)

    a_m = 6378137.0;
    flattening = 1/298.257223563;
    eccentricitySquared = flattening*(2-flattening);

    lat = deg2rad(lat_deg);
    lon = deg2rad(lon_deg);

    primeVerticalRadius = ...
        a_m/sqrt(1-eccentricitySquared*sin(lat)^2);

    positionECEF = [ ...
        (primeVerticalRadius+h_m)*cos(lat)*cos(lon), ...
        (primeVerticalRadius+h_m)*cos(lat)*sin(lon), ...
        (primeVerticalRadius*(1-eccentricitySquared)+h_m)*sin(lat)];
end

function [H,sigmaRho_m,visible,CN0_dBHz,azimuth_deg,elevation_deg] = ...
    epochGeometryAndQuality( ...
        positionECEF,architecture,userPositionECEF,cfg, ...
        archList,zenithGasLoss_dB,cn0Calibration_dB)

    N = size(positionECEF,1);
    userPositionECEF = userPositionECEF(:).';

    delta = positionECEF-userPositionECEF;
    slantRange_m = vecnorm(delta,2,2);

    if any(~isfinite(slantRange_m)) || any(slantRange_m <= eps)
        error("Invalid transmitter-user range in temporal analysis.");
    end

    unitReceiverToTransmitter = delta./slantRange_m;
    H = [-unitReceiverToTransmitter,ones(N,1)];

    [azimuth_deg,elevation_deg] = ecefAzEl( ...
        positionECEF,userPositionECEF, ...
        cfg.user.lat_deg,cfg.user.lon_deg);

    visible = false(N,1);
    frequency_GHz = nan(N,1);
    eirpPsd_dBW_MHz = nan(N,1);
    receiverFoM_dBK = nan(N,1);
    gasLoss_dB = nan(N,1);

    for i = 1:N
        arch = architecture(i);
        archIdx = find(archList == arch,1,'first');

        if isempty(archIdx)
            error("Unsupported architecture: %s",arch);
        end

        switch arch
            case "HAPS"
                mask_deg = cfg.mask.HAPS_deg;
                frequency_GHz(i) = cfg.rf.HAPS.Frequency_GHz;
                eirpPsd_dBW_MHz(i) = cfg.rf.HAPS.EIRP_PSD_dBW_MHz;
                receiverFoM_dBK(i) = cfg.rf.HAPS.ReceiverFoM_dBK;

            case "LEO"
                mask_deg = cfg.mask.LEO_deg;
                frequency_GHz(i) = cfg.rf.LEO.Frequency_GHz;
                receiverFoM_dBK(i) = cfg.rf.LEO.ReceiverFoM_dBK;

                elevationForPower_deg = min(max( ...
                    elevation_deg(i), ...
                    min(cfg.rf.LEO.ElevationControl_deg)), ...
                    max(cfg.rf.LEO.ElevationControl_deg));

                eirpPsd_dBW_MHz(i) = interp1( ...
                    cfg.rf.LEO.ElevationControl_deg, ...
                    cfg.rf.LEO.EIRP_Control_dBW_MHz, ...
                    elevationForPower_deg,'linear');

            case "MEO"
                mask_deg = cfg.mask.MEO_deg;
                frequency_GHz(i) = cfg.rf.MEO.Frequency_GHz;
                eirpPsd_dBW_MHz(i) = ...
                    cfg.rf.MEO.EIRP_total_dBW-10*log10(cfg.rf.MEO.BW_MHz);
                receiverFoM_dBK(i) = cfg.rf.MEO.ReceiverFoM_dBK;

            case "GEO"
                mask_deg = cfg.mask.GEO_deg;
                frequency_GHz(i) = cfg.rf.GEO.Frequency_GHz;
                eirpPsd_dBW_MHz(i) = ...
                    cfg.rf.GEO.EIRP_total_dBW-10*log10(cfg.rf.GEO.BW_MHz);
                receiverFoM_dBK(i) = ...
                    cfg.rf.GEO.G_rx_dBi-10*log10(cfg.rf.GEO.Tsys_K);
        end

        visible(i) = elevation_deg(i) >= mask_deg;

        elevationForGas_deg = max(elevation_deg(i),0.1);
        gasLoss_dB(i) = ...
            zenithGasLoss_dB(archIdx)/sind(elevationForGas_deg);
    end

    FSPL_dB = 20*log10( ...
        4*pi.*slantRange_m.*frequency_GHz*1e9/cfg.ranging.c_mps);

    EIRPpositioning_dBW = ...
        eirpPsd_dBW_MHz ...
        + 10*log10(cfg.ranging.beta_Hz/1e6) ...
        + 10*log10(cfg.ranging.etaPos);

    CN0_dBHz = ...
        EIRPpositioning_dBW ...
        - FSPL_dB ...
        - gasLoss_dB ...
        + receiverFoM_dBK ...
        + 228.6 ...
        + cn0Calibration_dB;

    CN0_linear = 10.^(CN0_dBHz/10);

    sigmaRho_m = ...
        cfg.ranging.c_mps./ ...
        (2*pi*cfg.ranging.beta_Hz.* ...
        sqrt(CN0_linear*cfg.ranging.Tcoh_s));

    if any(~isfinite(sigmaRho_m(visible))) || ...
       any(sigmaRho_m(visible) <= 0)
        error("Invalid temporal pseudorange standard deviation.");
    end
end

function [azimuth_deg,elevation_deg] = ...
    ecefAzEl(positionECEF,userECEF,userLat_deg,userLon_deg)

    lat = deg2rad(userLat_deg);
    lon = deg2rad(userLon_deg);

    rotationECEFToENU = [ ...
        -sin(lon),             cos(lon),            0
        -sin(lat)*cos(lon),   -sin(lat)*sin(lon),   cos(lat)
         cos(lat)*cos(lon),    cos(lat)*sin(lon),   sin(lat)];

    delta = positionECEF-userECEF(:).';
    enu = (rotationECEFToENU*delta.').';

    azimuth_deg = mod(atan2d(enu(:,1),enu(:,2)),360);
    elevation_deg = atan2d(enu(:,3),hypot(enu(:,1),enu(:,2)));
end

function [selected,finalBound_m,targetReached] = backwardElimination( ...
    initialIdx,Hfull,sigmaRho_m,target_m,targetTol)

    selected = initialIdx(:).';
    [currentBound_m,~,currentRank] = ...
        subsetMetrics(selected,Hfull,sigmaRho_m);

    if currentRank ~= 4 || currentBound_m > target_m+targetTol
        finalBound_m = currentBound_m;
        targetReached = false;
        return;
    end

    while numel(selected) > 4
        bestRemoved = NaN;
        bestTrialBound_m = inf;

        for k = 1:numel(selected)
            candidateToRemove = selected(k);
            trialIdx = selected(selected ~= candidateToRemove);

            [trialBound_m,~,trialRank] = ...
                subsetMetrics(trialIdx,Hfull,sigmaRho_m);

            if trialRank ~= 4 || trialBound_m > target_m+targetTol
                continue;
            end

            if trialBound_m < bestTrialBound_m
                bestTrialBound_m = trialBound_m;
                bestRemoved = candidateToRemove;
            elseif abs(trialBound_m-bestTrialBound_m) <= 1e-14
                if isnan(bestRemoved) || candidateToRemove < bestRemoved
                    bestTrialBound_m = trialBound_m;
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

    finalBound_m = currentBound_m;
    targetReached = finalBound_m <= target_m+targetTol;
end

function [bound_m,PDOP,rankH] = subsetMetrics(idx,Hfull,sigmaRho_m)

    idx = idx(:);
    H = Hfull(idx,:);
    sigma = sigmaRho_m(idx);

    rankH = rank(H);
    bound_m = inf;
    PDOP = inf;

    if rankH < 4
        return;
    end

    normalMatrix = H.'*H;

    if rcond(normalMatrix) > 1e-14
        covarianceGeometry = normalMatrix\eye(4);
        traceGeometry = trace(covarianceGeometry(1:3,1:3));

        if traceGeometry > 0 && isfinite(traceGeometry)
            PDOP = sqrt(traceGeometry);
        end
    end

    J = H.'*((diag(sigma.^2))\H);

    if rcond(J) <= 1e-14
        return;
    end

    covarianceBound = J\eye(4);
    tracePosition = trace(covarianceBound(1:3,1:3));

    if tracePosition > 0 && isfinite(tracePosition)
        bound_m = sqrt(tracePosition);
    end
end

function fig = paperFigure(figureName,width_in,height_in)
    fig = figure('Name',figureName,'NumberTitle','off', ...
        'WindowStyle','normal','Units','inches', ...
        'Position',[1 1 width_in height_in],'Color','w', ...
        'PaperUnits','inches','PaperSize',[width_in height_in], ...
        'PaperPosition',[0 0 width_in height_in], ...
        'PaperPositionMode','manual','InvertHardcopy','off');
    % Preserve the intended print dimensions if the popup is later resized.
    setappdata(fig,'PaperSizeInches',[width_in height_in]);
end

function ax = paperAxes(fig,axesPosition,fontName,fontSize)
    ax = axes('Parent',fig,'Units','normalized','Position',axesPosition, ...
        'PositionConstraint','innerposition', ...
        'FontName',fontName,'FontUnits','points','FontSize',fontSize, ...
        'LabelFontSizeMultiplier',1,'TitleFontSizeMultiplier',1, ...
        'TitleFontWeight','normal','LineWidth',0.6,'Box','on', ...
        'GridAlpha',0.18);
    grid(ax,'on');
    % Preserve fonts/positions when high-level plot/bar commands run.
    hold(ax,'on');
end

function exportFigurePair(fig,outDir,baseName)
    if ~isfolder(outDir)
        [ok,msg] = mkdir(outDir);
        if ~ok
            error('Could not create output folder: %s',msg);
        end
    end

    pageSize = getappdata(fig,'PaperSizeInches');
    set(fig,'PaperUnits','inches','PaperSize',pageSize, ...
        'PaperPosition',[0 0 pageSize],'PaperPositionMode','manual');
    drawnow;

    % Explicit page dimensions preserve the intended column width in PDF
    % and PNG; exportgraphics tight-cropping is deliberately not used here.
    % Short temporary paths also avoid image-library failures on long paths.
    % MATLAB documentation: https://www.mathworks.com/help/matlab/ref/print.html
    tempPdf = [tempname '.pdf'];
    tempPng = [tempname '.png'];
    cleanup = onCleanup(@() removeTemporaryExports(tempPdf,tempPng));
    print(fig,tempPdf,'-dpdf','-painters');
    print(fig,tempPng,'-dpng','-r300');

    pdfFile = fullfile(outDir,baseName+".pdf");
    pngFile = fullfile(outDir,baseName+".png");
    [ok,msg] = copyfile(tempPdf,char(pdfFile),'f');
    if ~ok
        error('Could not save PDF to %s: %s',char(pdfFile),msg);
    end
    [ok,msg] = copyfile(tempPng,char(pngFile),'f');
    if ~ok
        error('Could not save PNG to %s: %s',char(pngFile),msg);
    end
    savefig(fig,char(fullfile(outDir,baseName+".fig")));
    fprintf('Paper figure saved: %s\n',char(baseName));
end

function removeTemporaryExports(varargin)
    for k = 1:nargin
        if isfile(varargin{k})
            delete(varargin{k});
        end
    end
end
