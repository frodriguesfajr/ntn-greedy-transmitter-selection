clear;
clc;
format long;

%% ============================================================
% NO-HAPS LEO AVAILABILITY AND COMPENSATION
%
% Objective:
%   Evaluate whether increasing the number of available LEO candidates
%   can compensate for the absence of HAPS while preserving the
%   transmitter-selection target.
%
% Cases:
%   8, 12, 16, and 20 available LEO candidates
%
% Fixed in every case:
%   7 MEO
%   7 GEO
%   no HAPS
%
% The LEO sets are nested according to the predefined experiment ordering.
% Orbital geometry is obtained from the archived TLEs using the included
% Vallado SGP4 implementation with WGS-72 constants.
%
% The atmospheric gaseous-loss values for the fixed experiment are read
% from:
%
%   data/atmospheric_gas_loss_reference.csv
%
% For each case, the script evaluates the full candidate set and then
% applies Forward Adding (FA) from the best full-rank four-transmitter
% initialization whenever the 0.6 m positional-bound target is feasible.
%
% Outputs:
%   results/LEO_availability_results.mat
%   results/LEO_availability_summary.csv
%% ============================================================

rootDir = fileparts(mfilename('fullpath'));

if isempty(rootDir)
    rootDir = pwd;
end

cd(rootDir);
setup_paths;

fprintf('\n============================================================\n');
fprintf('NO-HAPS LEO AVAILABILITY\n');
fprintf('============================================================\n');

%% ============================================================
% Input files
%% ============================================================

catalogFile = fullfile( ...
    rootDir, ...
    'results', ...
    'orbital_catalog.mat');

gasReferenceFile = fullfile( ...
    rootDir, ...
    'data', ...
    'atmospheric_gas_loss_reference.csv');

assert(isfile(catalogFile), ...
    ['Orbital catalog not found:\n%s\n' ...
     'Run build_orbital_catalog.m first.'], ...
    catalogFile);

assert(isfile(gasReferenceFile), ...
    'Atmospheric gas-loss reference not found:\n%s', ...
    gasReferenceFile);

%% ============================================================
% Load orbital catalog
%% ============================================================

S = load(catalogFile);

assert(isfield(S,'TOrbitalCatalog'), ...
    'TOrbitalCatalog not found in orbital_catalog.mat.');

assert(isfield(S,'UserPositionECEF'), ...
    'UserPositionECEF not found in orbital_catalog.mat.');

T = S.TOrbitalCatalog;

UserPositionECEF = S.UserPositionECEF(:).';

if isfield(S,'analysisTimeUTC')
    analysisTimeUTC = S.analysisTimeUTC;
else
    analysisTimeUTC = datetime( ...
        2026,8,1,12,0,0, ...
        'TimeZone','UTC');
end

%% ============================================================
% Load atmospheric gas-loss reference
%% ============================================================

Tgas = readtable( ...
    gasReferenceFile, ...
    'TextType','string');

requiredGasVars = { ...
    'Architecture', ...
    'NORAD_CAT_ID', ...
    'AtmosphericGasLoss_dB'};

for k = 1:numel(requiredGasVars)

    assert(ismember( ...
        requiredGasVars{k}, ...
        Tgas.Properties.VariableNames), ...
        'Missing gas-reference variable: %s', ...
        requiredGasVars{k});

end

Tgas.Architecture = ...
    string(Tgas.Architecture);

Tgas.NORAD_CAT_ID = ...
    string(Tgas.NORAD_CAT_ID);

%% ============================================================
% Orbital ordering used in the experiment
%% ============================================================

leoOrder = [ ...
    "65175"
    "66864"
    "65348"
    "59340"
    "67736"
    "59960"
    "65279"
    "48482"
    "59544"
    "57470"
    "67011"
    "66267"
    "65388"
    "64047"
    "64037"
    "56537"
    "60271"
    "58736"
    "67734"
    "67950"];

meoOrder = [ ...
    "62362"
    "40081"
    "44115"
    "43231"
    "43234"
    "40080"
    "40079"];

geoOrder = [ ...
    "41589"
    "41904"
    "42692"
    "43562"
    "43175"
    "43228"
    "38087"];

TLEO20 = reorderArchitecture( ...
    T, ...
    "LEO", ...
    leoOrder, ...
    "LEO");

TMEO7 = reorderArchitecture( ...
    T, ...
    "MEO", ...
    meoOrder, ...
    "MEO");

TGEO7 = reorderArchitecture( ...
    T, ...
    "GEO", ...
    geoOrder, ...
    "GEO");

%% LEO nested rank

TLEO20.LEORank = (1:20).';

%% ============================================================
% Add link-budget quantities
%% ============================================================

TLEO20 = addLinkBudgetFromReference( ...
    TLEO20,Tgas);

TMEO7 = addLinkBudgetFromReference( ...
    TMEO7,Tgas);

TGEO7 = addLinkBudgetFromReference( ...
    TGEO7,Tgas);

%% ============================================================
% Cases
%% ============================================================

target_m = 0.6;

leoCases = [8 12 16 20];

nCases = numel(leoCases);

%% ============================================================
% Output arrays
%% ============================================================

LEOAvailable = leoCases(:);

CandidatePoolSize = zeros(nCases,1);

AllCandidateBound_m = nan(nCases,1);
AllCandidatePDOP    = nan(nCases,1);

Best4Bound_m = nan(nCases,1);

TargetReached = false(nCases,1);

Kfinal       = nan(nCases,1);
FinalBound_m = nan(nCases,1);
FinalPDOP    = nan(nCases,1);

LEOSelected = zeros(nCases,1);
MEOSelected = zeros(nCases,1);
GEOSelected = zeros(nCases,1);

LEOCN0Min_dBHz = nan(nCases,1);
LEOCN0Median_dBHz = nan(nCases,1);
LEOCN0Max_dBHz = nan(nCases,1);

LEOSigmaMin_m = nan(nCases,1);
LEOSigmaMedian_m = nan(nCases,1);
LEOSigmaMax_m = nan(nCases,1);

LEOGasLossMin_dB = nan(nCases,1);
LEOGasLossMedian_dB = nan(nCases,1);
LEOGasLossMax_dB = nan(nCases,1);

selectedPerCase = cell(nCases,1);

%% ============================================================
% Main loop
%% ============================================================

for caseIdx = 1:nCases

    nLEO = leoCases(caseIdx);

    fprintf('\n\n############################################################\n');
    fprintf('CASE: %d LEO + 7 MEO + 7 GEO | NO HAPS\n',nLEO);
    fprintf('############################################################\n');

    %% --------------------------------------------------------
    % Nested LEO subset
    %% --------------------------------------------------------

    TL = TLEO20(1:nLEO,:);

    % LEORank exists only in TLEO20.
    % Remove it from the temporary table before concatenation.
    TLcase = removevars(TL,'LEORank');

    Tcase = [ ...
        TLcase;
        TMEO7;
        TGEO7];

    Npool = height(Tcase);

    CandidatePoolSize(caseIdx) = Npool;

    Architecture = ...
        string(Tcase.Architecture);

    Label = ...
        string(Tcase.Label);

    PositionECEF = [ ...
        Tcase.X_ECEF_m, ...
        Tcase.Y_ECEF_m, ...
        Tcase.Z_ECEF_m];

    sigmaRho_m = ...
        Tcase.SigmaRho_m(:);

    %% --------------------------------------------------------
    % Geometry
    %% --------------------------------------------------------

    r = ...
        PositionECEF - UserPositionECEF;

    d = ...
        sqrt(sum(r.^2,2));

    assert(all(isfinite(d)));
    assert(all(d > 0));

    u = r ./ d;

    Hfull = ...
        [-u,ones(Npool,1)];

    %% --------------------------------------------------------
    % Full available set
    %% --------------------------------------------------------

    allIdx = 1:Npool;

    [allBound,allPDOP,allRank] = ...
        subsetMetrics( ...
        allIdx, ...
        Hfull, ...
        sigmaRho_m);

    AllCandidateBound_m(caseIdx) = ...
        allBound;

    AllCandidatePDOP(caseIdx) = ...
        allPDOP;

    fprintf('\nFull candidate set\n');
    fprintf('N             : %d\n',Npool);
    fprintf('Bound         : %.12f m\n',allBound);
    fprintf('PDOP          : %.12f\n',allPDOP);
    fprintf('rank(H)       : %d\n',allRank);

    fprintf('Target met    : %s\n', ...
        string(allBound <= target_m));

    %% --------------------------------------------------------
    % LEO measurement-quality statistics
    %% --------------------------------------------------------

    LEOCN0Min_dBHz(caseIdx) = ...
        min(TL.CN0_dBHz);

    LEOCN0Median_dBHz(caseIdx) = ...
        median(TL.CN0_dBHz);

    LEOCN0Max_dBHz(caseIdx) = ...
        max(TL.CN0_dBHz);

    LEOSigmaMin_m(caseIdx) = ...
        min(TL.SigmaRho_m);

    LEOSigmaMedian_m(caseIdx) = ...
        median(TL.SigmaRho_m);

    LEOSigmaMax_m(caseIdx) = ...
        max(TL.SigmaRho_m);

    LEOGasLossMin_dB(caseIdx) = ...
        min(TL.AtmosphericGasLoss_dB);

    LEOGasLossMedian_dB(caseIdx) = ...
        median(TL.AtmosphericGasLoss_dB);

    LEOGasLossMax_dB(caseIdx) = ...
        max(TL.AtmosphericGasLoss_dB);

    %% --------------------------------------------------------
    % Exhaustive best K=4 initialization
    %% --------------------------------------------------------

    comb4 = ...
        nchoosek(1:Npool,4);

    best4Bound = inf;
    best4 = [];

    for k = 1:size(comb4,1)

        idx = comb4(k,:);

        [b,~,rankH] = ...
            subsetMetrics( ...
            idx, ...
            Hfull, ...
            sigmaRho_m);

        if rankH == 4 && ...
           b < best4Bound

            best4Bound = b;
            best4 = idx;

        end
    end

    if isempty(best4)

        error( ...
            'No valid initial K=4 subset for %d LEO.', ...
            nLEO);

    end

    Best4Bound_m(caseIdx) = ...
        best4Bound;

    fprintf('\nBest K=4\n');

    fprintf('Bound         : %.12f m\n', ...
        best4Bound);

    fprintf('Indices       : ');
    fprintf('%d ',best4);
    fprintf('\n');

    fprintf('Labels        : ');
    fprintf('%s ',Label(best4));
    fprintf('\n');

    %% --------------------------------------------------------
    % If even full set is infeasible
    %% --------------------------------------------------------

    if allBound > target_m

        selected = allIdx;

        finalBound = allBound;
        finalPDOP  = allPDOP;

        reached = false;

        fprintf(['\nTarget is infeasible with all ' ...
                 'available transmitters.\n']);

    else

        %% ----------------------------------------------------
        % Forward Adding
        %% ----------------------------------------------------

        selected = best4(:).';

        currentBound = best4Bound;

        fprintf('\nForward Adding\n');

        fprintf( ...
            'K=%d | bound=%.12f m | initialization\n', ...
            numel(selected), ...
            currentBound);

        while currentBound > target_m && ...
              numel(selected) < Npool

            remaining = ...
                setdiff( ...
                1:Npool, ...
                selected, ...
                'stable');

            bestCandidate = NaN;
            bestCandidateBound = inf;

            for k = 1:numel(remaining)

                candidate = ...
                    remaining(k);

                idxTrial = ...
                    [selected,candidate];

                [b,~,rankH] = ...
                    subsetMetrics( ...
                    idxTrial, ...
                    Hfull, ...
                    sigmaRho_m);

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

            fprintf( ...
                'K=%d | add %-6s (%s) | bound=%.12f m\n', ...
                numel(selected), ...
                Label(bestCandidate), ...
                Architecture(bestCandidate), ...
                currentBound);

        end

        [finalBound,finalPDOP,~] = ...
            subsetMetrics( ...
            selected, ...
            Hfull, ...
            sigmaRho_m);

        reached = ...
            finalBound <= target_m;

    end

    %% --------------------------------------------------------
    % Final result
    %% --------------------------------------------------------

    Kfinal(caseIdx) = ...
        numel(selected);

    FinalBound_m(caseIdx) = ...
        finalBound;

    FinalPDOP(caseIdx) = ...
        finalPDOP;

    TargetReached(caseIdx) = ...
        reached;

    LEOSelected(caseIdx) = ...
        sum(Architecture(selected)=="LEO");

    MEOSelected(caseIdx) = ...
        sum(Architecture(selected)=="MEO");

    GEOSelected(caseIdx) = ...
        sum(Architecture(selected)=="GEO");

    selectedPerCase{caseIdx} = ...
        selected;

    fprintf('\nRESULT\n');

    fprintf('Final K       : %d\n', ...
        numel(selected));

    fprintf('Final bound   : %.12f m\n', ...
        finalBound);

    fprintf('PDOP          : %.12f\n', ...
        finalPDOP);

    fprintf('Target met    : %s\n', ...
        string(reached));

    fprintf( ...
        'Composition   : LEO=%d | MEO=%d | GEO=%d\n', ...
        LEOSelected(caseIdx), ...
        MEOSelected(caseIdx), ...
        GEOSelected(caseIdx));

    fprintf('Selected      : ');
    fprintf('%d ',selected);
    fprintf('\n');

end

%% ============================================================
% Summary
%% ============================================================

Tcomparison = table( ...
    LEOAvailable, ...
    CandidatePoolSize, ...
    AllCandidateBound_m, ...
    AllCandidatePDOP, ...
    Best4Bound_m, ...
    TargetReached, ...
    Kfinal, ...
    FinalBound_m, ...
    FinalPDOP, ...
    LEOSelected, ...
    MEOSelected, ...
    GEOSelected, ...
    LEOCN0Min_dBHz, ...
    LEOCN0Median_dBHz, ...
    LEOCN0Max_dBHz, ...
    LEOSigmaMin_m, ...
    LEOSigmaMedian_m, ...
    LEOSigmaMax_m, ...
    LEOGasLossMin_dB, ...
    LEOGasLossMedian_dB, ...
    LEOGasLossMax_dB);

fprintf('\n============================================================\n');
fprintf('LEO-AVAILABILITY SUMMARY\n');
fprintf('============================================================\n');

disp(Tcomparison);

%% ============================================================
% Save
%% ============================================================

resultsDir = ...
    fullfile(rootDir,'results');

if ~exist(resultsDir,'dir')
    mkdir(resultsDir);
end

matFile = fullfile( ...
    resultsDir, ...
    'LEO_availability_results.mat');

csvFile = fullfile( ...
    resultsDir, ...
    'LEO_availability_summary.csv');

save(matFile, ...
    'TLEO20', ...
    'TMEO7', ...
    'TGEO7', ...
    'Tcomparison', ...
    'selectedPerCase', ...
    'leoCases', ...
    'target_m', ...
    'UserPositionECEF', ...
    'analysisTimeUTC');

writetable( ...
    Tcomparison, ...
    csvFile);

fprintf('\nSaved:\n%s\n%s\n', ...
    matFile,csvFile);

%% ============================================================
% Reorder architecture according to the experiment catalog
%% ============================================================

function Tout = reorderArchitecture( ...
    T,arch,noradOrder,labelPrefix)

    Tout = T([],:);

    for k = 1:numel(noradOrder)

        idx = ...
            string(T.Architecture)==arch & ...
            string(T.NORAD_CAT_ID)==noradOrder(k);

        if sum(idx) ~= 1

            error( ...
                'Could not uniquely locate %s NORAD %s.', ...
                arch,noradOrder(k));

        end

        Tout = [Tout;T(idx,:)]; %#ok<AGROW>

    end

    Tout.Label = ...
        labelPrefix + string((1:height(Tout)).');

end


%% ============================================================
% Add link-budget quantities from fixed atmospheric-loss table
%% ============================================================

function T = addLinkBudgetFromReference(T,Tgas)

    N = height(T);

    Architecture = ...
        string(T.Architecture);

    NORAD_CAT_ID = ...
        string(T.NORAD_CAT_ID);

    %% --------------------------------------------------------
    % Atmospheric gas loss
    %% --------------------------------------------------------

    AtmosphericGasLoss_dB = ...
        nan(N,1);

    for i = 1:N

        idx = ...
            Tgas.Architecture == Architecture(i) & ...
            Tgas.NORAD_CAT_ID == NORAD_CAT_ID(i);

        if sum(idx) ~= 1

            error([ ...
                'Atmospheric loss not uniquely found for ' ...
                '%s NORAD %s.'], ...
                Architecture(i), ...
                NORAD_CAT_ID(i));

        end

        AtmosphericGasLoss_dB(i) = ...
            Tgas.AtmosphericGasLoss_dB(idx);

    end

    %% --------------------------------------------------------
    % Architecture RF parameters
    %% --------------------------------------------------------

    Frequency_GHz = nan(N,1);

    EIRP_PSD_dBW_MHz = nan(N,1);

    ReceiverFoM_dBK = nan(N,1);

    %% LEO

    idx = Architecture=="LEO";

    Frequency_GHz(idx) = 20.0;

    ReceiverFoM_dBK(idx) = 10.5;

    controlElevation_deg = [10;25;90];

    controlEIRP_dBW_MHz = [ ...
        -2.7
        17.9
        11.3];

    elev = T.Elevation_deg(idx);

    elevClamped = ...
        min( ...
        max(elev,min(controlElevation_deg)), ...
        max(controlElevation_deg));

    EIRP_PSD_dBW_MHz(idx) = ...
        interp1( ...
        controlElevation_deg, ...
        controlEIRP_dBW_MHz, ...
        elevClamped, ...
        'linear');

    %% MEO

    idx = Architecture=="MEO";

    Frequency_GHz(idx) = 20.0;

    EIRP_PSD_dBW_MHz(idx) = ...
        49.7 - 10*log10(100);

    ReceiverFoM_dBK(idx) = ...
        18.4;

    %% GEO

    idx = Architecture=="GEO";

    Frequency_GHz(idx) = ...
        19.95;

    EIRP_PSD_dBW_MHz(idx) = ...
        57.0 - 10*log10(81);

    ReceiverFoM_dBK(idx) = ...
        44.5 - 10*log10(250);

    assert(all(isfinite(Frequency_GHz)));
    assert(all(isfinite(EIRP_PSD_dBW_MHz)));
    assert(all(isfinite(ReceiverFoM_dBK)));

    %% --------------------------------------------------------
    % FSPL
    %% --------------------------------------------------------

    FSPL_dB = ...
        92.45 ...
        + 20*log10(T.SlantRange_km) ...
        + 20*log10(Frequency_GHz);

    %% --------------------------------------------------------
    % Positioning EIRP
    %% --------------------------------------------------------

    beta = 1.023e6;

    beta_MHz = beta/1e6;

    etaPos = 0.01;

    EIRP_positioning_dBW = ...
        EIRP_PSD_dBW_MHz ...
        + 10*log10(beta_MHz) ...
        + 10*log10(etaPos);

    %% --------------------------------------------------------
    % C/N0
    %% --------------------------------------------------------

    CN0_dBHz = ...
        EIRP_positioning_dBW ...
        - FSPL_dB ...
        - AtmosphericGasLoss_dB ...
        + ReceiverFoM_dBK ...
        + 228.6;

    CN0_linear = ...
        10.^(CN0_dBHz/10);

    %% --------------------------------------------------------
    % Pseudorange uncertainty
    %% --------------------------------------------------------

    c = 299792458;

    Tcoh = 20e-3;

    SigmaRho_m = ...
        c ./ ...
        (2*pi*beta .* ...
        sqrt(CN0_linear*Tcoh));

    VarianceRho_m2 = ...
        SigmaRho_m.^2;

    %% --------------------------------------------------------
    % Store
    %% --------------------------------------------------------

    T.Frequency_GHz = ...
        Frequency_GHz;

    T.EIRP_PSD_dBW_MHz = ...
        EIRP_PSD_dBW_MHz;

    T.ReceiverFoM_dBK = ...
        ReceiverFoM_dBK;

    T.AtmosphericGasLoss_dB = ...
        AtmosphericGasLoss_dB;

    T.FSPL_dB = ...
        FSPL_dB;

    T.EIRP_positioning_dBW = ...
        EIRP_positioning_dBW;

    T.CN0_dBHz = ...
        CN0_dBHz;

    T.SigmaRho_m = ...
        SigmaRho_m;

    T.VarianceRho_m2 = ...
        VarianceRho_m2;

end


%% ============================================================
% Weighted positional metric
%% ============================================================

function [bound_m,PDOP,rankH] = ...
    subsetMetrics(idx,Hfull,sigmaRho_m)

    idx = idx(:);

    H = ...
        Hfull(idx,:);

    sigma = ...
        sigmaRho_m(idx);

    rankH = ...
        rank(H);

    bound_m = inf;
    PDOP = inf;

    if rankH < 4
        return;
    end

    %% Geometric PDOP

    A = H.'*H;

    if rcond(A) > 1e-14

        Cgeom = ...
            A\eye(4);

        trGeom = ...
            trace(Cgeom(1:3,1:3));

        if trGeom > 0 && ...
           isfinite(trGeom)

            PDOP = sqrt(trGeom);

        end
    end

    %% Weighted positional bound

    R = ...
        diag(sigma.^2);

    J = ...
        H.'*(R\H);

    if rcond(J) <= 1e-14
        return;
    end

    C = ...
        J\eye(4);

    trPos = ...
        trace(C(1:3,1:3));

    if trPos > 0 && ...
       isfinite(trPos)

        bound_m = sqrt(trPos);

    end

end