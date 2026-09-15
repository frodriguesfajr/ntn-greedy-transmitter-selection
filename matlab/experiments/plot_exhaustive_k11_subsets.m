%% Plot saved exhaustive K=11 results (no subset recalculation).
% Save in the repository root or experiments folder and run this file.
% Outputs: PDF (vector), PNG (300 dpi), and editable MATLAB FIG.
% Requires the completed MAT file produced by run_exhaustive_k11_subsets.
% FA/BE are displayed as FTS/BTE; change methodLabel if desired.
clear; clc;

methodLabel = 'FTS/BTE';
nBins = 180;
nBestToShow = 20;
fontSize = 11;

relativeFile = fullfile('results','experiment2_selection_robustness', ...
    'experiment2_exhaustive_K11_all_subsets.mat');
rootDir = fileparts(mfilename('fullpath'));
assert(~isempty(rootDir),'Save this script before running it.');
while ~isfile(fullfile(rootDir,relativeFile))
    parentDir = fileparts(rootDir);
    if strcmp(parentDir,rootDir)
        error('Place this script inside the project containing %s.',relativeFile);
    end
    rootDir = parentDir;
end
inputFile = fullfile(rootDir,relativeFile);
outDir = fullfile(rootDir,'results', ...
    'experiment2_selection_robustness','figures_exhaustive_K11');
if ~isfolder(outDir), mkdir(outDir); end

% Load only bounds and metadata, not the 7.7 million index rows.
D = load(inputFile,'AlphaCRB_m','Tsummary','K', ...
    'totalCombinations','target_m','targetTol','BestIndices');
assert(isfield(D,'Tsummary'), ...
    'No final summary found. Complete the exhaustive run first.');
assert(numel(D.AlphaCRB_m)==D.totalCombinations);
assert(~any(isnan(D.AlphaCRB_m)), 'Unwritten or invalid NaN entries.');
assert(D.K==11,'Expected K=11.');
if ~isfield(D,'targetTol'), D.targetTol=1e-12; end
alpha = D.AlphaCRB_m(:);
finiteAlpha = alpha(isfinite(alpha));
assert(~isempty(finiteAlpha),'No finite bounds.');
nFeasible = sum(alpha<=D.target_m+D.targetTol);
assert(nFeasible==D.Tsummary.NTargetReached, ...
    'Summary and saved bounds disagree.');
fa = D.Tsummary.FAAlphaCRB_m;
be = D.Tsummary.BEAlphaCRB_m;
assert(abs(fa-be)<=D.targetTol, ...
    'FA and BE differ: use separate reference labels for this data.');

blue = [0.12 0.36 0.65];
red = [0.75 0.16 0.12];
green = [0.05 0.48 0.30];

%% Figure 1: full population, including the complete right tail.
% Log counts keep rare bins visible. Empty bins remain empty.
fig1 = figure('Color','w','Units','centimeters', ...
    'Position',[2 2 16 10]);
ax1 = axes(fig1); hold(ax1,'on');
edges = linspace(min(finiteAlpha),max(finiteAlpha),nBins+1);
h = histogram(ax1,finiteAlpha,edges,'FaceColor',blue, ...
    'EdgeColor','none','FaceAlpha',0.85);
set(ax1,'YScale','log','FontName','Times New Roman', ...
    'FontSize',fontSize,'LineWidth',0.8,'Box','on', ...
    'TickLabelInterpreter','latex');
t = xline(ax1,D.target_m,'--','Color',red,'LineWidth',1.5);
f = xline(ax1,fa,':','Color',green,'LineWidth',1.6);
maxCount = max(h.Values);
ylim(ax1,[0.8 max(10,2*maxCount)]);
xlim(ax1,[min(finiteAlpha)-0.03 max(finiteAlpha)+0.03]);
grid(ax1,'on'); ax1.GridAlpha=0.15;
xlabel(ax1,'$\alpha_{\mathrm{CRB}}$ (m)','Interpreter','latex');
ylabel(ax1,'Number of subsets (log scale)','Interpreter','latex');
title(ax1,sprintf('All %d subsets, K = %d', ...
    D.totalCombinations,D.K),'FontWeight','normal');
legend(ax1,[h t f],{'Exhaustive population', ...
    sprintf('Target: %.3f m',D.target_m), ...
    sprintf('%s: %.6f m',methodLabel,fa)}, ...
    'Location','northeast','Interpreter','none','Box','off');
% Near-overlapping method/target lines are deliberately not offset.
text(ax1,0.97,0.58,sprintf('%d feasible / %d subsets', ...
    nFeasible,D.totalCombinations),'Units','normalized', ...
    'HorizontalAlignment','right','FontSize',fontSize-1);
exportFigure(fig1,outDir,'fig_exhaustive_K11_distribution');

%% Figure 2: exact best order statistics, not a random sample.
% mink avoids sorting or plotting all 7.7 million points.
nShow = min(nBestToShow,numel(finiteAlpha));
bestAlpha = mink(finiteAlpha,nShow);
fig2 = figure('Color','w','Units','centimeters', ...
    'Position',[3 3 16 10]);
ax2 = axes(fig2); hold(ax2,'on');
p = plot(ax2,1:nShow,bestAlpha,'o-','Color',blue, ...
    'MarkerFaceColor','w','MarkerSize',5,'LineWidth',1);
t2 = yline(ax2,D.target_m,'--','Color',red,'LineWidth',1.5);
f2 = plot(ax2,1,fa,'p','Color',green,'MarkerFaceColor',green, ...
    'MarkerSize',11,'LineWidth',1);
% Only identify rank 1 as FTS/BTE if the saved result supports it.
assert(abs(bestAlpha(1)-fa)<=D.targetTol, ...
    'FTS/BTE is not the minimum for these data.');
set(ax2,'FontName','Times New Roman','FontSize',fontSize, ...
    'LineWidth',0.8,'Box','on','TickLabelInterpreter','latex');
xlim(ax2,[0.5 nShow+0.5]);
xticks(ax2,unique([1 5:5:nShow nShow]));
span = max([bestAlpha(end)-bestAlpha(1), ...
    abs(D.target_m-bestAlpha(1)),1e-4]);
ylim(ax2,[min(bestAlpha(1),D.target_m)-0.15*span, ...
    max(bestAlpha(end),D.target_m)+0.15*span]);
grid(ax2,'on'); ax2.GridAlpha=0.15;
xlabel(ax2,'Subset rank (increasing positional bound)', ...
    'Interpreter','latex');
ylabel(ax2,'$\alpha_{\mathrm{CRB}}$ (m)','Interpreter','latex');
title(ax2,sprintf('Best %d subsets, K = %d',nShow,D.K), ...
    'FontWeight','normal');
legend(ax2,[p t2 f2],{'Exhaustive order statistics', ...
    sprintf('Target: %.3f m',D.target_m),methodLabel}, ...
    'Location','southeast','Interpreter','none','Box','off');
exportFigure(fig2,outDir,'fig_exhaustive_K11_best20');

fprintf('\nTotal subsets: %d\n',D.totalCombinations);
fprintf('Nonfinite bounds: %d\n',sum(~isfinite(alpha)));
fprintf('Feasible (target + tolerance): %d\n',nFeasible);
fprintf('Minimum: %.15f m\n',bestAlpha(1));
if nShow>=2
    fprintf('Second best: %.15f m\n',bestAlpha(2));
    fprintf('Second best minus target: %.15f m\n', ...
        bestAlpha(2)-D.target_m);
end
fprintf('Figures saved in:\n%s\n',outDir);

function exportFigure(fig,outDir,baseName)
    % These named figure exports are replaced if this script is rerun.
    exportgraphics(fig,fullfile(outDir,[baseName '.pdf']), ...
        'ContentType','vector','BackgroundColor','white');
    exportgraphics(fig,fullfile(outDir,[baseName '.png']), ...
        'Resolution',300,'BackgroundColor','white');
    savefig(fig,fullfile(outDir,[baseName '.fig']));
end