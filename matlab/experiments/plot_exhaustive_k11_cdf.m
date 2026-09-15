%% CDF of ALL saved K=11 subsets -- no new bound evaluations.
% Put this file in experiments and run plot_exhaustive_k11_cdf.
% Uses base MATLAB only. All finite observations are plotted (no sampling).
% PDF uses a high-resolution image to avoid millions of vector segments.
clear; clc;

methodLabel = 'FTS/BTE';
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
outDir = fullfile(rootDir,'results', ...
    'experiment2_selection_robustness','figures_exhaustive_K11');
if ~isfolder(outDir), mkdir(outDir); end

fprintf('Loading saved bounds (no exhaustive recalculation)...\n');
D = load(fullfile(rootDir,relativeFile),'AlphaCRB_m','Tsummary', ...
    'totalCombinations','K','target_m','targetTol');
assert(isfield(D,'Tsummary'),'Complete the exhaustive run first.');
N = D.totalCombinations;
target = D.target_m;
assert(numel(D.AlphaCRB_m)==N);
assert(~any(isnan(D.AlphaCRB_m)),'NaN entries: check the saved run.');
assert(D.K==11,'Expected K=11.');
if ~isfield(D,'targetTol'), D.targetTol=1e-12; end
alpha = D.AlphaCRB_m(:);
fa = D.Tsummary.FAAlphaCRB_m;
be = D.Tsummary.BEAlphaCRB_m;
assert(isfinite(fa) && abs(fa-be)<=D.targetTol, ...
    'This combined marker requires coincident finite FA/BE bounds.');

% Strict CDF: fraction of ALL subsets with bound <= x.
% An Inf entry remains in the denominator but is not plotted at finite x.
nTarget = sum(alpha<=D.target_m);
nTargetTol = sum(alpha<=D.target_m+D.targetTol);
assert(nTargetTol==D.Tsummary.NTargetReached);
fTarget = nTarget/N;
fMethod = sum(alpha<=fa)/N;
fprintf('Sorting all finite bounds...\n');
a = sort(alpha(isfinite(alpha)));
assert(~isempty(a));
assert(abs(a(1)-fa)<=D.targetTol,'Method is not the global K=11 minimum.');
clear alpha D;

% Keep the LAST rank at each tied value: correct right-continuous CDF.
lastOfTie = [find(diff(a)~=0);numel(a)];
x = a(lastOfTie);
F = double(lastOfTie)/N;
clear a lastOfTie;

fig = figure('Color','w','Units','centimeters', ...
    'Position',[2 2 14 9]);
ax = axes(fig,'Position',[0.15 0.19 0.81 0.76]);
hold(ax,'on');
set(ax,'YScale','log','FontName','Times New Roman', ...
    'FontSize',11,'LineWidth',0.8,'Box','on', ...
    'TickLabelInterpreter','latex','YMinorGrid','off', ...
    'XMinorGrid','off','YMinorTick','off');
blue = [0.10 0.35 0.65];
red = [0.75 0.15 0.12];
green = [0.05 0.48 0.30];

fprintf('Drawing the complete CDF...\n');
c = stairs(ax,x,F,'Color',blue,'LineWidth',1.3);
t = xline(ax,target,'--','Color',red,'LineWidth',1.2);
assert(target==0.6,'Change the annotation for this dataset.');
m = plot(ax,fa,fMethod,'p','MarkerSize',10, ...
    'MarkerFaceColor',green,'MarkerEdgeColor',green);
xlim(ax,[min(x(1),0.6)-0.04 x(end)+0.04]);
ylim(ax,[0.65/N 1.25]);
yticks(ax,10.^(-7:0));
grid(ax,'on'); ax.GridAlpha=0.16;
xlabel(ax,'$\alpha_{\mathrm{CRB}}$ (m)','Interpreter','latex');
ylabel(ax,'CDF (log scale)','Interpreter','latex');
legend(ax,[c t m],{'All K = 11 subsets','Target: 0.6 m',methodLabel}, ...
    'Location','southeast','Box','off','Interpreter','none', ...
    'FontSize',10);
if nTarget==1
    note = sprintf('Only one feasible subset\nF(0.6) = 1 / %d',N);
else
    note = sprintf('%d feasible subsets\nF(0.6) = %.4g',nTarget,fTarget);
end
text(ax,0.36,0.38,note,'Units','normalized', ...
    'FontSize',10,'FontName','Times New Roman', ...
    'VerticalAlignment','middle','BackgroundColor','w');
drawnow;

% Existing exports with these names are replaced when rerunning.
base = fullfile(outDir,'fig_exhaustive_K11_cdf');
fprintf('Exporting PNG and PDF...\n');
exportgraphics(fig,[base '.png'],'Resolution',600, ...
    'BackgroundColor','white');
exportgraphics(fig,[base '.pdf'],'ContentType','image', ...
    'Resolution',600,'BackgroundColor','white');
% Optional editable FIG (large, since it contains all CDF points):
% savefig(fig,[base '.fig']);
fprintf('FTS/BTE bound: %.15f m\n',fa);
fprintf('F(target), strict <=: %.15g (%d / %d)\n',fTarget,nTarget,N);
fprintf('Feasible with stored tolerance: %d\n',nTargetTol);
fprintf('CDF COMPLETE. Outputs:\n%s.png\n%s.pdf\n',base,base);

%% Formato final para o paper
xlim(ax,[0.58 1.20]);
xticks(ax,0.6:0.1:1.2);
ylim(ax,[5e-8 1.25]);

set(ax,'Position',[0.17 0.20 0.79 0.75], ...
    'FontSize',9);

ax.XLabel.FontSize = 10;
ax.YLabel.FontSize = 10;

lgd = legend(ax);
lgd.FontSize = 8;
lgd.Location = 'southeast';

set(findobj(ax,'Type','text'),'FontSize',8);

set(fig,'PaperUnits','centimeters', ...
    'PaperPositionMode','manual', ...
    'PaperPosition',[0 0 8.8 6.5], ...
    'PaperSize',[8.8 6.5]);

drawnow;

print(fig,[base '_compact.png'],'-dpng','-r600');
print(fig,[base '_compact.pdf'],'-dpdf','-r600');