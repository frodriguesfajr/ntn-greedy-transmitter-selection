clear; clc;

rootDir = fileparts(fileparts(mfilename('fullpath')));
cd(rootDir);
setup_paths;

oldBase = '/MATLAB Drive/master-thesis-ntn-positioning/matlab_code';
catalogFile = fullfile(oldBase, ...
    'results_scenario4_candidate_pool', ...
    'scenario4_orbital_catalog.mat');

O = load(catalogFile);

epochUTC = O.analysisTimeUTC;
if isempty(epochUTC.TimeZone)
    epochUTC.TimeZone = 'UTC';
end

eop = paper_epoch_eop();

tables = {O.TLEO20, O.TMEO7, O.TGEO7};
archNames = {'LEO','MEO','GEO'};

Architecture = strings(0,1);
NORAD = zeros(0,1);
Name = strings(0,1);
dX_m = zeros(0,1);
dY_m = zeros(0,1);
dZ_m = zeros(0,1);
PositionError3D_m = zeros(0,1);

for a = 1:numel(tables)

    Tref = tables{a};
    arch = archNames{a};

    for k = 1:height(Tref)

        idText = char(string(Tref.NORAD_CAT_ID(k)));
        noradID = str2double(idText);

        tleFile = fullfile(rootDir,'tle',arch,idText,[idText '.tle']);

        v = propagate_tle_vallado(tleFile,epochUTC,eop);

        rRef_m = [ ...
            Tref.X_ECEF_m(k); ...
            Tref.Y_ECEF_m(k); ...
            Tref.Z_ECEF_m(k)];

        dr = v.rECEF_m - rRef_m;

        Architecture(end+1,1) = string(arch);
        NORAD(end+1,1) = noradID;
        Name(end+1,1) = string(Tref.ObjectName(k));
        dX_m(end+1,1) = dr(1);
        dY_m(end+1,1) = dr(2);
        dZ_m(end+1,1) = dr(3);
        PositionError3D_m(end+1,1) = norm(dr);

    end
end

R = table(Architecture,NORAD,Name,dX_m,dY_m,dZ_m,PositionError3D_m);

fprintf('\n============================================================\n');
fprintf('VALLADO vs PAPER ORBITAL CATALOG\n');
fprintf('Epoch: 01-Aug-2026 12:00:00 UTC\n');
fprintf('============================================================\n');

for a = 1:numel(archNames)
    arch = archNames{a};
    idx = R.Architecture == string(arch);

    fprintf('\n%s (%d satellites)\n',arch,sum(idx));
    fprintf('Position error: mean   = %.6f m\n', ...
        mean(R.PositionError3D_m(idx)));
    fprintf('                median = %.6f m\n', ...
        median(R.PositionError3D_m(idx)));
    fprintf('                max    = %.6f m\n', ...
        max(R.PositionError3D_m(idx)));
end

fprintf('\nALL ORBITAL SATELLITES (%d)\n',height(R));
fprintf('Position error: mean   = %.6f m\n',mean(R.PositionError3D_m));
fprintf('                median = %.6f m\n',median(R.PositionError3D_m));
fprintf('                max    = %.6f m\n',max(R.PositionError3D_m));

[~,idxWorst] = sort(R.PositionError3D_m,'descend');
nShow = min(5,height(R));

fprintf('\nWorst %d cases:\n',nShow);
for k = 1:nShow
    j = idxWorst(k);
    fprintf('%s  NORAD %d  %-24s  %.3f m\n', ...
        R.Architecture(j),R.NORAD(j),R.Name(j),R.PositionError3D_m(j));
end

resultsDir = fullfile(rootDir,'results');
if ~exist(resultsDir,'dir')
    mkdir(resultsDir);
end

outFile = fullfile(resultsDir, ...
    'validation_vallado_vs_paper_baseline.csv');
writetable(R,outFile);

fprintf('\nSaved to:\n%s\n',outFile);
