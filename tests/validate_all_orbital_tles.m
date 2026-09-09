clear; clc;

rootDir = fileparts(fileparts(mfilename('fullpath')));
cd(rootDir);
setup_paths;

epochUTC = datetime(2026,8,1,12,0,0,'TimeZone','UTC');
eop = paper_epoch_eop();

architectures = {'LEO','MEO','GEO'};

Architecture = strings(0,1);
NORAD = zeros(0,1);
Name = strings(0,1);
dX_m = zeros(0,1);
dY_m = zeros(0,1);
dZ_m = zeros(0,1);
PositionError3D_m = zeros(0,1);
VelocityError3D_mps = zeros(0,1);

for a = 1:numel(architectures)

    arch = architectures{a};
    satDirs = dir(fullfile(rootDir,'tle',arch));
    satDirs = satDirs([satDirs.isdir]);
    satDirs = satDirs(~ismember({satDirs.name},{'.','..'}));

    for k = 1:numel(satDirs)

        noradText = satDirs(k).name;

        tleList = dir(fullfile(rootDir,'tle',arch,noradText,'*.tle'));
        ephList = dir(fullfile(rootDir,'tle','stk_ephemerides',arch,noradText,'*.e'));

        if numel(tleList) ~= 1
            error('Expected exactly one TLE for %s %s.',arch,noradText);
        end

        if numel(ephList) ~= 1
            error('Expected exactly one STK ephemeris for %s %s.',arch,noradText);
        end

        tleFile = fullfile(tleList(1).folder,tleList(1).name);
        ephFile = fullfile(ephList(1).folder,ephList(1).name);

        v = propagate_tle_vallado(tleFile,epochUTC,eop);
        s = read_stk_ephemeris_at_epoch(ephFile,epochUTC);

        dr = v.rECEF_m - s.rECEF_m;
        dv = v.vECEF_mps - s.vECEF_mps;

        Architecture(end+1,1) = string(arch);
        NORAD(end+1,1) = v.norad;
        Name(end+1,1) = string(v.name);
        dX_m(end+1,1) = dr(1);
        dY_m(end+1,1) = dr(2);
        dZ_m(end+1,1) = dr(3);
        PositionError3D_m(end+1,1) = norm(dr);
        VelocityError3D_mps(end+1,1) = norm(dv);

    end
end

T = table(Architecture,NORAD,Name,dX_m,dY_m,dZ_m, ...
    PositionError3D_m,VelocityError3D_mps);

T = sortrows(T,{'Architecture','NORAD'});

fprintf('\n============================================================\n');
fprintf('VALLADO SGP4 vs STK VALIDATION\n');
fprintf('Epoch: 01-Aug-2026 12:00:00 UTC\n');
fprintf('============================================================\n\n');

disp(T)

for a = 1:numel(architectures)
    arch = architectures{a};
    idx = T.Architecture == string(arch);

    fprintf('\n%s (%d satellites)\n',arch,sum(idx));
    fprintf('Position error: mean   = %.6f m\n',mean(T.PositionError3D_m(idx)));
    fprintf('                median = %.6f m\n',median(T.PositionError3D_m(idx)));
    fprintf('                max    = %.6f m\n',max(T.PositionError3D_m(idx)));
    fprintf('Velocity error: max    = %.9f m/s\n',max(T.VelocityError3D_mps(idx)));
end

fprintf('\nALL ORBITAL SATELLITES (%d)\n',height(T));
fprintf('Position error: mean   = %.6f m\n',mean(T.PositionError3D_m));
fprintf('                median = %.6f m\n',median(T.PositionError3D_m));
fprintf('                max    = %.6f m\n',max(T.PositionError3D_m));
fprintf('Velocity error: max    = %.9f m/s\n',max(T.VelocityError3D_mps));

resultsDir = fullfile(rootDir,'results');
if ~exist(resultsDir,'dir')
    mkdir(resultsDir);
end

outFile = fullfile(resultsDir,'validation_vallado_vs_stk.csv');
writetable(T,outFile);

fprintf('\nValidation table saved to:\n%s\n',outFile);
