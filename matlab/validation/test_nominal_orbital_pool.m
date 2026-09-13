clear;
clc;

%% Repository paths
thisFile = mfilename('fullpath');
propDir  = fileparts(thisFile);
matlabDir = fileparts(propDir);
repoRoot = fileparts(matlabDir);

addpath(fullfile(repoRoot,'external','sgp4'));
addpath(propDir);

%% Reference epoch
epochUTC = datetime(2026,8,1,12,0,0, ...
    'TimeZone','UTC');

%% Nominal orbital pool
leo = [ ...
    65175 66864 65348 59340 ...
    67736 59960 65279 48482 ];

meo = [ ...
    62362 40081 44115 43231 ...
    43234 40080 40079 ];

geo = [ ...
    41589 41904 42692 43562 ...
    43175 43228 38087 ];

idsByType = {leo, meo, geo};
types     = {'LEO','MEO','GEO'};

fprintf('\n');
fprintf('============================================================\n');
fprintf('NOMINAL ORBITAL POOL - SGP4 VALIDATION\n');
fprintf('Epoch: 01-Aug-2026 12:00:00 UTC\n');
fprintf('============================================================\n');

nOK = 0;

for g = 1:numel(types)

    ids = idsByType{g};

    fprintf('\n%s\n',types{g});
    fprintf('------------------------------------------------------------\n');

    for k = 1:numel(ids)

        norad = ids(k);

        tleFile = fullfile( ...
            repoRoot,'tle',lower(types{g}), ...
            sprintf('%d',norad), ...
            sprintf('%d.tle',norad));

        assert(isfile(tleFile), ...
            'Missing TLE for NORAD %d: %s',norad,tleFile);

        out = propagate_tle_to_ecef(tleFile,epochUTC);

        fprintf('%5d  %-18s  |r| = %10.3f km   OK\n', ...
            out.norad,out.name,norm(out.rECEF_km));

        nOK = nOK + 1;
    end
end

fprintf('\n============================================================\n');
fprintf('Successfully propagated: %d / 22 satellites\n',nOK);
fprintf('============================================================\n');