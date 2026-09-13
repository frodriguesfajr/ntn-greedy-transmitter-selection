clear;
clc;

%% Repository paths

thisFile  = mfilename('fullpath');
atmDir    = fileparts(thisFile);
matlabDir = fileparts(atmDir);
repoRoot  = fileparts(matlabDir);

p2145Dir = fullfile( ...
    repoRoot, ...
    'data','itu','p2145','annual');

%% Files

files = struct;

files.P    = fullfile(p2145Dir,'P_Annual','P_50.TXT');
files.PSCH = fullfile(p2145Dir,'P_Annual','PSCH.TXT');
files.ZP   = fullfile(p2145Dir,'P_Annual','Z_ground.TXT');

files.T    = fullfile(p2145Dir,'T_Annual','T_50.TXT');
files.TSCH = fullfile(p2145Dir,'T_Annual','TSCH.TXT');
files.ZT   = fullfile(p2145Dir,'T_Annual','Z_ground.TXT');

files.RHO  = fullfile(p2145Dir,'RHO_Annual','RHO_50.TXT');
files.VSCH = fullfile(p2145Dir,'RHO_Annual','VSCH.TXT');
files.ZRHO = fullfile(p2145Dir,'RHO_Annual','Z_ground.TXT');

files.V    = fullfile(p2145Dir,'V_Annual','V_50.TXT');
files.VV   = fullfile(p2145Dir,'V_Annual','VSCH.TXT');
files.ZV   = fullfile(p2145Dir,'V_Annual','Z_ground.TXT');

%% Check files

names = fieldnames(files);

fprintf('\n');
fprintf('============================================================\n');
fprintf('ITU-R P.2145 MAP CHECK\n');
fprintf('============================================================\n');

for k = 1:numel(names)

    name = names{k};
    file = files.(name);

    assert(isfile(file), ...
        'Missing file: %s',file);

    A = readmatrix(file);

    fprintf('\n%-5s\n',name);
    fprintf('  Size : %d x %d\n',size(A,1),size(A,2));
    fprintf('  Min  : %.8f\n',min(A(:),[],'omitnan'));
    fprintf('  Max  : %.8f\n',max(A(:),[],'omitnan'));

end

fprintf('\n============================================================\n');