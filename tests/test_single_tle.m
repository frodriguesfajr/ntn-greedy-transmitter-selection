clear; clc;

rootDir = fileparts(fileparts(mfilename('fullpath')));
cd(rootDir);
setup_paths;

epochUTC = datetime(2026,8,1,12,0,0);
eop = paper_epoch_eop();

tleFile = fullfile(rootDir,'tle','LEO','48482','48482.tle');

out = propagate_tle_vallado(tleFile,epochUTC,eop);

fprintf('\nSatellite : %s\n',out.name);
fprintf('NORAD ID  : %d\n',out.norad);
fprintf('tsince    : %.6f min\n',out.tsince_min);
fprintf('SGP4 error: %d\n',out.sgp4_error);

fprintf('\nECEF position [m]\n');
fprintf('X = %.3f\n',out.rECEF_m(1));
fprintf('Y = %.3f\n',out.rECEF_m(2));
fprintf('Z = %.3f\n',out.rECEF_m(3));

% Independent STK reference at the same epoch
rSTK_m = [ ...
     3903419.237863; ...
    -4983805.074331; ...
    -2623393.702025];

dr_m = out.rECEF_m - rSTK_m;

fprintf('\nVallado vs STK\n');
fprintf('dX = %.6f m\n',dr_m(1));
fprintf('dY = %.6f m\n',dr_m(2));
fprintf('dZ = %.6f m\n',dr_m(3));
fprintf('3-D position difference = %.6f m\n',norm(dr_m));
