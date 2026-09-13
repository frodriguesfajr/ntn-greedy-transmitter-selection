clear;
clc;

%% Repository paths

thisFile  = mfilename('fullpath');
atmDir    = fileparts(thisFile);
matlabDir = fileparts(atmDir);
repoRoot  = fileparts(matlabDir);

addpath(atmDir);

p2145Dir = fullfile( ...
    repoRoot, ...
    'data','itu','p2145','annual');

%% Receiver used in paper simulations

lat_deg  = -22.8596582;
lon_deg  = -43.2303236;
height_m = 10;

%% P.2145 annual 50%

met = itu2145_annual50_at_location( ...
    p2145Dir, ...
    lat_deg, ...
    lon_deg, ...
    height_m);

%% Report

fprintf('\n');
fprintf('============================================================\n');
fprintf('ITU-R P.2145 ANNUAL 50%% - RIO DE JANEIRO\n');
fprintf('============================================================\n');

fprintf('Latitude          : %.7f deg\n',met.latitude_deg);
fprintf('Longitude         : %.7f deg\n',met.longitude_deg);
fprintf('Station height    : %.3f km\n',met.height_m/1000);

fprintf('\n');
fprintf('Pressure          : %.6f hPa\n',met.P_hPa);
fprintf('Temperature       : %.6f K\n',met.T_K);
fprintf('Water-vapour rho  : %.6f g/m^3\n',met.rho_gm3);
fprintf('Integrated vapour : %.6f kg/m^2\n',met.V_kgm2);

fprintf('\n');
fprintf('Pressure scale    : %.6f km\n',met.PSCH_km);
fprintf('Temperature scale : %.6f K/km\n',met.TSCH_Kkm);
fprintf('Vapour scale      : %.6f km\n',met.VSCH_km);
fprintf('Ground height     : %.6f km\n',met.Zground_km);

fprintf('\n');
fprintf('Latitude bracket  : %.2f ... %.2f deg\n', ...
    met.gridLatitude_deg);

fprintf('Longitude bracket : %.2f ... %.2f deg\n', ...
    met.gridLongitude_deg);

fprintf('============================================================\n');