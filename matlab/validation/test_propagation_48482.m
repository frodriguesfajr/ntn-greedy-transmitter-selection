clear;
clc;

%% Repository paths
thisFile = mfilename('fullpath');
propagationDir = fileparts(thisFile);
matlabDir = fileparts(propagationDir);
repoRoot = fileparts(matlabDir);

sgp4Dir = fullfile(repoRoot, 'external', 'sgp4');
addpath(sgp4Dir);

%% TLE
norad = 48482;
tleFile = fullfile(repoRoot, 'tle', 'leo', ...
    sprintf('%d', norad), sprintf('%d.tle', norad));

assert(isfile(tleFile), 'TLE file not found: %s', tleFile);

lines = readlines(tleFile);
lines = strip(lines);
lines(lines == "") = [];

assert(numel(lines) >= 3, ...
    'Expected name + two TLE lines in %s.', tleFile);

satName  = char(lines(1));
longstr1 = char(lines(2));
longstr2 = char(lines(3));

%% Vallado SGP4 initialization
% WGS-72 constants are conventionally used with TLE/SGP4.
typerun   = 'c';
typeinput = 'e';
opsmode   = 'i';
whichconst = 72;

[~, ~, ~, satrec] = twoline2rv( ...
    longstr1, longstr2, ...
    typerun, typeinput, opsmode, whichconst);

%% Target epoch: 01-Aug-2026 12:00:00 UTC
year   = 2026;
month  = 8;
day    = 1;
hour   = 12;
minute = 0;
second = 0;

[jdUTC, jdUTCfrac] = jday( ...
    year, month, day, hour, minute, second);

jdTarget = jdUTC + jdUTCfrac;

%% Recover TLE epoch
if ~isfield(satrec, 'jdsatepoch')
    error('satrec.jdsatepoch was not created by twoline2rv.');
end

jdEpoch = satrec.jdsatepoch;

% Vallado versions may store the fractional part separately.
if isfield(satrec, 'jdsatepochF')
    jdEpoch = jdEpoch + satrec.jdsatepochF;
elseif isfield(satrec, 'jdsatepochf')
    jdEpoch = jdEpoch + satrec.jdsatepochf;
end

%% Propagation interval [min]
tsince = (jdTarget - jdEpoch) * 1440.0;

%% SGP4 propagation: TEME
[satrec, rTEME, vTEME] = sgp4(satrec, tsince);

if satrec.error ~= 0
    error('SGP4 propagation failed with error code %d.', satrec.error);
end

rTEME = rTEME(:);   % km
vTEME = vTEME(:);   % km/s

%% TEME -> ECEF
%
% Initial portable approximation:
%   UT1 - UTC = 0
%   polar motion xp = yp = 0
%   LOD = 0
%
% TT - UTC = 69.184 s (TAI-UTC = 37 s).
%
% These assumptions are sufficient for the first propagation
% validation and can later be replaced by Earth-orientation data.

TTminusUTC = 69.184;              % s
jdTT = jdTarget + TTminusUTC/86400;

ttt   = (jdTT - 2451545.0) / 36525.0;
jdUT1 = jdTarget;

lod = 0.0;     % s
xp  = 0.0;     % rad
yp  = 0.0;     % rad

aTEME = zeros(3,1);

[rECEF, vECEF, ~] = teme2ecef( ...
    rTEME, vTEME, aTEME, ...
    ttt, jdUT1, lod, xp, yp);

%% Report
fprintf('\n');
fprintf('============================================================\n');
fprintf('VALLADO SGP4 PROPAGATION TEST\n');
fprintf('============================================================\n');
fprintf('Satellite      : %s\n', satName);
fprintf('NORAD ID       : %d\n', norad);
fprintf('Target epoch   : 01-Aug-2026 12:00:00 UTC\n');
fprintf('Minutes from TLE epoch: %.6f\n', tsince);

fprintf('\nTEME position [km]\n');
fprintf('X = %14.6f\n', rTEME(1));
fprintf('Y = %14.6f\n', rTEME(2));
fprintf('Z = %14.6f\n', rTEME(3));
fprintf('|r| = %14.6f km\n', norm(rTEME));

fprintf('\nECEF position [km]\n');
fprintf('X = %14.6f\n', rECEF(1));
fprintf('Y = %14.6f\n', rECEF(2));
fprintf('Z = %14.6f\n', rECEF(3));
fprintf('|r| = %14.6f km\n', norm(rECEF));

fprintf('\nSGP4 status    : OK\n');
fprintf('============================================================\n');