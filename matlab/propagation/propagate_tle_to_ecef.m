function out = propagate_tle_to_ecef(tleFile, epochUTC)
%PROPAGATE_TLE_TO_ECEF Propagate one TLE with Vallado SGP4.
%
% out = propagate_tle_to_ecef(tleFile, epochUTC)
%
% INPUTS
%   tleFile  : path to TLE file
%   epochUTC : datetime in UTC
%
% OUTPUT
%   out.name
%   out.norad
%   out.epochUTC
%   out.rTEME_km
%   out.vTEME_kmps
%   out.rECEF_km
%   out.vECEF_kmps
%   out.tsince_min
%   out.sgp4Error

arguments
    tleFile (1,:) char
    epochUTC (1,1) datetime
end

%% Read TLE
lines = readlines(tleFile);
lines = strip(lines);
lines(lines == "") = [];

if numel(lines) < 3
    error('Expected satellite name and two TLE lines in: %s', tleFile);
end

satName  = char(lines(1));
longstr1 = char(lines(2));
longstr2 = char(lines(3));

norad = str2double(strtrim(longstr1(3:7)));

%% Initialize SGP4
typerun    = 'c';
typeinput  = 'e';
opsmode    = 'i';
whichconst = 72;     % WGS-72

[~, ~, ~, satrec] = twoline2rv( ...
    longstr1, longstr2, ...
    typerun, typeinput, opsmode, whichconst);

%% Target Julian date
epochUTC.TimeZone = 'UTC';

yr = year(epochUTC);
mo = month(epochUTC);
dy = day(epochUTC);
hr = hour(epochUTC);
mn = minute(epochUTC);
sc = second(epochUTC);

[jd, jdfrac] = jday(yr, mo, dy, hr, mn, sc);
jdTarget = jd + jdfrac;

%% TLE epoch
jdEpoch = satrec.jdsatepoch;

if isfield(satrec, 'jdsatepochF')
    jdEpoch = jdEpoch + satrec.jdsatepochF;
elseif isfield(satrec, 'jdsatepochf')
    jdEpoch = jdEpoch + satrec.jdsatepochf;
end

%% Propagate with SGP4
tsince = (jdTarget - jdEpoch) * 1440.0;

[satrec, rTEME, vTEME] = sgp4(satrec, tsince);

if satrec.error ~= 0
    error('SGP4 error %d for NORAD %d.', satrec.error, norad);
end

rTEME = rTEME(:);
vTEME = vTEME(:);

%% TEME -> ECEF
% Portable Earth-orientation approximation:
% UT1-UTC = 0, polar motion = 0, LOD = 0.

TTminusUTC = 69.184;     % s
jdTT = jdTarget + TTminusUTC/86400;

ttt   = (jdTT - 2451545.0) / 36525.0;
jdUT1 = jdTarget;

lod = 0.0;
xp  = 0.0;
yp  = 0.0;

aTEME = zeros(3,1);

[rECEF, vECEF, ~] = teme2ecef( ...
    rTEME, vTEME, aTEME, ...
    ttt, jdUT1, lod, xp, yp);

%% Output
out.name        = satName;
out.norad       = norad;
out.epochUTC    = epochUTC;
out.rTEME_km    = rTEME;
out.vTEME_kmps  = vTEME;
out.rECEF_km    = rECEF(:);
out.vECEF_kmps  = vECEF(:);
out.tsince_min  = tsince;
out.sgp4Error   = satrec.error;

end