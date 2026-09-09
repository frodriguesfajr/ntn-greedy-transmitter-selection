function out = propagate_tle_vallado(tleFile, epochUTC, eop)
%PROPAGATE_TLE_VALLADO Propagate a TLE using Vallado SGP4.
%
% Inputs:
%   tleFile  - path to a three-line TLE file
%   epochUTC - MATLAB datetime interpreted as UTC
%   eop      - Earth orientation parameter structure
%
% Output fields include TEME and ECEF position/velocity.

if nargin < 3 || isempty(eop)
    eop = paper_epoch_eop();
end

if ~isfile(tleFile)
    error('TLE file not found: %s',tleFile);
end

txt = fileread(tleFile);
lines = regexp(strtrim(txt), '\r\n|\n|\r', 'split');

if numel(lines) < 3
    error('Expected a name line followed by two TLE lines.');
end

name     = strtrim(lines{1});
longstr1 = char(lines{2});
longstr2 = char(lines{3});

% Standard SGP4 configuration for TLE propagation
typerun    = 'c';
typeinput  = 'e';
opsmode    = 'a';
whichconst = 72;

[~,~,~,satrec] = twoline2rv( ...
    longstr1,longstr2,typerun,typeinput,opsmode,whichconst);

% Convert requested UTC epoch to Julian date
[yy,mo,dd,hh,mi,ss] = datevec(epochUTC);
[jdUTC,jdUTCF] = jday(yy,mo,dd,hh,mi,ss);
jdUTCtotal = jdUTC + jdUTCF;

% Minutes from the TLE epoch
tsince_min = ((jdUTCtotal) - ...
    (satrec.jdsatepoch + satrec.jdsatepochf))*1440.0;

% SGP4 propagation: TEME frame
[satrec,rTEME_km,vTEME_kms] = sgp4(satrec,tsince_min);

if satrec.error ~= 0
    error('SGP4 propagation failed with error code %d.',satrec.error);
end

% UTC -> UT1
jdut1 = jdUTCtotal + eop.dut1_s/86400.0;

% UTC -> TT
tt_utc_s = eop.tai_utc_s + 32.184;
jdtt = jdUTCtotal + tt_utc_s/86400.0;
ttt = (jdtt - 2451545.0)/36525.0;

% TEME -> Earth-fixed frame
ateme = zeros(3,1);

[rECEF_km,vECEF_kms,~] = teme2ecef( ...
    rTEME_km(:),vTEME_kms(:),ateme,ttt,jdut1, ...
    eop.lod_s,eop.xp_rad,eop.yp_rad);

out.name        = name;
out.norad       = satrec.satnum;
out.epochUTC    = epochUTC;
out.tsince_min  = tsince_min;
out.sgp4_error  = satrec.error;
out.rTEME_km    = rTEME_km(:);
out.vTEME_kms   = vTEME_kms(:);
out.rECEF_m     = 1000.0*rECEF_km(:);
out.vECEF_mps   = 1000.0*vECEF_kms(:);

end
