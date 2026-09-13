function met = itu2145_annual50_at_location( ...
    p2145Dir, lat_deg, lon_deg, height_m)
%ITU2145_ANNUAL50_AT_LOCATION
% Annual 50%-exceedance meteorological parameters from ITU-R P.2145-0.
%
% INPUTS
%   p2145Dir : directory containing
%              P_Annual/
%              T_Annual/
%              RHO_Annual/
%              V_Annual/
%
%   lat_deg  : geodetic latitude [deg], -90 <= lat <= 90
%   lon_deg  : longitude [deg], -180 <= lon <= 180
%   height_m : station height above mean sea level [m]
%
% OUTPUT
%   met.P_hPa
%   met.T_K
%   met.rho_gm3
%   met.V_kgm2
%   met.PSCH_km
%   met.TSCH_Kkm
%   met.VSCH_km
%   met.Zground_km
%
% Method:
%   ITU-R P.2145-0:
%   1) obtain the four surrounding 0.25-deg grid points;
%   2) correct each grid-point value to the desired station height;
%   3) perform bilinear interpolation of the corrected values.
%
% Annual exceedance probability used here: 50%.

arguments
    p2145Dir (1,:) char
    lat_deg (1,1) double
    lon_deg (1,1) double
    height_m (1,1) double
end

assert(lat_deg >= -90 && lat_deg <= 90, ...
    'Latitude must be between -90 and 90 deg.');

assert(lon_deg >= -180 && lon_deg <= 180, ...
    'Longitude must be between -180 and 180 deg.');

alt_km = height_m/1000;


%% ============================================================
% FILES
% ============================================================

pDir   = fullfile(p2145Dir,'P_Annual');
tDir   = fullfile(p2145Dir,'T_Annual');
rhoDir = fullfile(p2145Dir,'RHO_Annual');
vDir   = fullfile(p2145Dir,'V_Annual');

files.P    = fullfile(pDir,'P_50.TXT');
files.PSCH = fullfile(pDir,'PSCH.TXT');
files.ZP   = fullfile(pDir,'Z_ground.TXT');

files.T    = fullfile(tDir,'T_50.TXT');
files.TSCH = fullfile(tDir,'TSCH.TXT');
files.ZT   = fullfile(tDir,'Z_ground.TXT');

files.RHO  = fullfile(rhoDir,'RHO_50.TXT');
files.VSCH_RHO = fullfile(rhoDir,'VSCH.TXT');
files.ZRHO = fullfile(rhoDir,'Z_ground.TXT');

files.V    = fullfile(vDir,'V_50.TXT');
files.VSCH_V = fullfile(vDir,'VSCH.TXT');
files.ZV   = fullfile(vDir,'Z_ground.TXT');

names = fieldnames(files);

for k = 1:numel(names)
    assert(isfile(files.(names{k})), ...
        'Missing P.2145 file: %s',files.(names{k}));
end


%% ============================================================
% LOAD MAPS
% ============================================================

P       = readmatrix(files.P);
PSCH    = readmatrix(files.PSCH);
ZP      = readmatrix(files.ZP);

T       = readmatrix(files.T);
TSCH    = readmatrix(files.TSCH);
ZT      = readmatrix(files.ZT);

RHO     = readmatrix(files.RHO);
VSCHRHO = readmatrix(files.VSCH_RHO);
ZRHO    = readmatrix(files.ZRHO);

V       = readmatrix(files.V);
VSCHV   = readmatrix(files.VSCH_V);
ZV      = readmatrix(files.ZV);


%% ============================================================
% GRID VALIDATION
% ============================================================

expectedSize = [721 1441];

assert(isequal(size(P),expectedSize));
assert(isequal(size(T),expectedSize));
assert(isequal(size(RHO),expectedSize));
assert(isequal(size(V),expectedSize));

assert(isequal(size(PSCH),expectedSize));
assert(isequal(size(TSCH),expectedSize));
assert(isequal(size(VSCHRHO),expectedSize));
assert(isequal(size(VSCHV),expectedSize));


%% ============================================================
% FOUR SURROUNDING GRID POINTS
%
% P.2145 grid:
%
% latitude  : -90 : 0.25 : +90
% longitude : -180 : 0.25 : +180
%
% MATLAB indices start at 1.
% ============================================================

dGrid = 0.25;

rowFloat = (lat_deg + 90)/dGrid + 1;
colFloat = (lon_deg + 180)/dGrid + 1;

r0 = floor(rowFloat);
c0 = floor(colFloat);

r0 = max(1,min(r0,720));
c0 = max(1,min(c0,1440));

r1 = r0 + 1;
c1 = c0 + 1;

fy = rowFloat-r0;
fx = colFloat-c0;

fy = min(max(fy,0),1);
fx = min(max(fx,0),1);


%% ============================================================
% EXTRACT FOUR CORNERS
% ============================================================

Pc    = corners(P,r0,r1,c0,c1);
PSCHc = corners(PSCH,r0,r1,c0,c1);
ZPc   = corners(ZP,r0,r1,c0,c1);

Tc    = corners(T,r0,r1,c0,c1);
TSCHc = corners(TSCH,r0,r1,c0,c1);
ZTc   = corners(ZT,r0,r1,c0,c1);

RHOc  = corners(RHO,r0,r1,c0,c1);
VSCHRHOc = corners(VSCHRHO,r0,r1,c0,c1);
ZRHOc = corners(ZRHO,r0,r1,c0,c1);

Vc    = corners(V,r0,r1,c0,c1);
VSCHVc = corners(VSCHV,r0,r1,c0,c1);
ZVc   = corners(ZV,r0,r1,c0,c1);


%% ============================================================
% HEIGHT CORRECTION AT EACH GRID POINT
%
% ITU-R P.2145:
%
% P    = P'   exp[-(alt-alt_i)/PSCH]
% T    = T' + TSCH (alt-alt_i)
% rho  = rho' exp[-(alt-alt_i)/VSCH]
% V    = V'   exp[-(alt-alt_i)/VSCH]
%
% alt and scale heights are in km.
% ============================================================

Pcorrected = ...
    Pc .* exp(-(alt_km-ZPc)./PSCHc);

Tcorrected = ...
    Tc + TSCHc.*(alt_km-ZTc);

RHOcorrected = ...
    RHOc .* exp(-(alt_km-ZRHOc)./VSCHRHOc);

Vcorrected = ...
    Vc .* exp(-(alt_km-ZVc)./VSCHVc);


%% ============================================================
% BILINEAR INTERPOLATION
% ============================================================

met.P_hPa = bilinear4(Pcorrected,fx,fy);

met.T_K = bilinear4(Tcorrected,fx,fy);

met.rho_gm3 = bilinear4(RHOcorrected,fx,fy);

met.V_kgm2 = bilinear4(Vcorrected,fx,fy);


%% ============================================================
% AUXILIARY INTERPOLATED PARAMETERS
% ============================================================

met.PSCH_km = ...
    bilinear4(PSCHc,fx,fy);

met.TSCH_Kkm = ...
    bilinear4(TSCHc,fx,fy);

met.VSCH_km = ...
    bilinear4(VSCHRHOc,fx,fy);

met.Zground_km = ...
    bilinear4(ZPc,fx,fy);


%% ============================================================
% METADATA
% ============================================================

met.latitude_deg  = lat_deg;
met.longitude_deg = lon_deg;
met.height_m      = height_m;

met.exceedance_pct = 50;

met.gridRows = [r0 r1];
met.gridCols = [c0 c1];

met.gridLatitude_deg = [ ...
    -90+(r0-1)*dGrid ...
    -90+(r1-1)*dGrid ];

met.gridLongitude_deg = [ ...
    -180+(c0-1)*dGrid ...
    -180+(c1-1)*dGrid ];

end


%% ============================================================
% LOCAL FUNCTIONS
% ============================================================

function x = corners(A,r0,r1,c0,c1)

x = [ ...
    A(r0,c0)
    A(r1,c0)
    A(r0,c1)
    A(r1,c1) ];

end


function y = bilinear4(x,fx,fy)
% x order:
%   [lat0/lon0
%    lat1/lon0
%    lat0/lon1
%    lat1/lon1]

y = ...
    (1-fx)*(1-fy)*x(1) + ...
    (1-fx)*fy    *x(2) + ...
    fx    *(1-fy)*x(3) + ...
    fx    *fy    *x(4);

end