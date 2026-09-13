function [T_K,P_hPa,rho_gm3] = ...
    itu835_seasonal_profile(h_km,lat_deg,season)
%ITU835_SEASONAL_PROFILE
% Seasonal atmosphere according to ITU-R P.835-7.
%
% For 15 < |latitude| < 45 deg, P.835-7 specifies linear
% interpolation between the low-latitude annual atmosphere
% and the appropriate mid-latitude seasonal atmosphere.
%
% This implementation currently supports the winter case used
% for Rio de Janeiro in August.
%
% INPUTS
%   h_km    geometric height [km]
%   lat_deg geodetic latitude [deg]
%   season  "winter"
%
% OUTPUTS
%   T_K      temperature [K]
%   P_hPa    total pressure [hPa]
%   rho_gm3  water-vapour density [g/m^3]

arguments
    h_km
    lat_deg (1,1) double
    season (1,1) string = "winter"
end

assert(all(h_km(:) >= 0 & h_km(:) <= 100), ...
    'Height must be between 0 and 100 km.');

if season ~= "winter"
    error('This implementation currently supports winter only.');
end

absLat = abs(lat_deg);

%% Low-latitude annual atmosphere
[Tlow,Plow,rholow] = lowLatitudeAnnual(h_km);

if absLat <= 15

    T_K     = Tlow;
    P_hPa   = Plow;
    rho_gm3 = rholow;
    return;

end

%% Mid-latitude winter atmosphere
[Tmid,Pmid,rhomid] = midLatitudeWinter(h_km);

if absLat >= 45

    T_K     = Tmid;
    P_hPa   = Pmid;
    rho_gm3 = rhomid;
    return;

end

%% P.835-7 latitude interpolation
w = (absLat - 15)/(45 - 15);

T_K     = (1-w).*Tlow     + w.*Tmid;
P_hPa   = (1-w).*Plow     + w.*Pmid;
rho_gm3 = (1-w).*rholow   + w.*rhomid;

end


%% ============================================================
% LOW-LATITUDE ANNUAL PROFILE
% ITU-R P.835-7, Annex 2
% ============================================================

function [T,P,rho] = lowLatitudeAnnual(h)

T   = zeros(size(h));
P   = zeros(size(h));
rho = zeros(size(h));

%% Temperature

idx = h < 17;
T(idx) = ...
    300.4222 ...
    - 6.3533*h(idx) ...
    + 0.005886*h(idx).^2;

idx = h >= 17 & h < 47;
T(idx) = 194 + 2.533*(h(idx)-17);

idx = h >= 47 & h < 52;
T(idx) = 270;

idx = h >= 52 & h < 80;
T(idx) = 270 - 3.0714*(h(idx)-52);

idx = h >= 80;
T(idx) = 184;

%% Pressure

P10 = ...
    1012.0306 ...
    - 109.0338*10 ...
    + 3.6316*10^2;

P72 = P10*exp(-0.147*(72-10));

idx = h <= 10;
P(idx) = ...
    1012.0306 ...
    - 109.0338*h(idx) ...
    + 3.6316*h(idx).^2;

idx = h > 10 & h <= 72;
P(idx) = ...
    P10 .* exp(-0.147*(h(idx)-10));

idx = h > 72;
P(idx) = ...
    P72 .* exp(-0.165*(h(idx)-72));

%% Water-vapour density

idx = h <= 15;

rho(idx) = ...
    19.6542 .* exp( ...
    -0.2313*h(idx) ...
    -0.1122*h(idx).^2 ...
    +0.01351*h(idx).^3 ...
    -0.0005923*h(idx).^4);

rho(h > 15) = 0;

end


%% ============================================================
% MID-LATITUDE WINTER PROFILE
% ITU-R P.835-7, Annex 2
% ============================================================

function [T,P,rho] = midLatitudeWinter(h)

T   = zeros(size(h));
P   = zeros(size(h));
rho = zeros(size(h));

%% Temperature

idx = h < 10;
T(idx) = ...
    272.7241 ...
    - 3.6217*h(idx) ...
    - 0.1759*h(idx).^2;

idx = h >= 10 & h < 33;
T(idx) = 218;

idx = h >= 33 & h < 47;
T(idx) = 218 + 3.3571*(h(idx)-33);

idx = h >= 47 & h < 53;
T(idx) = 265;

idx = h >= 53 & h < 80;
T(idx) = 265 - 2.0370*(h(idx)-53);

idx = h >= 80;
T(idx) = 210;

%% Pressure

P10 = ...
    1018.8627 ...
    - 124.2954*10 ...
    + 4.8307*10^2;

P72 = P10*exp(-0.147*(72-10));

idx = h <= 10;
P(idx) = ...
    1018.8627 ...
    - 124.2954*h(idx) ...
    + 4.8307*h(idx).^2;

idx = h > 10 & h <= 72;
P(idx) = ...
    P10 .* exp(-0.147*(h(idx)-10));

idx = h > 72;
P(idx) = ...
    P72 .* exp(-0.155*(h(idx)-72));

%% Water-vapour density

idx = h <= 10;

rho(idx) = ...
    3.4742 .* exp( ...
    -0.2697*h(idx) ...
    -0.03604*h(idx).^2 ...
    +0.0004489*h(idx).^3);

rho(h > 10) = 0;

end