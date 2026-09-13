function [T_K, P_hPa, rho_gm3] = itu835_midlatitude_winter(h_km)
%ITU835_MIDLATITUDE_WINTER
% ITU-R P.835-7 mid-latitude winter reference atmosphere.
%
% INPUT
%   h_km : geometric height [km], 0 <= h <= 100
%
% OUTPUTS
%   T_K      : temperature [K]
%   P_hPa    : total barometric pressure [hPa]
%   rho_gm3  : water-vapour density [g/m^3]
%
% Source:
%   ITU-R Recommendation P.835-7, Annex 2,
%   mid-latitude winter atmosphere.

h = h_km;

assert(all(h(:) >= 0 & h(:) <= 100), ...
    'Height must be between 0 and 100 km.');

T_K     = zeros(size(h));
P_hPa   = zeros(size(h));
rho_gm3 = zeros(size(h));

%% Temperature

idx = h < 10;
T_K(idx) = ...
    272.7241 ...
    - 3.6217*h(idx) ...
    - 0.1759*h(idx).^2;

idx = h >= 10 & h < 33;
T_K(idx) = 218;

idx = h >= 33 & h < 47;
T_K(idx) = 218 + 3.3571*(h(idx)-33);

idx = h >= 47 & h < 53;
T_K(idx) = 265;

idx = h >= 53 & h < 80;
T_K(idx) = 265 - 2.0370*(h(idx)-53);

idx = h >= 80;
T_K(idx) = 210;

%% Pressure

P10 = 1018.8627 ...
    - 124.2954*10 ...
    + 4.8307*10^2;

P72 = P10 * exp(-0.147*(72-10));

idx = h <= 10;
P_hPa(idx) = ...
    1018.8627 ...
    - 124.2954*h(idx) ...
    + 4.8307*h(idx).^2;

idx = h > 10 & h <= 72;
P_hPa(idx) = ...
    P10 .* exp(-0.147*(h(idx)-10));

idx = h > 72;
P_hPa(idx) = ...
    P72 .* exp(-0.155*(h(idx)-72));

%% Water-vapour density

idx = h <= 10;
rho_gm3(idx) = ...
    3.4742 .* exp( ...
    -0.2697*h(idx) ...
    -0.03604*h(idx).^2 ...
    +0.0004489*h(idx).^3);

rho_gm3(h > 10) = 0;

end