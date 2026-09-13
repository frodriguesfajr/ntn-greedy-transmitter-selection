function A_dB = itu676_slant_path( ...
    f_GHz, elevation_deg, stationHeight_km, topHeight_km)
%ITU676_SLANT_PATH
% Slant-path gaseous attenuation using ITU-R P.676-13 Annex 1
% with the ITU-R P.835-7 mid-latitude winter atmosphere.
%
% INPUTS
%   f_GHz            frequency [GHz]
%   elevation_deg    local elevation angle [deg]
%   stationHeight_km receiver height above mean sea level [km]
%   topHeight_km     upper path height [km]
%
% OUTPUT
%   A_dB             integrated gaseous attenuation [dB]

assert(elevation_deg > 0);
assert(topHeight_km > stationHeight_km);
assert(topHeight_km <= 100);

RE_km = 6371;

%% P.676 exponentially increasing layers

i = (1:922).';

delta_km = 1e-4 .* exp((i-1)/100);

bottom_km = [0; cumsum(delta_km(1:end-1))];
top_km    = bottom_km + delta_km;

%% Clip layers to requested path

lo = max(bottom_km,stationHeight_km);
hi = min(top_km,topHeight_km);

valid = hi > lo;

lo = lo(valid);
hi = hi(valid);

dh = hi-lo;
hm = 0.5*(lo+hi);

%% Atmospheric state at station

% [T0,P0,rho0] = ...
    % itu835_midlatitude_winter(stationHeight_km);
[T0,P0,rho0] =  ...
    itu835_seasonal_profile(stationHeight_km,-22.8596582,"winter");

e0 = rho0*T0/216.7;
p0 = P0-e0;

N0 = ...
    77.6*p0/T0 + ...
    72*e0/T0 + ...
    3.75e5*e0/T0^2;

n0 = 1 + N0*1e-6;

r0 = RE_km + stationHeight_km;

%% Integrate

A_dB = 0;

for k = 1:numel(hm)

    h = hm(k);

    % [T,P,rho] = itu835_midlatitude_winter(h);
    [T,P,rho] = itu835_seasonal_profile(h,-22.8596582,"winter");

    [gamma,~,~] = ...
        itu676_specific_attenuation( ...
        f_GHz,P,T,rho);

    %% Refractive index: ITU-R P.453 formulation

    e = rho*T/216.7;
    p = P-e;

    N = ...
        77.6*p/T + ...
        72*e/T + ...
        3.75e5*e/T^2;

    n = 1 + N*1e-6;

    r = RE_km+h;

    %% P.676 spherical Snell-law geometry

    cosPhi = ...
        (r0*n0)/(r*n) * cosd(elevation_deg);

    cosPhi = min(max(cosPhi,-1),1);

    sinPhi = sqrt(max(1-cosPhi^2,eps));

    ds_km = dh(k)/sinPhi;

    A_dB = A_dB + gamma*ds_km;
end

end