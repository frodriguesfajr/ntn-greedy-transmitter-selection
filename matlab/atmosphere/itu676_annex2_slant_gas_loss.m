function A_dB = itu676_annex2_slant_gas_loss( ...
    f_GHz, elevation_deg, stationHeight_km, met)
%ITU676_ANNEX2_SLANT_GAS_LOSS
% Statistical slant-path gaseous attenuation according to
% ITU-R P.676-13 Annex 2.
%
% INPUTS
%   f_GHz            : frequency [GHz]
%   elevation_deg    : elevation angle [deg]
%   stationHeight_km : station height above mean sea level [km]
%
%   met.P_hPa        : total surface pressure [hPa]
%   met.T_K          : surface temperature [K]
%   met.rho_gm3      : surface water-vapour density [g/m^3]
%   met.V_kgm2       : integrated water-vapour content [kg/m^2]
%
% OUTPUT
%   A_dB             : slant-path gaseous attenuation [dB]
%
% The implementation below is specialized for the frequencies
% used in the paper: 19.95, 20.0, and 31.15 GHz.

assert(elevation_deg >= 5 && elevation_deg <= 90, ...
    'P.676 Annex 2 requires elevation between 5 and 90 deg.');

P   = met.P_hPa;
T   = met.T_K;
rho = met.rho_gm3;
V   = met.V_kgm2;

%% ============================================================
% SURFACE SPECIFIC ATTENUATION
% ============================================================

[~,gammaO_dBkm,~] = ...
    itu676_specific_attenuation( ...
        f_GHz,P,T,rho);


%% ============================================================
% OXYGEN EQUIVALENT HEIGHT
%
% P.676-13 Annex 2, Eq. (31):
%
% h0 = a0(f) + b0(f) T + c0(f) Ps + d0(f) rho
%
% The coefficient points below are the P.676-13 Part 1 values
% surrounding the three frequencies used in the paper.
% ============================================================

if f_GHz >= 19.5 && f_GHz <= 20.0

    fTab = [19.5 20.0];

    a0Tab = [ ...
        -2.353538 ...
        -2.357536 ];

    b0Tab = [ ...
         2.766938e-2 ...
         2.768932e-2 ];

    c0Tab = [ ...
        -6.005149e-4 ...
        -6.014528e-4 ];

    d0Tab = [ ...
        -7.529316e-4 ...
        -7.629822e-4 ];

elseif f_GHz >= 31.0 && f_GHz <= 31.5

    fTab = [31.0 31.5];

    a0Tab = [ ...
        -2.461645 ...
        -2.466936 ];

    b0Tab = [ ...
         2.819933e-2 ...
         2.822499e-2 ];

    c0Tab = [ ...
        -6.230776e-4 ...
        -6.241492e-4 ];

    d0Tab = [ ...
        -1.035226e-3 ...
        -1.049731e-3 ];

else
    error(['Frequency %.3f GHz is outside the coefficient ' ...
           'ranges implemented for this paper.'],f_GHz);
end

a0 = interp1(fTab,a0Tab,f_GHz,'linear');
b0 = interp1(fTab,b0Tab,f_GHz,'linear');
c0 = interp1(fTab,c0Tab,f_GHz,'linear');
d0 = interp1(fTab,d0Tab,f_GHz,'linear');

h0_km = ...
    a0 ...
    + b0*T ...
    + c0*P ...
    + d0*rho;

A0_zenith_dB = ...
    gammaO_dBkm*h0_km;


%% ============================================================
% WATER-VAPOUR ZENITH ATTENUATION
%
% P.676-13 Annex 2 statistical method based on integrated
% water-vapour content V.
% ============================================================

fRef_GHz = 20.6;
pRefDry_hPa = 845;

rhoRef_gm3 = V/2.38;

tRef_C = ...
    14*log(0.22*V/2.38) + 3;

TRef_K = tRef_C + 273.15;

% itu676_specific_attenuation expects total pressure.
eRef_hPa = rhoRef_gm3*TRef_K/216.7;

PRefTotal_hPa = ...
    pRefDry_hPa + eRef_hPa;

[~,~,gammaW_f] = ...
    itu676_specific_attenuation( ...
        f_GHz, ...
        PRefTotal_hPa, ...
        TRef_K, ...
        rhoRef_gm3);

[~,~,gammaW_ref] = ...
    itu676_specific_attenuation( ...
        fRef_GHz, ...
        PRefTotal_hPa, ...
        TRef_K, ...
        rhoRef_gm3);

Aw_zenith_dB = ...
    0.0176*V*gammaW_f/gammaW_ref;


%% ============================================================
% HEIGHT CORRECTION FOR WATER VAPOUR
% ============================================================

if f_GHz > 20

    a = ...
        0.2048*exp(-((f_GHz-22.43)/3.097)^2) ...
        + 0.2326*exp(-((f_GHz-183.5)/4.096)^2) ...
        + 0.2073*exp(-((f_GHz-325)/3.651)^2) ...
        - 0.1113;

    b = ...
        8.741e4*exp(-0.587*f_GHz) ...
        + 312.2*f_GHz^(-2.38) ...
        + 0.723;

    h = min(max(stationHeight_km,0),4);

    Aw_zenith_dB = ...
        Aw_zenith_dB*(1 + a*h^b);
end


%% ============================================================
% ZENITH -> SLANT PATH
% ============================================================

A_dB = ...
    (A0_zenith_dB + Aw_zenith_dB) ...
    / sind(elevation_deg);

end