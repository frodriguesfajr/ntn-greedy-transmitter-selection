clear;
clc;

%% ============================================================
%  DIAGNOSE GASEOUS-ATTENUATION DIFFERENCES - FINAL
%
%  Validation only.
%
%  Purpose
%  -------
%  This script separates two effects that must not be mixed:
%
%    1) meteorological-input convention used by MATLAB R2023b;
%    2) ITU-R P.676 revision used by MATLAB R2023b versus the
%       portable implementation in this repository.
%
%  MATLAB R2023b p618PropagationLosses documents:
%    - Pressure input as DRY-AIR pressure;
%    - ITU-R P.676 (08/2019), i.e. P.676-12.
%
%  The public portable implementation uses:
%    - P.2145 surface meteorology;
%    - ITU-R P.676-13 Annex 2.
%
%  Therefore, direct equality between MATLAB R2023b and the
%  portable P.676-13 implementation is NOT expected.
%
%  This script compares:
%
%    1) Original/reference MATLAB gaseous losses
%    2) MATLAB with internally obtained meteorology
%    3) MATLAB with P.2145 meteorology, using the documented
%       R2023b DRY-AIR pressure convention
%    4) Portable P.2145 + P.676-13 implementation
%    5) A validation-only compatibility reconstruction that
%       replaces ONLY the P.676-13 oxygen equivalent height h0
%       by the P.676-12 Annex 2 expression used by R2023b
%
%  IMPORTANT
%  ---------
%  The compatibility reconstruction is for validation only.
%  It must NOT replace the public P.676-13 implementation.
%
%  Satellite Communications Toolbox is required only for this
%  validation script. The public simulation pipeline does not
%  require it.
%  ============================================================


%% Repository paths

thisFile = mfilename('fullpath');
validationDir = fileparts(thisFile);
matlabDir = fileparts(validationDir);
repoRoot = fileparts(matlabDir);

addpath(fullfile(repoRoot,'matlab','atmosphere'));


%% MATLAB information

fprintf('\n');
fprintf('============================================================\n');
fprintf('MATLAB / TOOLBOX INFORMATION\n');
fprintf('============================================================\n');

fprintf('MATLAB: %s\n',version);
releaseTag = version('-release');
fprintf('Release: %s\n',releaseTag);
isR2023b = strcmpi(releaseTag,'2023b') || strcmpi(releaseTag,'R2023b');

v = ver;

idxSatcom = strcmpi({v.Name}, ...
    'Satellite Communications Toolbox');

if any(idxSatcom)
    fprintf('Satellite Communications Toolbox: %s\n', ...
        v(find(idxSatcom,1)).Version);
else
    error(['Satellite Communications Toolbox is required only ' ...
           'for this validation script.']);
end

if ~isR2023b
    warning([ ...
        'This diagnostic is calibrated to MATLAB R2023b. ' ...
        'Other releases can use different ITU-R revisions or ' ...
        'input conventions.']);
end


%% Candidate pool

P = load(fullfile( ...
    repoRoot,'data','generated', ...
    'nominal_candidate_pool.mat'));

candidatePool = P.candidatePool;
config = P.config;

N = height(candidatePool);


%% P.2145 meteorology from the portable implementation

p2145Dir = fullfile( ...
    repoRoot,'data','itu','p2145','annual');

met = itu2145_annual50_at_location( ...
    p2145Dir, ...
    config.receiver.lat_deg, ...
    config.receiver.lon_deg, ...
    config.receiver.h_m);

% met.P_hPa is the TOTAL barometric surface pressure used by
% the portable P.676-13 implementation.
Ptotal_hPa = met.P_hPa;
T_K        = met.T_K;
rho_gm3    = met.rho_gm3;
V_kgm2     = met.V_kgm2;

% Water-vapour partial pressure, P.676 convention.
e_hPa = rho_gm3*T_K/216.7;

% MATLAB R2023b documents the p618PropagationLosses Pressure
% argument as DRY-AIR pressure.
Pdry_hPa = Ptotal_hPa-e_hPa;

assert(Pdry_hPa > 0, ...
    'Computed dry-air pressure must be positive.');


%% Original MATLAB reference gaseous losses
%
% Historical values retained exactly as originally generated.

reference_dB = [ ...
    0.464872051720006
    0.915619201385095
    0.915619201385093
    1.083269801941140
    0.579644683798001
    0.746529809423433
    0.823309107604283
    1.011147561592230
    1.014773052147900
    1.223282587959560
    1.307988823149580
    1.451125653190190
    0.736578957477030
    0.740967662299461
    0.879755276729633
    1.095074353889510
    2.077132816727250
    2.077451994238410
    4.049861773076530
    5.131094883426020
    0.924362770672161
    0.784988441756104
    0.682794565160070
    0.630299567781078
    0.651756984573601
    0.692143534373514 ];

assert(numel(reference_dB) == N, ...
    'Reference vector and candidate pool have different sizes.');


%% Allocate

matlabAuto_dB          = nan(N,1);
matlabExplicitTotal_dB = nan(N,1);  % diagnostic only: wrong R2023b convention
matlabExplicitDry_dB   = nan(N,1);  % documented R2023b convention
portableP67613_dB      = nan(N,1);
portableCompat12_dB    = nan(N,1);

h_km = config.receiver.h_m/1000;


%% Evaluate all links

for i = 1:N

    fGHz = candidatePool.Frequency_GHz(i);
    elDeg = candidatePool.Elevation_deg(i);

    %% MATLAB P.618 configuration

    cfg = p618Config;

    cfg.Frequency = fGHz*1e9;
    cfg.ElevationAngle = elDeg;
    cfg.Latitude = config.receiver.lat_deg;
    cfg.Longitude = config.receiver.lon_deg;
    cfg.GasAnnualExceedance = 50;


    %% --------------------------------------------------------
    % 1) MATLAB with internally obtained meteorology
    % ---------------------------------------------------------

    plAuto = p618PropagationLosses( ...
        cfg, ...
        'StationHeight',h_km);

    matlabAuto_dB(i) = plAuto.Ag;


    %% --------------------------------------------------------
    % 2) MATLAB explicit meteorology using TOTAL pressure
    %
    % Retained only to document the original diagnostic mismatch.
    % This is NOT the documented R2023b Pressure convention.
    % ---------------------------------------------------------

    plTotal = p618PropagationLosses( ...
        cfg, ...
        'StationHeight',h_km, ...
        'Pressure',Ptotal_hPa, ...
        'Temperature',T_K, ...
        'WaterVaporDensity',rho_gm3, ...
        'IntegratedWaterVaporContent',V_kgm2);

    matlabExplicitTotal_dB(i) = plTotal.Ag;


    %% --------------------------------------------------------
    % 3) MATLAB explicit P.2145 meteorology using DRY pressure
    %
    % This follows the R2023b documented Pressure convention.
    % ---------------------------------------------------------

    plDry = p618PropagationLosses( ...
        cfg, ...
        'StationHeight',h_km, ...
        'Pressure',Pdry_hPa, ...
        'Temperature',T_K, ...
        'WaterVaporDensity',rho_gm3, ...
        'IntegratedWaterVaporContent',V_kgm2);

    matlabExplicitDry_dB(i) = plDry.Ag;


    %% --------------------------------------------------------
    % 4) Public portable implementation: P.2145 + P.676-13
    % ---------------------------------------------------------

    portable13 = ...
        itu676_annex2_slant_gas_loss( ...
            fGHz, ...
            elDeg, ...
            h_km, ...
            met);

    portableP67613_dB(i) = portable13;


    %% --------------------------------------------------------
    % 5) Validation-only P.676-12 compatibility reconstruction
    %
    % The specific oxygen attenuation gammaO and the complete
    % water-vapour contribution remain exactly those of the
    % portable implementation. Only the oxygen equivalent
    % height h0 is changed from P.676-13 to P.676-12.
    % ---------------------------------------------------------

    [~,gammaO_dBkm,~] = ...
        itu676_specific_attenuation( ...
            fGHz, ...
            Ptotal_hPa, ...
            T_K, ...
            rho_gm3);

    h0_13_km = ...
        oxygenEquivalentHeightP676_13PaperBands( ...
            fGHz, ...
            Ptotal_hPa, ...
            T_K, ...
            rho_gm3);

    h0_12_km = ...
        oxygenEquivalentHeightP676_12( ...
            fGHz, ...
            Pdry_hPa, ...
            T_K, ...
            rho_gm3);

    Ao13_slant_dB = ...
        gammaO_dBkm*h0_13_km/sind(elDeg);

    % The remainder is the complete portable water-vapour term,
    % including the station-height correction when applicable.
    AwPortable_slant_dB = ...
        portable13-Ao13_slant_dB;

    Ao12_slant_dB = ...
        gammaO_dBkm*h0_12_km/sind(elDeg);

    portableCompat12_dB(i) = ...
        Ao12_slant_dB+AwPortable_slant_dB;
end


%% Differences

dAutoRef = ...
    matlabAuto_dB-reference_dB;

dTotalVsDry = ...
    matlabExplicitTotal_dB-matlabExplicitDry_dB;

dPortable13VsMAT12 = ...
    portableP67613_dB-matlabExplicitDry_dB;

dCompat12VsMAT12 = ...
    portableCompat12_dB-matlabExplicitDry_dB;


%% Comparison table

comparison = table( ...
    candidatePool.Label, ...
    candidatePool.Architecture, ...
    candidatePool.Frequency_GHz, ...
    candidatePool.Elevation_deg, ...
    reference_dB, ...
    matlabAuto_dB, ...
    matlabExplicitTotal_dB, ...
    matlabExplicitDry_dB, ...
    portableP67613_dB, ...
    portableCompat12_dB, ...
    dAutoRef, ...
    dTotalVsDry, ...
    dPortable13VsMAT12, ...
    dCompat12VsMAT12, ...
    'VariableNames',{ ...
    'Label', ...
    'Architecture', ...
    'Frequency_GHz', ...
    'Elevation_deg', ...
    'Reference_dB', ...
    'MATLAB_Auto_dB', ...
    'MATLAB_ExplicitTotalP_dB', ...
    'MATLAB_ExplicitDryP_dB', ...
    'Portable_P676_13_dB', ...
    'Portable_Compat_P676_12_dB', ...
    'Auto_minus_Ref_dB', ...
    'MATLAB_TotalP_minus_DryP_dB', ...
    'Portable13_minus_MATLAB12_dB', ...
    'Compat12_minus_MATLAB12_dB'});

disp(comparison);


%% Summary

fprintf('\n');
fprintf('============================================================\n');
fprintf('FINAL DIAGNOSTIC SUMMARY\n');
fprintf('============================================================\n');

fprintf('\nA) MATLAB auto vs original historical reference\n');
fprintf('Mean abs difference : %.12g dB\n', ...
    mean(abs(dAutoRef)));
fprintf('Max  abs difference : %.12g dB\n', ...
    max(abs(dAutoRef)));

fprintf('\nB) Effect of pressure convention inside R2023b\n');
fprintf('   MATLAB(total-P input) - MATLAB(dry-P input)\n');
fprintf('Mean abs difference : %.12g dB\n', ...
    mean(abs(dTotalVsDry)));
fprintf('Max  abs difference : %.12g dB\n', ...
    max(abs(dTotalVsDry)));

fprintf('\nC) Portable P.676-13 vs MATLAB R2023b P.676-12\n');
fprintf('   Both use the same explicit surface meteorology, with\n');
fprintf('   pressure converted to each implementation''s convention.\n');
fprintf('Mean abs difference : %.12g dB\n', ...
    mean(abs(dPortable13VsMAT12)));
fprintf('Max  abs difference : %.12g dB\n', ...
    max(abs(dPortable13VsMAT12)));
fprintf('Interpretation       : expected ITU-revision difference.\n');

fprintf('\nD) Compatibility validation: P.676-12 h0 reconstruction\n');
fprintf('Mean abs difference : %.12g dB\n', ...
    mean(abs(dCompat12VsMAT12)));
fprintf('Max  abs difference : %.12g dB\n', ...
    max(abs(dCompat12VsMAT12)));

fprintf('\nP.2145 meteorological inputs\n');
fprintf('Total pressure P     : %.12f hPa\n',Ptotal_hPa);
fprintf('Water-vapour e       : %.12f hPa\n',e_hPa);
fprintf('Dry-air pressure p   : %.12f hPa\n',Pdry_hPa);
fprintf('Temperature T        : %.12f K\n',T_K);
fprintf('Water-vapour rho     : %.12f g/m^3\n',rho_gm3);
fprintf('Integrated vapour V  : %.12f kg/m^2\n',V_kgm2);

compatTol_dB = 1e-6;
maxCompat_dB = max(abs(dCompat12VsMAT12));

fprintf('\nValidation tolerance : %.3e dB\n',compatTol_dB);

if isR2023b

    if maxCompat_dB <= compatTol_dB
        fprintf('Compatibility check  : PASS\n');
    else
        fprintf('Compatibility check  : FAIL\n');
        error([ ...
            'P.676-12 compatibility reconstruction does not ' ...
            'match MATLAB R2023b within tolerance.']);
    end

else
    fprintf(['Compatibility check  : NOT ENFORCED ' ...
             '(script calibrated to R2023b)\n']);
end

fprintf('\nConclusion\n');
fprintf(['The portable P.676-13 implementation must NOT be changed ' ...
         'to force equality with MATLAB R2023b.\n']);
fprintf(['The residual direct difference is explained by the ' ...
         'different P.676 revisions after the documented ' ...
         'pressure convention is respected.\n']);

fprintf('============================================================\n');


%% ============================================================
% Local function - P.676-13 oxygen equivalent height
%
% Specialized to the three carrier-frequency regions used by
% the paper: 19.95, 20.0 and 31.15 GHz.
%
% This reproduces the h0 branch used by
% itu676_annex2_slant_gas_loss.m.
% =============================================================

function h0_km = ...
    oxygenEquivalentHeightP676_13PaperBands( ...
        f_GHz,Ptotal_hPa,T_K,rho_gm3)

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
        error([ ...
            'Frequency %.3f GHz is outside the P.676-13 ' ...
            'coefficient ranges used by this paper.'], ...
            f_GHz);
    end

    a0 = interp1(fTab,a0Tab,f_GHz,'linear');
    b0 = interp1(fTab,b0Tab,f_GHz,'linear');
    c0 = interp1(fTab,c0Tab,f_GHz,'linear');
    d0 = interp1(fTab,d0Tab,f_GHz,'linear');

    h0_km = ...
        a0 ...
        + b0*T_K ...
        + c0*Ptotal_hPa ...
        + d0*rho_gm3;
end


%% ============================================================
% Local function - P.676-12 Annex 2 oxygen equivalent height
%
% Recommendation ITU-R P.676-12 (08/2019), Annex 2:
%   equations (30)-(35a).
%
% Input pressure p is DRY-AIR pressure [hPa].
% =============================================================

function h0_km = ...
    oxygenEquivalentHeightP676_12( ...
        f_GHz,Pdry_hPa,T_K,rho_gm3)

    e_hPa = rho_gm3*T_K/216.7;

    rp = ...
        (Pdry_hPa+e_hPa)/1013.25;

    A = ...
        0.7832 ...
        + 0.00709*(T_K-273.15);

    t1 = ...
        5.1040/(1+0.066*rp^(-2.3)) ...
        * exp( ...
        -( ...
        (f_GHz-59.7) ...
        /(2.87+12.4*exp(-7.9*rp)) ...
        )^2);

    ci = [ ...
        0.1597
        0.1066
        0.1325
        0.1242
        0.0938
        0.1448
        0.1374 ];

    fi = [ ...
        118.750334
        368.498246
        424.763020
        487.249273
        715.392902
        773.839490
        834.145546 ];

    t2 = 0;

    for m = 1:numel(ci)

        t2 = t2 + ...
            ci(m)*exp(2.12*rp) ...
            /( ...
            (f_GHz-fi(m))^2 ...
            + 0.025*exp(2.2*rp) ...
            );
    end

    t3 = ...
        (0.0114*f_GHz)/(1+0.14*rp^(-2.6)) ...
        * ( ...
        15.02*f_GHz^2 ...
        - 1353*f_GHz ...
        + 5.333e4 ...
        ) ...
        /( ...
        f_GHz^3 ...
        - 151.3*f_GHz^2 ...
        + 9629*f_GHz ...
        - 6803 ...
        );

    h0_km = ...
        (6.1*A)/(1+0.17*rp^(-1.1)) ...
        * (1+t1+t2+t3);

    if f_GHz < 70
        h0_km = ...
            min(h0_km,10.7*rp^0.3);
    end
end
