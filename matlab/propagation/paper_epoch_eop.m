function eop = paper_epoch_eop()
%PAPER_EPOCH_EOP Earth orientation parameters for the paper epoch.
%
% Epoch: 01-Aug-2026 12:00:00 UTC
% Values are linearly interpolated between the surrounding IERS
% daily values at 00:00 UTC.

eop.xp_arcsec = 0.2219685;
eop.yp_arcsec = 0.3648005;
eop.dut1_s    = 0.0124859;
eop.lod_s     = 0.00042555;
eop.tai_utc_s = 37.0;
eop.eqeterms  = 2;

arcsec2rad = pi/(180*3600);
eop.xp_rad = eop.xp_arcsec * arcsec2rad;
eop.yp_rad = eop.yp_arcsec * arcsec2rad;

end
