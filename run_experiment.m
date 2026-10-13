% RUN_EXPERIMENT User geometry: +X/+Y 0.03m stage shifts, 1deg clockwise.
% Run in this file's folder. Files must share fixed detector coordinates.
clear; here=fileparts(mfilename('fullpath')); cd(here);
files={'wafer_original.xyz','wafer_xshift.xyz','wafer_yshift.xyz','wafer_rotation.xyz'};
cfg.GridSize=[875 875];
cfg.PixelPitch=0.000362; % confirmed m/pixel
cfg.K=81; cfg.IndexBase=1; cfg.HeaderLines=0;
cfg.CenterPixel=[438 438]; cfg.RadiusPixels=437; % confirm physical aperture
cfg.YAxisSign=1; % use -1 if pixel y points down and physical +Y points up
% Cartesian +Y-up viewed from measurement side; clockwise stage=-1 deg,
% inverse coordinate rotation=+1 deg. If viewing side differs, update sign.
shiftXArgument_m=-0.03;
shiftYArgument_m=-0.03;
rotationArgument_deg=1;
assert(all(isfinite([shiftXArgument_m shiftYArgument_m rotationArgument_deg])), ...
    'Fill confirmed signed geometry and check PixelPitch before experimental use.');
cfg.ShiftPixels=[shiftXArgument_m shiftYArgument_m]/cfg.PixelPitch;
cfg.RotationDeg=rotationArgument_deg;
cfg.IncludePiston=true; cfg.FitStride=1; cfg.BlockSize=10000;
% cfg.Truth='independent_reference_wafer.xyz'; % optional for true error RMS
[wafer,out]=fun_ShiftRotation(files{:},cfg);
save('wafer_absolute_result.mat','wafer','out','-v7.3');
[px,py]=meshgrid(1:875,1:875); fid=fopen('wafer_absolute.xyz','w');
assert(fid~=-1,'Cannot create output file.');
fprintf(fid,'%d\t%d\t%.17g\n',[px(:) py(:) wafer(:)]'); fclose(fid);
fprintf('Difference-equation residual RMS = %.6g nm (NOT absolute accuracy)\n',out.allResidualRMS*1e9);
fprintf('Truth-aligned error RMS = %.6g nm (NaN without truth)\n',out.errorRMSAligned*1e9);
