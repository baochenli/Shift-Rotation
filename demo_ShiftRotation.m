function result = demo_ShiftRotation(N)
%DEMO_SHIFTROTATION Generate four synthetic XYZ files and solve G*x=R.
% demo_ShiftRotation defaults to the requested full 875 x 875, 81 modes.
% Geometry matches user stage movements under a Cartesian +Y-up convention.
if nargin<1, N=875; end
here=fileparts(mfilename('fullpath')); folder=fullfile(here,'simulation');
rng(20261001,'twister'); tic;
[files,cfg,truth,a]=sr_make_simulation(folder,N,81,0.1e-9,true);
[wafer,out]=fun_ShiftRotation(files{:},cfg);
elapsed=toc;
save(fullfile(folder,'reconstruction.mat'),'wafer','out','truth','a','-v7.3');
[x,y]=meshgrid(1:N); fid=fopen(fullfile(folder,'wafer_reconstructed.xyz'),'w');
assert(fid~=-1,'Cannot write reconstruction.');
fprintf(fid,'%d\t%d\t%.17g\n',[x(:) y(:) wafer(:)]'); fclose(fid);
result=struct('simulationOnly',true,'gridSize',[N N],'K',81, ...
    'shiftPixels',cfg.ShiftPixels,'rotationDeg',cfg.RotationDeg, ...
    'noiseStd_nm',0.1,'errorRMSAligned_nm',out.errorRMSAligned*1e9, ...
    'errorRMSRaw_nm',out.errorRMSRaw*1e9,'residualRMS_nm',out.allResidualRMS*1e9, ...
    'surfaceRMS_nm',out.surfaceRMS*1e9,'rank',out.rank,'systemSize',out.systemSize, ...
    'scaledConditionNumber',out.scaledConditionNumber,'seconds',elapsed);
fid=fopen(fullfile(folder,'validation_summary.json'),'w');
fprintf(fid,'%s',jsonencode(result,'PrettyPrint',true)); fclose(fid);
fig=figure('Visible','off','Position',[100 100 1200 360]);
maps={truth,wafer,out.errorMapAligned};
titles={'Synthetic truth (nm)','Recovered wafer (nm)','Aligned error (nm)'};
for j=1:3
    subplot(1,3,j); imagesc(maps{j}*1e9,'AlphaData',isfinite(maps{j}));
    axis image; set(gca,'YDir','normal'); colorbar; title(titles{j});
    xlabel('x pixel'); ylabel('y pixel');
end
exportgraphics(fig,fullfile(folder,'validation.png'),'Resolution',180); close(fig);
disp(result);
end
