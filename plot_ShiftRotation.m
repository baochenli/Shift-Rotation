function figFile = plot_ShiftRotation(folder)
%PLOT_SHIFTROTATION Save editable full-resolution surface maps in MATLAB FIG.
% Uses existing simulation XYZ and reconstruction.mat; does not rerun solver.
% All plots: X/Y in mm, height in nm, missing pixels transparent.
if nargin<1
    folder=fullfile(fileparts(mfilename('fullpath')),'simulation');
end
data=load(fullfile(folder,'reconstruction.mat'),'wafer','truth','out');
cfg=data.out.cfg;
names={'sim_original.xyz','sim_xshift.xyz','sim_yshift.xyz','sim_rotation.xyz'};
maps=cell(1,6);
for k=1:4, maps{k}=sr_read_xyz(fullfile(folder,names{k}),cfg); end
maps{5}=data.truth; maps{6}=data.wafer;
labels={'原位置测量面形 W1','沿 +x 平移 30 mm 后测量 W2', ...
    '沿 +y 平移 30 mm 后测量 W3','顺时针旋转 1° 后测量 W4', ...
    'Wafer 原始真值（模拟）','计算的 Wafer 绝对面形'};
x=data.out.x_m(1,:)*1e3; y=data.out.y_m(:,1)*1e3;
limits=zeros(2,2);
for group=1:2
    limits(group,:)=[Inf -Inf];
    if group==1, indices=1:4; else, indices=5:6; end
    for k=indices
        z=maps{k}(isfinite(maps{k}))*1e9;
        assert(~isempty(z),'Surface has no valid pixels.');
        limits(group,1)=min(limits(group,1),min(z));
        limits(group,2)=max(limits(group,2),max(z));
    end
    if limits(group,1)==limits(group,2), limits(group,:)=limits(group,:)+[-1 1]; end
end
fig=figure('Visible','off','Color','w','Position',[60 60 1550 930], ...
    'Name','Shift-Rotation 面形对比（模拟数据）','NumberTitle','off');
cleanup=onCleanup(@()close(fig));
layout=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
title(layout,{'Shift-Rotation 面形对比（模拟数据）', ...
    '四次测量含参考面及姿态项；绝对面形采用去活塞、倾斜和 power 的基准'}, ...
    'FontName','Microsoft YaHei','FontSize',16,'FontWeight','bold');
for k=1:6
    ax=nexttile(layout); z=maps{k}*1e9;
    imagesc(ax,x,y,z,'AlphaData',isfinite(z),'Tag',sprintf('surface_%d',k));
    axis(ax,'image'); set(ax,'YDir','normal','Color',[.93 .93 .93], ...
        'FontName','Microsoft YaHei','FontSize',11);
    xlabel(ax,'x (mm)'); ylabel(ax,'y (mm)');
    if k<=4, group=1; else, group=2; end
    caxis(ax,limits(group,:)); colormap(ax,jet(256));
    cb=colorbar(ax); cb.Label.String='高度 (nm)';
    cb.FontName='Microsoft YaHei'; cb.Label.FontName='Microsoft YaHei';
    v=z(isfinite(z)); surfaceRMS=sqrt(mean((v-mean(v)).^2));
    title(ax,{labels{k},sprintf('面形 RMS = %.4f nm',surfaceRMS)}, ...
        'FontName','Microsoft YaHei','FontSize',12,'FontWeight','normal');
end
% Preserve the numerical arrays and conventions in editable FIG metadata.
fig.UserData=struct('simulationOnly',true,'coordinateUnit','mm','heightUnit','nm', ...
    'measurementColorLimits_nm',limits(1,:),'waferColorLimits_nm',limits(2,:), ...
    'note','W1-W4 share one color scale; truth and recovered wafer share another.');
outdir=fullfile(folder,'figures'); if ~isfolder(outdir), mkdir(outdir); end
figFile=fullfile(outdir,'ShiftRotation_surface_maps.fig');
drawnow;
savefig(fig,figFile,'compact');
exportgraphics(fig,fullfile(outdir,'ShiftRotation_surface_maps.png'),'Resolution',160);
% Reopen the saved file and verify all six full-resolution data objects.
check=openfig(figFile,'new','invisible');
checkCleanup=onCleanup(@()close(check));
for k=1:6
    obj=findobj(check,'Tag',sprintf('surface_%d',k));
    assert(isscalar(obj)&&isequaln(obj.CData,maps{k}*1e9),'Saved FIG data mismatch.');
end
fprintf('Saved and verified: %s\n',figFile);
end
