function [files,cfg,truth,a] = sr_make_simulation(folder,N,K,noiseStd,includePiston)
%SR_MAKE_SIMULATION Synthetic registered detector maps, NEVER experiment data.
% N=875 default use in demo; arbitrary reference includes non-Zernike content.
% K=15 uses independently written Cartesian polynomials as a forward oracle.
% Larger K adds modes of the same basis for an exact-model recovery test.
if ~isfolder(folder), mkdir(folder); end
cfg=struct('GridSize',[N N],'PixelPitch',0.000362,'IndexBase',1, ...
    'ShiftPixels',[-0.03 -0.03]/0.000362,'RotationDeg',1, ...
    'RadiusPixels',(N-1)/2,'CenterPixel',[(N+1)/2 (N+1)/2], ...
    'YAxisSign',1,'K',K,'IncludePiston',includePiston,'HeaderLines',0);
% For small unit tests use the SAME normalized shifts as the 875-grid example.
cfg.ShiftPixels=cfg.ShiftPixels*(N-1)/874;
[px,py]=meshgrid(1:N,1:N); x=(px-(N+1)/2)/cfg.RadiusPixels;
y=(py-(N+1)/2)/cfg.RadiusPixels; r2=x.^2+y.^2; disk=r2<=1;
a=zeros(K,1); a(5)=2.1e-9/sqrt(6); a(6)=-1.4e-9/sqrt(6);
a(7)=.8e-9/sqrt(8); a(8)=.3e-9/sqrt(8); a(9)=.35e-9/sqrt(8);
a(11)=.45e-9/sqrt(5);
if K>15, a(16:K)=0.15e-9*sin((16:K)') ./ sqrt((16:K)'-14); end
reference=4e-9*(sin(2.7*x).*cos(3.1*y)+0.2*x.*y)+1e-9*r2;
truth=forward(x,y,a); truth(~disk)=NaN;
pose=1e-9*[.3 -.7 .8 -.2; 1 -.6 .4 .2; -.5 .2 -.7 .6; .6 -.2 .1 -.4];
if ~includePiston, pose(1,:)=0; end
files=cell(1,4); names={'original','xshift','yshift','rotation'};
for j=1:4
    u=x; v=y;
    if j==2, u=x+cfg.ShiftPixels(1)/cfg.RadiusPixels; end
    if j==3, v=y+cfg.ShiftPixels(2)/cfg.RadiusPixels; end
    if j==4
        u=cosd(cfg.RotationDeg)*x-sind(cfg.RotationDeg)*y;
        v=sind(cfg.RotationDeg)*x+cosd(cfg.RotationDeg)*y;
    end
    low=pose(1,j)+2*pose(2,j)*x+2*pose(3,j)*y+sqrt(3)*pose(4,j)*(2*r2-1);
    measurement=reference+forward(u,v,a)+low+noiseStd*randn(N);
    valid=disk & hypot(u,v)<=1+1e-12; measurement(~valid)=NaN;
    files{j}=fullfile(folder,['sim_' names{j} '.xyz']);
    write_xyz(files{j},[px(:) py(:) measurement(:)]);
end
cfg.Truth=truth;
save(fullfile(folder,'simulation_truth.mat'),'truth','a','cfg','pose','noiseStd');
end

function z=forward(x,y,a)
r2=x.^2+y.^2;
z=1e-9*(2.1*(x.^2-y.^2)-1.4*(2*x.*y)+ ...
    .8*(3*r2-2).*x+.3*(3*r2-2).*y+ ...
    .35*(x.^3-3*x.*y.^2)+.45*(6*r2.^2-6*r2+1));
if numel(a)>15
    % Chunk the synthesis so the 875x875x81 basis is never retained.
    for first=1:10000:numel(x)
        ii=first:min(first+9999,numel(x)); Z=sr_zernike(x(ii),y(ii),numel(a));
        extra=Z(:,16:end)*a(16:end); z(ii)=z(ii)+reshape(extra,size(z(ii)));
    end
end
end

function write_xyz(file,A)
fid=fopen(file,'w'); assert(fid~=-1,'Cannot write XYZ.');
c=onCleanup(@()fclose(fid));
fprintf(fid,'%d\t%d\t%.17g\n',A');
end
