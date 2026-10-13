function [wafer,out] = fun_ShiftRotation(file1,file2,file3,file4,cfg)
%FUN_SHIFTROTATION Absolute wafer figure from original/X/Y/rotation XYZ maps.
% [wafer,out]=fun_ShiftRotation(f1,f2,f3,f4,cfg)
% cfg.ShiftPixels=[dx dy], cfg.RotationDeg=angle are REQUIRED, in the
% coordinate-argument convention of paper Eqs.(3)-(5), NOT stage motion.
% In Cartesian +Y-up coordinates: T2=(x+dx,y), T3=(x,y+dy),
% T4=(cos(a)*x-sin(a)*y,sin(a)*x+cos(a)*y).
% A physical positive stage displacement/CCW rotation normally needs the
% opposite sign here (inverse coordinate map). Keep measurements registered
% to the FIXED detector/reference; do not rotate/translate measured W maps.
%
% Defaults: GridSize=[875 875], PixelPitch=0.000362 [m/pixel], IndexBase=1,
% CenterPixel=[438 438] (INTERNAL one-based coordinates), RadiusPixels=437,
% YAxisSign=1 (pixel Y increases in Cartesian +Y), K=81, HeaderLines=0,
% IncludePiston=true, FitStride=1, BlockSize=10000, RankTolerance=1e-11.
% Optional cfg.Truth: same-size surface [m], or XYZ path, for error RMS.
% Optional cfg.ValidMasks: ny-by-nx-by-4 logical mask; for zero-filled voids.
% Optional cfg.ReturnSystem=true: stores sampled G,R for small diagnostics.
%
% Output wafer(y,x), height in METERS, NaN outside original valid aperture.
% Gauge: a1:a4=0 (piston, two tilts, power); those terms are NOT recovered.
% IncludePiston=false reproduces paper's 9 nuisance terms; true adds three
% independent differential offsets, robust to acquisition piston changes.
% No normal equations, explicit inverse, regularization or silent pinv.
if nargin<5, cfg=struct; end
if ~isfield(cfg,'ShiftPixels') || ~isfield(cfg,'RotationDeg')
    error('ShiftRotation:Geometry','Supply cfg.ShiftPixels and cfg.RotationDeg.');
end
cfg=defaults(cfg);
validateattributes(cfg.GridSize,{'numeric'},{'vector','numel',2,'integer','>=',3});
validateattributes(cfg.ShiftPixels,{'numeric'},{'vector','numel',2,'finite','real'});
validateattributes(cfg.RotationDeg,{'numeric'},{'scalar','finite','real'});
validateattributes(cfg.PixelPitch,{'numeric'},{'scalar','positive','finite'});
validateattributes(cfg.RadiusPixels,{'numeric'},{'scalar','positive','finite'});
validateattributes(cfg.CenterPixel,{'numeric'},{'vector','numel',2,'finite'});
validateattributes(cfg.K,{'numeric'},{'scalar','integer','>=',6,'<=',153});
validateattributes(cfg.FitStride,{'numeric'},{'scalar','integer','positive'});
validateattributes(cfg.BlockSize,{'numeric'},{'scalar','integer','positive'});
validateattributes(cfg.HeaderLines,{'numeric'},{'scalar','integer','nonnegative'});
validateattributes(cfg.RankTolerance,{'numeric'},{'scalar','positive','<',1});
assert(ismember(cfg.IndexBase,[0 1]),'IndexBase must be zero or one.');
assert(ismember(cfg.YAxisSign,[-1 1]),'YAxisSign must be +1 or -1.');
assert(isscalar(cfg.IncludePiston)&&ismember(cfg.IncludePiston,[0 1]),'Invalid IncludePiston.');
assert(isscalar(cfg.ReturnSystem)&&ismember(cfg.ReturnSystem,[0 1]),'Invalid ReturnSystem.');
files={file1,file2,file3,file4}; W=cell(1,4); readInfo=cell(1,4);
for j=1:4, [W{j},readInfo{j}]=sr_read_xyz(files{j},cfg); end
[px,py]=meshgrid(1:cfg.GridSize(2),1:cfg.GridSize(1));
x=(px-cfg.CenterPixel(1))/cfg.RadiusPixels;
y=cfg.YAxisSign*(py-cfg.CenterPixel(2))/cfg.RadiusPixels;
disk=hypot(x,y)<=1+1e-12;
if isfield(cfg,'ValidMasks') && ~isempty(cfg.ValidMasks)
    assert(islogical(cfg.ValidMasks)&&isequal(size(cfg.ValidMasks),[cfg.GridSize 4]), ...
        'ValidMasks must be a logical ny-by-nx-by-4 array.');
    for j=1:4, W{j}(~cfg.ValidMasks(:,:,j))=NaN; end
end
base=disk & isfinite(W{1});
sample=mod(px-1,cfg.FitStride)==0 & mod(py-1,cfg.FitStride)==0;
active=5:cfg.K; ka=numel(active); nuisance=2:4;
if cfg.IncludePiston, nuisance=1:4; end
q=numel(nuisance); p=ka+3*q;
masks=cell(1,3); ids=cell(1,3); count=zeros(1,3);
for j=1:3
    [u,v]=transform(x,y,j,cfg);
    masks{j}=base & isfinite(W{j+1}) & hypot(u,v)<=1+1e-12;
    ids{j}=find(masks{j}&sample); count(j)=numel(ids{j});
end
if any(count==0)||sum(count)<p
    error('ShiftRotation:Data','Insufficient overlap/sample points: [%d %d %d].',count);
end
if cfg.ReturnSystem && sum(count)*p>2e7
    error('ShiftRotation:Memory','ReturnSystem is only for small diagnostics (<20M entries).');
end
T=zeros(0,p+1); Gsaved=cell(0,1); Rsaved=cell(0,1);
% Streaming QR of augmented [G R]. Orthogonal transformations preserve
% ||G*c-R||, including residual energy in the final augmented column.
for j=1:3
    for first=1:cfg.BlockSize:numel(ids{j})
        ii=ids{j}(first:min(first+cfg.BlockSize-1,end));
        A=design(x(ii),y(ii),j,cfg,active,nuisance);
        b=W{j+1}(ii)-W{1}(ii);
        T=qr([T; A b],0);
        % One-output qr(X,0) is economy R, no huge Q allocated.
        if cfg.ReturnSystem
            Gsaved{end+1,1}=A; Rsaved{end+1,1}=b; %#ok<AGROW>
        end
    end
end
C=T(:,1:p); d=T(:,p+1); scale=sqrt(sum(C.^2,1));
if any(scale==0), error('ShiftRotation:Rank','Zero design column: degenerate geometry.'); end
[U,S,V]=svd(C./scale,'econ'); s=diag(S);
numericalRank=sum(s>cfg.RankTolerance*s(1));
if numericalRank<p
    error('ShiftRotation:Rank', ...
        'Rank %d < %d after gauge removal. Change shifts/angle/mask/K.',numericalRank,p);
end
condition=s(1)/s(end);
if condition>1e8
    warning('ShiftRotation:Condition','Scaled condition number is %.3g; results may amplify noise.',condition);
end
solution=(V*((U'*d)./s))./scale';
a=zeros(cfg.K,1); a(active)=solution(1:ka);
wafer=nan(cfg.GridSize); ii=find(base);
for first=1:cfg.BlockSize:numel(ii)
    jj=ii(first:min(first+cfg.BlockSize-1,end));
    wafer(jj)=sr_zernike(x(jj),y(jj),cfg.K)*a;
end
residual=cell(1,3); residualRMS=zeros(1,3); sumsq=0; nall=0;
% Validate on ALL available pixels even when a fitting stride was requested.
for j=1:3
    ii=find(masks{j}); e=nan(cfg.GridSize); ss=0;
    for first=1:cfg.BlockSize:numel(ii)
        jj=ii(first:min(first+cfg.BlockSize-1,end));
        A=design(x(jj),y(jj),j,cfg,active,nuisance);
        v=A*solution-(W{j+1}(jj)-W{1}(jj));
        e(jj)=v; ss=ss+sum(v.^2);
    end
    residual{j}=e; residualRMS(j)=sqrt(ss/numel(ii));
    sumsq=sumsq+ss; nall=nall+numel(ii);
end
[~,modes]=sr_zernike(0,0,cfg.K);
out=struct;
out.coefficients=a; out.modeTable=modes; % cols: n, |m|, 0=radial/1=cos/2=sin
out.solution=solution; out.nuisance=reshape(solution(ka+1:end),q,3);
out.nuisanceModeIndices=nuisance; out.activeModeIndices=active;
out.cfg=cfg; out.readInfo=readInfo;
out.x_m=(px-cfg.CenterPixel(1))*cfg.PixelPitch;
out.y_m=cfg.YAxisSign*(py-cfg.CenterPixel(2))*cfg.PixelPitch;
out.apertureMask=base; out.differenceMasks=masks;
out.fitRowCounts=count; out.systemSize=[sum(count) p];
out.rank=numericalRank; out.scaledConditionNumber=condition;
out.scaledSingularValues=s;
out.fitResidualRMS=sqrt(sum((C*solution-d).^2)/sum(count));
out.residualMaps=residual; out.residualRMS=residualRMS;
out.allResidualRMS=sqrt(sumsq/nall);
v=wafer(base); out.surfaceRMS=sqrt(mean((v-mean(v)).^2));
out.surfacePV=max(v)-min(v);
out.gauge='a(1:4)=0: piston, x/y tilt, power NOT calibrated';
out.referenceWithOriginalAlignment=W{1}-wafer;
% This last map contains unknown original pose errors; not pure reference.
out.errorRMSRaw=NaN; out.errorRMSAligned=NaN; out.truthCount=0;
if isfield(cfg,'Truth') && ~isempty(cfg.Truth)
    truth=cfg.Truth;
    if ischar(truth)||isstring(truth), truth=sr_read_xyz(truth,cfg); end
    assert(isequal(size(truth),cfg.GridSize),'Truth must match GridSize.');
    good=base & isfinite(truth); jj=find(good);
    if numel(jj)<4, error('ShiftRotation:Truth','Insufficient truth overlap.'); end
    E=wafer(jj)-truth(jj); B=sr_zernike(x(jj),y(jj),4);
    low=B\E; aligned=E-B*low;
    out.errorRMSRaw=sqrt(mean(E.^2));
    out.errorRMSAligned=sqrt(mean(aligned.^2));
    out.truthAlignmentCoefficients=low; out.truthCount=numel(jj);
    out.errorMapAligned=nan(cfg.GridSize); out.errorMapAligned(jj)=aligned;
end
if cfg.ReturnSystem, out.G=vertcat(Gsaved{:}); out.R=vertcat(Rsaved{:}); end
end

function cfg=defaults(cfg)
d=struct('GridSize',[875 875],'PixelPitch',0.000362,'IndexBase',1, ...
    'YAxisSign',1,'K',81,'HeaderLines',0,'IncludePiston',true, ...
    'FitStride',1,'BlockSize',10000,'RankTolerance',1e-11,'ReturnSystem',false);
f=fieldnames(d);
for j=1:numel(f), if ~isfield(cfg,f{j}), cfg.(f{j})=d.(f{j}); end, end
if ~isfield(cfg,'CenterPixel'), cfg.CenterPixel=[(cfg.GridSize(2)+1)/2 (cfg.GridSize(1)+1)/2]; end
if ~isfield(cfg,'RadiusPixels'), cfg.RadiusPixels=(min(cfg.GridSize)-1)/2; end
end

function [u,v]=transform(x,y,j,c)
if j==1
    u=x+c.ShiftPixels(1)/c.RadiusPixels; v=y;
elseif j==2
    u=x; v=y+c.ShiftPixels(2)/c.RadiusPixels;
else
    ca=cosd(c.RotationDeg); sa=sind(c.RotationDeg);
    u=ca*x-sa*y; v=sa*x+ca*y;
end
end

function A=design(x,y,j,c,active,nuisance)
[u,v]=transform(x,y,j,c); Z=sr_zernike(x,y,c.K);
D=sr_zernike(u,v,c.K)-Z; q=numel(nuisance); ka=numel(active);
A=zeros(numel(x),ka+3*q); A(:,1:ka)=D(:,active);
A(:,ka+(j-1)*q+(1:q))=Z(:,nuisance);
end
