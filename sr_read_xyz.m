function [W,info] = sr_read_xyz(file,cfg)
%SR_READ_XYZ Read unsorted [x_pixel,y_pixel,z_m] into W(y,x).
% Invalid Z is NaN/Inf. Missing pixels are allowed. Zero is VALID height.
% cfg: GridSize=[ny nx], IndexBase=0 or 1, HeaderLines>=0.
if ~isfile(file), error('ShiftRotation:File','Missing file: %s',file); end
A=readmatrix(file,'FileType','text','NumHeaderLines',cfg.HeaderLines);
if size(A,2)~=3 || isempty(A)
    error('ShiftRotation:Format','Expected exactly three numeric columns: %s',file);
end
x=A(:,1); y=A(:,2); z=A(:,3); ny=cfg.GridSize(1); nx=cfg.GridSize(2);
if any(~isfinite(x)|~isfinite(y)|x~=round(x)|y~=round(y))
    error('ShiftRotation:Index','X/Y must be finite integer pixel indices: %s',file);
end
x=x-cfg.IndexBase+1; y=y-cfg.IndexBase+1;
if any(x<1|x>nx|y<1|y>ny)
    error('ShiftRotation:Index','Indices outside grid; check IndexBase/GridSize: %s',file);
end
idx=sub2ind([ny nx],y,x);
if numel(unique(idx))~=numel(idx)
    error('ShiftRotation:Duplicate','Duplicate (x,y) coordinates: %s',file);
end
z(~isfinite(z))=NaN; W=nan(ny,nx); W(idx)=z;
info=struct('file',char(file),'rows',numel(z),'finiteHeights',sum(isfinite(z)));
end
