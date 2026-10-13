function results = test_ShiftRotation
%TEST_SHIFTROTATION Independent polynomial oracle + rank/IO/full-basis tests.
here=fileparts(mfilename('fullpath')); root=fullfile(here,'test_output');
if ~isfolder(root), mkdir(root); end
rng(20261001,'twister');
[f,c,t,a]=sr_make_simulation(fullfile(root,'polynomial'),65,15,0,true);
c.ReturnSystem=true; c.BlockSize=173;
[w,o]=fun_ShiftRotation(f{:},c);
assert(o.errorRMSAligned<1e-19,'Independent polynomial recovery failed.');
assert(norm(o.coefficients-a)<1e-18,'Coefficient recovery failed.');
direct=o.G\o.R;
assert(norm(direct-o.solution)<1e-18,'Streaming QR differs from direct solve.');
assert(isequal(size(w),[65 65])&&all(isnan(w(~o.apertureMask))),'Output mask failed.');
results.polynomialErrorRMS=o.errorRMSAligned;
results.streamingVsDirect=norm(direct-o.solution);
% Zero-based, unsorted rows, one-line header, different NaN/void masks.
for j=1:4
    A=readmatrix(f{j},'FileType','text'); A(:,1:2)=A(:,1:2)-1;
    A=A(randperm(size(A,1)),:);
    fid=fopen(f{j},'w'); fprintf(fid,'x_pixel y_pixel z_m\n');
    fprintf(fid,'%d\t%d\t%.17g\n',A'); fclose(fid);
end
c.IndexBase=0; c.HeaderLines=1; c.FitStride=2;
c.ValidMasks=true(65,65,4); c.ValidMasks(23:27,25:31,2)=false;
[~,o]=fun_ShiftRotation(f{:},c);
assert(o.errorRMSAligned<1e-18,'Header/index/mask/stride test failed.');
results.ioMaskErrorRMS=o.errorRMSAligned;
% Exactly paper's 9 nuisance terms (without differential piston).
[f,c]=sr_make_simulation(fullfile(root,'paper9'),65,15,0,false);
[~,o]=fun_ShiftRotation(f{:},c);
assert(o.systemSize(2)==15-4+9 && o.errorRMSAligned<1e-18,'Paper mode failed.');
results.paper9ErrorRMS=o.errorRMSAligned;
% Additional low-order truth is unobservable: aligned error remains zero.
[xx,yy]=meshgrid(linspace(-1,1,65));
c.Truth=t+7e-9+3e-9*xx-2e-9*yy+5e-9*(xx.^2+yy.^2);
[~,o]=fun_ShiftRotation(f{:},c);
assert(o.errorRMSRaw>1e-9 && o.errorRMSAligned<1e-18,'Gauge alignment failed.');
results.gaugeAlignedErrorRMS=o.errorRMSAligned;
% Zero or 180-degree rotation cannot identify astigmatism with free tilts.
for angle=[0 180]
    bad=c; bad.RotationDeg=angle; caught=false;
    try
        fun_ShiftRotation(f{:},bad);
    catch ME, caught=strcmp(ME.identifier,'ShiftRotation:Rank'); end
    assert(caught,'Degenerate angle must be rejected.');
end
% Full 81-mode exact/noisy reconstruction at manageable grid size.
[f,c,~,a]=sr_make_simulation(fullfile(root,'basis81'),101,81,0,true);
[~,o]=fun_ShiftRotation(f{:},c);
assert(o.errorRMSAligned<1e-17 && norm(o.coefficients-a)<1e-17,'81-mode recovery failed.');
results.basis81ErrorRMS=o.errorRMSAligned;
results.basis81Condition=o.scaledConditionNumber;
[f,c]=sr_make_simulation(fullfile(root,'noise81'),101,81,0.1e-9,true);
[~,o]=fun_ShiftRotation(f{:},c);
assert(isfinite(o.errorRMSAligned) && o.errorRMSAligned<0.2e-9,'Noise recovery unexpected.');
results.noise81ErrorRMS=o.errorRMSAligned;
results.noise81ResidualRMS=o.allResidualRMS;
save(fullfile(root,'test_results.mat'),'results');
fid=fopen(fullfile(root,'test_results.json'),'w');
fprintf(fid,'%s',jsonencode(results,'PrettyPrint',true)); fclose(fid);
disp(results); fprintf('All ShiftRotation tests passed.\n');
end
