function cat_vol_defs(job)
% Apply deformations to images. In contrast to spm_defs images are saved
% in the original directory.
%_______________________________________________________________________
% Christian Gaser
% $Id: cat_vol_defs.m 1020 2016-09-26 12:46:28Z dahnke $
%
% based on spm_deformations.m

many_images = 0;

try
  PU = job.field;
catch
  PU = job.field1;
  many_images = 1;
end

def.verb = cat_get_defaults('extopts.verb'); 
job      = cat_io_checkinopt(job,def);

PI   = job.images;
intp = job.interp;

for i=1:numel(PU),

  [pth,nam,ext] = spm_fileparts(PU{i});
  PU{i} = fullfile(pth,[nam ext]);

  [Def,mat] = get_def(PU{i});
 
  for m=1:numel(PI)
    
    if many_images % many images
      PIi = char(PI{m}); 
    else % many subjects
      PIi = char(PI{m}{i}); 
    end
    PIri = apply_def(Def,mat,PIi,intp,job.modulate);
    
    if job.verb
      fprintf('Display resampled %s\n',spm_file(PIri,'link','spm_image(''Display'',''%s'')'));
    end
  end
end
return

%_______________________________________________________________________
function [Def,mat] = get_def(job)
% Load a deformation field saved as an image
Nii = nifti(job);
Def = single(Nii.dat(:,:,:,1,:));
d   = size(Def);
if d(4)~=1 || d(5)~=3, error('Deformation field is wrong!'); end
Def = reshape(Def,[d(1:3) d(5)]);
mat = Nii.mat;

%_______________________________________________________________________
function fname = apply_def(Def,mat,fnames,intrp,modulate)
% Warp an image or series of images according to a deformation field

intrp = [intrp*[1 1 1], 0 0 0];

for i=1:size(fnames,1),

    % Generate headers etc for output images
    %----------------------------------------------------------------------
    [pth,nam,ext,num] = spm_fileparts(deblank(fnames(i,:))); ext = '.nii';  %#ok<ASGLU>
    NI = nifti(fullfile(pth,[nam ext]));
    j_range = 1:size(NI.dat,4);
    k_range = 1:size(NI.dat,5);
    l_range = 1:size(NI.dat,6);
    if ~isempty(num),
        num = sscanf(num,',%d');
        if numel(num)>=1, j_range = num(1); end
        if numel(num)>=2, k_range = num(2); end
        if numel(num)>=3, l_range = num(3); end
    end

    NO = NI;
    ext = '.nii'; 

    % use float for modulated images
    if modulate
        NO.dat.scl_slope = 1.0;
        NO.dat.scl_inter = 0.0;
        NO.dat.dtype     = 'float32-le';
    end

    dim            = size(Def);
    dim            = dim(1:3);
    NO.dat.dim     = [dim NI.dat.dim(4:end)];
    NO.dat.offset  = 0; % For situations where input .nii images have an extension.
    NO.mat         = mat;
    NO.mat0        = mat;
    NO.mat_intent  = 'Aligned';
    NO.mat0_intent = 'Aligned';

    switch modulate
    case 0
        NO.dat.fname = fullfile(pth,['w',nam,ext]);
        NO.descrip   = sprintf('Warped');
    case 1
        NO.dat.fname = fullfile(pth,['mw',nam,ext]);
        NO.descrip   = sprintf('Warped & Jac scaled');
    case 2
        NO.dat.fname = fullfile(pth,['m0w',nam,ext]);
        NO.descrip   = sprintf('Warped & Jac scaled (nonlinear only)');
    end
    fname = NO.dat.fname; 
    
    NO.extras      = [];
    create(NO);

    if modulate
      dt = spm_diffeo('def2det',Def)/det(mat(1:3,1:3));
      dt(:,:,[1 end]) = NaN;
      dt(:,[1 end],:) = NaN;
      dt([1 end],:,:) = NaN;

      % for modulation of non-linear parts only we have to remove the affine part
      % of the jacobian determinant
      if modulate == 2
        [x1,x2,x3] = ndgrid(single(1:size(Def,1)),single(1:size(Def,2)),single(1:size(Def,3)));
        X = cat(4,x1,x2,x3);
        Ma = spm_get_closest_affine(X,Def);
        M3 = Ma\mat;
        dt = dt*abs(det(M3));
      end
    
    end
    
    for j=j_range

        M0 = NI.mat;
        if ~isempty(NI.extras) && isstruct(NI.extras) && isfield(NI.extras,'mat')
            M1 = NI.extras.mat;
            if size(M1,3) >= j && sum(sum(M1(:,:,j).^2)) ~=0
                M0 = M1(:,:,j);
            end
        end
        M  = inv(M0);
        % Generate new deformation (if needed)
        Y     = affine(Def,M);
        % Write the warped data for this time point
        %------------------------------------------------------------------
        for k=k_range
            for l=l_range
                C   = spm_diffeo('bsplinc',single(NI.dat(:,:,:,j,k,l)),intrp);
                dat = spm_diffeo('bsplins',C,Y,intrp);
                if modulate
                  dat = dat.*double(dt);
                end
                NO.dat(:,:,:,j,k,l) = dat;
            end
        end
    end

end;
return;

%==========================================================================
% function Def = affine(y,M)
%==========================================================================
function Def = affine(y,M)
Def          = zeros(size(y),'single');
Def(:,:,:,1) = y(:,:,:,1)*M(1,1) + y(:,:,:,2)*M(1,2) + y(:,:,:,3)*M(1,3) + M(1,4);
Def(:,:,:,2) = y(:,:,:,1)*M(2,1) + y(:,:,:,2)*M(2,2) + y(:,:,:,3)*M(2,3) + M(2,4);
Def(:,:,:,3) = y(:,:,:,1)*M(3,1) + y(:,:,:,2)*M(3,2) + y(:,:,:,3)*M(3,3) + M(3,4);
