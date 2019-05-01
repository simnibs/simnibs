function out = cat_vol_groupwise_ls(Nii, output, prec, b_settings, ord)
% Groupwise registration via least squares
% FORMAT out = spm_groupwise_ls(Nii, output, prec, b_settings, ord)
% Nii    - a nifti object for two or more image volumes.
% output - a cell array of output options (as scharacter strings).
%          'wimg   - write realigned images to disk
%          'avg'   - return average in out.avg
%          'wavg'  - write average to disk, and return filename in out.avg
%
% prec       - reciprocal of noise variance on images.
% b_settings - regularisation settings for nonuniformity field.
% ord        - degree of B-spline interpolation used for sampling images.
%
%_______________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging
%
% modified version (added masked registration) of
% John Ashburner
% spm_groupwise_ls.m 6008 2014-05-22 12:08:01Z john
%
% $Id cat_vol_groupwise_ls.m $

% Use brainmask to obtain better registration
%-----------------------------------------------------------------------
use_brainmask = 1;

% Get handles to NIfTI data
%-----------------------------------------------------------------------
if ~isa(Nii,'nifti')
    if isa(Nii,'char')
        Nii = nifti(Nii);
    else
        error('Unrecognised NIfTI images');
    end
end

% Specify default settings
%-----------------------------------------------------------------------
if nargin<3, prec       = NaN; end
if nargin<4, b_settings = [0 0 1e6]; end
if nargin<5, ord        = [3 3 3 0 0 0]; end

% If settings are not subject-specific, then generate
%-----------------------------------------------------------------------
if size(b_settings,1)==1, b_settings = repmat(b_settings,numel(Nii),1); end
if numel(prec)       ==1, prec       = repmat(prec,1,numel(Nii));       end

% Determine noise estimates when unknown
for i=find(~isfinite(prec)),
    prec0 = spm_noise_estimate(Nii(i));
    if isfinite(prec0)
      prec(i) = 1/prec0.^2;
    end
end

% set all values to constant if NaN values were found in noise estimation
prec(find(~isfinite(prec))) = 1;

% Basis functions for algebra of rigid-body transform
%-----------------------------------------------------------------------
B = se3_basis;

% Set boundary conditions 
%-----------------------------------------------------------------------
spm_field('boundary',1); % Bias correction - Neumann
spm_diffeo('boundary',0);     % Diffeomorphism  - circulant

% Computations for figuring out how many grid levels are likely to work
%-----------------------------------------------------------------------
d = [0 0 0];
for i=1:numel(Nii),
    dm = [size(Nii(i).dat) 1];
    d  = max(d, dm(1:3));
end
d  = min(d);

% Specify highest resolution data
%-----------------------------------------------------------------------
clear pyramid
pyramid(max(ceil(log2(d)-log2(4)),1)) = struct('d',[1 1 1],'mat',eye(4),'img',[]);
for i=numel(Nii):-1:1,
    pyramid(1).img(i).f   = single(Nii(i).dat(:,:,:,1,1));
    pyramid(1).img(i).mat = Nii(i).mat;
end

% Generate sucessively lower resolution versions
%-----------------------------------------------------------------------
for level = 2:numel(pyramid),
    for i=numel(Nii):-1:1,
        pyramid(level).img(i).f   = spm_diffeo('restrict',pyramid(level-1).img(i).f);
        pyramid(level).img(i).f(~isfinite(pyramid(level).img(i).f)) = 0;
        s1 = [size(pyramid(level-1).img(i).f) 1];
        s2 = [size(pyramid(level  ).img(i).f) 1];
        s  = s1(1:3)./s2(1:3);
        pyramid(level).img(i).mat = pyramid(level-1).img(i).mat*[diag(s), (1-s(:))*0.5; 0 0 0 1];
        clear s1 s2
    end
end

% Convert all image data into B-spline coefficients (for interpolation)
%-----------------------------------------------------------------------
for level=1:numel(pyramid),
    for i=1:numel(Nii)
        pyramid(level).img(i).f = spm_diffeo('bsplinc',pyramid(level).img(i).f,ord);
    end
end

% Adjust precision for number of subjects
%-----------------------------------------------------------------------
%nscan = numel(pyramid(1).img);
%prec  = prec*(nscan-1)/nscan;

% Stuff for figuring out the orientation, dimensions etc of the highest resolution template
%-----------------------------------------------------------------------
Mat0 = cat(3,pyramid(1).img.mat);
dims = zeros(numel(Nii),3);
for i=1:size(dims,1),
    dims(i,:) = Nii(i).dat.dim(1:3);
end
[pyramid(1).mat,pyramid(1).d] = compute_avg_mat(Mat0,dims);
pyramid(1).sc   = abs(det(pyramid(1).mat(1:3,1:3)));
pyramid(1).prec = prec;

% Figure out template info for each sucessively lower resolution version
%-----------------------------------------------------------------------
for level=2:numel(pyramid),
    pyramid(level).d    = ceil(pyramid(level-1).d/2);
    s                   = pyramid(level-1).d./pyramid(level).d;
    pyramid(level).mat  = pyramid(level-1).mat*[diag(s), (1-s(:))*0.5; 0 0 0 1];

    % Relative scaling of regularisation
    pyramid(level).sc   = abs(det(pyramid(level).mat(1:3,1:3)));
    pyramid(level).prec = prec*sqrt(pyramid(level).sc/pyramid(1).sc); % Note that the sqrt is ad hoc
end


nlevels = numel(pyramid);

for level=nlevels:-1:1, % Loop over resolutions, starting with the lowest

    % Collect data
    %-----------------------------------------------------------------------
    img       = pyramid(level).img;
    M_avg     = pyramid(level).mat;
    d         = pyramid(level).d;
    vx        = sqrt(sum(pyramid(level).mat(1:3,1:3).^2));
    sc        = pyramid(level).sc;
    prec      = pyramid(level).prec;

    if level==nlevels,
        % If lowest resolution, initialise parameter estimates to zero
        %-----------------------------------------------------------------------
        clear param
        bias_est = zeros(numel(Nii),1);
        for i=numel(Nii):-1:1,
            bias_est(i) = log(mean(mean(mean(img(i).f))));
        end
        bias_est = bias_est - mean(bias_est);
        for i=numel(Nii):-1:1,
            param(i) = struct('R',   eye(4), 'r', zeros(6,1),...
                              'bias',[],     'eb',0,...
                              'v0',  [],     'ev',0, 'y',[], 'J',[],...
                              's2',  1,      'ss',1);

            if all(isfinite(b_settings(i,:))),
                param(i).bias = zeros(size(img(i).f),'single')+bias_est(i);
            end

        end
    else
        % Initialise parameter estimates by prolongation of previous lower resolution versions.
        %-----------------------------------------------------------------------
        for i=1:numel(Nii),

            if all(isfinite(b_settings(i,:))),
                vxi           = sqrt(sum(img(i).mat(1:3,1:3).^2));
                spm_diffeo('boundary',1);
                param(i).bias = spm_diffeo('resize',param(i).bias,size(img(i).f));
                spm_diffeo('boundary',0);
                bmom          = spm_field('vel2mom', param(i).bias, [vxi b_settings(i,:)*sc]);
                param(i).eb   = sum(bmom(:).*param(i).bias(:));
                clear bmom
            else
                param(i).bias = [];
                param(i).eb   = 0;
            end

            param(i).v0 = [];
            param(i).ev = 0;
            param(i).y  = [];
            param(i).J  = [];
        end

        % Remove lower resolution versions that are no longer needed
        pyramid = pyramid(1:(end-1));
    end

    spm_plot_convergence('Clear');
    spm_plot_convergence('Init',['Optimising (level ' num2str(level) ')'],'Objective Function','Step');
    for iter=1:(2*2^(level-1)+1), % Use more iterations at lower resolutions (its faster, so may as well)


        % Compute deformations from initial velocities
        %-----------------------------------------------------------------------

        if true,
            % Rigid-body
            %=======================================================================
            % Recompute template data (with gradients)
            %-----------------------------------------------------------------------
            [mu,ss,nvox,D] = compute_mean(pyramid(level), param, ord);
            % for i=1:numel(param), fprintf('  %12.5g %12.5g %12.5g', prec(i)*ss(i), param(i).eb, param(i).ev); end; fprintf('  0\n');
            
            % create mask at final level to obtain masked registration
            if (level == 1) && (iter == 1) && use_brainmask
                PG = fullfile(spm('Dir'),'toolbox','FieldMap','T1.nii');
                PB = fullfile(spm('Dir'),'toolbox','FieldMap','brainmask.nii');

                [pth,nam]   = fileparts(Nii(1).dat.fname);
                nam         = fullfile(pth,['avg_' nam '.nii']);
                Nio         = nifti;
                Nio.dat     = file_array(nam,size(mu),'int16-be',0,max(max(mu(:))/32767,-min(mu(:))/32768),0);
                Nio.mat     = M_avg;
                Nio.mat0    = Nio.mat;
                Nio.mat_intent  = 'Aligned';
                Nio.mat0_intent = Nio.mat_intent;
                Nio.descrip = sprintf('Average of %d', numel(param));
                create(Nio);
                Nio.dat(:,:,:) = mu;
                PF = nam;

                [Ym, brainmask] = cat_long_APP(PF,PG,PB);

                spm_plot_convergence('Clear');
                spm_plot_convergence('Init',['Optimising (level ' num2str(level) ')'],'Objective Function','Step');
                clear Ym
            else 
              brainmask = [];
            end

            % Compute objective function (approximately)
            %-----------------------------------------------------------------------
            ll = 0;
            for i=1:numel(param),
                param(i).ss = ss(i);
                ll          = ll - 0.5*prec(i)*param(i).ss - 0.5*param(i).eb - 0.5*param(i).ev;
            end
            spm_plot_convergence('set',ll);

            for i=1:numel(img),
                % Gauss-Newton update of logs of rigid-body matrices
                %-----------------------------------------------------------------------
                [R,dR]    = spm_dexpm(param(i).r,B);
                M         = img(i).mat\R*M_avg;
                [x1a,x2a] = ndgrid(1:d(1),1:d(2));

                Hess = zeros(12);
                gra  = zeros(12,1);
                for m=1:d(3)
                    dt    = ones(d(1:2),'single')*abs(det(M(1:3,1:3)));
                    y     = zeros([d(1:2) 1 3],'single');
                    [y(:,:,1),y(:,:,2),y(:,:,3)] = ndgrid(single(1:d(1)),single(1:d(2)),single(m));
                    y     = transform_warp(M,y);

                    f     = spm_diffeo('bsplins',img(i).f,y,ord);

                    if all(isfinite(b_settings(i,:))),
                        ebias = exp(spm_diffeo('samp',param(i).bias,y));
                    else
                        ebias = ones(size(f),'single');
                    end

                    b     = f-mu(:,:,m).*ebias;

                    if (level == 1) && use_brainmask && ~isempty(brainmask)
                      msk = isfinite(b) & (brainmask(:,:,m) > 0.25);
                    else
                      msk = isfinite(b);
                    end
                    ebias = ebias(msk);
                    b     = b(msk);
                    dt    = dt(msk);
                    x1    = x1a(msk);
                    x2    = x2a(msk);
                    d1    = D{1}(:,:,m);d1 = d1(msk).*ebias;
                    d2    = D{2}(:,:,m);d2 = d2(msk).*ebias;
                    d3    = D{3}(:,:,m);d3 = d3(msk).*ebias;

                    A     = [x1(:).*d1(:) x1(:).*d2(:) x1(:).*d3(:) ...
                             x2(:).*d1(:) x2(:).*d2(:) x2(:).*d3(:) ...
                                  m*d1(:)      m*d2(:)      m*d3(:) ...
                                    d1(:)        d2(:)        d3(:)];

                    Hess  = Hess + double(A'*bsxfun(@times,A,dt));
                    gra   = gra  + double(A'*(dt.*b));

                    clear dt y f ebias b msk x1 x2 d1 d2 d3 A 
                end

                dA = zeros(12,6);
                for m=1:6,
                    tmp     = (R*M_avg)\dR(:,:,m)*M_avg;
                    dA(:,m) = reshape(tmp(1:3,:),12,1);
                end

                Hess       = dA'*Hess*dA*prec(i);
                gra        = dA'*gra*prec(i);
                param(i).r = param(i).r - Hess\gra;

                clear R dR M x1a x2a dA tmp Hess gra
            end
            clear mu D

            % Mean correct the rigid-body transforms and compute exponentials
            % Note that this gives us a Karcher mean.
            %-----------------------------------------------------------------------
            r_avg = mean(cat(2,param.r),2);
            for i=1:numel(param),
                param(i).r = param(i).r-r_avg;
                param(i).R = spm_dexpm(param(i).r,B);
            end
            clear r_avg
        end

        if any(all(isfinite(b_settings),2)),
            % Bias field
            %=======================================================================
            % Recompute template data
            %-----------------------------------------------------------------------
            [mu,ss] = compute_mean(pyramid(level), param, ord);
            % for i=1:numel(param), fprintf('  %12.5g %12.5g %12.5g', prec(i)*ss(i), param(i).eb, param(i).ev); end; fprintf('  1\n');

            % Compute objective function (approximately)
            %-----------------------------------------------------------------------
            ll = 0;
            for i=1:numel(param),
                param(i).ss = ss(i);
                ll          = ll - 0.5*prec(i)*param(i).ss - 0.5*param(i).eb - 0.5*param(i).ev;
            end
            spm_plot_convergence('set',ll);

            for i=1:numel(img),
                if all(isfinite(b_settings(i,:))),
                    % Gauss-Newton update of logs of bias field.
                    % Note that 1st and second derivatives are computed in template space
                    % and subsequently pushed back to native space for re-estimation.
                    %-----------------------------------------------------------------------
                    M    = img(i).mat\param(i).R*M_avg;
                    gra  = zeros(d,'single');
                    Hess = zeros(d,'single');

                    for m=1:d(3)
                        dt    = ones(d(1:2),'single')*abs(det(M(1:3,1:3)));
                        y     = zeros([d(1:2) 1 3],'single');
                        [y(:,:,1),y(:,:,2),y(:,:,3)] = ndgrid(single(1:d(1)),single(1:d(2)),m);
                        y     = transform_warp(M,y);

                        f           = spm_diffeo('bsplins',img(i).f,y,ord);
                        ebias       = exp(spm_diffeo('samp',param(i).bias,y));

                        msk         = isfinite(f) & isfinite(ebias);
                        smu         = mu(:,:,m).*ebias;
                        f(~msk)     = 0;
                        smu(~msk)   = 0;
                        gra(:,:,m)  = smu.*(smu-f).*dt*prec(i);
                        Hess(:,:,m) = smu.*smu.*dt*prec(i);

                        clear dt y f ebias msk smu
                    end

                    % Push derivatives to native space
                    %-----------------------------------------------------------------------
                    y    = transform_warp(M,identity(d));
                    gra  = spm_diffeo('push',gra,y,size(param(i).bias));
                    Hess = spm_diffeo('push',Hess,y,size(param(i).bias));
                    clear y

                    vxi           = sqrt(sum(img(i).mat(1:3,1:3).^2));
                    gra           = gra + spm_field('vel2mom', param(i).bias, [vxi b_settings(i,:)*sc]);
                    param(i).bias = param(i).bias - spm_field(Hess,gra,[vxi b_settings(i,:)*sc 2 2]); % Gauss-Newton update
                    clear M gra Hess

                    % Compute part of objective function
                    %-----------------------------------------------------------------------
                    bmom          = spm_field('vel2mom', param(i).bias, [vxi b_settings(i,:)*sc]);
                    param(i).eb   = sum(bmom(:).*param(i).bias(:));
                    clear bmom vxi
                end
            end
            clear mu
        end

    end
end

% Figure out what needs to be saved
%-----------------------------------------------------------------------
need_avg  = false;
need_wimg = false;

if any(strcmp('avg',output)) || any(strcmp('wavg',output)),  need_avg  = true; end
if any(strcmp('wimg',output)),                               need_wimg = true; end

if need_avg,
    mu = compute_mean(pyramid(1), param, ord);

    if any(strcmp('wavg',output)),
        [pth,nam]   = fileparts(Nii(1).dat.fname);
        nam         = fullfile(pth,['avg_' nam '.nii']);
        Nio         = nifti;
        Nio.dat     = file_array(nam,size(mu),'int16-be',0,max(max(mu(:))/32767,-min(mu(:))/32768),0);
        Nio.mat     = M_avg;
        Nio.mat0    = Nio.mat;
        Nio.mat_intent  = 'Aligned';
        Nio.mat0_intent = Nio.mat_intent;
        Nio.descrip = sprintf('Average of %d', numel(param));
        create(Nio);
        Nio.dat(:,:,:) = mu;
        out.avg       = nam;
    else
        out.avg       = mu;
    end
end

vol = get_transformed_images(pyramid(1), param, ord);

% do some final bias correction between scans which is more effective and also masked
bias_nits = 8;
bias_fwhm = 60;
bias_reg = 1e-6;
bias_lmreg = 1e-6;

for i=1:numel(img)
  vol(:,:,:,i) = bias_correction(mu,vol(:,:,:,i),brainmask,pyramid(1),bias_nits,bias_fwhm,bias_reg,bias_lmreg);
end

if need_wimg,
    for i=1:numel(param),
        img = vol(:,:,:,i);
        [pth,nam]   = fileparts(Nii(i).dat.fname);
        nam         = fullfile(pth,['r' nam '.nii']);
        Nio         = nifti;
        Nio.dat     = file_array(nam,size(img),'int16-be',0,max(max(img(:))/32767,-min(img(:))/32768),0);
        Nio.mat     = M_avg;
        Nio.mat0    = Nio.mat;
        Nio.mat_intent  = 'Aligned';
        Nio.mat0_intent = Nio.mat_intent;
        Nio.descrip = sprintf('Realigned %d', numel(param));
        create(Nio);
        Nio.dat(:,:,:) = img;
        out.rimg{i}    = nam;
    end
end

clear mu vol;

spm_plot_convergence('Clear');
return;
%_______________________________________________________________________


%_______________________________________________________________________
function vol = get_transformed_images(data, param, ord)
d     = data.d;
M_avg = data.mat;
img   = data.img;
prec  = data.prec;

vol   = zeros([d numel(img)],'single');

for m=1:d(3),
    F  = cell(1,numel(img));
    Dt = cell(1,numel(img));
    Bf = cell(1,numel(img));

    mum = zeros(d(1:2),'single');

    for i=1:numel(img),
        M = img(i).mat\param(i).R*M_avg;
        if ~isempty(param(i).y),
            y     = transform_warp(M,param(i).y(:,:,m,:));
            Dt{i} = spm_diffeo('det',param(i).J(:,:,m,:,:))*abs(det(M(1:3,1:3)));
        else
            Dt{i} = ones(d(1:2),'single')*abs(det(M(1:3,1:3)));
            y     = zeros([d(1:2) 1 3],'single');
            [y(:,:,1),y(:,:,2),y(:,:,3)] = ndgrid(single(1:d(1)),single(1:d(2)),m);
            y     = transform_warp(M,y);
        end

        F{i}  = spm_diffeo('bsplins',img(i).f,y,ord);
        if ~isempty(param(i).bias),
            Bf{i} = exp(spm_diffeo('bsplins',param(i).bias,y,[1 1 1 ord(4:end)])); % Trilinear
        else
            Bf{i} = ones(d(1:2),'single');
        end

        f        = F{i};
        ebias    = Bf{i};
        dt       = Dt{i};
        scal     = Bf{i}.*Dt{i}*prec(i);
        vol(:,:,m,i) = F{i}.*scal;

    end
end

return;
%_______________________________________________________________________

%_______________________________________________________________________
function [mu,ss,nvox,D] = compute_mean(data, param, ord)
d     = data.d;
M_avg = data.mat;
img   = data.img;
prec  = data.prec;

mu   = zeros(d,'single');
nvox = zeros(numel(img),1);
ss   = zeros(numel(img),1);
if nargout>=4, % Compute gradients of template
    D  = {zeros(d,'single'),zeros(d,'single'),zeros(d,'single')};
end

for m=1:d(3),
    if nargout>=4,
        Dm1  = {zeros(d(1:2),'single'),zeros(d(1:2),'single'),zeros(d(1:2),'single')};
        Dm2  = {zeros(d(1:2),'single'),zeros(d(1:2),'single'),zeros(d(1:2),'single')};
        Df   = cell(3,1);
        Db   = cell(3,1);
    end
    F  = cell(1,numel(img));
    Dt = cell(1,numel(img));
    Bf = cell(1,numel(img));
    Msk= cell(1,numel(img));

    mum = zeros(d(1:2),'single');
    mgm = zeros(d(1:2),'single');

    for i=1:numel(img),
        M = img(i).mat\param(i).R*M_avg;
        if ~isempty(param(i).y),
            y     = transform_warp(M,param(i).y(:,:,m,:));
            Dt{i} = spm_diffeo('det',param(i).J(:,:,m,:,:))*abs(det(M(1:3,1:3)));
        else
            Dt{i} = ones(d(1:2),'single')*abs(det(M(1:3,1:3)));
            y     = zeros([d(1:2) 1 3],'single');
            [y(:,:,1),y(:,:,2),y(:,:,3)] = ndgrid(single(1:d(1)),single(1:d(2)),m);
            y     = transform_warp(M,y);
        end

        if nargout>=4,
            % Sample image and bias field, along with their gradients.  Gradients are
            % then transformed by multiplying with the transpose of the Jacobain matrices
            % of the deformation.
            if ~isempty(param(i).J),
                Jm = reshape(param(i).J(:,:,m,:,:),[d(1)*d(2),3,3]);
                Jm = reshape(reshape(permute(Jm,[1 2 3]),d(1)*d(2)*3,3)*M(1:3,1:3),[d(1) d(2) 3 3]);
            else
                Jm = repmat(reshape(single(M(1:3,1:3)),[1 1 3 3]),[d(1) d(2) 1 1]);
            end

            [F{i} ,d1,d2,d3]  = spm_diffeo('bsplins',img(i).f,y,ord); 
            Df{1} = Jm(:,:,1,1).*d1 + Jm(:,:,2,1).*d2 + Jm(:,:,3,1).*d3;
            Df{2} = Jm(:,:,1,2).*d1 + Jm(:,:,2,2).*d2 + Jm(:,:,3,2).*d3;
            Df{3} = Jm(:,:,1,3).*d1 + Jm(:,:,2,3).*d2 + Jm(:,:,3,3).*d3;

            if ~isempty(param(i).bias),
                [Bf{i},d1,d2,d3]  = spm_diffeo('bsplins',param(i).bias,y,[1 1 1 ord(4:end)]); % Trilinear
                Bf{i} = exp(Bf{i});
                Db{1} = Jm(:,:,1,1).*d1 + Jm(:,:,2,1).*d2 + Jm(:,:,3,1).*d3;
                Db{2} = Jm(:,:,1,2).*d1 + Jm(:,:,2,2).*d2 + Jm(:,:,3,2).*d3;
                Db{3} = Jm(:,:,1,3).*d1 + Jm(:,:,2,3).*d2 + Jm(:,:,3,3).*d3;
            else
                Bf{i} =  ones(d(1:2),'single');
                Db{1} = zeros(d(1:2),'single');
                Db{2} = zeros(d(1:2),'single');
                Db{3} = zeros(d(1:2),'single');
            end
            clear d1 d2 d3
        else
            F{i}  = spm_diffeo('bsplins',img(i).f,y,ord);
            if ~isempty(param(i).bias),
                Bf{i} = exp(spm_diffeo('bsplins',param(i).bias,y,[1 1 1 ord(4:end)])); % Trilinear
            else
                Bf{i} = ones(d(1:2),'single');
            end
        end

        msk      = isfinite(F{i}) & isfinite(Bf{i});
        Msk{i}   = msk;
        f        = F{i}(msk);
        ebias    = Bf{i}(msk);
        dt       = Dt{i}(msk);
        scal     = ebias.*dt*prec(i);
        mum(msk) = mum(msk) + f.*scal;
        mgm(msk) = mgm(msk) + ebias.*scal;

        if nargout>=4
            % For computing gradients
            Dm1{1}(msk) = Dm1{1}(msk) + (Df{1}(msk) + f.*Db{1}(msk)).*scal;
            Dm1{2}(msk) = Dm1{2}(msk) + (Df{2}(msk) + f.*Db{2}(msk)).*scal;
            Dm1{3}(msk) = Dm1{3}(msk) + (Df{3}(msk) + f.*Db{3}(msk)).*scal;

            scal        = ebias.*scal;
            Dm2{1}(msk) = Dm2{1}(msk) + Db{1}(msk).*scal;
            Dm2{2}(msk) = Dm2{2}(msk) + Db{2}(msk).*scal;
            Dm2{3}(msk) = Dm2{3}(msk) + Db{3}(msk).*scal;
        end
    end
    mgm       = mgm + eps;
    mu(:,:,m) = mum./mgm; % Weighted mean

    if nargout>=2,
        if nargout>=4,
            % Compute "gradients of template (mu)".  Note that the true gradients
            % would incorporate the gradients of the Jacobians, but we do not want
            % these to be part of the "template gradients".
            wt          = 2*mum./(mgm.*mgm);
            D{1}(:,:,m) = Dm1{1}./mgm - Dm2{1}.*wt;
            D{2}(:,:,m) = Dm1{2}./mgm - Dm2{2}.*wt;
            D{3}(:,:,m) = Dm1{3}./mgm - Dm2{3}.*wt;
        end

        % Compute matching term
        for i=1:numel(img),
            msk      = Msk{i};
            f        = F{i}(msk);
            ebias    = Bf{i}(msk);
            dt       = Dt{i}(msk);
            mum      = mu(:,:,m);
            mum      = mum(msk);
            nvox(i)  = nvox(i) + sum(dt);
            ss(i)    = ss(i) + sum((f-mum.*ebias).^2.*dt);
        end
    end
end

for i=1:numel(img)
    ss(i) = ss(i)/nvox(i)*numel(img(i).f);
end
return;
%_______________________________________________________________________

%_______________________________________________________________________
function y1 = transform_warp(M,y)
% Affine transformation of a deformation
d  = size(y);
y1 = reshape(bsxfun(@plus,reshape(y,[prod(d(1:3)),3])*single(M(1:3,1:3)'),single(M(1:3,4)')),d);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function y = identity(d)
% Generate an identity transform of size d(1) x d(2) x d(3)
y = zeros([d(1:3) 3],'single');
[y(:,:,:,1),y(:,:,:,2),y(:,:,:,3)] = ndgrid(single(1:d(1)),single(1:d(2)),single(1:d(3)));
%_______________________________________________________________________

%_______________________________________________________________________
function [M_avg,d] = compute_avg_mat(Mat0,dims)
% Compute an average voxel-to-world mapping and suitable dimensions
% FORMAT [M_avg,d] = compute_avg_mat(Mat0,dims)
% Mat0  - array of matrices (4x4xN)
% dims  - image dimensions (Nx3)
% M_avg - voxel-to-world mapping
% d     - dimensions for average image
%

% Rigid-body matrices computed from exp(p(1)*B(:,:,1)+p(2)+B(:,:,2)...)
%-----------------------------------------------------------------------
B = se3_basis;

% Find combination of 90 degree rotations and flips that brings all
% the matrices closest to axial
%-----------------------------------------------------------------------
Matrices = Mat0;
pmatrix  = [1,2,3; 2,1,3; 3,1,2; 3,2,1; 1,3,2; 2,3,1];
for i=1:size(Matrices,3)
    vx    = sqrt(sum(Matrices(1:3,1:3,i).^2));
    tmp   = Matrices(:,:,i)/diag([vx 1]);
    R     = tmp(1:3,1:3);
    minss = Inf;
    minR  = eye(3);
    for i1=1:6,
        R1 = zeros(3);
        R1(pmatrix(i1,1),1)=1;
        R1(pmatrix(i1,2),2)=1;
        R1(pmatrix(i1,3),3)=1;
        for i2=0:7,
            F  = diag([bitand(i2,1)*2-1, bitand(i2,2)-1, bitand(i2,4)/2-1]);
            R2 = F*R1;
            ss = sum(sum((R/R2-eye(3)).^2));
            if ss<minss,
                minss = ss;
                minR  = R2;
            end
        end
    end
    rdim = abs(minR*dims(i,:)');
    R2   = inv(minR);
    minR = [R2 R2*((sum(R2,1)'-1)/2.*(rdim+1)); 0 0 0 1];
    Matrices(:,:,i) = Matrices(:,:,i)*minR;
end

% Average of these matrices
%-----------------------------------------------------------------------
M_avg = spm_meanm(Matrices);

% If average involves shears, then find the closest matrix that does not
% require them
%-----------------------------------------------------------------------
p = spm_imatrix(M_avg);
if sum(p(10:12).^2)>1e-8,

    % Zooms computed from exp(p(7)*B2(:,:,1)+p(8)*B2(:,:,2)+p(9)*B2(:,:,3))
    %-----------------------------------------------------------------------
    B2        = zeros(4,4,3);
    B2(1,1,1) = 1;
    B2(2,2,2) = 1;
    B2(3,3,3) = 1;

    p      = zeros(9,1); % Parameters
    for it=1:10000,
        [R,dR] = spm_dexpm(p(1:6),B);  % Rotations + Translations
        [Z,dZ] = spm_dexpm(p(7:9),B2); % Zooms

        M  = R*Z; % Voxel-to-world estimate
        dM = zeros(4,4,6);
        for i=1:6, dM(:,:,i)   = dR(:,:,i)*Z; end
        for i=1:3, dM(:,:,i+6) = R*dZ(:,:,i); end
        dM = reshape(dM,[16,9]);

        d   = M(:)-M_avg(:); % Difference
        gr  = dM'*d;         % Gradient
        Hes = dM'*dM;        % Hessian
        p   = p - Hes\gr;    % Gauss-Newton update
        if sum(gr.^2)<1e-8, break; end
    end
    M_avg = M;
end

% Ensure that the FoV covers all images, with a few voxels to spare
%-----------------------------------------------------------------------
mn    =  Inf*ones(3,1);
mx    = -Inf*ones(3,1);
for i=1:size(Mat0,3),
    dm      = [dims(i,:) 1 1];
    corners = [
        1 dm(1)    1  dm(1)   1  dm(1)    1  dm(1)
        1    1  dm(2) dm(2)   1     1  dm(2) dm(2)
        1    1     1     1 dm(3) dm(3) dm(3) dm(3)
        1    1     1     1    1     1     1     1];
    M  = M_avg\Mat0(:,:,i);
    vx = M(1:3,:)*corners;
    mx = max(mx,max(vx,[],2));
    mn = min(mn,min(vx,[],2));
end
mx    = ceil(mx);
mn    = floor(mn);
d     = (mx-mn+7)';
M_avg = M_avg * [eye(3) mn-4; 0 0 0 1];

% check whether entry for x-voxel-size is positive and correct it
[mx, ind] = max(abs(M_avg(1,1:3)));
if (M_avg(1,ind) > 0)
  M_avg(1,:) = -M_avg(1,:);
end
return;
%_______________________________________________________________________

%_______________________________________________________________________
function B = se3_basis
% Basis functions for the lie algebra of the special Eucliden group
% (SE(3)).
B        = zeros(4,4,6);
B(1,4,1) = 1;
B(2,4,2) = 1;
B(3,4,3) = 1;
B([1,2],[1,2],4) = [0 1;-1 0];
B([3,1],[3,1],5) = [0 1;-1 0];
B([2,3],[2,3],6) = [0 1;-1 0];
return
%_______________________________________________________________________

%_______________________________________________________________________
function t = transf(B1,B2,B3,T)
d2 = [size(T) 1];
t1 = reshape(reshape(T, d2(1)*d2(2),d2(3))*B3', d2(1), d2(2));
t  = B1*t1*B2';
return;
%_______________________________________________________________________

%_______________________________________________________________________
function dat = bias_correction(volG,volF,brainmask,pyramid,nits,fwhm,reg1,reg2)
% This function is intended for doing bias correction between scans inside a defined mask
% A version of the second image is returned out, that has the
% same bias as that of the first image.  

vx      = sqrt(sum(pyramid.mat(1:3,1:3).^2));
d       = size(volG);
sd      = vx(1)*size(volG,1)/fwhm; d3(1) = ceil(sd*2); krn_x   = exp(-(0:(d3(1)-1)).^2/sd.^2)/sqrt(vx(1));
sd      = vx(2)*size(volG,2)/fwhm; d3(2) = ceil(sd*2); krn_y   = exp(-(0:(d3(2)-1)).^2/sd.^2)/sqrt(vx(2));
sd      = vx(3)*size(volG,3)/fwhm; d3(3) = ceil(sd*2); krn_z   = exp(-(0:(d3(3)-1)).^2/sd.^2)/sqrt(vx(3));
Cbias   = kron(krn_z,kron(krn_y,krn_x)).^(-2)*prod(d)*reg1;
Cbias   = sparse(1:length(Cbias),1:length(Cbias),Cbias,length(Cbias),length(Cbias));
B3bias  = spm_dctmtx(d(3),d3(3));
B2bias  = spm_dctmtx(d(2),d3(2));
B1bias  = spm_dctmtx(d(1),d3(1));
lmRb    = speye(size(Cbias))*prod(d)*reg2;
Tbias   = zeros(d3);

% correct global scaling
thG = mean(volG(isfinite(volG)))/8; thG = mean(volG(volG>thG));
thF = mean(volF(isfinite(volF)))/8; thF = mean(volF(volF>thF));
volF = volF*thG/thF;
thF = mean(volF(isfinite(volF)))/8; thF = mean(volF(volF>thF));

ll = Inf;
try
    spm_plot_convergence('Init','Bias Correction','- Log-likelihood','Iteration');
catch
    spm_chi2_plot('Init','Bias Correction','- Log-likelihood','Iteration');
end
for subit=1:nits,

    % Compute objective function and its 1st and second derivatives
    Alpha = zeros(prod(d3),prod(d3)); % Second derivatives
    Beta  = zeros(prod(d3),1); % First derivatives
    oll   = ll;
    ll    = 0.5*Tbias(:)'*Cbias*Tbias(:);

    for z=1:size(volG,3),
        f1o = volF(:,:,z);
        f2o = volG(:,:,z);
        f1o(~isfinite(f1o)) = 0;
        f2o(~isfinite(f2o)) = 0;
        if ~isempty(brainmask)
          msk = (f1o==0) & (f2o==0) & (brainmask(:,:,z) < 0.25);
        else
          msk = (f1o==0) & (f2o==0);
        end
        f1o(msk) = 0;
        f2o(msk) = 0;
        ro       = transf(B1bias,B2bias,B3bias(z,:),Tbias);
        if ~isempty(brainmask)
          msk      = (abs(ro)>0.01) & (brainmask(:,:,z) > 0.25); 
        else
          msk      = abs(ro)>0.01; 
        end

        % Use the form based on an integral for bias that is
        % far from uniform.
        f1  = f1o(msk);
        f2  = f2o(msk);
        r   = ro(msk);
        e   = exp(r);
        t1  = (f2.*e-f1);
        t2  = (f1./e-f2);
        ll  = ll + 1/4*sum(sum((t1.^2-t2.^2)./r));
        wt1 = zeros(size(f1o));
        wt2 = zeros(size(f1o));
        wt1(msk) = (2*(t1.*f2.*e+t2.*f1./e)./r + (t2.^2-t1.^2)./r.^2)/4;
        wt2(msk) = ((f2.^2.*e.^2-f1.^2./e.^2+t1.*f2.*e-t2.*f1./e)./r/2 ...
                 - (t1.*f2.*e+t2.*f1./e)./r.^2 + (t1.^2-t2.^2)./r.^3/2);

        % Use the simple symmetric form for bias close to uniform
        f1  = f1o(~msk);
        f2  = f2o(~msk);
        r   = ro(~msk);
        e   = exp(r);
        t1  = (f2.*e-f1);
        t2  = (f1./e-f2);
        ll  = ll + (sum(t1.^2)+sum(t2.^2))/4;
        wt1(~msk) = (t1.*f2.*e-t2.*f1./e)/2;
        wt2(~msk) = ((f2.*e).^2+t1.*f2.*e + (f1./e).^2+t2.*f1./e)/2;

        b3    = B3bias(z,:)';
        Beta  = Beta  + kron(b3,spm_krutil(wt1,B1bias,B2bias,0));
        Alpha = Alpha + kron(b3*b3',spm_krutil(wt2,B1bias,B2bias,1));
    end;
    try
        spm_plot_convergence('Set',ll/prod(d));
    catch
        spm_chi2_plot('Set',ll/prod(d));
    end

    if subit > 1 && ll>oll,
        % Hasn't improved, so go back to previous solution
        Tbias = oTbias;
        ll    = oll;
        lmRb  = lmRb*10;
    else
        % Accept new solution
        oTbias = Tbias;
        Tbias  = Tbias(:);
        Tbias  = Tbias - (Alpha + Cbias + lmRb)\(Beta + Cbias*Tbias);
        Tbias  = reshape(Tbias,d3);
    end;
end;

dat = zeros(size(volG));

for z=1:size(volG,3),
    tmp = volF(:,:,z);
    tmp(~isfinite(tmp)) = 0;
    dat(:,:,z) = tmp;
    r  = transf(B1bias,B2bias,B3bias(z,:),Tbias);
    dat(:,:,z) = dat(:,:,z)./exp(r);
end;

return;
%_______________________________________________________________________

