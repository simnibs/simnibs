function Ycls = cat_main(res,tpm,job)
% Write out CAT preprocessed data
%
% FORMAT Ycls = cat_main(res,tpm,job)
%
% based on John Ashburners version of
% spm_preproc_write8.m 2531 2008-12-05 18:59:26Z john $
%
% ______________________________________________________________________
% Christian Gaser
% $Id: cat_main.m 1274 2018-02-12 08:36:30Z gaser $

%#ok<*ASGLU>

% if there is a breakpoint in this file set debug=1 and do not clear temporary variables 
dbs   = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end


global cat_err_res; % for CAT error report

tc = [cat(1,job.tissue(:).native) cat(1,job.tissue(:).warped)]; 
clsint = @(x) round( sum(res.mn(res.lkp==x) .* res.mg(res.lkp==x)') * 10^5)/10^5;

%% complete job structure
defr.ppe = struct(); 
res = cat_io_checkinopt(res,defr);


def.cati            = 0;
def.color.error     = [0.8 0.0 0.0];
def.color.warning   = [0.0 0.0 1.0];
def.color.warning   = [0.8 0.9 0.3];
def.color.highlight = [0.2 0.2 0.8];
job = cat_io_checkinopt(job,def);

% definition of subfolders
if job.extopts.subfolders
  mrifolder     = 'mri';
  reportfolder  = 'report';
  labelfolder   = 'label';
else
  mrifolder     = '';
  reportfolder  = '';
  labelfolder   = '';
end

M1   = tpm.M;

% Sort out bounding box etc
[bb1,vx1] = spm_get_bbox(tpm.V(1), 'old');
bb = job.extopts.bb;
vx = job.extopts.vox(1);
bb(~isfinite(bb)) = bb1(~isfinite(bb));
if ~isfinite(vx), vx = abs(prod(vx1))^(1/3); end; 
bb(1,:) = vx.*round(bb(1,:)./vx);
bb(2,:) = vx.*round(bb(2,:)./vx);
res.bb = bb; 
clear vx vx1 bb1   

if isfield(res,'mg'),
    lkp = res.lkp;
    Kb  = max(lkp);
else
    Kb  = size(res.intensity(1).lik,2);
end

N = numel(res.image);
if N > 1
  warning('CAT12:noMultiChannel',...
    'CAT12 does not support multiple channels. Only the first channel will be used.');
end

do_dartel = 1 + (job.extopts.regstr(1)~=0);      % always use dartel (do_dartel=1) or shooting (do_dartel=2) normalization
if do_dartel
  need_dartel = any(job.output.warps) || ...
    job.output.bias.warped || ... job.output.bias.dartel || ...
    job.output.label.warped || ... job.output.label.dartel || ...
    any(any(tc(:,[4 5 6]))) || job.output.jacobian.warped || ...
    ... job.output.surface || ...
    job.output.ROI || ...
    any([job.output.atlas.warped]) || ...
    numel(job.extopts.regstr)>1 || ...
    numel(job.extopts.vox)>1;
  if ~need_dartel
    %fprintf('Option for Dartel output was deselected because no normalized images need to be saved.\n');  
    do_dartel = 0;
  end
end
if do_dartel<2, job.extopts.templates = job.extopts.darteltpms; else job.extopts.templates = job.extopts.shootingtpms; end % for LAS
res.do_dartel = do_dartel;


stime = cat_io_cmd('SPM preprocessing 2 (write)');

% remove noise/interpolation prefix
VT  = res.image(1);  % denoised/interpolated n*.nii
VT0 = res.image0(1); % original 
fname0 = VT0.fname;
[pth,nam] = spm_fileparts(VT0.fname); 

% delete old xml file 
oldxml = fullfile(pth,reportfolder,['cat_' nam '.xml']);  
if exist(oldxml,'file'), delete(oldxml); end
clear oldxml

d = VT.dim(1:3);
[x1,x2,o] = ndgrid(1:d(1),1:d(2),1);
x3  = 1:d(3);



chan(N) = struct('B1',[],'B2',[],'B3',[],'T',[],'Nc',[],'Nf',[],'ind',[]);
for n=1:N,
    d3         = [size(res.Tbias{n}) 1];
    chan(n).B3 = spm_dctmtx(d(3),d3(3),x3);
    chan(n).B2 = spm_dctmtx(d(2),d3(2),x2(1,:)');
    chan(n).B1 = spm_dctmtx(d(1),d3(1),x1(:,1));
    chan(n).T  = res.Tbias{n};

    [pth1,nam1] = spm_fileparts(res.image(n).fname);
    
    %[pth2,mri2] = spm_fileparts(pth1); if strcmp(mri2,'mri'), pth1=pth2; end; clear pth2 mri2; % mri subdir
    
    if job.extopts.NCstr~=0
      nam1 = nam1(2:end);
    end
    chan(n).ind      = res.image(n).n;

    if job.output.bias.native,
        chan(n).Nc      = nifti;
        chan(n).Nc.dat  = file_array(fullfile(pth,mrifolder,['m', nam1, '.nii']),...
                                 res.image(n).dim(1:3),...
                                 [spm_type('float32') spm_platform('bigend')],...
                                 0,1,0);
        chan(n).Nc.mat  = res.image(n).mat;
        chan(n).Nc.mat0 = res.image(n).mat;
        chan(n).Nc.descrip = 'Bias corrected';
        create(chan(n).Nc);
    end
    clear pth1 
end
clear d3

prm     = [3 3 3 0 0 0];
Coef    = cell(1,3);
Coef{1} = spm_bsplinc(res.Twarp(:,:,:,1),prm);
Coef{2} = spm_bsplinc(res.Twarp(:,:,:,2),prm);
Coef{3} = spm_bsplinc(res.Twarp(:,:,:,3),prm);

if job.output.warps(2),
    Ndef      = nifti;
    Ndef.dat  = file_array(fullfile(pth,mrifolder,['iy_', nam1, '.nii']),...
                           [VT.dim(1:3),1,3],...
                           [spm_type('float32') spm_platform('bigend')],...
                           0,1,0);
    Ndef.mat  = VT.mat;
    Ndef.mat0 = VT.mat;
    Ndef.descrip = 'Inverse SPM Deformation';
    create(Ndef);
end

spm_progress_bar('init',length(x3),['Working on ' nam],'Planes completed');
M = M1\res.Affine*VT.mat;

Q = zeros([d(1:3),Kb],'single');

for z=1:length(x3),

    %% Bias corrected image
    cr = cell(1,N);
    for n=1:N,
        f = spm_sample_vol(res.image(n),x1,x2,o*x3(z),0);
        bf1 = exp(transf(chan(n).B1,chan(n).B2,chan(n).B3(z,:),chan(n).T));
        bf1(bf1>100) = 100;
        cr{n} = bf1.*f;
        
        % Write a plane of bias corrected data
        if job.output.bias.native,
            chan(n).Nc.dat(:,:,z,chan(n).ind(1),chan(n).ind(2)) = cr{n};
        end
        if ~isempty(chan(n).Nf),
            % Write a plane of bias field
            chan(n).Nf.dat(:,:,z,chan(n).ind(1),chan(n).ind(2)) = bf1;
        end;
    end

    % Compute the deformation (mapping voxels in image to voxels in TPM)
    [t1,t2,t3] = defs(Coef,z,res.MT,prm,x1,x2,x3,M);

    if exist('Ndef','var'),
        % Write out the deformation to file, adjusting it so mapping is
        % to voxels (voxels in image to mm in TPM)
        tmp = M1(1,1)*t1 + M1(1,2)*t2 + M1(1,3)*t3 + M1(1,4);
        Ndef.dat(:,:,z,1,1) = tmp;
        tmp = M1(2,1)*t1 + M1(2,2)*t2 + M1(2,3)*t3 + M1(2,4);
        Ndef.dat(:,:,z,1,2) = tmp;
        tmp = M1(3,1)*t1 + M1(3,2)*t2 + M1(3,3)*t3 + M1(3,4);
        Ndef.dat(:,:,z,1,3) = tmp;
    end

    if isfield(res,'mg'),
        q   = zeros([d(1:2) Kb]);
        q1  = likelihoods(cr,[],res.mg,res.mn,res.vr);
        q1  = reshape(q1,[d(1:2),numel(res.mg)]);
        b   = spm_sample_priors8(tpm,t1,t2,t3);
        wp  = res.wp;
        s   = zeros(size(b{1}));
        for k1 = 1:Kb,
            b{k1} = wp(k1)*b{k1};
            s     = s + b{k1};
        end
        for k1=1:Kb,
            q(:,:,k1) = sum(q1(:,:,lkp==k1),3).*(b{k1}./s);
        end
    else
        % Nonparametric representation of intensity distributions
        q   = spm_sample_priors8(tpm,t1,t2,t3);
        wp  = res.wp;
        s   = zeros(size(q{1}));
        for k1 = 1:Kb,
            q{k1} = wp(k1)*q{k1};
            s     = s + q{k1};
        end
        for k1 = 1:Kb,
            q{k1} = q{k1}./s;
        end
        q   = cat(3,q{:});

        for n=1:N,
            tmp = round(cr{n}*res.intensity(n).interscal(2) + res.intensity(n).interscal(1));
            tmp = min(max(tmp,1),size(res.intensity(n).lik,1));
            for k1=1:Kb,
                likelihood = res.intensity(n).lik(:,k1);
                q(:,:,k1)  = q(:,:,k1).*likelihood(tmp);
            end
        end
    end
    Q(:,:,z,:) = reshape(q,[d(1:2),1,Kb]);


    % initialize Yy only at first slice
    if z==1
        Yy = zeros([VT.dim(1:3),3],'single');
    end
    Yy(:,:,z,1) = t1;
    Yy(:,:,z,2) = t2;
    Yy(:,:,z,3) = t3;

    spm_progress_bar('set',z);
end


if ~isfield(res,'spmpp')
  % cleanup with brain mask - required for ngaus [1 1 2 4 3 2] and R1/MP2Rage like data 
  YbA = single(spm_sample_vol(tpm.V(1),double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),1)) + ...
        single(spm_sample_vol(tpm.V(2),double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),1)) + ...
        single(spm_sample_vol(tpm.V(3),double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),1));
  YbA = reshape(YbA,VT.dim(1:3));   
  YbA = cat_vol_smooth3X(cat_vol_smooth3X(YbA,2)>0.1,2); % dilate + smooth 
  for i=1:3, Q(:,:,:,i) = Q(:,:,:,i) .* YbA; end
  
  clear YbA
  
  % sum up classes 
  sQ = (sum(Q,4)+eps)/255; P = zeros([d(1:3),Kb],'uint8');
  for k1=1:size(Q,4)
    P(:,:,:,k1) = cat_vol_ctype(round(Q(:,:,:,k1)./sQ));
  end
  clear sQ Qspm_progress_bar('clear');

  vx_vol  = sqrt(sum(VT.mat(1:3,1:3).^2));    % voxel size of the processed image
  vx_volr = sqrt(sum(VT0.mat(1:3,1:3).^2));   % voxel size of the original image 
  vx_volp = prod(vx_vol)/1000;
  voli    = @(v) (v ./ (pi * 4./3)).^(1/3);   % volume > radius

  % cleanup for high resolution data
  if max(vx_vol)<1.5 && mean(vx_vol)<1.3
    P = clean_gwc(P,1);
  end

  % load bias corrected image
  % restrict bias field to maximum of 3 and a minimum of 0.1
  % (sometimes artefacts at the borders can cause huge values in bias field)
  Ybf  = zeros(VT.dim(1:3),'single');
  Ysrc = zeros(VT.dim(1:3),'single');
  Ysrc(isnan(Ysrc)) = 0; % prevent NaN
  for z=1:length(x3),
      f = spm_sample_vol(VT,x1,x2,o*x3(z),0);
      bf1 = exp(transf(chan(1).B1,chan(1).B2,chan(1).B3(z,:),chan(1).T));
      Ysrc(:,:,z) = single(bf1 .* f); 
      Ybf(:,:,z)  = single(bf1 .* ones(size(f)));
  end
  clear chan o x1 x2 x3 bf1 f z

  Ycls = {zeros(d,'uint8') zeros(d,'uint8') zeros(d,'uint8') ...
          zeros(d,'uint8') zeros(d,'uint8') zeros(d,'uint8')};

  

  %% create a new brainmask

  % median in case of WMHs!
  WMth = double(max( clsint(2),...
          cat_stat_nanmedian(cat_stat_nanmedian(cat_stat_nanmedian(Ysrc(P(:,:,:,2)>192)))))); 
  if clsint(3)>clsint(2) % invers
    CMth = clsint(3); 
  else
    CMth = min( [  clsint(1) - diff([clsint(1),WMth]) , clsint(3) ]);
  end
  T3th = [ CMth , clsint(1) , WMth];


  %    ds('l2','',vx_vol,Ysrc./WMth,Yp0>0.3,Ysrc./WMth,Yp0,80)
  Yp0  = single(P(:,:,:,3))/255/3 + single(P(:,:,:,1))/255*2/3 + single(P(:,:,:,2))/255;;
  if sum(Yp0(:)>0.3)<100 
    % this error often depends on a failed affine registration, where SPM
    % have to find the brain in the head or background
    BGth  = min(cat_stat_nanmean(Ysrc( P(:,:,:,end)>128 )),clsint(6));
    HDHth = clsint(5);
    HDLth = clsint(4);
    clsvol = nan(1,size(P,4)); for ci=1:size(P,4), Yct = P(:,:,:,ci)>128; clsvol(ci) = sum(Yct(:))*vx_volp; end; clear Yct; 
    if size(P,4)==6
        error('CAT:cat_main:SPMpreprocessing:emptySegmentation', ...
         sprintf(['Empty Segmentation: \n ' ...
          'Possibly the affine registration failed. Pleace check image orientation.\n' ...
          ' Tissue class:           %10s%10s%10s%10s%10s%10s\n' ...
          ' Rel. to image volume:   %10.2f%10.2f%10.2f%10.2f%10.2f%10.2f\n' ...
          ' Rel. to brain volume:   %10.2f%10.2f%10.2f%10.2f%10.2f%10.2f\n' ...
          ' Tissue intensity:       %10.2f%10.2f%10.2f%10.2f%10.2f%10.2f'],...
          'BG','CSF','GM','WM','HDH','HDL', ...
          [ clsvol([6 3 1 2 4 5])/cat_stat_nansum(clsvol)*100, clsvol([6 3 1 2 4 5])/cat_stat_nansum(clsvol(1:3))*100, BGth,T3th,HDHth,HDLth]));  %#ok<SPERR>
    elseif size(P,4)==4 % skull-stripped
        error('CAT:cat_main:SPMpreprocessing:emptySegmentation', ...
         sprintf(['Empty Segmentation: \n ' ...
          'Possibly the affine registration failed. Pleace check image orientation.\n' ...
          ' Tissue class:           %10s%10s%10s%10s\n' ...
          ' Rel. to image volume:   %10.2f%10.2f%10.2f%10.2f\n' ...
          ' Rel. to brain volume:   %10.2f%10.2f%10.2f%10.2f\n' ...
          ' Tissue intensity:       %10.2f%10.2f%10.2f%10.2f'],...
          'BG','CSF','GM','WM', ...
          [ clsvol([4 3 1 2])/cat_stat_nansum(clsvol)*100, clsvol([4 3 1 2])/cat_stat_nansum(clsvol(1:3))*100, BGth,T3th]));  %#ok<SPERR>
    else
        error('CAT:cat_main:SPMpreprocessing:emptySegmentation', ...
         sprintf(['Empty Segmentation: \n ' ...
         'Possibly the affine registration failed. Pleace check image orientation.\n'])); 
    end
  end
  
  %%
  Yp0(smooth3(cat_vol_morph(Yp0>0.3,'lo'))<0.5)=0; % not 1/6 because some ADNI scans have large "CSF" areas in the background 
  Yp0     = Yp0 .* cat_vol_morph(Yp0 & (Ysrc>WMth*0.05),'lc',2);
  Yp0toC  = @(Yp0,c) 1-min(1,abs(Yp0-c));

  % values are only used if errors occur
  cat_err_res.init.T3th = T3th; 
  cat_err_res.init.subjectmeasures.vol_abs_CGW = [prod(vx_vol)/1000 .* sum(Yp0toC(Yp0(:),1)), ... CSF
                                                  prod(vx_vol)/1000 .* sum(Yp0toC(Yp0(:),2)), ... GM 
                                                  prod(vx_vol)/1000 .* sum(Yp0toC(Yp0(:),3)), ... WM
                                                  prod(vx_vol)/1000 .* sum(Yp0toC(Yp0(:),4))];  % WMH
  cat_err_res.init.subjectmeasures.vol_TIV     =  sum(cat_err_res.init.subjectmeasures.vol_abs_CGW); 
  cat_err_res.init.subjectmeasures.vol_rel_CGW =  cat_err_res.init.subjectmeasures.vol_abs_CGW ./ ...
                                                  cat_err_res.init.subjectmeasures.vol_TIV;
  [cat_err_res.init.Yp0,cat_err_res.init.BB] = cat_vol_resize(Yp0,'reduceBrain',vx_vol,2,Yp0>0.5); 
  cat_err_res.init.Yp0 = cat_vol_ctype(cat_err_res.init.Yp0/3*255);


  if isfield(res,'msk')
    Ybg = ~res.msk.dat; 
    P4  = cat_vol_ctype( single(P(:,:,:,6)) .* (Ysrc<T3th(2))  .* (Ybg==0) + single(P(:,:,:,4)) .* (Ybg<1) ); % remove air in head
    P5  = cat_vol_ctype( single(P(:,:,:,6)) .* (Ysrc>=T3th(2)) .* (Ybg==0) + single(P(:,:,:,5)) .* (Ybg<1) ); % remove air in head
    P6  = cat_vol_ctype( single(sum(P(:,:,:,4:5),4)) .* (Ybg==1) + single(P(:,:,:,6)) .* (Ybg>0) ); % add objects/artifacts to background
    P(:,:,:,4) = P4;
    P(:,:,:,5) = P5;
    P(:,:,:,6) = P6;
    clear P4 P5 P6; 
  end
  
  %%
  skullstripped=max(res.lkp)==4; 
  if skullstripped % skull-stripped
    Yp0  = single(P(:,:,:,3))/255/3 + single(P(:,:,:,1))/255*2/3 + single(P(:,:,:,2))/255;
    Yb   = Yp0>0.5/3; 
    Ybb  = Yb; 

    P(:,:,:,6) = P(:,:,:,4); 
    P(:,:,:,4) = zeros(size(Yp0),'uint8');
    P(:,:,:,5) = zeros(size(Yp0),'uint8'); 
    res.lkp = [res.lkp 5 6];
    res.mn  = [res.mn(1:end-1),0,0,0];
    res.mg  = [res.mg(1:end-1);1;1;1];
    res.vr(1,1,numel(res.lkp)-1:numel(res.lkp)) = 0;
     
    [Ysrcb,Yp0,BB] = cat_vol_resize({Ysrc,Yp0},'reduceBrain',vx_vol,round(6/mean(vx_vol)),Yp0>1/3);
    Yg   = cat_vol_grad(Ysrcb/T3th(3),vx_vol);
    Ydiv = cat_vol_div(Ysrcb/T3th(3),vx_vol);
    Yg   = cat_vol_resize(Yg ,'dereduceBrain',BB);
    Ydiv = cat_vol_resize(Ydiv ,'dereduceBrain',BB);
  elseif (job.extopts.experimental || (job.extopts.INV && any(sign(diff(T3th))==-1))) && job.extopts.gcutstr>0
   % (sum( abs( (Ysrc(:)==0) - (Yp0(:)<0.5) ) ) / sum(Ysrc(:)==0)) < 0.1  || ...
  % use gcut2

    %   brad = voli(sum(Yp0(:)>0).*prod(vx_vol)/1000); 
    %   noise = nanstd(Ysrc(cat_vol_morph(Yp0>0.8,'o') & Yp0>0.99)/diff(T3th(2:3))); 

    Vl1 = spm_vol(job.extopts.cat12atlas{1});
    Yl1 = cat_vol_ctype(spm_sample_vol(Vl1,double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),0));
    Yl1 = reshape(Yl1,size(Ysrc)); [D,I] = cat_vbdist(single(Yl1>0)); Yl1 = Yl1(I);   

    %%
    clear D I
    Ybg = P(:,:,:,6)>128; 
    T3ths = [min(min(min(single(Ysrc(P(:,:,:,6)>128))))),...
             min( cat_stat_nanmean(Ysrc(Ybg(:))) + 2*cat_stat_nanstd(Ysrc(Ybg(:))) , ...
              mean([cat_stat_nanmean(Ysrc(Ybg(:))),min(T3th)])), ...
             T3th, T3th(3) + cat_stat_nanmean(diff(T3th))];
    T3thx = [0,0.05,1,2,3,4];
    if T3th(1)>T3th(3), T3thx = [0,0.05,1,2,3,2]; T3ths(end) = T3ths(2); end; 
    [T3ths,si] = sort(T3ths);
    T3thx      = T3thx(si);
    Ym = Ysrc+0; 
    for i=numel(T3ths):-1:2
      M = Ysrc>T3ths(i-1) & Ysrc<=T3ths(i);
      Ym(M(:)) = T3thx(i-1) + (Ysrc(M(:)) - T3ths(i-1))/diff(T3ths(i-1:i))*diff(T3thx(i-1:i));
    end
    M  = Ysrc>=T3ths(end); 
    Ym(M(:)) = numel(T3ths)/6 + (Ysrc(M(:)) - T3ths(i))/diff(T3ths(end-1:end))*diff(T3thx(i-1:i));    
    Ym = Ym / 3; 

    for k1=1:3
        Ycls{k1} = P(:,:,:,k1);
    end

    %%
    Yb = cat_main_gcut(Ym,Yp0>0.1,Ycls,Yl1,false(size(Ym)),vx_vol,...
      struct('gcutstr',0.1,'verb',0,'LAB',job.extopts.LAB,'LASstr',0));
    Ybb  = cat_vol_smooth3X(Yb,2); 
    
    [Ysrcb,Yp0,BB] = cat_vol_resize({Ysrc,Yp0},'reduceBrain',vx_vol,round(6/mean(vx_vol)),Yp0>1/3);
    Ysrcb = max(0,min(Ysrcb,max(T3th)*2));
    Yg   = cat_vol_grad(Ysrcb/T3th(3),vx_vol);
    Ydiv = cat_vol_div(Ysrcb/T3th(3),vx_vol);

    Yb   = smooth3(Yb)>0.5; 
    Yg   = cat_vol_resize(Yg ,'dereduceBrain',BB);
    Ydiv = cat_vol_resize(Ydiv ,'dereduceBrain',BB);
  elseif job.extopts.gcutstr==0
    % brain mask
    Yp0  = single(P(:,:,:,3))/255/3 + single(P(:,:,:,1))/255*2/3 + single(P(:,:,:,2))/255;
    Yb   = Yp0>0.5/3; 
    Ybb  = cat_vol_smooth3X(Yb,2); 
    
    [Ysrcb,Yp0,BB] = cat_vol_resize({Ysrc,Yp0},'reduceBrain',vx_vol,round(6/mean(vx_vol)),Yp0>1/3);
    Yg   = cat_vol_grad(Ysrcb/T3th(3),vx_vol);
    Ydiv = cat_vol_div(Ysrcb/T3th(3),vx_vol);
    Yg   = cat_vol_resize(Yg ,'dereduceBrain',BB);
    Ydiv = cat_vol_resize(Ydiv ,'dereduceBrain',BB);
  else
    % old skull-stripping
    brad = voli(sum(Yp0(:)>0.5).*prod(vx_vol)/1000); 
    Yp0  = single(P(:,:,:,3))/255/3 + single(P(:,:,:,1))/255*2/3 + single(P(:,:,:,2))/255;
    [Ysrcb,Yp0,BB] = cat_vol_resize({Ysrc,Yp0},'reduceBrain',vx_vol,round(6/mean(vx_vol)),Yp0>1/3);
    %Ysrcb = max(0,min(Ysrcb,max(T3th)*2));
    BGth = min(cat_stat_nanmean(Ysrc( P(:,:,:,6)>128 )),clsint(6));
    Yg   = cat_vol_grad((Ysrcb-BGth)/diff([BGth,T3th(3)]),vx_vol);
    Ydiv = cat_vol_div((Ysrcb-BGth)/diff([BGth,T3th(3)]),vx_vol);
    Ybo  = cat_vol_morph(cat_vol_morph(Yp0>0.3,'lc',2),'d',brad/2/mean(vx_vol)); 
    BVth = diff(T3th(1:2:3))/abs(T3th(3))*1.5; 
    RGth = diff(T3th(2:3))/abs(T3th(3))*0.1; 
    Yb   = single(cat_vol_morph((Yp0>1.9/3) | (Ybo & Ysrcb>mean(T3th(2)) & Ysrcb<T3th(3)*1.5 & Yg<0.5),'lo',max(0,0.6/mean(vx_vol)))); 
    %% region-growing GM 1
    Yb(~Yb & (~Ybo | Ysrcb<cat_stat_nanmean(T3th(2)) | Ysrcb>cat_stat_nanmean(T3th(3)*1.2) | Yg>BVth))=nan;
    [Yb1,YD] = cat_vol_downcut(Yb,Ysrcb/T3th(3),RGth); Yb(isnan(Yb))=0; Yb(YD<400/mean(vx_vol))=1; clear Yb1; 
    Yb(smooth3(Yb)<0.5)=0; Yb = single(Yb | (Ysrcb>T3th(1) & Ysrcb<1.2*T3th(3) & cat_vol_morph(Yb,'lc',4)));
    %% region-growing GM 2
    Yb(~Yb & (~Ybo | Ysrcb<T3th(1) | Ysrcb>cat_stat_nanmean(T3th(3)*1.2) | Yg>BVth))=nan;
    [Yb1,YD] = cat_vol_downcut(Yb,Ysrcb/T3th(3),RGth/2); Yb(isnan(Yb))=0; Yb(YD<400/mean(vx_vol))=1; clear Yb1; 
    Yb(smooth3(Yb)<0.5)=0; Yb = single(Yb | (Ysrcb>T3th(1) & Ysrcb<1.2*T3th(3) & cat_vol_morph(Yb,'lc',4)));
    %% region-growing GM 3
    Yb(~Yb & (~Ybo | Ysrcb<mean([BGth,T3th(1)]) | Ysrcb>cat_stat_nanmean(T3th(3)*1.2) | Yg>BVth))=nan;
    [Yb1,YD] = cat_vol_downcut(Yb,Ysrcb/T3th(3),RGth/10); Yb(isnan(Yb))=0; Yb(YD<400/mean(vx_vol))=1; clear Yb1; 
    Yb(smooth3(Yb)<0.5)=0; Yb(Yp0toC(Yp0*3,1)>0.9 & Yg<0.3 & Ysrcb>BGth & Ysrcb<T3th(2)) = 1; 
    %% ventrile closing
    [Ybr,Ymr,resT2] = cat_vol_resize({Yb>0,Ysrcb/T3th(3)},'reduceV',vx_vol,2,32); clear Ysrcb
    Ybr = Ybr | (Ymr<0.8 & cat_vol_morph(Ybr,'lc',6)); % large ventricle closing
    Ybr = cat_vol_morph(Ybr,'lc',2);                 % standard closing
    Yb  = Yb | cat_vol_resize(cat_vol_smooth3X(Ybr,2),'dereduceV',resT2)>0.7; 
    Yb  = smooth3(Yb)>0.5; 
    Ybb = cat_vol_smooth3X(Yb,2); 
    Yb   = cat_vol_resize(Yb ,'dereduceBrain',BB);
    Ybb  = cat_vol_resize(Ybb,'dereduceBrain',BB);
    Yg   = cat_vol_resize(Yg ,'dereduceBrain',BB);
    Ydiv = cat_vol_resize(Ydiv ,'dereduceBrain',BB);
    clear Ybo;
  end
	clear Ysrcb

  
  if ~(job.extopts.INV && any(sign(diff(T3th))==-1))
    %% Update probability maps
    % background vs. head - important for noisy backgrounds such as in MT weighting
    if job.extopts.gcutstr==0
      Ybg = ~Yb; 
    else
      if sum(sum(sum(P(:,:,:,6)>240 & Ysrc<cat_stat_nanmean(T3th(1:2)))))>10000
        Ybg = P(:,:,:,6); 
        [Ybgr,Ysrcr,resT2] = cat_vol_resize({Ybg,Ysrc},'reduceV',vx_vol,2,32); 
        Ybgrth = max(cat_stat_nanmean(Ysrcr(Ybgr(:)>128)) + 2*std(Ysrcr(Ybgr(:)>128)),T3th(1));
        Ybgr = cat_vol_morph(cat_vol_morph(cat_vol_morph(Ybgr>128,'d') & Ysrcr<Ybgrth,'lo',1),'lc',1);
        Ybg  = cat_vol_resize(cat_vol_smooth3X(Ybgr,1),'dereduceV',resT2); 
        clear Ysrcr Ybgr; 
      else
        Ybg = ~Yb;
      end
    end
    P4   = cat_vol_ctype( single(P(:,:,:,6)) .* (Ysrc<T3th(2))  .* (Ybg==0) + single(P(:,:,:,4)) .* (Ybg<1) ); % remove air in head
    P5   = cat_vol_ctype( single(P(:,:,:,6)) .* (Ysrc>=T3th(2)) .* (Ybg==0) + single(P(:,:,:,5)) .* (Ybg<1) ); % remove air in head
    P6   = cat_vol_ctype( single(sum(P(:,:,:,4:5),4)) .* (Ybg==1) + single(P(:,:,:,6)) .* (Ybg>0) ); % add objects/artifacts to background
    P(:,:,:,4) = P4;
    P(:,:,:,5) = P5;
    P(:,:,:,6) = P6;
    clear P4 P5 P6;

    %% correct probability maps to 100% 
    sumP = cat_vol_ctype(255 - sum(P(:,:,:,1:6),4));
    P(:,:,:,1) = P(:,:,:,1) + sumP .* uint8( Ybg<0.5  &  Yb & Ysrc>cat_stat_nanmean(T3th(1:2)) & Ysrc<cat_stat_nanmean(T3th(2:3)));
    P(:,:,:,2) = P(:,:,:,2) + sumP .* uint8( Ybg<0.5  &  Yb & Ysrc>=cat_stat_nanmean(T3th(2:3)));
    P(:,:,:,3) = P(:,:,:,3) + sumP .* uint8( Ybg<0.5  &  Yb & Ysrc<=cat_stat_nanmean(T3th(1:2)));
    P(:,:,:,4) = P(:,:,:,4) + sumP .* uint8( Ybg<0.5  & ~Yb & Ysrc<T3th(2));
    P(:,:,:,5) = P(:,:,:,5) + sumP .* uint8( Ybg<0.5  & ~Yb & Ysrc>=T3th(2));
    P(:,:,:,6) = P(:,:,:,6) + sumP .* uint8( Ybg>=0.5 & ~Yb );
    clear Ybg;

    %% head to WM 
    % Undercorrection of strong inhomogeneities in high field scans 
    % (>1.5T) can cause missalignments of the template and therefore 
    % miss classifications of the tissues that finally avoid further 
    % corrections in by LAS. 
    % Typically the alginment failed in this cases because the high 
    % intensities next to the head that were counted as head and not
    % corrected by SPM.
    % e.g. HR075, Magdeburg7T, SRS_SRS_Jena_DaRo81_T1_20150320-191509_MPR-08mm-G2-bw330-nbc.nii, ...
    Ywm = single(P(:,:,:,2)>128 & Yg<0.3 & Ydiv<0.03); Ywm(Ybb<0.5 | (P(:,:,:,1)>128 & abs(Ysrc/T3th(3)-2/3)<1/3) | Ydiv>0.03) = nan;
    [Ywm1,YD] = cat_vol_downcut(Ywm,1-Ysrc/T3th(3),0.02); Yb(isnan(Yb))=0; Ywm(YD<300)=1; Ywm(isnan(Ywm))=0; clear Ywm1 YD;
    Ywmc = uint8(smooth3(Ywm)>0.7);
    Ygmc = uint8(cat_vol_morph(Ywmc,'d',2) & ~Ywmc & Ydiv>0 & Yb & cat_vol_smooth3X(Yb,8)<0.9);
    P(:,:,:,[1,3:6]) = P(:,:,:,[1,3:6]) .* repmat(1-Ywmc,[1,1,1,5]);
    P(:,:,:,2:6)     = P(:,:,:,2:6)     .* repmat(1-Ygmc,[1,1,1,5]);
    P(:,:,:,1)       = max(P(:,:,:,1),255*Ygmc);
    P(:,:,:,2)       = max(P(:,:,:,2),255*Ywmc);
    Yp0  = single(P(:,:,:,3))/255/3 + single(P(:,:,:,1))/255*2/3 + single(P(:,:,:,2))/255;
    clear Ygmc Ywmc Yg Ydiv;

    % head to GM
    Ygm = uint8(cat_vol_morph(Ywm>0.5,'d',2) & Ywm<0.9 & (Ysrc>cat_stat_nanmean(T3th(2:3))) & Yp0<2/3);
    P(:,:,:,5) = P(:,:,:,5) .* (1-Ygm);
    P(:,:,:,3) = P(:,:,:,3) .* (1-Ygm);
    P(:,:,:,2) = P(:,:,:,2) .* (1-Ygm);
    P(:,:,:,1) = cat_vol_ctype(single(P(:,:,:,1)) + 255*single(Ygm));
    Yp0  = single(P(:,:,:,3))/255/3 + single(P(:,:,:,1))/255*2/3 + single(P(:,:,:,2))/255;
    clear Ywm Ygm;

    %% remove brain tissues outside the brainmask ...
    % tissues > skull (within the brainmask)
    Yhdc = uint8(smooth3( Ysrc/T3th(3).*(Ybb>0.2) - Yp0 )>0.5); 
    sumP = sum(P(:,:,:,1:3),4); 
    P(:,:,:,4)   =  cat_vol_ctype( single(P(:,:,:,4)) + sumP .* ((Ybb<=0.05) | Yhdc ) .* (Ysrc<T3th(2)));
    P(:,:,:,5)   =  cat_vol_ctype( single(P(:,:,:,5)) + sumP .* ((Ybb<=0.05) | Yhdc ) .* (Ysrc>=T3th(2)));
    P(:,:,:,1:3) =  P(:,:,:,1:3) .* repmat(uint8(~(Ybb<=0.05) | Yhdc ),[1,1,1,3]);
    clear sumP Yp0 Yhdc; 
  end
  clear Ybb;
  

  %% MRF
  % Used spm_mrf help and tested the probability TPM map for Q without good results.         
  nmrf_its = 0; % 10 interations better to get full probability in thin GM areas 
  spm_progress_bar('init',nmrf_its,['MRF: Working on ' nam],'Iterations completed');
  G   = ones([Kb,1],'single');
  vx2 = single(sum(VT.mat(1:3,1:3).^2));
  % P = zeros([d(1:3),Kb],'uint8');
  % P = spm_mrf(P,Q,G,vx2); % init: transfer data from Q to P 
  if 0
    %% use TPM as Q
    Q = zeros(size(P),'uint8');
    for di=1:6
      vol = cat_vol_ctype(spm_sample_vol(tpm.V(di),...
        double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),0)*255,'uint8');
      Q(:,:,:,di) = reshape(vol,d);
    end
  end
  for iter=1:nmrf_its,
      P = spm_mrf(P,single(P),G,vx2); % spm_mrf(P,Q,G,vx2);
      spm_progress_bar('set',iter);
  end
  % cleanup
  %P = clean_gwc(P,1);
  spm_progress_bar('clear');
  for k1=1:size(P,4)
      Ycls{k1} = P(:,:,:,k1);
  end
  clear Q P q q1 Coef b cr s t1 t2 t3 N lkp n wp M k1


  if job.extopts.verb>2
    % save information for debugging and OS test
    % input variables + bias corrected, bias field, class image
    % strong differences in bias fields can be the result of different 
    % registration > check 'res.image.mat' and 'res.Affine'
    [pth,nam] = spm_fileparts(res.image0(1).fname); 
    tpmci  = 1;
    tmpmat = fullfile(pth,reportfolder,sprintf('%s_%s%02d%s.mat',nam,'write',tpmci,'postbias'));
    save(tmpmat,'res','tpm','job','Ysrc','Ybf','Ycls');
  end

  clear Ybf
  fprintf('%5.0fs\n',etime(clock,stime));



  %% ---------------------------------------------------------------------
  %  check if the brain hits a image boundary
  %  ---------------------------------------------------------------------

  %Yib1 = cat_vol_morph(true(size(Yb)),'e',1);
  %Yib3 = cat_vol_morph(true(size(Yb)),'e',2);
  %Ybbs = sum(Yb & Yib);



  %% ---------------------------------------------------------------------
  %  Global (and local) intensity normalization and partioning 
  %  ---------------------------------------------------------------------
  %  Global and local intensity corrections are the basis of most of the 
  %  following functions. The global normalization based on the SPM tissue
  %  thresholds (res.mn) and were used anyway. For strong differences 
  %  (mostly by the CSF) the median will used, because it is e.g. more 
  %  stable. This will cause a warning by the cat_main_gintnorm.
  %
  %  The local adaptive segmentation include a further bias correction  
  %  and a global  and local intensity correction. The local intensity 
  %  correction refines the tissue maps to aproximate the local tissue 
  %  peaks of WM (maximum-based), GM, and CSF. 
  %
  %  If you want to see intermediate steps of the processing use the "ds"
  %  function:
  %    ds('l2','',vx_vol,Ym,Yb,Ym,Yp0,80)
  %  that display 4 images (unterlay, overlay, image1, image2) for one 
  %  slice. The images were scaled in a range of 0 to 1. The overlay 
  %  allows up to 20 colors
  %  
  %  ---------------------------------------------------------------------
  stime = cat_io_cmd('Global intensity correction');
  if debug
    Ym   = Ysrc / T3th(3); %#ok<NASGU> % only WM scaling
    Yp0  = (single(Ycls{1})/255*2 + single(Ycls{2})/255*3 + single(Ycls{3})/255)/3; %#ok<NASGU> % label map
  end
  if any( min(vx_vol*2,1.4)./vx_vol >= 2 )
    %%
    Ysrcr = cat_vol_resize(Ysrc,'reduceV',vx_vol,min(vx_vol*2,1.4),32,'meanm');
    Ybr   = cat_vol_resize(single(Yb),'reduceV',vx_vol,min(vx_vol*2,1.4),32,'meanm')>0.5;
    Yclsr = cell(size(Ycls)); for i=1:6, Yclsr{i} = cat_vol_resize(Ycls{i},'reduceV',vx_vol,min(vx_vol*2,1.4),32); end
    [Ymr,Ybr,T3th,Tth,job.inv_weighting,noise,cat_warnings] = cat_main_gintnorm(Ysrcr,Yclsr,Ybr,vx_vol,res,Yy,job.extopts);
    clear Ymr Ybr Ysrcr Yclsr; 
    Ym = cat_main_gintnorm(Ysrc,Tth); 
  else
    [Ym,Yb,T3th,Tth,job.inv_weighting,noise,cat_warnings] = cat_main_gintnorm(Ysrc,Ycls,Yb,vx_vol,res,Yy,job.extopts);;
  end

  % update in inverse case ... required for LAS
  if job.inv_weighting
%    Ysrc = Ym * Tth.T3th(5); Tth.T3th = Tth.T3thx * Tth.T3th(5);
    if T3th(1)>T3th(3) && T3th(2)<T3th(3) && T3th(1)>T3th(2)
      Yp0  = single(Ycls{3})/255/3 + single(Ycls{1})/255*2/3 + single(Ycls{2})/255;
      %ds('l2','',vx_vol,Ym,Yb,Ym,Yp0,90)
      Yb2  = cat_vol_morph(Yp0>0.5,'lc',2); 
      prob = cat(4,cat_vol_ctype(Yb2.*Yp0toC(Ym*3,2)*255),...
                   cat_vol_ctype(Yb2.*Yp0toC(Ym*3,3)*255),...
                   cat_vol_ctype(Yb2.*Yp0toC(min(3,Ym*3),1)*255)); 
      prob = clean_gwc(prob);
      for ci=1:3, Ycls{ci} = prob(:,:,:,ci); end; 
      %job.extopts.mrf = 0.3; 
      clear prob;  
    end
 
    Ysrc = Ym; Tth.T3thx(3:5) = 1/3:1/3:1; Tth.T3th = Tth.T3thx; T3th = 1/3:1/3:1;
    
  end
  fprintf('%5.0fs\n',etime(clock,stime));

  if job.extopts.verb>2
    tpmci  = tpmci + 1;
    tmpmat = fullfile(pth,reportfolder,sprintf('%s_%s%02d%s.mat',nam,'write',tpmci,'postgintnorm'));
    save(tmpmat,'Ysrc','Ycls','Ym','Yb','T3th','vx_vol');
  end


  %% After the intensity scaling and with correct information about the
  % variance of the tissue, a further harder noise correction is meaningful.
  % Finally, a stronger NLM-filter is better than a strong MRF filter!

  if job.extopts.NCstr~=0  

    % use a boundary box (BB) of the brain to speed up denoising
    useBB = 1; %job.extopts.output.bias.native; 

    if useBB
      [Yms,Ybr,BB] = cat_vol_resize({Ym,Yb},'reduceBrain',vx_vol,round(2/cat_stat_nanmean(vx_vol)),Yb); Ybr = Ybr>0.5; 


      % apply NLM filter
      if job.extopts.NCstr==2 || job.extopts.NCstr==3
        cat_io_cmd(sprintf('ISARNLM noise correction (NCstr=%d)',job.extopts.NCstr));
        if job.extopts.verb>1, fprintf('\n'); end
        if job.extopts.NCstr==2, NCstr=-inf; else NCstr=1; end 
        Yms = cat_vol_isarnlm(Yms,res.image,job.extopts.verb>1,NCstr); 

        Ym(BB.BB(1):BB.BB(2),BB.BB(3):BB.BB(4),BB.BB(5):BB.BB(6)) = ...
          Ym(BB.BB(1):BB.BB(2),BB.BB(3):BB.BB(4),BB.BB(5):BB.BB(6)) .* (1-Ybr) + ...
          Yms .* Ybr;
      else
        %%
        Ymo = Yms + 0;
        if isinf(job.extopts.NCstr) 
          stime = cat_io_cmd(sprintf('SANLM noise correction'));
        else
          stime = cat_io_cmd(sprintf('SANLM noise correction (NCstr=%0.2f)',job.extopts.NCstr));
        end

        % filter
        cat_sanlm(Yms,3,1,0);

        % merging
        if isinf(job.extopts.NCstr) || sign(job.extopts.NCstr)==-1
          job.extopts.NCstr = min(1,max(0,cat_stat_nanmean(abs(Yms(Ybr(:)) - Ymo(Ybr(:)))) * 15 * min(1,max(0,abs(job.extopts.NCstr))) )); 
          NC     = min(2,abs(Yms - Ymo) ./ max(eps,Yms) * 15 * 2 * min(1,max(0,abs(job.extopts.NCstr)))); 
          NCs    = NC + 0; spm_smooth(NCs,NCs,2); NCs = NCs .* cat_stat_nanmean(NCs(Ybr(:))) / cat_stat_nanmean(NC(Ybr(:))); clear NC;
          NCs  = max(0,min(1,NCs)); 
        end
        fprintf('%5.0fs\n',etime(clock,stime));  

        % mix original and noise corrected image and go back to original resolution
        if exist('NCs','var');
          Ym(BB.BB(1):BB.BB(2),BB.BB(3):BB.BB(4),BB.BB(5):BB.BB(6)) = ...
            Ym(BB.BB(1):BB.BB(2),BB.BB(3):BB.BB(4),BB.BB(5):BB.BB(6)) .* (1-Ybr) + ...
            (1-NCs) .* Ym(BB.BB(1):BB.BB(2),BB.BB(3):BB.BB(4),BB.BB(5):BB.BB(6)) .* Ybr + ...
            NCs .* Yms .* Ybr;
        else
          Ym(BB.BB(1):BB.BB(2),BB.BB(3):BB.BB(4),BB.BB(5):BB.BB(6)) = ...
            Ym(BB.BB(1):BB.BB(2),BB.BB(3):BB.BB(4),BB.BB(5):BB.BB(6)) .* (1-Ybr) + ...
            (1-job.extopts.NCstr) .* Ym(BB.BB(1):BB.BB(2),BB.BB(3):BB.BB(4),BB.BB(5):BB.BB(6)) .* Ybr + ...
            job.extopts.NCstr .* Yms .* Ybr;
        end    
        clear NCs; 

      end
      clear Yms Ybr BB;
    else
      if job.extopts.NCstr==2 || job.extopts.NCstr==3
        cat_io_cmd(sprintf('ISARNLM noise correction (NCstr=%0.2f)',job.extopts.NCstr));
        if job.extopts.verb>1, fprintf('\n'); end
        if job.extopts.NCstr==2, NCstr=-inf; else NCstr=1; end 
        Ym = cat_vol_isarnlm(Ym,res.image,job.extopts.verb>1,NCstr); 
      else
        Ymo = Ym + 0;

        if isinf(job.extopts.NCstr) 
          stime = cat_io_cmd(sprintf('SANLM noise correction'));
        else
          stime = cat_io_cmd(sprintf('SANLM noise correction (NCstr=%0.2f)',job.extopts.NCstr));
        end

        % filter
        cat_sanlm(Ym,3,1,0);

        % merging
        if isinf(job.extopts.NCstr) || sign(job.extopts.NCstr)==-1
          NCstr2 = min(1,max(0,cat_stat_nanmean(abs(Ym(Ybr(:)) - Ymo(Ybr(:)))) * 15 * min(1,max(0,abs(job.extopts.NCstr))) )); 
          NC     = min(2,abs(Ym - Ymo) ./ max(eps,Ym) * 15 * 2 * min(1,max(0,abs(NCstr2)))); 
          NCs    = NC + 0; spm_smooth(NCs,NCs,2); NCs = NCs .* cat_stat_nanmean(NCs(Ybr(:))) / cat_stat_nanmean(NC(Ybr(:)));
          NCs  = max(0,min(1,NCs));
        end
        fprintf('%5.0fs\n',etime(clock,stime));  

        % mix original and noise corrected image and go back to original resolution
        if exist('NCs','var');
          Ym = (1-NCs) .* Ymo  +  NCs .* Ym ;
        else
          Ym = (1-NCstr2) .* Ymo + NCstr2 .* Ym;
        end    
      end
    end

    if job.inv_weighting
      Ysrc = Ym;
    else
      Ysrc = cat_main_gintnormi(Ym,Tth);
    end

  end
  
  
  
  if 0
    % prepared for improved partitioning - RD20170320
    job2=job; job2.extopts.regstr=eps; res2=res; res2.do_dartel=2; 
    [trans,res.ppe.reginitp] = cat_main_registration(job2,res2,Ycls,Yy,tpm.M);
    Yy = trans.warped.y; clear trans job2 res2; 
  end
  
  
  


  %% Local Intensity Correction 
  Ymo = Ym;
  if job.extopts.LASstr>0
    if job.extopts.LASstr>1 
      extoptsLAS2 = job.extopts;
      extoptsLAS2.LASstr = extoptsLAS2.LASstr-1; 
      stime = cat_io_cmd(sprintf('Local adaptive segmentation 2 (LASstr=%0.2f)',extoptsLAS2.LASstr));
      [Ymi,Ym,Yclsi] = cat_main_LASs(Ysrc,Ycls,Ym,Yb,Yy,Tth,res,vx_vol,extoptsLAS2); % use Yclsi after cat_vol_partvol
    else
      stime = cat_io_cmd(sprintf('Local adaptive segmentation (LASstr=%0.2f)',job.extopts.LASstr)); 
      [Ymi,Ym] = cat_main_LAS(Ysrc,Ycls,Ym,Yb,Yy,T3th,res,vx_vol,job.extopts,Tth); 
    end

    %Ymioc = Ymi+0; 
    if job.extopts.NCstr>0 
      %% noise correction of the local normalized image Ymi, whereas only small changes are expected in Ym by the WM bias correction
      [Ymis,Ymior,BB]  = cat_vol_resize({Ymi,Ymo},'reduceBrain',vx_vol,round(2/mean(vx_vol)),Yb);
      Yc = abs(Ymis - Ymior); Yc = Yc * 6 * min(2,max(0,abs(job.extopts.NCstr))); 
      spm_smooth(Yc,Yc,2./vx_vol); Yc = max(0,min(1,Yc)); clear Ymior; 

      if job.extopts.NCstr==2 || job.extopts.NCstr==3
        if job.extopts.verb
          cat_io_cmd(sprintf('  ISARNLM noise correction for LAS'),'g5');
        end
        if job.extopts.verb>1, fprintf('\n'); end
        if job.extopts.NCstr==2, NCstr=-inf; else NCstr=1; end 
        Ymis = cat_vol_isarnlm(Ymis,res.image,(job.extopts.verb>1)*2,NCstr);
      else
        if job.extopts.verb>1
          stime2 = cat_io_cmd(sprintf('  SANLM noise correction for LAS'),'g5'); 
        end
        cat_sanlm(Ymis,3,1,0);
        if job.extopts.verb>1
          cat_io_cmd(' ','','',job.extopts.verb,stime2); 
        end 
        fprintf('%5.0fs\n',etime(clock,stime));
      end
      % mix original and noise corrected image and go back to original resolution
      Ybr = Yb(BB.BB(1):BB.BB(2),BB.BB(3):BB.BB(4),BB.BB(5):BB.BB(6));
      Ymi(BB.BB(1):BB.BB(2),BB.BB(3):BB.BB(4),BB.BB(5):BB.BB(6)) = ...
        Ymi(BB.BB(1):BB.BB(2),BB.BB(3):BB.BB(4),BB.BB(5):BB.BB(6)) .* (1-Ybr) + ...
        (1-Yc) .* Ymi(BB.BB(1):BB.BB(2),BB.BB(3):BB.BB(4),BB.BB(5):BB.BB(6)) .* Ybr + ...
        Yc .* Ymis .* Ybr;

      % extreme background denoising to remove wholes?
      Ymis = cat_vol_median3(Ymi,Ymi>0 & Ymi<0.4,Ymi<0.4); Ymi = Ymi.*max(0.1,Ymi>0.4) + Ymis.*min(0.9,Ymi<=0.4);
      Ymis = cat_vol_median3(Ym,Ym>0 & Ym<0.4,Ym<0.4); Ym = Ym.*max(0.1,Ym>0.4) + Ymis.*min(0.9,Ym<=0.4);
      
      clear Ymis;
    end
    %clear Ymioc; 
    
  else
    Ymi = Ym; 
  end
  if ~debug; clear Ysrc ; end


  if job.extopts.verb>2
    tpmci=tpmci+1; tmpmat = fullfile(pth,reportfolder,sprintf('%s_%s%02d%s.mat',nam,'write',tpmci,'postLAS'));
    save(tmpmat,'Ysrc','Ycls','Ymi','Yb','T3th','vx_vol');
  end
  
  %  ---------------------------------------------------------------------
  %  Partitioning: 
  %  --------------------------------------------------------------------- 
  %  For most of the following adaptions further knowledge of special 
  %  regions is helpfull. Also Ymi is maybe still a little bit inhomogen 
  %  the alignment should work. Only strong inhomogenities can cause 
  %  problems, especially for the blood vessel detection. 
  %  But for bias correction the ROIs are important too, to avoid over
  %  corrections in special regions like the cerbellum and subcortex. 
  %  ---------------------------------------------------------------------
  NS = @(Ys,s) Ys==s | Ys==s+1; 
  stime = cat_io_cmd('ROI segmentation (partitioning)');
  [Yl1,Ycls,YBG,YMF] = cat_vol_partvol(Ymi,Ycls,Yb,Yy,vx_vol,job.extopts,tpm.V,noise);
  if exist('Yclsi','var'), Ycls = Yclsi; clear Yclsi; end % new Ycls from LAS
  fprintf('%5.0fs\n',etime(clock,stime));

  if ~debug; clear YBG Ycr; end


  %  ---------------------------------------------------------------------
  %  Blood Vessel Correction 
  %  ---------------------------------------------------------------------
  %  Blood vessel correction has to be done before the segmentation to 
  %  remove high frequency strutures and avoid missclassifications.
  %  Problems can occure for strong biased images, because the partioning 
  %  has to be done before bias correction.
  %  Of course we only want to do this for highres T1 data!
  %  ---------------------------------------------------------------------
  if job.extopts.BVCstr && ~job.inv_weighting && all(vx_vol<2); 
    stime = cat_io_cmd(sprintf('Blood vessel correction (BVCstr=%0.2f)',job.extopts.BVCstr));

    Ybv  = cat_vol_smooth3X(cat_vol_smooth3X( ...
      NS(Yl1,7) .* (Ymi*3 - (1.5-job.extopts.BVCstr)),0.3).^4,0.1)/3;

    % correct src images
    Ymi   = max(0,Ymi - Ybv*2/3); 
    Ymi   = cat_vol_median3(Ymi,cat_vol_morph(Ybv>0.5,'dilate')); 
    Ymis  = cat_vol_smooth3X(Ymi); Ymi(Ybv>0.5) = Ymis(Ybv>0.5); clear Ymis;

    % update classes
    Ycls{1} = min(Ycls{1},cat_vol_ctype(255 - Ybv*127)); 
    Ycls{2} = min(Ycls{2},cat_vol_ctype(255 - Ybv*127)); 
    Ycls{3} = max(Ycls{3},cat_vol_ctype(127*Ybv)); 

    fprintf('%5.0fs\n',etime(clock,stime));
    clear Ybv p0; 
  end





  %% ---------------------------------------------------------------------
  %  Segmentation part
  %  ---------------------------------------------------------------------
  %  Now, it is time for skull-stripping (gcut,morph), AMAP tissue 
  %  segmentation, and further tissue corrections (cleanup,LAS,finalmask).
  %  ---------------------------------------------------------------------


  %  -------------------------------------------------------------------
  %  skull-stipping
  %  -------------------------------------------------------------------
  %  For skull-stripping gcut is used in general, but a simple and very 
  %  old function is still available as backup solution.
  %  Futhermore, both parts prepare the initial segmentation map for the 
  %  AMAP function.
  %  -------------------------------------------------------------------
  if job.extopts.gcutstr>0
    %  -----------------------------------------------------------------
    %  gcut+: skull-stripping using graph-cut
    %  -----------------------------------------------------------------
    try 
      stime = cat_io_cmd(sprintf('Skull-stripping using graph-cut (gcutstr=%0.2f)',job.extopts.gcutstr));
      [Yb,Yl1] = cat_main_gcut(Ymo,Yb,Ycls,Yl1,YMF,vx_vol,job.extopts);
      fprintf('%5.0fs\n',etime(clock,stime));
      if 0
        %% just for manual debuging / development of gcut and gcutstr > remove this in 201709?
        job.extopts.gcutstr=0.5; [Yb05,Yl105] = cat_main_gcut(Ym,Yb,Ycls,Yl1,YMF,vx_vol,job.extopts); 
        job.extopts.gcutstr=0.1; [Yb01,Yl101] = cat_main_gcut(Ym,Yb,Ycls,Yl1,YMF,vx_vol,job.extopts); 
        job.extopts.gcutstr=0.9; [Yb09,Yl109] = cat_main_gcut(Ym,Yb,Ycls,Yl1,YMF,vx_vol,job.extopts);
        ds('d2','',vx_vol,(Yb01 + Yb05+Yb09)/3,Ymi.*(0.2+0.8*Yb01),Ymi.*(0.2+0.8*Yb05),Ymi.*(0.2+0.8*Yb09),50)
      end
    catch %#ok<CTCH>
      fprintf('\n'); cat_warnings = cat_io_addwarning(cat_warnings,'CAT:cat_main_gcut:err99','Unknown error in cat_main_gcut. Use old brainmask.'); fprintf('\n');
      job.extopts.gcutstr = 99;
    end
  end
  %if ~debug; clear Ymo; end




  %% -------------------------------------------------------------------
  %  AMAP segmentation
  %  -------------------------------------------------------------------
  %  Most corrections were done before and the AMAP routine is used with 
  %  a low level of iterations and no further bias correction, because
  %  some images get tile artifacts. 
  %  -------------------------------------------------------------------

  % correct for harder brain mask to avoid meninges in the segmentation
  Ymib = Ymi; Ymib(~Yb) = 0; 
  rf = 10^4; Ymib = round(Ymib*rf)/rf;

  %  prepare data for segmentation
  if 1
    %% classic approach, consider the WMH!
    Kb2 = 4;
    cls2 = zeros([d(1:2) Kb2]);
    Yp0  = zeros(d,'uint8');
    for i=1:d(3)
        for k1 = 1:Kb2, cls2(:,:,k1) = Ycls{k1}(:,:,i); end
        % find maximum for reordered segmentations
        [maxi,maxind] = max(cls2(:,:,[3,1,2,4:Kb2]),[],3);
        k1ind = [1 2 3 1 0 0]; 
        for k1 = 1:Kb2
          Yp0(:,:,i) = Yp0(:,:,i) + cat_vol_ctype((maxind == k1) .* (maxi~=0) * k1ind(k1) .* Yb(:,:,i)); 
        end
    end
    %clear maxi maxind Kb k1 cls2;
  else
    % more direct method ... a little bit more WM, less CSF
    % Yp0 = uint8(max(Yb,min(3,round(Ymi*3)))); Yp0(~Yb) = 0;
  end  


  % use index to speed up and save memory
  sz = size(Yb);
  [indx, indy, indz] = ind2sub(sz,find(Yb>0));
  indx = max((min(indx) - 1),1):min((max(indx) + 1),sz(1));
  indy = max((min(indy) - 1),1):min((max(indy) + 1),sz(2));
  indz = max((min(indz) - 1),1):min((max(indz) + 1),sz(3));

  % Yb source image because Amap needs a skull stripped image
  % set Yp0b and source inside outside Yb to 0
  Yp0b = Yp0(indx,indy,indz); %clear Yp0
  Ymib = Ymib(indx,indy,indz); 


  % remove non-brain tissue with a smooth mask and set values inside the
  % brain at least to CSF to avoid wholes for images with CSF==BG.
  if job.extopts.LASstr>0 
    Ywmstd = cat_vol_localstat(single(Ymib),Yp0b==3,1,4); 
    CSFnoise(1) = cat_stat_nanmean(Ywmstd(Ywmstd(:)>0))/mean(vx_vol); 
    Ywmstd = cat_vol_localstat(cat_vol_resize(single(Ymib),'reduceV',vx_vol,vx_vol*2,16,'meanm'),...
      cat_vol_resize(Yp0==3,'reduceV',vx_vol,vx_vol*2,16,'meanm')>0.5,1,4); 
    CSFnoise(2) = cat_stat_nanmean(Ywmstd(Ywmstd(:)>0))/mean(vx_vol); 
    Ycsf = double(0.33 * Yb(indx,indy,indz)); spm_smooth(Ycsf,Ycsf,0.6*vx_vol);
    Ycsf = Ycsf + cat_vol_smooth3X(randn(size(Ycsf)),0.5) * max(0.005,min(0.2,CSFnoise(1)/4)); % high-frequency noise
    Ycsf = Ycsf + cat_vol_smooth3X(randn(size(Ycsf)),1.0) * max(0.005,min(0.2,CSFnoise(2)*1)); % high-frequency noise
    Ymib = max(Ycsf*0.8 .* cat_vol_smooth3X(Ycsf>0,2),Ymib); 
    clear Ycsf; 
    % Yb is needed for surface reconstruction
  %     if ~job.output.surface, clear Yb; end
  end
  
  % Amap parameters  - default sub=16 caused errors with highres data!
  Ymib = double(Ymib); n_iters = 16; sub = round(32/min(vx_vol)); n_classes = 3; pve = 5; bias_fwhm = 60; init_kmeans = 0;  %#ok<NASGU>
  if job.extopts.mrf~=0, iters_icm = 50; else iters_icm = 0; end %#ok<NASGU>


  % adaptive mrf noise 
  if job.extopts.mrf>=1 || job.extopts.mrf<0; 
    % estimate noise
    [Yw,Yg] = cat_vol_resize({Ymi.*(Ycls{1}>240),Ymi.*(Ycls{2}>240)},'reduceV',vx_vol,3,32,'meanm');
    Yn = max(cat(4,cat_vol_localstat(Yw,Yw>0,2,4),cat_vol_localstat(Yg,Yg>0,2,4)),[],4);
    job.extopts.mrf = double(min(0.15,3*cat_stat_nanmean(Yn(Yn(:)>0)))) * 0.5; 
    %clear Yn Ycls1 Ycls2 Yg;
  end

  % display something
  stime = cat_io_cmd(sprintf('Amap using initial SPM12 segmentations (MRF filter strength %0.2f)',job.extopts.mrf));       

  % do segmentation  
  amapres = evalc(['prob = cat_amap(Ymib, Yp0b, n_classes, n_iters, sub, pve, init_kmeans, ' ...
    'job.extopts.mrf, vx_vol, iters_icm, bias_fwhm);']);
  fprintf('%5.0fs\n',etime(clock,stime));
  
  % analyse segmentation ... the input Ym is normalized an the tissue peaks should be around [1/3 2/3 3/3]
  amapres = textscan(amapres,'%s'); amapres = amapres{1}; 
  th{1}   = cell2mat(textscan(amapres{11},'%f*%f')); 
  th{2}   = cell2mat(textscan(amapres{12},'%f*%f')); 
  th{3}   = cell2mat(textscan(amapres{13},'%f*%f')); 
  
  
  if job.extopts.verb>1 
    fprintf('    AMAP peaks: [CSF,GM,WM] = [%0.2f%s%0.2f,%0.2f%s%0.2f,%0.2f%s%0.2f]\n',...
      th{1}(1),char(177),th{1}(2),th{2}(1),char(177),th{2}(2),th{3}(1),char(177),th{3}(2));
  end
  if th{1}(1)<0 || th{1}(1)>0.6 || th{2}(1)<0.5 || th{2}(1)>0.9 || th{3}(1)<0.95-th{3}(2) || th{3}(1)>1.1
    error('cat_main:amap',['AMAP estimated untypical tissue peaks that point to an \n' ...
                           'error in the preprocessing bevor the AMAP segmentation. ']);
  end
  % reorder probability maps according to spm order
  clear Yp0b Ymib; 
  prob = prob(:,:,:,[2 3 1]); 
  clear vol %Ymib
  %fprintf(sprintf('%s',repmat('\b',1,94+4)));
  



  %% -------------------------------------------------------------------
  %  final cleanup
  %  There is one major parameter to controll the strength of the cleanup.
  %  As far as the cleanup has a strong relation to the skull-stripping, 
  %  cleanupstr is controlled by the gcutstr. 
  %     Yp0ox = single(prob(:,:,:,1))/255*2 + single(prob(:,:,:,2))/255*3 + single(prob(:,:,:,3))/255; Yp0o = zeros(d,'single'); Yp0o(indx,indy,indz) = Yp0ox; 
  %     Yp0   = zeros(d,'uint8'); Yp0(indx,indy,indz) = Yp0b; 
  %  -------------------------------------------------------------------
  if job.extopts.cleanupstr>0 && job.extopts.cleanupstr<=1 
    %%
    prob = clean_gwc(prob,round(job.extopts.cleanupstr*2)); % old cleanup
    
    %for i=1:3
    %   Ycls{i}(:) = 0; Ycls{i}(indx,indy,indz) = prob(:,:,:,i);
    %end
    %Yp0b = Yb(indx,indy,indz); 
    
    [Ycls,Yp0b] = cat_main_cleanup(Ycls,prob,Yl1(indx,indy,indz),Ymo(indx,indy,indz),job.extopts,job.inv_weighting,vx_volr,indx,indy,indz);
  else
    if 0% job.extopts.cleanupstr == 2 % old cleanup for tests
      stime = cat_io_cmd('Old cleanup');
      prob = clean_gwc(prob,1);
      fprintf('%5.0fs\n',etime(clock,stime));
    end
    for i=1:3
       Ycls{i}(:) = 0; Ycls{i}(indx,indy,indz) = prob(:,:,:,i);
    end
    Yp0b = Yb(indx,indy,indz); 
  end;
  if ~debug; clear Ymo; end
  clear prob




  %% -------------------------------------------------------------------
  %  Correction of WM hyperintensities
  %  -------------------------------------------------------------------
  %  The correction of WMH should be important for a correct normalization.
  %  It is only important to close the mayor WMH structures, and further
  %  closing can lead to problems with small gyri. So keep it simple here 
  %  and maybe add further refinements in the partitioning function.
  %  -------------------------------------------------------------------
  LAB  = job.extopts.LAB;
  vxv  = 1/max(vx_vol);

  Yp0 = zeros(d,'uint8'); Yp0(indx,indy,indz) = Yp0b; 
  Ywmhrel = NS(Yl1,23);
  qa.subjectmeasures.WMH_abs    = sum(Ywmhrel(:));                                            % absolute WMH volume without PVE
  qa.subjectmeasures.WMH_rel    = 100*qa.subjectmeasures.WMH_abs / sum(Yp0(:)>(0.5/3*255));   % relative WMH volume to TIV without PVE
  qa.subjectmeasures.WMH_WM_rel = 100*qa.subjectmeasures.WMH_abs / sum(Yp0(:)>(2.5/3*255));   % relative WMH volume to WM without PVE
  qa.subjectmeasures.WMH_abs    = prod(vx_vol)/1000 * qa.subjectmeasures.WMH_abs;             % absolute WMH volume without PVE in cm^3
  clear Ywmhrel Yp0

  Yp0b = cat_vol_ctype(single(Ycls{1})*2/3 + single(Ycls{2}) + single(Ycls{3})*1/3,'uint8');
  Yp0b = Yp0b(indx,indy,indz); 


  % correction for normalization [and final segmentation]
  if job.extopts.WMHC && job.extopts.WMHCstr>0 && ~job.inv_weighting; 

    if job.extopts.WMHC==1
      stime = cat_io_cmd(sprintf('Internal WMH correction for spatial normalization (WMHCstr=%0.2f)',job.extopts.WMHCstr));
    elseif job.extopts.WMHC>1
      stime = cat_io_cmd(sprintf('Permanent WMH correction (WMHCstr=%0.2f)',job.extopts.WMHCstr));
    end

    % setting of furter WMHC that can now be detected by the further
    % evalution of the segmentation. 
    % estimation of WMHC is important for LAS (do I use it?)
    %  ... code ...

    % prepare correction map
    Ywmh = cat_vol_morph(NS(Yl1,LAB.HI),'d'); 
    Yp0  = zeros(d,'uint8'); Yp0(indx,indy,indz) = Yp0b; 
    Ywmh = Ywmh | (cat_vol_morph(Ycls{2}>64 | Yp0>2.25,'e') & cat_vol_morph(Yl1==1 | Yl1==2,'e',2) & Ymi<2.9 & Ymi>2.2); 
    Ywmh = Ywmh .* (1-Ymi)*3.1; clear Yp0; 
    % only the peaks WMHs
    Ywmh2 = nan(size(Ywmh),'single'); Ywmh2(Ywmh>0)=0; 
    Ywmh2(smooth3(Ywmh>0.5 & Ymi>0.5 & ~cat_vol_morph(NS(Yl1,LAB.VT),'d',2))>0.5)=1; 
    Ywmh2 = cat_vol_downcut(Ywmh2,Ywmh,-0.01);
    Ywmh2(smooth3(Ywmh2)<0.5)=0;
    %
    Ywmh = single(max(0,min(1,Ywmh.*Ywmh2 - smooth3(cat_vol_morph(NS(Yl1,LAB.VT),'d',2) & Ymi<0.66) ))*255);

    %% WMH as separate class 
    Yclso = Ycls;
    Yclssum = max(eps,single(Ycls{1})+single(Ycls{2})+single(Ycls{3}));
    Ycls{1} = cat_vol_ctype(single(Ycls{1}) - Ywmh .* (single(Ycls{1})./Yclssum));
    Ycls{2} = cat_vol_ctype(single(Ycls{2}) - Ywmh .* (single(Ycls{2})./Yclssum));
    Ycls{3} = cat_vol_ctype(single(Ycls{3}) - Ywmh .* (single(Ycls{3})./Yclssum));
    Ywmh    = cat_vol_ctype(Yclssum - single(Ycls{1}) - single(Ycls{2}) - single(Ycls{3}));
    clear Yclssum Ywmh2;  

    %% if the segmentation should be corrected later...
    if job.extopts.WMHC>1
      Yclso = Ycls;
    end

    % upate of the actual segmentation only for Dartel
    Ycls{2} = cat_vol_ctype(single(Ycls{2}) + single(Ywmh));

    if job.extopts.WMHC>1
      % update of Yp0b for WM/CSF PVE ROI

      Yp0  = single(Ycls{1})/255*2 + single(Ycls{2})/255*3 + single(Ycls{3})/255;
      Yp0  = Yp0(indx,indy,indz);
      Yl1b = Yl1(indx,indy,indz);
      Ymib = Ymi(indx,indy,indz);

      Ybs  = NS(Yl1b,LAB.BS) & Ymib>2/3; clear Ymib; 
      YpveVB = cat_vol_morph(NS(Yl1b,LAB.VT) | Ybs,'d',2);                % ventricle and brainstem
      YpveCC = cat_vol_morph(Yl1b==1,'d',3*vxv) & cat_vol_morph(Yl1b==2,'d',3*vxv) & ...
               cat_vol_morph(NS(Yl1b,LAB.VT),'d',2);                      % corpus callosum
      Ynpve  = smooth3(NS(Yl1b,LAB.BG) | NS(Yl1b,LAB.TH))>0.3;            % no subcortical structure 
      Yroi = (YpveVB | YpveCC) & ~Ynpve & ...
             cat_vol_morph(Yp0==3,'d',2) & cat_vol_morph(Yp0==1,'d',2) & ...
             Yp0<3 & Yp0>1 & ...
             smooth3((Yp0<3 & Yp0>1) & ~cat_vol_morph(Yp0<3 & Yp0>1,'o',1))>0.1;
      clear YpveVB YpveCC Ybs Ynpve ;         
      Yncm = (3-Yp0)/2.*Yroi; 

      Ycls{1}(indx,indy,indz) = min(Ycls{1}(indx,indy,indz),uint8(~Yroi*255));
      Ycls{2}(indx,indy,indz) = cat_vol_ctype(single(Ycls{2}(indx,indy,indz)).*~Yroi + (Yroi - Yncm)*255,'uint8');
      Ycls{3}(indx,indy,indz) = cat_vol_ctype(single(Ycls{3}(indx,indy,indz)).*~Yroi + Yncm*255,'uint8');
      clear Yp0 Yroi;

      Yp0b = cat_vol_ctype(single(Ycls{1})*2/3 + single(Ycls{2}) + single(Ycls{3})*1/3,'uint8');
      Yp0b = Yp0b(indx,indy,indz); 
    else
      if qa.subjectmeasures.WMH_rel>3 || qa.subjectmeasures.WMH_WM_rel>5 % #% of the TIV or the WM are affected
        cat_warnings = cat_io_addwarning(cat_warnings,...
          'MATLAB:SPM:CAT:cat_main:uncorrectedWMH',...
          sprintf('Uncorrected WM lesions (%2.2f%%%%%%%% of the WM)!',qa.subjectmeasures.WMH_WM_rel),1);
        fprintf('\n'); cat_io_cmd(' ','','',1);
      end
    end
    fprintf('%5.0fs\n',etime(clock,stime));
  else
    if qa.subjectmeasures.WMH_rel>3 || qa.subjectmeasures.WMH_WM_rel>5 % #% of the TIV or the WM are affected
      cat_warnings = cat_io_addwarning(cat_warnings,...
        'MATLAB:SPM:CAT:cat_main:uncorrectedWMH',...
        sprintf('Uncorrected WM lesions greater (%2.2f%%%%%%%% of the WM)!\\n',qa.subjectmeasures.WMH_rel));
    end
  end
  clear Yclsb;

  if job.extopts.verb>2
    Yp0  = single(Ycls{1})/255*2 + single(Ycls{2})/255*3 + single(Ycls{3})/255; %#ok<NASGU>
    tpmci=tpmci+1; tmpmat = fullfile(pth,reportfolder,sprintf('%s_%s%02d%s.mat',nam,'write',tpmci,'preDartel'));
    save(tmpmat,'Yp0','Ycls','Ymi','T3th','vx_vol','Yl1');
    clear Yp0;
  end


  % clear last 3 tissue classes to save memory
  % please do not try to write out these segmentations because class 4-6 are form SPM12
  % and class 1-3 from CAT12 and these are completely different segmentation approaches
  if all(cell2mat(struct2cell(job.output.TPMC)')==0)
    for i=6:-1:4, Ycls(i)=[]; end   
  end

  %% ---------------------------------------------------------------------
  %  Affine registration:
  %  Use the segmentation result (p0 label image) for the final affine 
  %  registration to a pseudo T1 image of the TPM. The dartel tempalte 
  %  can not be used here, due to special TPM image for chrildren or other
  %  species!
  %  ---------------------------------------------------------------------
  if 0
    % parameter
    aflags = struct('sep',job.opts.samp,'regtype',job.opts.affreg,'WG',[],'WF',[],'globnorm',0);

    % VG template
    cid = [2 3 1]; 
    VG = tpm.V(1); VG.dat = zeros(VG.dim,'uint8'); VG.dt = [2 0]; VG.pinfo(3) = 0;
    for ci=1:3 
      Yt = spm_read_vols(tpm.V(ci)); 
      VG.dat = VG.dat + cat_vol_ctype(Yt * 255 * cid(ci)/3,'uint8'); clear Yt 
    end

    % VF image
    VF = VT; VF.fname = [tempname '.nii']; VF = rmfield(VF,'private'); VF.pinfo(3)=0;
    VF.dat = zeros(d,'uint8'); VF.dt = [2 0]; 
    VF.dat(indx,indy,indz) = Yp0b; 

    % just brain mask - this improves affine registration but lead to worse
    % resultes in dartel normalization .. retry later without CSF
    % VG.dat = uint8(VG.dat>85)*255; VF.dat = uint8(VF.dat>85)*255;

    % smooth source with job.opts.samp mm
    VF = spm_smoothto8bit(VF,job.opts.samp/2); % smoothing, because of the TPM smoothness!

    %% affreg registration for one tissue segment

    % deactivated on 20160405 because something failed in the logitudinal process 
    warning off;  
    Affine = spm_affreg(VG, VF, aflags, res.Affine); 
    warning on; 
    rf=10^6; Affine    = fix(Affine * rf)/rf;
  else
    Affine = res.Affine; 
  end

  res.Affine0        = res.Affine;
  res.Affine         = Affine;

  clear VG VF cid %tpm 

else
%% SPM segmentation input  
%  ------------------------------------------------------------------------
%  Here, DARTEL and PBT processing is prepared. 
%  We simply use the SPM segmentation as it is without further modelling of
%  a PVE or other refinements. 
%  ------------------------------------------------------------------------
  
  % here we need the c1 filename
  VT0 = res.imagec(1);
  [pth,nam] = spm_fileparts(VT0.fname); 
  
  vx_vol              = sqrt(sum(VT.mat(1:3,1:3).^2));          % voxel size of the processed image
  cat_warnings        = struct('identifier',{},'message',{});   % warning structure from cat_main_gintnorm 
  NS                  = @(Ys,s) Ys==s | Ys==s+1;                % for side independent atlas labels
  
  % QA WMH values required by cat_vol_qa later
  qa.subjectmeasures.WMH_abs    = nan;  % absolute WMH volume without PVE
  qa.subjectmeasures.WMH_rel    = nan;  % relative WMH volume to TIV without PVE
  qa.subjectmeasures.WMH_WM_rel = nan;  % relative WMH volume to WM without PVE
  qa.subjectmeasures.WMH_abs    = nan;  % absolute WMH volume without PVE in cm^3
  
  % load SPM segments
  [pp,ff,ee] = spm_fileparts(res.image0(1).fname);
  Ycls{1} = uint8(spm_read_vols(spm_vol(fullfile(pp,['c1' ff ee])))*255); 
  Ycls{2} = uint8(spm_read_vols(spm_vol(fullfile(pp,['c2' ff ee])))*255); 
  Ycls{3} = uint8(spm_read_vols(spm_vol(fullfile(pp,['c3' ff ee])))*255); 

  % create (resized) label map and brainmask
  Yp0  = single(Ycls{3})/3 + single(Ycls{1})/3*2 + single(Ycls{2})/3*3;
  Yb   = Yp0>0.5;
  
  % load original images and get tissue thresholds
  Ysrc = spm_read_vols(spm_vol(fullfile(pp,[ff ee])));
  WMth = double(max(clsint(2),...
           cat_stat_nanmedian(cat_stat_nanmedian(cat_stat_nanmedian(Ysrc(Ycls{2}>192)))))); 
  T3th = [ min([  clsint(1) - diff([clsint(1),WMth]) ,clsint(3)]) , clsint(2) , WMth];
  if T3th(3)<T3th(2) % inverse weighting allowed 
    job.inv_weighting   = 1;                                     
  else
    job.inv_weighting   = 0; 
  end
  if ~debug, clear Ysrc; end
  
  
  % the intensity normalized images are here represented by the segmentation 
  Ym   = Yp0/255;
  Ymi  = Yp0/255; 
  
 
  % low resolution Yp0b
  sz = size(Yb);
  [indx, indy, indz] = ind2sub(sz,find(Yb>0));
  indx = max((min(indx) - 1),1):min((max(indx) + 1),sz(1));
  indy = max((min(indy) - 1),1):min((max(indy) + 1),sz(2));
  indz = max((min(indz) - 1),1):min((max(indz) + 1),sz(3));
  Yp0b = Yp0(indx,indy,indz);
  if ~debug, clear Yp0 Yb; end
  
  % load atlas map and prepare filling mask YMF
  % compared to CAT default processing, we have here the DARTEL mapping, but no individual refinement 
  Vl1 = spm_vol(job.extopts.cat12atlas{1});
  Yl1 = cat_vol_ctype(spm_sample_vol(Vl1,double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),0));
  Yl1 = reshape(Yl1,size(Ym)); [D,I] = cat_vbdist(single(Yl1>0)); Yl1 = Yl1(I);   
  YMF = NS(Yl1,job.extopts.LAB.VT) | NS(Yl1,job.extopts.LAB.BG) | NS(Yl1,job.extopts.LAB.BG); 
  
  fprintf('%5.0fs\n',etime(clock,stime));  
end



%% ---------------------------------------------------------------------
%  Spatial Registration with Dartel or Shooting
%  ---------------------------------------------------------------------
  [trans,res.ppe.reg] = cat_main_registration(job,res,Ycls,Yy,tpm.M);; 


%%  --------------------------------------------------------------------
%  write results
%  ---------------------------------------------------------------------

% reset WMHC
if (job.extopts.WMHC==1 || job.extopts.WMHC==3) && job.extopts.WMHCstr>0 && ~job.inv_weighting && exist('Yclso','var')
  Ycls = Yclso; 
end
clear Yclso;


stime = cat_io_cmd('Write result maps');
Yp0 = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)*3/255; 

% bias, noise and global corrected without masking for subject space and with masking for other spaces 
cat_io_writenii(VT0,Ym,mrifolder,'m', ...Dartel
  'bias and noise corrected, global intensity normalized','uint16',[0,0.0001], ... 
  min([1 0 2],[job.output.bias.native job.output.bias.warped job.output.bias.dartel]),trans);
cat_io_writenii(VT0,Ym.*(Yp0>0.1),mrifolder,'m', ... 
  'bias and noise corrected, global intensity normalized (masked due to normalization)','uint16',[0,0.0001], ...
  min([0 1 0],[job.output.bias.native job.output.bias.warped job.output.bias.dartel]),trans);

% bias, noise and local intensity corrected without masking for subject space and with masking for other spaces 
cat_io_writenii(VT0,Ymi,mrifolder,'mi', ...
  'bias and noise corrected, local intensity normalized','uint16',[0,0.0001], ... 
  min([1 0 2],[job.output.las.native job.output.las.warped job.output.las.dartel]),trans);
cat_io_writenii(VT0,Ymi.*(Yp0>0.1),mrifolder,'mi', ... 
  'bias and noise corrected, local intensity normalized (masked due to normalization)','uint16',[0,0.0001], ...
  min([0 1 0],[job.output.las.native job.output.las.warped job.output.las.dartel]),trans);
  
% Yp0b maps
if job.extopts.WMHC==3 && job.extopts.WMHCstr>0 && ~job.inv_weighting; 
  Yp0 = Yp0 + single(Ywmh)/255; 
end

cat_io_writenii(VT0,Yp0,mrifolder,'p0','Yp0b map','uint8',[0,4/255],job.output.label,trans);
clear Yp0; 

% partitioning
cat_io_writenii(VT0,Yl1,mrifolder,'a0','brain atlas map for major structures and sides',...
  'uint8',[0,1],job.output.atlas,trans);

% class maps 1-3
fn = {'GM','WM','CSF','head','head','background'};
for clsi=1:3
  cat_io_writenii(VT0,single(Ycls{clsi})/255,mrifolder,sprintf('p%d',clsi),...
    sprintf('%s tissue map',fn{clsi}),'uint8',[0,1/255],...
    min([1 1 0 3],[job.output.(fn{clsi}).native job.output.(fn{clsi}).warped ...
    job.output.(fn{clsi}).mod job.output.(fn{clsi}).dartel]),trans);
  cat_io_writenii(VT0,single(Ycls{clsi})/255,mrifolder,sprintf('p%d',clsi),...
    sprintf('%s tissue map',fn{clsi}),'uint16',[0,1/255],...
    min([0 0 3 0],[job.output.(fn{clsi}).native job.output.(fn{clsi}).warped ...
    job.output.(fn{clsi}).mod job.output.(fn{clsi}).dartel]),trans);
end

% write WMH class maps
if job.extopts.WMHC==3 && job.extopts.WMHCstr>0 && ~job.inv_weighting;
  cat_io_writenii(VT0,single(Ywmh)/255,mrifolder,'p7','WMH tissue map','uint8',[0,1/255],...
    min([1 1 0 2],[job.output.WMH.native job.output.WMH.warped ...
    job.output.WMH.mod job.output.WMH.dartel]),trans); % 1 0 0 0
  cat_io_writenii(VT0,single(Ywmh)/255,mrifolder,'p7','WMH tissue map','uint16',[0,1/255],...
    min([0 0 3 0],[job.output.WMH.native job.output.WMH.warped ...
    job.output.WMH.mod job.output.WMH.dartel]),trans); % 0 1 2 2
end 


% developer - intensity scaled tissue classe maps
% ----------------------------------------------------------------------
% The strong normalization of the T1 data can directly be used as tissue
% segmentation. The Ymi images is scaled to get similar maps for each 
% tissue class, with good visible differences in the sulci.
job.output.intsegments = job.extopts.experimental;
if job.output.intsegments
  if (any(tc(:)) || job.extopts.WMHC==3 && job.extopts.WMHCstr>0 && ~job.inv_weighting); 

    % intensity scaled tissue maps
    Yclsi = cell(1,3);
    for clsi=1:3
      clsid = [2 3 1];
      Yp0 = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)*3/255; 
      Yclsi{clsi} = Yp0toC(max(1/3,min(1,Ymi)) * 3 .* (Yp0>0),clsid(clsi));
      switch clsi
        case 1
          Yclsi{clsi} = Yclsi{clsi} .* (Ycls{clsi}>0) + ...
            Yp0toC(max(1/3,min(1,Ymi)) * 3 .* (Yp0>0),3) .* (Ycls{2}==0 & (Ycls{1}>0));
        case 2 
          Yclsi{clsi} = Yclsi{clsi} .* (Ycls{clsi}>0) +  ...
            Yp0toC(max(1/3,min(1,Ymi)) * 3 .* (Yp0>0),2) .* (Ycls{2}==255 & ~cat_vol_morph(Ycls{2}>192,'e'));
          Ywmhp = cat_vol_ctype( Yp0toC(max(1/3,min(1,Ymi)) * 3 .* (Yp0>0),2) .* cat_vol_morph(Ycls{2}>192,'e') * 255);
        case 3 
          Yclsi{clsi} = Yclsi{clsi} + ...
            (Yp0toC(max(1/3,min(1,Ymi)) * 3 .* (Yp0>0),2)) .* (Ycls{1}==0 & (Ycls{3}>0)) + ...
            (Yp0toC(max(1/3,min(1,Ymi)) * 3 .* (Yp0>0),3)) .* (Ycls{2}==0 & (Ycls{3}>0));
      end
      Yclsi{clsi} = cat_vol_ctype(Yclsi{clsi} * 255);
    end
    clear Yp0; 
    
    % class maps 1-3
    % Yclss = single(Yclsi{1} + Yclsi{2} + Yclsi{3} + Ywmhp + Ywmh)/255; % + single(Ywmh)/255;
    for clsi=1:3
      cat_io_writenii(VT0,single(Yclsi{clsi})/255,mrifolder,sprintf('pi%d',clsi),...
        sprintf('%s tissue map',fn{clsi}),'uint8',[0,1/255],...
        min([1 1 0 3],[job.output.(fn{clsi}).native job.output.(fn{clsi}).warped ...
        job.output.(fn{clsi}).mod job.output.(fn{clsi}).dartel]),trans);
      cat_io_writenii(VT0,single(Yclsi{clsi})/255,mrifolder,sprintf('pi%d',clsi),...
        sprintf('%s tissue map',fn{clsi}),'uint16',[0,1/255],...
        min([0 0 2 0],[job.output.(fn{clsi}).native job.output.(fn{clsi}).warped ...
        job.output.(fn{clsi}).mod job.output.(fn{clsi}).dartel]),trans);
    end
    clear Yclsi; 
    
    % write WMH class maps
    if job.extopts.WMHC==3 && job.extopts.WMHCstr>0 && ~job.inv_weighting;
      cat_io_writenii(VT0,(single(Ywmhp) + single(Ywmh))/255,mrifolder,...
        'pi7','WMH tissue map','uint8',[0,1/255],...
        min([1 1 0 2],[job.output.WMH.native job.output.WMH.warped ...
        job.output.WMH.mod job.output.WMH.dartel]),trans); % 1 0 0 0
      cat_io_writenii(VT0,(single(Ywmhp) + single(Ywmh))/255,mrifolder,...
        'pi7','WMH tissue map','uint16',[0,1/255],...
        min([0 0 2 0],[job.output.WMH.native job.output.WMH.warped ...
        job.output.WMH.mod job.output.WMH.dartel]),trans); % 0 1 2 2
    end 
    clear Ywmhp;
  end
end
% ----------------------------------------------------------------------


% classe maps 4-6 (for full TPM/template creation, e.g. for apes)
if any(cell2mat(struct2cell(job.output.TPMC)'))
  for clsi=4:6
    cat_io_writenii(VT0,single(Ycls{clsi})/255,mrifolder,sprintf('p%d',clsi),...
      sprintf('%s tissue map',fn{clsi}),'uint8',[0,1/255],...
      min([1 1 0 3],[job.output.TPMC.native job.output.TPMC.warped ...
      job.output.TPMC.mod job.output.TPMC.dartel]),trans);
    cat_io_writenii(VT0,single(Ycls{clsi})/255,mrifolder,sprintf('p%d',clsi),...
      sprintf('%s tissue map',fn{clsi}),'uint16',[0,1/255],...
      min([0 0 3 0],[job.output.TPMC.native job.output.TPMC.warped ...
      job.output.TPMC.mod job.output.TPMC.dartel]),trans);
  end
end
%clear cls clsi fn Ycls; % we need this maps later for the ROIs

% write jacobian determinant
if job.output.jacobian.warped
  %%
  if do_dartel==2 % shooting
    %[y0, dt] = spm_dartel_integrate(reshape(trans.jc.u,[trans.jc.odim(1:3) 1 3]),[1 0], 6);
    dt = trans.jc.dt2; 
    dx = 10; % smaller values are more accurate, but large look better; 
    [D,I] = cat_vbdist(single(~(isnan(dt) | dt<0 | dt>100) )); D=min(1,D/min(dx,max(D(:)))); 
    dt = dt(I); dt = dt .* (1-D) + D; dt(isnan(dt))=1; 
    %dt = 1/max(eps,dt); 
    clear y0 D I
  else %dartel
    [y0, dt] = spm_dartel_integrate(reshape(trans.jc.u,[trans.warped.odim(1:3) 1 3]),[1 0], 6);
    clear y0
  end
  N      = nifti;
  N.dat  = file_array(fullfile(pth,mrifolder,['wj_', nam, '.nii']),trans.warped.odim(1:3),...
             [spm_type('float32') spm_platform('bigend')],0,1,0);
  N.mat  = trans.warped.M1;
  N.mat0 = trans.warped.M1;
  N.descrip = ['Jacobian' VT0.descrip];
  create(N);
  N.dat(:,:,:) = dt;
end

% deformations y - dartel > subject
if job.output.warps(1)
    Yy        = spm_diffeo('invdef',trans.warped.yx,trans.warped.odim,eye(4),trans.warped.M0);
    N         = nifti;
    N.dat     = file_array(fullfile(pth,mrifolder,['y_', nam1, '.nii']),[trans.warped.odim(1:3),1,3],'float32',0,1,0);
    N.mat     = trans.warped.M1;
    N.mat0    = trans.warped.M1;
    N.descrip = 'Deformation';
    create(N);
    N.dat(:,:,:,:,:) = reshape(Yy,[trans.warped.odim,1,3]);
end

if ~debug, clear Yy; end

% deformation iy - subject > dartel
if job.output.warps(2)
  % transformation from voxel to mm space
  yn = numel(trans.warped.yx); 
  p  = ones([4,yn/3],'single'); 
  p(1,:) = trans.warped.yx(1:yn/3);
  p(2,:) = trans.warped.yx(yn/3+1:yn/3*2);
  p(3,:) = trans.warped.yx(yn/3*2+1:yn);
  p      = M1(1:3,:) * p;
  
  if any(trans.native.Vo.dim~=trans.native.Vi.dim)
    %% update cat_stat_nanmedian for interpolated images
    Vdef = res.image0(1);
    Vdef.dt(1) = spm_type('float32');
    Vdef = rmfield(Vdef,'private');
    Vdef.dat = zeros(size(Vdef.dim),'single');
    Vdef.pinfo(3) = 0; 
    Vdef.fname = fullfile(pth,mrifolder,['iy2_r', nam1, '.nii']);
    Yy2 = zeros([trans.native.Vo.dim(1:3) 1 3],'double');
    Vyy = VT; Vyy.pinfo(3)=0; Vyy.dt=[16 0]; Vyy = rmfield(Vyy,'private');  
    Yy1 = zeros([res.image(1).dim(1:3),3],'single'); 
    Yy1(1:yn/3)        = p(1,:);
    Yy1(yn/3+1:yn/3*2) = p(2,:);
    Yy1(yn/3*2+1:yn)   = p(3,:);
    for i=1:3
      Vyy.dat=Yy1(:,:,:,i); 
      [Vt,Yy2(:,:,:,:,i)] = cat_vol_imcalc(Vyy,Vdef,'i1',struct('interp',6));
    end
    clear Vt Vdef Vyy
  else 
    yn = numel(trans.warped.y); 
    p  = ones([4,yn/3],'single'); 
    p(1,:) = trans.warped.y(1:yn/3);
    p(2,:) = trans.warped.y(yn/3+1:yn/3*2);
    p(3,:) = trans.warped.y(yn/3*2+1:yn);
    p      = M1(1:3,:) * p;

    Yy2 = zeros([res.image(1).dim(1:3),1,3],'single'); 
    Yy2(1:yn/3)        = p(1,:);
    Yy2(yn/3+1:yn/3*2) = p(2,:);
    Yy2(yn/3*2+1:yn)   = p(3,:);
  end
  clear p; 
  
  % f2 = spm_diffeo('resize', f1, dim)
  % write new output
  Ndef      = nifti;
  Ndef.dat  = file_array(fullfile(pth,mrifolder,['iy_', nam1, '.nii']),[res.image0(1).dim(1:3),1,3],...
              [spm_type('float32') spm_platform('bigend')],0,1,0);
  Ndef.mat  = res.image0(1).mat;
  Ndef.mat0 = res.image0(1).mat;
  Ndef.descrip = 'Inverse Deformation';
  create(Ndef);
  Ndef.dat(:,:,:,:,:) = Yy2;
  clear Yy2;
end
fprintf('%5.0fs\n',etime(clock,stime));



%% ---------------------------------------------------------------------
%  surface creation and thickness estimation
%  ---------------------------------------------------------------------
if job.output.surface
  stime = cat_io_cmd('Surface and thickness estimation');; 
  
  % specify WM/CSF width/depth/thickness estimation
  if job.output.surface>10
    job.output.surface=job.output.surface-10;
    WMT = 1; 
  else
    WMT = 0; 
  end
  if job.extopts.experimental || job.extopts.expertgui==2
    WMT = 1; 
  end
  % specify surface
  switch job.output.surface
    case 1, surf = {'lh','rh'};
    case 2, surf = {'lh','rh','lc','rc'};
    case 3, surf = {'lh'};
    case 4, surf = {'rh'};
    % fast surface reconstruction without simple spherical mapping     
    case 5, surf = {'lhfst','rhfst'};                   
    case 6, surf = {'lhfst','rhfst','lcfst','rcfst'};    
    % fast surface reconstruction with simple spherical mapping     
    case 7, surf = {'lhsfst','rhsfst'};                    
    case 8, surf = {'lhsfst','rhsfst','lcsfst','rcsfst'}; 
  end
  if job.output.surface>4 % fast 
    job.extopts.pbtres = max(0.8,min([((min(vx_vol)^3)/2)^(1/3) 1.0]));
  end
  
  % brain masking and correction of blood vessels 
  Yp0 = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)*3/255; 
  
  % surface creation and thickness estimation
  if 1
    %% using the Ymi map
    Ymix = Ymi .* (Yp0>0.5);
    if ~debug, clear Yp0; end 
    
    [Yth1,S,Psurf] = cat_surf_createCS(VT,Ymix,Yl1,YMF,...
      struct('interpV',job.extopts.pbtres,'Affine',res.Affine,'surf',{surf},'inv_weighting',job.inv_weighting,...
      'verb',job.extopts.verb,'WMT',WMT)); 
  else
    %% using the segmentation
    [Yth1,S,Psurf] = cat_surf_createCS(VT,Yp0/3,Yl1,YMF,...
      struct('interpV',job.extopts.pbtres,'Affine',res.Affine,'surf',{surf},'inv_weighting',job.inv_weighting,...
      'verb',job.extopts.verb,'WMT',WMT));
  end

  cat_io_cmd('Surface and thickness estimation');  
  fprintf('%5.0fs\n',etime(clock,stime));
  if ~debug; clear YMF; end
else
  if ~debug; clear Ymi; end
end




%% ---------------------------------------------------------------------
%  ROI Partitioning 
%  ---------------------------------------------------------------------
%  This part estimated indivudal measurements for different ROIs.
%  The ROIs are described in the CAT normalized space and there are to 
%  ways to estimate them - (1) in subject space, and (2) in normalized 
%  space. Estimation in normalized space is more direct an avoid further
%  transformations. The way over the subject space have the advantage 
%  that indivdiual anatomical refinients are possible, but the this has
%  to be done and evalutated for each atlas. 
%  ---------------------------------------------------------------------
if job.output.ROI
  stime = cat_io_cmd('ROI estimation');   
  if job.extopts.verb, fprintf('\n'); end; 
  
  % get atlases
  FAF = job.extopts.atlas; 
  FA  = {}; fai = 1;
  AN  = fieldnames(job.output.atlases);
  for ai = 1:numel(AN)
    fafi = find(cellfun('isempty',strfind(FAF(:,1),[AN{ai} '.']))==0);
    if ~isempty(fafi) && job.output.atlases.(AN{ai}), FA(fai,:) = FAF(fafi,:); fai = fai+1; end
  end
  
  % get atlas resolution 
  % we sort the atlases and to reduce data resampling
  VA = spm_vol(char(FA(1:end,1))); 
  for ai=1:numel(VA), VAvx_vol(ai,:) = sqrt(sum(VA(ai).mat(1:3,1:3).^2)); end   %#ok<AGROW>
  [VAs,VAi] = sortrows(VAvx_vol); 
  FA = FA(VAi,:); VA = VA(VAi,:); VAvx_vol = VAvx_vol(VAi,:); %clear VA; 
  
  %%
  for ai=1:size(FA,1)
    if ai==1 || any(VAvx_vol(ai,:)~=VAvx_vol(ai-1,:))
      % resampe data in atlas resolution for the first time or if the atlas resolution change 
      
      Yp0 = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)*3/255; 

      % map data to actual template space
      if ai==1
        stime2  = cat_io_cmd('  Data mapping to normalized space','g5','', job.extopts.verb-1); 
      else
        stime2  = cat_io_cmd('  Data mapping to normalized space','g5','', job.extopts.verb-1,stime2); 
      end  
      transw      = trans.warped;                     % dartel/shooting deformation data 
      transw.odim = VA(ai).dim;                       % adaption for atlas image size
      transw.ress = job.extopts.vox(1)./VAvx_vol(ai,:);  % adaption for atlas sampling resolution 
      
      wYp0     = cat_vol_ROInorm(Yp0,transw,1,0,job.extopts.atlas);
      wYcls    = cat_vol_ROInorm(Ycls,transw,1,1,job.extopts.atlas);

      if exist('Ywmh','var')
        wYcls{7} = cat_vol_ctype(cat_vol_ROInorm(single(Ywmh),transw,1,0,job.extopts.atlas));
      end
      for ci=1:numel(wYcls); wYcls{ci} = wYcls{ci} * prod(vx_vol); end      % volume
      %wYm      = cat_vol_ROInorm(Ym,transw,1,0,job.extopts.atlas);         % intensity
      Yp0toC   = @(Yp0,c) 1-min(1,abs(Yp0-c));
      if exist('Yth1','var')
        % ROI based thickness of all GM voxels per ROI
        Yth1x  = Yth1; Yth1x(Yp0toC(Yp0,2)<0.5)=nan;
        Ymim   = single(smooth3(Yp0>1.5 & Yp0<2.5 & (Yl1==1 | Yl1==2))>0.5);
        wYmim  = cat_vol_ROInorm(Ymim,transw,1,0,job.extopts.atlas)>0.5;
        wYth1  = cat_vol_ROInorm(Yth1x,transw,1,0,job.extopts.atlas);
        clear Ymim; 
      end
    end
    
   %% clear YMF
 
    tissue = FA{ai,3};
    [px,atlas] = fileparts(FA{ai,1}); 

    % ds('l2','',1.5,wYm,round(wYp0),wYm,single(wYa)/50 .* (wYp0<2.5),70)

    stime2 = cat_io_cmd(sprintf('  ROI estimation of ''%s'' atlas',atlas),'g5','', job.extopts.verb-1,stime2);

    % map atlas to actual template space 
    transa      = trans.warped; 
    transa.M1   = VA(ai).mat;
    transa.odim = transw.odim;
    wYa   = cat_vol_ROInorm([],transa,ai,0,FA);

    % write output
    if any(cell2mat(struct2cell(job.output.atlas)'))
      % map atlas in native space
      Vlai = spm_vol(FA{ai,1});
      Ylai = cat_vol_ctype(spm_sample_vol(Vlai,double(trans.warped.yx(:,:,:,1)),...
        double(trans.warped.yx(:,:,:,2)),double(trans.warped.yx(:,:,:,3)),0));
      Ylai = reshape(Ylai,size(Yp0)); 

      % write map (mri as tissue subforder and mri_atals as ROI subfolder)
      if isempty(mrifolder), amrifolder = ''; else amrifolder = 'mri_atlas'; end
      cat_io_writenii(VT0,Ylai,amrifolder,[atlas '_'],[atlas ' original'],...
        'uint8',[0,1],job.output.atlas,trans);
      clear Vlai Ylai;
    end

    % extract ROI data
    csv   = cat_vol_ROIestimate(wYp0,wYa,wYcls,ai,'V',[],tissue,FA);  % volume

    % thickness
    if exist('Yth1','var'),
    % for thickness we need special corrections to avoid values 
    % in poor ROIs that contribute to GM
      csv    = cat_vol_ROIestimate(wYp0,wYa,wYth1.*wYmim,ai,'T',csv,tissue,job.extopts.atlas);
      csvth1 = cat_vol_ROIestimate(wYp0,wYa,wYcls{2}.*wYmim,ai,'V',[] ,{''},job.extopts.atlas);
      corth1 = [csv{2:end,end}]; 
      try, corth1(corth1<mean(vx_vol)/2 | [csvth1{2:end,end}]<0.5)=nan; end
      csv(2:end,end) = num2cell(corth1);
    end

    % xml-export one file for all (this is a structure)
    ROI.(atlas) = csv;
    
  end 
  
  % write results
  catROI = cat_roi_fun('csvtab2xmlroi',ROI);
  cat_io_xml(fullfile(pth,labelfolder,['catROI_' nam '.xml']),catROI,'write'); 

  cat_io_cmd(' ','g5','',job.extopts.verb,stime2);
  cat_io_cmd('','n','',1,stime);
end
clear wYp0 wYcls wYv trans


%  estimate volumes and TIV
qa.subjectmeasures.vol_abs_CGW = [prod(vx_vol)/1000/255 .* sum(Ycls{1}(:)), ... CSF
                    prod(vx_vol)/1000/255 .* sum(Ycls{2}(:)), ... GM 
                    prod(vx_vol)/1000/255 .* sum(Ycls{3}(:)), ... WM
                    qa.subjectmeasures.WMH_abs]; 
qa.subjectmeasures.vol_TIV     =  sum(qa.subjectmeasures.vol_abs_CGW); 
qa.subjectmeasures.vol_rel_CGW =  qa.subjectmeasures.vol_abs_CGW ./ qa.subjectmeasures.vol_TIV;
clear Ycls

%% ---------------------------------------------------------------------
%  XML-report and Quality Assurance
%  ---------------------------------------------------------------------
stime = cat_io_cmd('Quality check'); job.stime = stime; 
Yp0   = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)/255*3; %qa2=qa;
qa    = cat_vol_qa('cat12',Yp0,fname0,Ym,res,cat_warnings,job.extopts.species, ...
          struct('write_csv',0,'write_xml',1,'method','cat12','job',job,'qa',qa));

% WMH updates? ... has to be done within cat_vol_qa?!
%qa.subjectmeasures.vol_abs_CGW(2) = qa.subjectmeasures.vol_abs_CGW(2) - qa2.subjectmeasures.WMH_abs;
%qa.subjectmeasures.vol_abs_CGW(4) = qa2.subjectmeasures.WMH_abs;
%qa.subjectmeasures.vol_rel_CGW    = qa.subjectmeasures.vol_abs_CGW ./ qa.subjectmeasures.vol_TIV;
if job.output.surface && exist('S','var')
  % metadata
  if isfield(S,'lh') && isfield(S.lh,'th1'), th=S.lh.th1; else th=[]; end;
  if isfield(S,'rh') && isfield(S.rh,'th1'), th=[th; S.rh.th1]; end
  qa.subjectmeasures.dist_thickness{1} = [cat_stat_nanmean(th(:)) cat_stat_nanstd(th(:))]; clear th; 
  if job.extopts.expertgui>1
    if isfield(S,'lh') && isfield(S.lh,'th2'), th=S.lh.th2; else th=[]; end; 
    if isfield(S,'rh') && isfield(S.lh,'th2'), th=[th; S.rh.th2]; end
    qa.subjectmeasures.dist_gyruswidth{1} = [cat_stat_nanmean(th(:)) cat_stat_nanstd(th(:))]; clear th; 
    if isfield(S,'lh') && isfield(S.lh,'th3'), th=S.lh.th3; else th=[]; end; 
    if isfield(S,'rh') && isfield(S.lh,'th3'), th=[th; S.rh.th3]; end
    qa.subjectmeasures.dist_sulcuswidth{1} = [cat_stat_nanmean(th(:)) cat_stat_nanstd(th(:))]; clear th; 
  end

  qam = cat_stat_marks('eval',job.cati,qa,'cat12');

  cat_io_xml(fullfile(pth,reportfolder,['cat_' nam '.xml']),struct(...
    ... 'subjectratings',qam.subjectmeasures, ... not ready
    'subjectmeasures',qa.subjectmeasures,'ppe',res.ppe),'write+'); % here we have to use the write+!

end  
clear Yo Yp0 qas;
fprintf('%5.0fs\n',etime(clock,stime));


if job.extopts.print
  %% display and print result if possible
  %  ---------------------------------------------------------------------
    QMC   = cat_io_colormaps('marks+',17);
    color = @(QMC,m) QMC(max(1,min(size(QMC,1),round(((m-1)*3)+1))),:);

    warning off; %#ok<WNOFF> % there is a div by 0 warning in spm_orthviews in linux


    %mark2str2 = @(mark,s,val) sprintf(sprintf('\\\\bf\\\\color[rgb]{%%0.2f %%0.2f %%0.2f}%s',s),color(QMC,mark),val);
    marks2str = @(mark,str) sprintf('\\bf\\color[rgb]{%0.2f %0.2f %0.2f}%s',color(QMC,mark),str);
    mark2rps    = @(mark) min(100,max(0,105 - mark*10));
    grades      = {'A+','A','A-','B+','B','B-','C+','C','C-','D+','D','D-','E+','E','E-','F'};
    mark2grad   = @(mark) grades{min(numel(grades),max(max(isnan(mark)*numel(grades),1),round((mark+2/3)*3-3)))};


    % CAT GUI parameter:
    % --------------------------------------------------------------------
    str = [];

    % 1 line: Matlab, SPM12, CAT12 version number and GUI and experimental mode 
    str = [str struct('name', 'Version: Matlab / SPM12 / CAT12:','value',...
      sprintf('%s / %s / %s',qa.software.version_matlab,qa.software.version_spm,qa.software.version_cat))];
    if     job.extopts.expertgui==1, str(end).value = [str(end).value '\bf\color[rgb]{0 0.2 1}e']; 
    elseif job.extopts.expertgui==2, str(end).value = [str(end).value '\bf\color[rgb]{0 0.2 1}d'];
    end  
    if job.extopts.experimental, str(end).value = [str(end).value '\bf\color[rgb]{0 0.2 1}x']; end  

    % 2 lines: TPM, Template, Normalization method with voxel size
    str = [str struct('name', 'Tissue Probability Map:','value',strrep(spm_str_manip(res.tpm(1).fname,'k40d'),'_','\_'))];
    if do_dartel
      if job.extopts.regstr==0 % Dartel
        str = [str struct('name', 'Dartel Registration to: ',...
                          'value',strrep(spm_str_manip(job.extopts.darteltpm{1},'k40d'),'_','\_'))];
      elseif job.extopts.regstr==4 % Dartel
        str = [str struct('name', 'Shooting Registration to: ',...
                          'value',strrep(spm_str_manip(job.extopts.shootingtpm{1},'k40d'),'_','\_'))];
      else
        if job.extopts.expertgui==0
          str = [str struct('name','Optimized Shooting Registration to:',...
                            'value',strrep(spm_str_manip(job.extopts.shootingtpm{1},'k40d'),'_','\_'))];
        else
          str = [str struct('name', sprintf('Optimized Shooting Registration (regstr:%s) to :',sprintf('%g ',job.extopts.regstr)),...
                            'value',strrep(spm_str_manip(job.extopts.shootingtpm{1},'k40d'),'_','\_'))];
        end
      end
    end
    
    % 1 line 1: Affreg
    str = [str struct('name', 'affreg:','value',sprintf('%s',job.opts.affreg))];
    % 1 line 2: APP
    APPstr = {'none','light','full','APPi','APPf','animal'}; APPstr{1071} = 'APP'; 
    str(end).name  = [str(end).name(1:end-1) ' / APP '];  
    str(end).value = [str(end).value sprintf(' / %s',APPstr{job.extopts.APP+1})];
    % 1 line 3: biasstr / biasreg+biasfwhm
    if job.opts.biasstr>0
      biasstr = {'ultralight','light','medium','strong','heavy'};
      str(end).name  = [str(end).name(1:end-1) ' / biasstr '];  
      str(end).value = [str(end).value sprintf(' / %s',biasstr{round(job.opts.biasstr*4)+1})];
      if job.extopts.expertgui % add the value
        str(end).value = [str(end).value sprintf('(%0.2f;breg:%0.2f;bfwhm:%0.2f)',job.opts.biasstr,job.opts.biasreg,job.opts.biasfwhm)]; 
      end
    else
      str(end).name  = [str(end).name(1:end-1) ' / biasreg / biasfwhm'];
      str(end).value = [str(end).value sprintf(' / %0.2f / %0.2f',job.opts.biasreg,job.opts.biasfwhm)]; 
    end


    % 1 line: adaptive noise parameter ( MRFstr + SANLM + NCstr )
    defstr  = {'none','ultralight','light','medium','strong','heavy',... sanlm vs. isarnlm
               'ultralight+','ultralight+','light+','medium+','strong+','heavy+'};
    defstrm = @(x) defstr{ round(max(0,min(2,x))*4) + 1 + (x>0) + (x>1)};
    str = [str struct('name', 'Noise reduction:','value','')]; 
    if job.extopts.NCstr==0 
      if job.extopts.mrf==0
        str(end).value = 'no noise correction';
      else
        if job.extopts.expertgui==0
          str(end).value = 'MRF'; 
        else
          str(end).value = sprintf('MRF(%0.2f)',job.extopts.mrf); 
        end  
      end
    elseif job.extopts.NCstr==2 
      str(end).value = 'ISARNLM(light)'; 
    elseif job.extopts.NCstr==3
      str(end).value = 'ISARNLM(full)'; 
    elseif job.extopts.NCstr==1
      str(end).value = 'SANLM(full)';
    else
      str(end).value = 'SANLM(light)';
    end
    if job.extopts.NCstr~=0 && job.extopts.mrf
      if job.extopts.expertgui==0
        str(end).value = ' + MRF'; 
      else
        str(end).value = sprintf(' + MRF(%0.2f)',job.extopts.mrf); 
      end 
    end


    % 1 line(s): LASstr / GCUTstr / CLEANUPstr
    str(end).name  = 'LASstr / GCUTstr / CLEANUPstr:';
    if ~job.extopts.expertgui
      str(end).value = sprintf('%s / %s / %s',defstrm(job.extopts.LASstr),defstrm(job.extopts.gcutstr),defstrm(job.extopts.cleanupstr)); 
    else
      str(end).value = sprintf('%s(%0.2f) / %s(%0.2f) / %s(%0.2f)',defstrm(job.extopts.LASstr),job.extopts.LASstr,...
        defstrm(job.extopts.gcutstr),job.extopts.gcutstr,defstrm(job.extopts.cleanupstr),job.extopts.cleanupstr); 
    end
    if job.extopts.expertgui 
      restype = char(fieldnames(job.extopts.restypes));
      str = [str struct('name', 'WMHC / WMHCstr / BVCstr / restype:','value',...
             sprintf('%d / %0.2f / %0.2f / %s',...
            job.extopts.WMHC,job.extopts.WMHCstr,job.extopts.BVCstr,restype))];
      if ~strcmp('native',restype)
        str(end).value = [str(end).value sprintf('(%0.2f %0.2f)',job.extopts.restypes.(restype))];
      end; 
    end

    % line 8: surfae parameter
    if job.output.surface
      str = [str struct('name', 'Voxel resolution (original > internal > PBT; vox):',...
             'value',sprintf('%4.2fx%4.2fx%4.2f mm%s > %4.2fx%4.2fx%4.2f mm%s > %4.2f mm%s; %4.2f mm  ', ...
             qa.qualitymeasures.res_vx_vol,char(179),qa.qualitymeasures.res_vx_voli,char(179),job.extopts.pbtres,char(179),job.extopts.vox(1)))];
    else
      str = [str struct('name', 'Voxel resolution (original > intern; vox):',...
             'value',sprintf('%4.2fx%4.2fx%4.2f mm%s > %4.2fx%4.2fx%4.2f mm%s; %4.2f mm', ...
             qa.qualitymeasures.res_vx_vol,char(179),qa.qualitymeasures.res_vx_voli,char(179),job.extopts.vox(1)))];
    end       
    % str = [str struct('name', 'Norm. voxel size:','value',sprintf('%0.2f mm',job.extopts.vox))]; % does not work yet 


    % Image Quality measures:
    % --------------------------------------------------------------------
    str2 =       struct('name', '\bfImage and Preprocessing Quality:','value',''); 
    str2 = [str2 struct('name',' Resolution:','value',marks2str(qa.qualityratings.res_RMS,...
      sprintf('%5.2f%% (%s)',mark2rps(qa.qualityratings.res_RMS),mark2grad(qa.qualityratings.res_RMS))))];
    str2 = [str2 struct('name',' Noise:','value',marks2str(qa.qualityratings.NCR,...
      sprintf('%5.2f%% (%s)',mark2rps(qa.qualityratings.NCR),mark2grad(qa.qualityratings.NCR))))];
    str2 = [str2 struct('name',' Bias:','value',marks2str(qa.qualityratings.ICR,...
      sprintf('%5.2f%% (%s)',mark2rps(qa.qualityratings.ICR),mark2grad(qa.qualityratings.ICR))))]; % not important and more confussing 
    str2 = [str2 struct('name','\bf Weighted average (IQR):','value',marks2str(qa.qualityratings.IQR,...
      sprintf('%5.2f%% (%s)',mark2rps(qa.qualityratings.IQR),mark2grad(qa.qualityratings.IQR))))];


    % Subject Measures
    % --------------------------------------------------------------------
    % Volume measures
    if job.extopts.WMHC>1
      str3 = struct('name', '\bfVolumes:','value',sprintf('%5s %5s %5s %5s%s','CSF','GM','WM','WMH')); 
      str3 = [str3 struct('name', ' Absolute volume:','value',sprintf(['%5.0f %5.0f %5.0f %5.0f cm' char(179)], ...
              qa.subjectmeasures.vol_abs_CGW(1:4)))];
      str3 = [str3 struct('name', ' Relative volume:','value',sprintf('%5.1f %5.1f %5.1f %5.1f %%', ...
              qa.subjectmeasures.vol_rel_CGW(1:4)*100))];
    else
      str3 = struct('name', '\bfVolumes:','value',sprintf('%5s %5s %5s %5s%s','CSF','GM','WM')); 
      str3 = [str3 struct('name', ' Absolute volume:','value',sprintf(['%5.0f %5.0f %5.0f cm' char(179)], ...
              qa.subjectmeasures.vol_abs_CGW(1:3)))];
      str3 = [str3 struct('name', ' Relative volume:','value',sprintf('%5.1f %5.1f %5.1f %%', ...
              qa.subjectmeasures.vol_rel_CGW(1:3)*100))];
      Ywmhrel = NS(Yl1,23);
      Yp0 = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)/255*3; 
      WMH_rel    = 100*sum(Ywmhrel(:)) / sum(Yp0(:)>0.5); 
      WMH_WM_rel = 100*sum(Ywmhrel(:)) / sum(Yp0(:)>2.5); 
      if WMH_rel>3 || WMH_WM_rel>5
        str3(end).value = [str3(end).value sprintf('\\bf\\color[rgb]{0.8 0 0} WMHs!')];   
        %str3(end).value = [str3(end).value sprintf('\\bf\\color[rgb]{0.8 0 0}~%0.0f%%WMHs!',WMH_WM_rel)];    
      end
      clear Yp0;
    end
    %
    str3 = [str3 struct('name', ' TIV:','value', sprintf(['%0.0f cm' char(179)],qa.subjectmeasures.vol_TIV))];  

    % Surface measures - Thickness, (Curvature, Depth, ...)
    if isfield(qa.subjectmeasures,'dist_thickness') && ~isempty(qa.subjectmeasures.dist_thickness)
      str3 = [str3 struct('name', '\bfThickness:','value',sprintf('%5.2f%s%5.2f mm', ...
             qa.subjectmeasures.dist_thickness{1}(1),177,qa.subjectmeasures.dist_thickness{1}(2)))];
      if isfield(qa.subjectmeasures,'dist_gyruswidth')
        str3 = [str3 struct('name', '\bfGyruswidth:','value',sprintf('%5.2f%s%5.2f mm', ...
               qa.subjectmeasures.dist_gyruswidth{1}(1),177,qa.subjectmeasures.dist_gyruswidth{1}(2)))];
      end
      if isfield(qa.subjectmeasures,'dist_sulcuswidth')
        str3 = [str3 struct('name', '\bfSulcuswidth:','value',sprintf('%5.2f%s%5.2f mm', ...
               qa.subjectmeasures.dist_sulcuswidth{1}(1),177,qa.subjectmeasures.dist_sulcuswidth{1}(2)))];
      end
    end

    % Preprocessing Time
    if job.extopts.experimental || job.extopts.expertgui>0 || 1
      str2 = [str2 struct('name','\bfProcessing time:','value',sprintf('%02.0f:%02.0f min', ...
      floor(round(etime(clock,res.stime))/60),mod(round(etime(clock,res.stime)),60)))]; 
    end

    % Warnings
    if numel(cat_warnings)>0 && job.extopts.expertgui>0
      str2 = [str2 struct('name', '','value','')]; 
      str2 = [str2 struct('name', '\bfWarnings:','value','')]; 
      for wi=1:numel(cat_warnings)
        shorter = cat_warnings(wi).identifier;
        % remove leading MATLAB, SPM or CAT elements
        dots    = max([min(strfind(shorter,'MATLAB')+7), ...
                       min(strfind(shorter,'SPM')+4), ...
                       min(strfind(shorter,'CAT')+4)]);
        if ~isempty(dots), shorter = shorter(dots:end); end
        % limit lenght of the string and replace critical character
        shorter = spm_str_manip(shorter,'l40');
        shorter = marks2str(4,shorter);
        shorter = strrep(shorter,'_','\_');
        str2    = [str2 struct('name',shorter,'value','')];  %#ok<AGROW>
      end
    end

    % adding one space for correct printing of bold fonts
    for si=1:numel(str)
      str(si).name   = [str(si).name  '  '];  str(si).value  = [str(si).value  '  '];
    end
    for si=1:numel(str2)
      str2(si).name  = [str2(si).name '  '];  str2(si).value = [str2(si).value '  '];
    end
    for si=1:numel(str3)
      str3(si).name  = [str3(si).name '  '];  str3(si).value = [str3(si).value '  '];
    end


    %
    fg = spm_figure('FindWin','Graphics'); 
    set(0,'CurrentFigure',fg)
    if isempty(fg)
      if job.nproc, fg = spm_figure('Create','Graphics','visible','off'); else fg = spm_figure('Create','Graphics'); end;
    else
      if job.nproc, set(fg,'visible','off'); end
    end
    set(fg,'windowstyle','normal'); 
    spm_figure('Clear','Graphics'); 
    switch computer
      case {'PCWIN','PCWIN64'}, fontsize = 8;
      case {'GLNXA','GLNXA64'}, fontsize = 8;
      case {'MACI','MACI64'},   fontsize = 9.5;
      otherwise,                fontsize = 9.5;
    end
    ax=axes('Position',[0.01 0.75 0.98 0.24],'Visible','off','Parent',fg);

    text(0,0.99,  ['Segmentation: ' spm_str_manip(res.image0(1).fname,'k60d') '       '],...
      'FontSize',fontsize+1,'FontWeight','Bold','Interpreter','none','Parent',ax);

    cm = job.extopts.colormap; 

    % check colormap name
    switch lower(cm)
      case {'jet','hsv','hot','cool','spring','summer','autumn','winter',...
          'gray','bone','copper','pink','bcgwhw','bcgwhn'}
      otherwise
        cat_io_cprintf(job.color.warning,'WARNING:Unknown Colormap - use default.\n'); 
        cm = 'gray';
    end

    % SPM_orthviews seams to allow only 60 values
    % It further requires a modified colormaps with lower values that the
    % colorscale and small adaption for the values. 
    surfcolors = 128; 
    switch lower(cm)
      case {'bcgwhw','bcgwhn'} % cat colormaps with larger range
        ytick       = [1,5:5:60];
        yticklabel  = {' BG',' ',' CSF',' CGM',' GM',' GWM',' WM',' ',' ',' ',' ',' ',' BV / HD '};
        yticklabelo = {' BG',' ','    ','    ','   ','     ',' avg WM  ',' ',' ',' ',' ',' ',' BV / HD '};
        yticklabeli = {' BG',' ','    ','    ','   ','  ','  ',' ',' ',' ',' ',' ',' BV / HD '};
        %colormap(cat_io_colormaps(cm,60));
        cmap = [cat_io_colormaps([cm 'ov'],60);flipud(cat_io_colormaps([cm 'ov'],60));jet(surfcolors)]; 
        cmmax = 2;
      case {'jet','hsv','hot','cool','spring','summer','autumn','winter','gray','bone','copper','pink'}
        ytick       = [1 20 40 60]; 
        yticklabel  = {' BG',' CSF',' GM',' WM'};
        yticklabelo = {' BG','    ','   ',' WM'};
        yticklabeli = {' BG','    ','   ','   '};
        cmap = [eval(sprintf('%s(60)',cm));flipud(eval(sprintf('%s(60)',cm)));jet(surfcolors)]; 
        cmmax = 1;
    end
    colormap(cmap);
    spm_orthviews('Redraw');

    htext = zeros(5,2,2);
    for i=1:size(str,2)   % main parameter
      htext(1,i,1) = text(0.01,0.98-(0.055*i), str(i).name  ,'FontSize',fontsize, 'Interpreter','tex','Parent',ax);
      htext(1,i,2) = text(0.51,0.98-(0.055*i), str(i).value ,'FontSize',fontsize, 'Interpreter','tex','Parent',ax);
    end
    for i=1:size(str2,2)  % qa-measurements
      htext(2,i,1) = text(0.01,0.40-(0.055*i), str2(i).name  ,'FontSize',fontsize, 'Interpreter','tex','Parent',ax);
      htext(2,i,2) = text(0.25,0.40-(0.055*i), str2(i).value ,'FontSize',fontsize, 'Interpreter','tex','Parent',ax);
    end
    % qa-scala
    %htext(5,1,1) = text(0.01,0.45-(0.055*(i+2)),str4(1).name,'FontSize',fontsize, 'Interpreter','tex','Parent',ax);
    for i=1:size(str3,2)  % subject-measurements
      htext(3,i,1) = text(0.51,0.40-(0.055*i), str3(i).name  ,'FontSize',fontsize, 'Interpreter','tex','Parent',ax);
      htext(3,i,2) = text(0.80,0.40-(0.055*i), str3(i).value ,'FontSize',fontsize, 'Interpreter','tex','Parent',ax);
    end



    pos = [0.01 0.38 0.48 0.36; 0.51 0.38 0.48 0.36; ...
           0.01 0.01 0.48 0.36; 0.51 0.01 0.48 0.36];
    spm_orthviews('Reset');




    % BB box is not optimal for all images
    disptype = 'affine'; 
    switch disptype
      case 'affine'
        dispmat = res.Affine; 
        spm_orthviews('BB', job.extopts.bb*0.95 );
      case 'ridid'
        % this does not work so good... AC has a little offset ...
        aff = spm_imatrix(res.Affine);  scale = aff(7:9); 
        spm_orthviews('BB', job.extopts.bb ./ mean(scale));
        dispmat = R; 
    end


    %% Yo - original image in original space
    % using of SPM peak values didn't work in some cases (5-10%), so we have to load the image and estimate the WM intensity 
    try %#ok<TRYNC>
      Yo  = single(VT.private.dat(:,:,:)); 
    end
    Yp0 = zeros(d,'single'); Yp0(indx,indy,indz) = single(Yp0b)*3/255; 
    if exist('Yo','var')
      if job.inv_weighting
        WMth = min([...
          cat_stat_nanmedian(Yo(Yp0(:)>0.8 & Yp0(:)<1.2))*2,...
          cat_stat_nanmedian(Yo(Yp0(:)>1.8 & Yp0(:)<2.2))*1.5,...
          ]);
        T1txt = '*.nii (Original PD/T2)'; 
      else
        WMth = cat_stat_nanmedian(Yo(Yp0(:)>2.8 & Yp0(:)<3.2)); clear Yo; 
        T1txt = '*.nii (Original T1)'; 
      end
      if ~debug, clear Yo; end

      if isfield(res,'spmpp')
        VT0x = res.image0(1); 
      else
        VT0x = VT0;
      end
      VT0x.mat = dispmat * VT0x.mat; 
      hho = spm_orthviews('Image',VT0x,pos(1,:)); 
      spm_orthviews('Caption',hho,{T1txt},'FontSize',fontsize,'FontWeight','Bold');
      spm_orthviews('window',hho,[0 WMth*cmmax]); caxis([0,2]);
      cc{1} = axes('Position',[pos(1,1) + 0.30 0.38 0.02 0.15],'Parent',fg); image((60:-1:1)');

      if job.inv_weighting
        set(cc{1},'YTick',ytick,'YTickLabel',fliplr(yticklabeli),'XTickLabel','','XTick',[],'TickLength',[0 0],...
          'FontSize',fontsize,'FontWeight','Bold','YAxisLocation','right');
      else  
        set(cc{1},'YTick',ytick,'YTickLabel',fliplr(yticklabelo),'XTickLabel','','XTick',[],'TickLength',[0 0],...
          'FontSize',fontsize,'FontWeight','Bold','YAxisLocation','right');
      end
    else
      cat_io_cprintf('warn','WARNING: Can''t display original file "%s"!\n',VT.fname); 
    end


    % Ym - normalized image in original space
    if ~isfield(res,'spmpp') 
      %%
      Vm        = res.image(1); 
      Vm.fname  = ''; 
      Vm.dt     = [spm_type('FLOAT32') spm_platform('bigend')];
      Vm.dat(:,:,:) = single(Ym); 
      Vm.pinfo  = repmat([1;0],1,size(Ym,3));
      Vm.mat    = dispmat * Vm.mat; 
      hhm = spm_orthviews('Image',Vm,pos(2,:));
      spm_orthviews('Caption',hhm,{'m*.nii (Int. Norm.)'},'FontSize',fontsize,'FontWeight','Bold');
      spm_orthviews('window',hhm,[0 cmmax]); caxis([0,2]);
      cc{2} = axes('Position',[pos(2,1) + 0.30 0.38 0.02 0.15],'Parent',fg); image((60:-1:1)');
      set(cc{2},'YTick',ytick,'YTickLabel',fliplr(yticklabel),'XTickLabel','','XTick',[],'TickLength',[0 0],...
        'FontSize',fontsize,'FontWeight','Bold','YAxisLocation','right');
    end

    % Yo - segmentation in original space
    VO        = res.image(1); 
    VO.fname  = ''; 
    VO.dt     = [spm_type('FLOAT32') spm_platform('bigend')];
    VO.dat(:,:,:) = single(Yp0/3); 
    VO.pinfo  = repmat([1;0],1,size(Yp0,3));
    VO.mat    = dispmat * VO.mat; 
    hhp0 = spm_orthviews('Image',VO,pos(3,:)); if ~debug, clear Yp0; end
    spm_orthviews('Caption',hhp0,'p0*.nii (Segmentation)','FontSize',fontsize,'FontWeight','Bold');
    spm_orthviews('window',hhp0,[0 cmmax]); caxis([0,2]);
    cc{3} = axes('Position',[pos(3,1) + 0.30 0.02 0.02 0.15],'Parent',fg); image((60:-1:1)');
    set(cc{3},'YTick',ytick,'YTickLabel',fliplr(yticklabel),'XTickLabel','','XTick',[],'TickLength',[0 0],...
      'FontSize',fontsize,'FontWeight','Bold','YAxisLocation','right');
    spm_orthviews('Reposition',[0 0 0]); 


    % surface
    if job.extopts.print>1
      if exist('Psurf','var')
        try
          hCS = subplot('Position',[0.50 0.05 0.55 0.30],'visible','off'); 
          hSD = cat_surf_display(struct('data',Psurf(1).Pthick,'readsurf',0,'expert',2,...
            'multisurf',job.output.surface,'view','s',...
            'parent',hCS,'verb',0,'caxis',[0 6],'imgprint',struct('do',0)));
          colormap(cmap);  set(hSD{1}.colourbar,'visible','off'); 
          cc{3} = axes('Position',[0.63 0.02 0.3 0.01],'Parent',fg); image((121:1:120+surfcolors));
          set(cc{3},'XTick',1:(surfcolors-1)/6:surfcolors,'XTickLabel',{'0','1','2','3','4','5','          6 mm'},...
            'YTickLabel','','YTick',[],'TickLength',[0 0],'FontSize',fontsize,'FontWeight','Bold');
        catch
          cat_io_cprintf('warn','WARNING: Can''t display surface!\n',VT.fname);   
        end
      end
    end



    %% print subject report file as standard PDF/PNG/... file
    job.imgprint.type  = 'pdf';
    job.imgprint.dpi   = 600;
    job.imgprint.fdpi  = @(x) ['-r' num2str(x)];
    job.imgprint.ftype = @(x) ['-d' num2str(x)];
    job.imgprint.fname     = fullfile(pth,reportfolder,['catreport_' nam '.' job.imgprint.type]); 

    fgold.PaperPositionMode = get(fg,'PaperPositionMode');
    fgold.PaperPosition     = get(fg,'PaperPosition');
    fgold.resize            = get(fg,'resize');

    % it is necessary to change some figure properties especialy the fontsizes 
    set(fg,'PaperPositionMode','auto','resize','on','PaperPosition',[0 0 1 1]);
    for hti = 1:numel(htext), if htext(hti)>0, set(htext(hti),'Fontsize',fontsize*0.8); end; end
    for hti = 1:numel(cc), set(cc{hti},'Fontsize',fontsize*0.8); end;
    print(fg, job.imgprint.ftype(job.imgprint.type), job.imgprint.fdpi(job.imgprint.dpi), job.imgprint.fname); 
    for hti = 1:numel(htext), if htext(hti)>0, set(htext(hti),'Fontsize',fontsize); end; end
    for hti = 1:numel(cc), set(cc{hti},'Fontsize',fontsize); end; 
    set(fg,'PaperPositionMode',fgold.PaperPositionMode,'resize',fgold.resize,'PaperPosition',fgold.PaperPosition);
    try
      fprintf('Print ''Graphics'' figure to: \n  %s\n',job.imgprint.fname);% windows error?
    end

    %% reset colormap to the simple SPM like gray60 colormap
    if exist('hSD','var')
      % if there is a surface than we have to use the gray colormap also here
      % because the colorbar change!
      try %#ok<TRYNC>
        cat_surf_render2('ColourMap',hSD{1}.axis,gray(128));
        cat_surf_render2('Clim',hSD{1}.axis,[0 6]);
        axes(cc{3}); image(0:60);
        set(cc{3},'XTick',max(1,0:10:60),'XTickLabel',{'0','1','2','3','4','5','          6 mm'},...
          'YTickLabel','','YTick',[],'TickLength',[0 0],'FontSize',fontsize,'FontWeight','Bold');
      end
    end

    WMfactor = 4/3; cmap = gray(60); colormap(cmap); caxis([0,numel(cmap)]); 

    % new colorscale
    if job.inv_weighting
      % if exist('hho' ,'var'), spm_orthviews('window',hho ,[0 T3th(3)*WMfactor]); end
    else
      if exist('hho' ,'var'), spm_orthviews('window',hho ,[0 T3th(3)*WMfactor]); end
    end
    if exist('hhm' ,'var'), spm_orthviews('window',hhm ,[0 WMfactor]); end
    if exist('hhp0','var'), spm_orthviews('window',hhp0,[0 WMfactor]); end

    warning on;  %#ok<WNON>


    %% command window output
    fprintf('\n%s',repmat('-',1,72));
    fprintf(1,'\nCAT preprocessing takes %0.0f minute(s) and %0.0f second(s).\n', ...
      floor(round(etime(clock,res.stime))/60),mod(round(etime(clock,res.stime)),60));
    cat_io_cprintf(color(QMC,qa.qualityratings.IQR), sprintf('Image Quality Rating (IQR):  %5.2f%%%% (%s)\n',...
      mark2rps(qa.qualityratings.IQR),mark2grad(qa.qualityratings.IQR)));

    % print subfolders
    if job.extopts.subfolders
      fprintf('Segmentations are saved in %s%s%s\n',pth,filesep,'mri');
      fprintf('Reports are saved in %s%s%s\n',pth,filesep,'report');
      if job.output.ROI
        fprintf('Labels are saved in %s%s%s\n',pth,filesep,'label');
      end
      if job.output.surface
        fprintf('Surface measurements are saved in %s%s%s\n',pth,filesep,'surf');
      end
    end

    fprintf('%s\n\n',repmat('-',1,72));
    
    % finish diary entry of "../report/cmdln_*.txt"
    % read diagry and add the command-line output to the *.xml and *.mat file
    diary off; 
    try
      fid=fopen(res.catlog);
      txt=fread(fid,200000,'uint8=>char');
      fclose(fid); 
      txt2=textscan(txt,'%s','Delimiter',''); 
      cat_io_xml(fullfile(pth,reportfolder,['cat_' nam '.xml']),struct(...
        'catlog',txt2),'write+'); % here we have to use the write+!
    end    
    

    clear C c Ymi Ymf Ym

end
%%
return;
%=======================================================================

%=======================================================================
function wYv = cat_vol_ROInorm(Yv,warped,ai,mod,FA)
% ----------------------------------------------------------------------
% normalized space:  
% ----------------------------------------------------------------------
% for normalized space no further adaptions are available, but 
% a masking based on the tissue map can be used
% ----------------------------------------------------------------------
 
  % load mask (and complete undefined parts)
  if isempty(Yv)
    % no input - load atlas
   
    if ~exist(FA{ai,1},'file')
      error('cat:cat_main:missAtlas','Miss cat atlas-file ''%s''!',FA{ai,1});
    end
    % try multiple times, because of read error in parallel processing
    for i=1:5
      try
        wVv = spm_vol(FA{ai,1});
        wYv = spm_read_vols(wVv);
        break
      catch 
        pause(0.5)
      end
    end
     
    % resample atlas, if the atlas resolution differs from the actual template resolution
    if wVv.mat(1) ~= warped.M1(1)
      wVv2 = wVv; wVv2.mat = warped.M1; wVv2.dim = warped.odim; 
      [t,wYv] = cat_vol_imcalc(wVv,wVv2,'i1',struct('interp',0,'verb',0));
    end
    wYv = cat_vol_ctype(wYv,wVv(1).private.dat.dtype);
  else
    % map image to atlas space
    for yi=1:numel(warped.ress), warped.y(:,:,:,yi) = warped.y(:,:,:,yi) * warped.ress(yi); end
    if mod==0
      [wYv,w] = spm_diffeo('push',Yv,warped.y,warped.odim(1:3)); spm_field('boundary',1);
      wYv = spm_field(w,wYv,[sqrt(sum(warped.M1(1:3,1:3).^2)) 1e-6 1e-4 0  3 2]);
    elseif mod==1 && iscell(Yv) % tissue case
      nicemapping = 1;
      if nicemapping
        % Modulation using spm_diffeo and push introduces aliasing artefacts,
        % thus we use the def2det function of the inverted deformations to obtain the old and 
        % in my view a more appropriate jacobian determinant 
        % The 2nd reason to use the old modulation is compatibility with cat_vol_defs.m
        Yy = spm_diffeo('invdef',warped.y,warped.odim,eye(4),warped.M0);
        w  = spm_diffeo('def2det',Yy)/det(warped.M0(1:3,1:3)); clear Yy;
        % ensure that jacobian det is positive (no clue why some times the sign is switched)
        if mean(w(~isnan(w))) < 0, w = -w; end 
        w(:,:,[1 end]) = NaN; w(:,[1 end],:) = NaN; w([1 end],:,:) = NaN;
      end
      
      wYv = cell(1,3);
      for i=1:3
        if nicemapping 
          [wYv{i},w2] = spm_diffeo('push',single(Yv{i})/255,warped.y,warped.odim(1:3)); 
          % divide by jacdet to get unmodulated data
          wYv{i} = wYv{i}./(w2+0.001); 
          wYv{i} = wYv{i} .* w;
        else
          % simple push
          [wYv{i},w] = spm_diffeo('push',single(Yv{i})/255,warped.y,warped.odim(1:3)); 
        end
      end
    else
      error('unknown case');
    end
    %spm_smooth(wYv,wYv,[1 1 1]); 
  end
  
  
return
%=======================================================================

%=======================================================================
function csv = cat_vol_ROIestimate(Yp0,Ya,Yv,ai,name,csv,tissue,FA)
% ----------------------------------------------------------------------
% estimate values
% ----------------------------------------------------------------------


% load atlas-csv-file

  [pp,ff] = fileparts(FA{ai,1});
  csvf = fullfile(pp,[ff '.csv']);

  if isempty(csv) 
    if exist(csvf,'file')
      csv = cat_io_csv(csvf,'','',struct('delimiter',';')); 
    else
      csv = [num2cell((1:max(Ya(:)))') ...
        cellstr([repmat('ROI',max(Ya(:)),1) num2str((1:max(Ya(:)))','%03d')]) ...
        cellstr([repmat('ROI',max(Ya(:)),1) num2str((1:max(Ya(:)))','%03d')])];
    end
    
    % remove empty rows and prepare structure names
    if size(csv,2)>2, csv(:,3:end)=[]; end
    for ri=size(csv,1):-1:1
      if isempty(csv{ri,1}) || isempty(csv{ri,2}); 
        csv(ri,:)=[];
      elseif csv{ri,1}==0
        csv(ri,:)=[];
      end       
    end
  end
  name = genvarname(strrep(strrep(name,'-','_'),' ','_'));
  
  
  %% volume case
  Yp0toC = @(Yp0,c) 1-min(1,abs(Yp0-c));   
  % other maps with masks
  for ti=1:numel(tissue)
    switch name(1)
      case 'V' % volume
        csv{1,end+1} = [name tissue{ti}];  %#ok<AGROW>
        for ri=2:size(csv,1)
          switch lower(tissue{ti})
            case 'csf',   Ymm=single(Yv{3}) .* single(Ya==csv{ri,1});
            case 'gm',    Ymm=single(Yv{1}) .* single(Ya==csv{ri,1});
            case 'wm',    Ymm=single(Yv{2}) .* single(Ya==csv{ri,1});
            case 'wmh',   Ymm=single(Yv{2}) .* single(Ya==csv{ri,1}); 
            case 'brain', Ymm=single(Yv{1} + Yv{2} + Yv{3}) .* single(Ya==csv{ri,1});
            case '',      Ymm=single(Ya==csv{ri,1});
          end
          csv{ri,end} = 1/1000 * cat_stat_nansum(Ymm(:));
        end
      case 'c'
        return
      otherwise % 
        csv{1,end+1} = [name tissue{ti}];  %#ok<AGROW>
        switch lower(tissue{ti})
          case 'csf',   Ymm=Yp0toC(Yp0,1); 
          case 'gm',    Ymm=Yp0toC(Yp0,2); 
          case 'wm',    Ymm=Yp0toC(Yp0,3); 
          case 'wmh',   Ymm=Yp0toC(Yp0,4); 
          case 'brain', Ymm=Yp0>0.5;
          case '',      Ymm=true(size(Yp0));
        end
        for ri=2:size(csv,1)
          csv{ri,end} = cat_stat_nanmean(Yv(Ya(:)==csv{ri,1} & Ymm(:)));
        end
    end
  end
  
return
%=======================================================================

%=======================================================================
function [x1,y1,z1] = defs(sol,z,MT,prm,x0,y0,z0,M)
iMT = inv(MT);
x1  = x0*iMT(1,1)+iMT(1,4);
y1  = y0*iMT(2,2)+iMT(2,4);
z1  = (z0(z)*iMT(3,3)+iMT(3,4))*ones(size(x1));
x1a = x0    + spm_bsplins(sol{1},x1,y1,z1,prm);
y1a = y0    + spm_bsplins(sol{2},x1,y1,z1,prm);
z1a = z0(z) + spm_bsplins(sol{3},x1,y1,z1,prm);
x1  = M(1,1)*x1a + M(1,2)*y1a + M(1,3)*z1a + M(1,4);
y1  = M(2,1)*x1a + M(2,2)*y1a + M(2,3)*z1a + M(2,4);
z1  = M(3,1)*x1a + M(3,2)*y1a + M(3,3)*z1a + M(3,4);
return;
%=======================================================================

%=======================================================================
function t = transf(B1,B2,B3,T)
if ~isempty(T)
    d2 = [size(T) 1];
    t1 = reshape(reshape(T, d2(1)*d2(2),d2(3))*B3', d2(1), d2(2));
    t  = B1*t1*B2';
else
    t = zeros(size(B1,1),size(B2,1),size(B3,1));
end;
return;
%=======================================================================

%=======================================================================
function p = likelihoods(f,bf,mg,mn,vr)
K  = numel(mg);
N  = numel(f);
M  = numel(f{1});
cr = zeros(M,N);
for n=1:N,
    if isempty(bf),
        cr(:,n) = double(f{n}(:));
    else
        cr(:,n) = double(f{n}(:).*bf{n}(:));
    end
end
p  = ones(numel(f{1}),K);
for k=1:K,
    amp    = mg(k)/sqrt((2*pi)^N * det(vr(:,:,k)));
    d      = bsxfun(@minus,cr,mn(:,k)')*inv(chol(vr(:,:,k))); %#ok<MINV>
    p(:,k) = amp*exp(-0.5*sum(d.*d,2)) + eps;
end
return;
%==========================================================================

%==========================================================================
% function [P] = clean_gwc(P,level)
%==========================================================================
function [P] = clean_gwc(P,level)
if nargin<4, level = 1; end

b    = P(:,:,:,2);

% Build a 3x3x3 seperable smoothing kernel
%--------------------------------------------------------------------------
kx=[0.75 1 0.75];
ky=[0.75 1 0.75];
kz=[0.75 1 0.75];
sm=sum(kron(kron(kz,ky),kx))^(1/3);
kx=kx/sm; ky=ky/sm; kz=kz/sm;

th1 = 0.15;
if level==2, th1 = 0.2; end
% Erosions and conditional dilations
%--------------------------------------------------------------------------
niter  = 32;
niter2 = 32;
spm_progress_bar('Init',niter+niter2,'Extracting Brain','Iterations completed');
for j=1:niter
    if j>2, th=th1; else th=0.6; end  % Dilate after two its of erosion
    for i=1:size(b,3)
        gp       = single(P(:,:,i,1));
        wp       = single(P(:,:,i,2));
        bp       = single(b(:,:,i))/255;
        bp       = (bp>th).*(wp+gp);
        b(:,:,i) = cat_vol_ctype(round(bp));
    end
    spm_conv_vol(b,b,kx,ky,kz,-[1 1 1]);
    spm_progress_bar('Set',j);
end

% Also clean up the CSF.
if niter2 > 0,
    c = b;
    for j=1:niter2
        for i=1:size(b,3)
            gp       = single(P(:,:,i,1));
            wp       = single(P(:,:,i,2));
            cp       = single(P(:,:,i,3));
            bp       = single(c(:,:,i))/255;
            bp       = (bp>th).*(wp+gp+cp);
            c(:,:,i) = cat_vol_ctype(round(bp));
        end
        spm_conv_vol(c,c,kx,ky,kz,-[1 1 1]);
        spm_progress_bar('Set',j+niter);
    end
end

th = 0.05;

for i=1:size(b,3)
    slices = cell(1,size(P,4));
    for k1=1:size(P,4),
        slices{k1} = single(P(:,:,i,k1))/255;
    end
    bp        = single(b(:,:,i))/255;
    bp        = ((bp>th).*(slices{1}+slices{2}))>th;
    slices{1} = slices{1}.*bp;
    slices{2} = slices{2}.*bp;

    if niter2>0,
        cp        = single(c(:,:,i))/255;
        cp        = ((cp>th).*(slices{1}+slices{2}+slices{3}))>th;
        slices{3} = slices{3}.*cp;
    end
    if numel(slices)>=5
      slices{5} = slices{5}+1e-4; % Add a little to the soft tissue class
    end
    tot       = zeros(size(bp))+eps;
    for k1=1:size(P,4),
        tot   = tot + slices{k1};
    end
    for k1=1:size(P,4),
        P(:,:,i,k1) = cat_vol_ctype(round(slices{k1}./tot*255));
    end 
end

spm_progress_bar('Clear');
return;
%==========================================================================

