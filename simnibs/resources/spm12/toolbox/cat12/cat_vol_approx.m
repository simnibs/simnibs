function TA=cat_vol_approx(T,method,vx_vol,res,opt)
% Approximation of missing values
% ______________________________________________________________________
% Approximation of missing values (nan's). First, a nearest neigbhor 
% approximation is used. After that all values within the convex hull 
% were corrected with a laplace filter. Depending on the 'method'
% variable the outside hull area is further improved for
% method=='nn'|'linear', but not for 'nh'. 
% Use a resolution 'res' similar to the voxel size for finer results 
% (i.e. 2 mm) or smaller 'res' for smoother images (default 4 mm).
% 
% TA = cat_vol_approx(T,method,vx_vol,res[,opt])
% 
% T       input image
% TA      output image (single)
% method  ['nh' | 'nn' | 'linear' | 'spm']
%         nh:     fast default method
%         nn:     fast improved method (additional update of the outside 
%                 hull area)
%         linear  slower improved method
% vx_vol  voxel resolution of T (default: 1 mm)
% res     voxel resolution for approximation (default: 4 mm)
% opt     further options for development and test 
%   .lfI  laplace filter stop criteria for the input  image
%   .lfI  laplace filter stop criteria for the output image
%   .hull use hull approximation (default: 1)
%
% Examples:
%   There is a cell mode test part in the file...%
%
% ______________________________________________________________________
% Robert Dahnke
% Structural Brain Mapping Group
% University Jena
%
% $Id: cat_vol_approx.m 1159 2017-08-04 14:08:14Z dahnke $

  if ~exist('res','var'); res=4; end
  if ~exist('vx_vol','var'); vx_vol=ones(1,3); end
  if ~exist('method','var'); method='nn'; end
  
  if ~exist('opt','var'), opt=struct(); end
  def.lfI  = 0.40;
  def.lfO  = 0.40;
  def.hull = 1;
  opt = cat_io_checkinopt(opt,def);
  opt.lfO = min(10,max(0.0001,opt.lfO));
  
  T(isnan(T) | isinf(T))=0; 
  maxT = max(T(T(:)<inf & ~isnan(T(:))));
  T = single(T/max(eps,maxT));
  
  if 0 % T(isnan(T))=0;
    % outlier removal ... not yet
    meanT = mean(T(T(:)>0));
    stdT  = std(T(T(:)>0));
    T(T>meanT + 2*stdT | T<meanT - 2*stdT)=0; 
  end

  [Tr,resTr]    = cat_vol_resize(T,'reduceV',vx_vol,res,16,'meanm');
  %strcmp(method,'linear') || 0 %
  if (opt.hull || strcmp(method,'linear')) && ~strcmp(method,'spm')
    [Brr,resTrr] = cat_vol_resize(Tr>0,'reduceV',resTr.vx_volr,16,16,'max');
    BMrr = cat_vol_morph(Brr>0,'distclose',20)>0;
    BMr  = cat_vol_resize(BMrr,'dereduceV',resTrr); 
  
    % inside hull approximation ...
    [MDr,MIr]  = cat_vbdist(single(Tr>0),Tr==0 | isnan(Tr),double(resTr.vx_volr)); 
    TAr=Tr(MIr); TAr(Tr>0) = Tr(Tr>0); 
    if opt.lfO >= 0.5
      meanTAr = cat_stat_nanmedian(TAr(Tr(:)>0));
      TAr     = TAr / meanTAr; 
      Ygr     = cat_vol_grad(TAr); 
      opt.lfO = min( 0.49 , max( 0.0001 , min(  mean(resTr.vx_volr)/10 , median(Ygr(Tr(:)>0)) /opt.lfO ))); 
      TAr     = cat_vol_laplace3R(TAr,true(size(TAr)),double(opt.lfO)) * meanTAr; 
    else
      TASr=cat_vol_smooth3X(TAr,2); TAr(~BMr)=TASr(~BMr); clear TASr; 
      opt.lfO = min(0.49,max(0.0001,opt.lfO));
      TAr = cat_vol_laplace3R(TAr,BMr & ~Tr,opt.lfO); TAr = cat_vol_median3(TAr); %,Tr>0,Tr>0,0.05); 
      %TAr = cat_vol_laplace3R(TAr,Tr>0,opt.lfI); 
      TAr = cat_vol_laplace3R(TAr,BMr & ~Tr,opt.lfO);
    end
  else
    TAr = Tr; 
    BMr = Tr>0; 
  end
  
  %ds('l2','',vx_vol,Tr,BMr,Tr/mean(Tr(Tr>0)),TAr/mean(Tr(Tr>0)),80)
  switch method
    case 'nh'
    case 'nn'
      TAr  = TAr .* (BMr | Tr);
      [MDr,MIr]  = cat_vbdist(single(TAr>0),TAr==0,double(resTr.vx_volr)); 
      TAr=TAr(MIr); TASr=cat_vol_smooth3X(TAr,4); TAr(~BMr)=TASr(~BMr);  clear TASr; 
      TAr = cat_vol_laplace3R(TAr,~BMr,double(opt.lfO)); TAr = cat_vol_median3(TAr,~BMr);
      TAr = cat_vol_laplace3R(TAr,~Tr,double(opt.lfO)); 
    case 'linear'
      TNr = TAr;
      Tr  = TAr .* BMr;
      % outside hull linear approximation ...
      vx_voln = resTr.vx_vol./mean(resTr.vx_vol);  
      [MDFr,EIFr] = cat_vbdist(single(cat_vol_morph(BMr>0,'disterode',max(3,8/res))),true(size(Tr)),vx_voln);  
      [MDNr,EINr] = cat_vbdist(single(cat_vol_morph(BMr>0,'disterode',max(1,6/res))),true(size(Tr)),vx_voln); 
      TAr = Tr; TAr(~Tr) = Tr(EINr(~Tr)) + ( (Tr(EINr(~Tr))-Tr(EIFr(~Tr))) ./ max(eps,( (MDFr(~Tr)-MDNr(~Tr))./MDFr(~Tr)) )); TAr(1)=TAr(2);
      % correction and smoothing
      TAr = min(max(TAr,TNr/2),TNr*2); % /2
      TAr = cat_vol_median3(TAr,~BMr); TAr=TAr(MIr); TASr=cat_vol_smooth3X(TAr,1); 
      TAr(~BMr)=TASr(~BMr); clear TASr; 
      TAr = cat_vol_laplace3R(TAr,~BMr,opt.lfO); TAr = cat_vol_median3(TAr,~BMr); 
      TAr = cat_vol_laplace3R(TAr,~BMr,opt.lfO); 
    case 'spm'
      fname     = fullfile(tempdir,'approxtst.nii');
      N         = nifti;
      N.dat     = file_array(fname,size(Tr),[spm_type('float32') spm_platform('bigend')],0,1,0);
      N.mat     = [eye(4,3), [size(Tr)'/2;1]]; N.mat([1,6,11])=resTr.vx_volr;
      N.mat0    = N.mat;
      N.descrip = 'approx test';
      create(N);
      N.dat(:,:,:) = double(Tr);

      fnameB    = fullfile(tempdir,'approxtstB.nii');
      B         = nifti;
      B.dat     = file_array(fnameB,size(Tr),[spm_type('float32') spm_platform('bigend')],0,1,0);
      B.mat     = [eye(4,3), [size(Tr)'/2;1]]; N.mat([1,6,11])=resTr.vx_volr;
      B.mat0    = B.mat;
      B.descrip = 'approx test B';
      create(B);
      B.dat(:,:,:) = ones(size(Tr));

      V  = spm_vol(fname); 
      VB = spm_vol(fnameB); 
    
    
      % SPM bias correction
      bT  = spm_bias_estimate(V,struct('nbins',256,'reg',0.001,'cutoff',30));
      VA  = spm_bias_apply(V ,bT); TA = spm_read_vols(VA); 
      VAr = spm_bias_apply(VB,bT); TAr = 1/max(eps,single(spm_read_vols(VAr)) * mean(TA(Tr(:)>0)));
%      ds('d2','',[1 1 1],PT{1},PT{2}/max(PT{2}(:)),TA,TAr,20)
      
      delete(fname,fnameB);
      
  end
  %{
  [TArr,BMrr,resTrr2] = cat_vol_resize({TAr,BMr},'reduceV',vx_vol,8,8,'linear');
  TArr = cat_vol_median3(TArr,~BMrr,~BMrr); 
  for i=1:round(8/bias),
    if mod(i,4)==0, TASrr=cat_vol_smooth3X(TArr,8/mean(resTr.vx_volr)); TArr(~BMrr)=TASrr(~BMrr); end;
  end
  TAr  = cat_vol_resize(TArr,'dereduceV',resTrr2); 
  %}
  %TAr = cat_vol_smooth3X(TAr,8/mean(resTr.vx_volr)*0.1/bias).^1.05; 
  
  TA  = cat_vol_resize(TAr,'dereduceV',resTr);
  TA  = TA*maxT;
end
function cat_tst_pre_approx
  %%
  PTsize  = repmat(64,1,3);
  PTrange = [0.1,2.4];
  vx_vol  = round(256./PTsize);
  
  % Bias type:
  % --------------------------------------------------------------------
  % linear bias
  PT{1} = zeros(PTsize,'single');
  for i=1:PTsize(2)
    PT{1}(:,i,:) = PTrange(1) + (i*diff(PTrange)/(size(PT{1},2))) * ones(PTsize(1),1,PTsize(3)); 
  end
  % circle bias
  PT{2} = zeros(PTsize,'single'); PT{2}(round(PTsize(1)*3/7),round(PTsize(2)*3/7),round(PTsize(3)*3/7)) = 1; 
  PT{2} = cat_vbdist(PT{2}); PT{2} = max(0,PTrange(2) - (diff(PTrange)*(PT{2}/max(PT{2}(PT{2}<inf)))));
  
  
  
  % Mask type:
  % --------------------------------------------------------------------
  PM{1} = zeros(PTsize,'single'); PM{1}(round(PTsize(1)/2),round(PTsize(2)/2),round(PTsize(3)/2)) = 1; 
  PM{1} = (cat_vbdist(PM{1},true(size(PM{1})),[1.5 1 1]) + ...
          cat_vol_smooth3X(rand(PTsize)*20,2) + cat_vol_smooth3X(randn(PTsize)*20)) <PTsize(1)/2; 
  
  ds('d2','',[1 1 1],PT{1},PT{2},PT{1}.*PM{1},PT{2}.*PM{1},32)

  
  %%
  % --------------------------------------------------------------------
  %method = {'val','nn','linear'} 
  
  for i=1:numel(PT)
    PN{i} = PT{i} + cat_vol_smooth3X(randn(PTsize)*0.4,2) + cat_vol_smooth3X(randn(PTsize)*0.2); 
  end
  
  
  %for 
  %=cat_vol_approx(T,method,vx_vol,res,bias)

  %%
    T = PT{1}.*PM{1}; %T(T==0)=nan;
  
    fname     = fullfile(tempdir,'approxtst.nii');
    N         = nifti;
    N.dat     = file_array(fname,PTsize,[spm_type('float32') spm_platform('bigend')],0,1,0);
    N.mat     = [eye(4,3)*vx_vol(1), [PTsize'/2;1]];
    N.mat0    = N.mat;
    N.descrip = 'approx test';
    create(N);
    N.dat(:,:,:) = double(T);

    fnameB    = fullfile(tempdir,'approxtstB.nii');
    B         = nifti;
    B.dat     = file_array(fnameB,PTsize,[spm_type('float32') spm_platform('bigend')],0,1,0);
    B.mat     = [eye(4,3)*vx_vol(1), [PTsize'/2;1]];
    B.mat0    = B.mat;
    B.descrip = 'approx test B';
    create(B);
    B.dat(:,:,:) = ones(PTsize);

    
    V  = spm_vol(fname); 
    VB = spm_vol(fnameB); 
    
 %% SPM spline interpolation  
    intp = 1; warp = 0;
    i = find(isnan(T)); [x,y,z]=ind2sub(PTsize,i);
    c = spm_bsplinc(V,repmat(intp,3,1));
    v = spm_bsplins(c,x,y,z,[repmat(intp,3,1) repmat(warp,3,1)]);
    TA=T; TA(i)=v;
    fprintf('done\n');
    
 %% SPM bias correction
    
    bT  = spm_bias_estimate(V,struct('nbins',256,'reg',0.001,'cutoff',30));
    VA  = spm_bias_apply(V ,bT); TA = spm_read_vols(VA); 
    VAB = spm_bias_apply(VB,bT); TAB = spm_read_vols(VAB);
    ds('d2','',[1 1 1],PT{1},PT{2}/max(PT{2}(:)),TA,TAB,20)
    
end


