function cat_vol_urqio(job)
% Ultrahigh Resolution Quantitative Image Optimization
% ______________________________________________________________________
%
% Skript to reduce inhomogeneity, noise, and blood vessels in ultra-high
% resolution quantitative MR data. The R1, PD (A), and R2s images are 
% required. The function use a simple tissue classification to apply a
% bias correction, a gobal intensity normalization, a blood vessel 
% correction, and finally a noise reduction. 
% 
% WARNING: This function is in an early development stage (alpha)
%
% cat_vol_urqio(job)
%
%  job.data .. structure with the filenames of the input images
%   .r1   .. R1 images
%   .pd   .. PD images (A)
%   .r2s  .. R2s images
%  job.ouput .. structure to control writing of output images
%   .pd  .. write corrected PD images
%   .r1  .. write corrected R1 images
%   .r2s .. write corrected R2s images
%   .t1  .. write synthetic T1 images (invertation of the PD)
%   .bv  .. write detected blood vessels
%  job.opts .. structure of option parameter
%   .prefix .. filename prefix (default = 1)
%   .verb   .. display processing notes
%   .bc     .. apply bias correction [0|1]
%   .in     .. apply global intensity normalisation [0|1]
%   .bvc    .. apply blood vessel correction [0|1]
%   .nc     .. apply noise correction
%
% ______________________________________________________________________
% Robert Dahnke
% Structural Brain Mapping Group
% University Jena
%
% $Id: cat_vol_urqio.m 1245 2017-12-12 10:35:00Z dahnke $


  dbs   = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end

  % default options
  if ~exist('job','var'), job = struct(); end
  % input data
  def.data.r1   = {};
  def.data.pd   = {};
  def.data.r2s  = {};
  % output data
  def.output.pd  = 1; 
  def.output.t1  = 1; 
  def.output.r1  = 1;
  def.output.r2s = 1; 
  def.output.bv  = 1; 
  % parameter
  def.opts.prefix = 'catsyn_';
  def.opts.verb = 2;
  def.opts.bc   = 1; 
  def.opts.in   = 1; 
  def.opts.bvc  = 1; 
  def.opts.nc   = 1; 
  def.opts.spm  = 1;  % use SPM preprocessing
  def.opts.spmres = 0.8;
  def.opts.ss   = 0;  % apply skull-stripping
  def.opts.SSbd = 20; % reduce the size of the image by removing parts with distance greater SSbd  
  % checkin
  job = cat_io_checkinopt(job,def); 
  
  % set empty data fields
  if isempty(job.data.r1) || isempty(job.data.r1{1})
    job.data.r1 = cellstr(spm_select(1,'image','Select R1 files'));
  end
  if isempty(job.data.pd) || isempty(job.data.pd{1})
    job.data.pd = cellstr(spm_select(1,'image','Select PD files'));
  end
  if  isempty(job.data.r2s) || isempty(job.data.r2s{1})
    job.data.r2s = cellstr(spm_select(1,'image','Select R2s files'));
  end
  job.data.r1   = cellstr(job.data.r1); 
  job.data.r2s  = cellstr(job.data.r2s); 
  job.data.pd   = cellstr(job.data.pd); 
  
  si = 1; %#ok<NASGU>
  
  
  
  %% main loop
  for si=1:numel(job.data.r1)
    stime2 = cat_io_cmd(sprintf('Preprocessing subject %d/%d:',si,numel(job.data.r1)),'','',job.opts.verb); fprintf('\n'); 
   
    % get header
    Vd  = spm_vol(job.data.r1{si});
    Vr  = spm_vol(job.data.r2s{si}); 
    Vp  = spm_vol(job.data.pd{si});

    % voxel size
    vx_vol = sqrt(sum(Vd.mat(1:3,1:3).^2));

    
    % SPM unified segmentation 
    if job.opts.spm
      %%
      stime = cat_io_cmd('  SPM preprocessing','g5');

      % SPM segmentation parameter
      ofname = spm_file(job.data.pd{si},'number',''); 
      nfname = spm_file(ofname,'prefix','n'); 
      preproc.channel.vols{1,1}   = nfname; 
      preproc.channel.biasreg     = 0.001;
      preproc.channel.biasfwhm    = 60;
      preproc.channel.write       = [0 1];
      for ti=1:6
        preproc.tissue(ti).tpm    = {[char(cat_get_defaults('opts.tpm')) ',' num2str(ti)]};
        preproc.tissue(ti).ngaus  = 3;
        preproc.tissue(ti).native = [0 0];                           % native, dartel
        preproc.tissue(ti).warped = [0 0];                           % unmod, mod
      end
      preproc.warp.mrf            = 1;                               % 
      preproc.warp.cleanup        = 1;                               % this is faster and give better results!
      preproc.warp.reg            = [0 0.001 0.5 0.05 0.2];
      preproc.warp.affreg         = 'none'; 
      preproc.warp.write          = [0 0];
      preproc.warp.samp           = 3; 
      preproc.Yclsout             = true(1,6); 
      
      %% lower resolution to improve SPM processing time and use smoothing for denoising
      copyfile(ofname,nfname);
      Vi        = spm_vol(nfname); 
      Vi        = rmfield(Vi,'private'); Vn = Vi; 
      imat      = spm_imatrix(Vi.mat); 
      Vi.dim    = round(Vi.dim .* vx_vol./repmat(job.opts.spmres,1,3));
      imat(7:9) = repmat(job.opts.spmres,1,3) .* sign(imat(7:9));
      Vi.mat    = spm_matrix(imat);

      spm_smooth(nfname,nfname,repmat(0.75,1,3)); % denoising
      [Vi,Yi] = cat_vol_imcalc(Vn,Vi,'i1',struct('interp',1,'verb',0));
      delete(nfname); spm_write_vol(Vi,Yi); 
      
      vout   = cat_spm_preproc_run(preproc,'run');
      %spmmat = open(fullfile(pp,mrifolder,['n' ff '_seg8.mat'])); 
      %AffineSPM = spmmat.Affine; 
      
      %% update Ycls
      if isfield(vout,'Ycls')
        for i=1:6
          Vn  = spm_vol(nfname); Vn = rmfield(Vn,'private');
          Vo  = spm_vol(ofname); Vo = rmfield(Vo,'private'); 
          Vn.pinfo = repmat([1;0],1,size(vout.Ycls{i},3));
          Vo.pinfo = repmat([1;0],1,Vo.dim(3));
          Vn.dt    = [2 0];
          Vo.dt    = [2 0]; 
          Vo.dat   = zeros(Vo.dim(1:3),'uint8');
          if ~isempty(vout.Ycls{i})
            Vn.dat   = vout.Ycls{i}; 
            [Vmx,vout.Ycls{i}]  = cat_vol_imcalc([Vo,Vn],Vo,'i2',struct('interp',1,'verb',0));
            vout.Ycls{i} = cat_vol_ctype(vout.Ycls{i}); 
          end
        end
      end
      %%
      Yp0 = (2*single(vout.Ycls{1})+3*single(vout.Ycls{2})+single(vout.Ycls{3}))/255;
      Ybm = cat_vol_morph(Yp0>0.5,'lc');
      Ybg = cat_vol_morph(~cat_vol_morph(vout.Ycls{6}<64,'lo',1),'lo',1);
            
      %if ~debug, clear vout; end

      %% load images
      Yd = single(spm_read_vols(Vd));
      Yr = single(spm_read_vols(Vr));
      Yp = single(spm_read_vols(Vp));
      
      %% correct errors in Yd 
      YM = Yd>max(Yd(:)*0.99) & Yp<mean(Yp(:)) & vout.Ycls{6}>64; 
      YM = cat_vol_morph(YM,'d',2); 
      Yd(YM) = 0; 
      
      %% display segmentation V2
      %ds('d2','a',vx_vol,Yd/dth,Yd/dth .* ~YM,Yr./rth,Yp/pth,120); colormap gray

    else
      stime  = cat_io_cmd('  Load images:','g5','',job.opts.verb-1);
      % load images
      Yd  = single(spm_read_vols(Vd));
      Yr  = single(spm_read_vols(Vr));
      Yp  = single(spm_read_vols(Vp));
      Ybm = ~(Yd==0 & Yr==0 & Yp==0); 
      Ybg = (Yd==0 & Yr==0 & Yp==0); 
    end
    
    if 0 % gcut>0 && job.opts.spm
      Yl1 = ones(size(Ybm)); 
      
      YMF = cat_vol_morph(Yp0==1 & cat_vol_smooth3X(Yp0,8)>1,'l',[0.4,2]); % ventricle
      cat_main_gcut(Yr,Ybm,vout.Ycls,Yl1,YMF,vx_vol,opt)
    end
    
    
    
    %% some kind of brain masking
    Ypg  = cat_vol_grad(Yp)./Yp;
    Ydg  = cat_vol_grad(Yd)./Yd;
    pgth = max( cat_stat_nanmean(Ypg(Ypg(:)>0.01 & Ypg(:)<0.7 & Ybm(:)) ) , ...
            abs(cat_stat_nanmean(Ypg(Ypg(:)>0.01 & Ypg(:)<0.7 & Ybm(:)) )) );
    pgth = min(0.5,max(0.25,pgth)); 
    if job.opts.spm
      Ytis = Ybm & Ypg>0 & Ypg<pgth*4 & Ydg<1 & Ydg>0;
      Ytis = cat_vol_morph(Ytis,'lc',2/mean(vx_vol)) & Ybm & Ypg>0 & Ypg<pgth*4 & Ydg<pgth*4 & Ydg>0;
    else
      Ytis = Ybm & Ypg>0 & Ypg<pgth*4 & Ydg<pgth*4 & Ydg>0;
      Ytis = cat_vol_morph(cat_vol_morph(Ytis,'lo',1),'lc',10/mean(vx_vol)) & Ypg>0 & Ypg<pgth*4 & Ydg<pgth*4 & Ydg>0;
    end
    
    %% some simple thresholds 
    dth = max( mean(Yd(Ytis(:))) , abs(mean(Yd(Ytis(:)))));
    rth = max( mean(Yr(Ytis(:))) , abs(mean(Yr(Ytis(:)))));
    pth = max( mean(Yp(Ytis(:))) , abs(mean(Yp(Ytis(:)))));
    if ~debug, clear Ytis; end
    
    %  gradient maps and thresholds
    %  only use PD contrast that have best SNR here
    Ypg  = cat_vol_grad(Yp);
    Ypd  = cat_vol_div(Yp);
    pgth = max( mean(Ypg(Ypg(:)>0)) , abs(mean(Ypg(Ypg(:)<0))));
    
    
        
    
    %% tissue classification
    %  --------------------------------------------------------------------
    %  Yw = white matter map
    %  Yg = gray matter map
    %  Yc = cerebrospinal fluid map
    %  --------------------------------------------------------------------
    stime  = cat_io_cmd('  Tissue classification:','g5','',job.opts.verb-1,stime);
    
    % WM
    Yw = Yp<pth*1 & Yp>pth*0.6  & Yd>dth & Yd<dth*3.5  &  Yr>rth/2 & Yr<rth*2.5  &  Ybm & ...
         Ypg<pgth*4 & Ypd>-pgth/2 & Ypd<pgth/2; 
    if job.opts.spm
      Yw = Yw | (Yp0>2.9  & Yp<pth*1 & Yp>pth*0.6 & Ypg<pgth*4 & Ypd>-pgth/4 & Ypd<pgth/4); 
      Yw = Yw | (Yp0>=2 & Yp<pth*1 & Yp>pth*0.6 & Ypg>pgth/8 & Ypd>0.01); 
    end 
    Yw = Yw & Yr<mean(Yr(Yw(:)))*1.5; % remove
    Yw(smooth3(Yw)<0.4)=0; % small open
    Yw = cat_vol_morph(Yw,'l');
    
    % GM
    Yg = Yp>pth & Yp<pth*2  &  Yd<dth*1.6 & Yd>dth*0.5  &  Yr>rth/2 &  Ybm &  ...  
         Ypg<pgth*2 & Ypd>-pgth/2 & Ypd<pgth/2  & ~Yw; 
    if job.opts.spm, Yg = Yg | (Yp>pth & Yp0>1.9 & Yp0<2.1 & Ypg<pgth*2 & Ypd>-pgth/2 & Ypd<pgth/2 & ~Yw); end        
    Yg(smooth3(Yg)<0.4)=0; % small open 
    
    % CSF ... have to use it ...
    Yc = Yp>pth*1.5 & Yp<pth*3  &  Yd<dth*0.5 & Yd>dth*0.01  &  Yr<rth/2 &  Ybm &  ...
         Ypg<pgth*8 & Ypd>-pgth*2  & ~Yw & ~Yg; 
    if job.opts.spm
      Yc = Yc | (Yp>pth*1.5 & Yp<pth*3  &  Yd<dth*0.5 & Yd>dth*0.01 & ...
        Yr<rth/2 & Yr<rth*2.5 & ~Yw & ~Yg & ...
        Yp0>0.9 & Yp0<1.1 & Ypg<pgth*4 & Ypd>-pgth/42 & Ypd<pgth/4);
    end        
    Yc(smooth3(Yc)<0.4)=0; % small open 
    
    %
    if job.opts.spm
      Yh = Yp>pth/8 & Yp<pth*2  &  Yr>rth & Yr<rth*12  & ~Ybm & Ypg<pgth*4 & ...  
          cat_vol_smooth3X(Ybm,4/mean(vx_vol))<0.1 & Ypd>-pgth & Ypd<pgth & vout.Ycls{5}>64; 
      Yh(smooth3(Yh)<0.5)=0;
      Yh = cat_vol_morph(Yh,'l',[0.02 50]); 
    end
    
    
    
    %% super WM (higher intensity in R1 and R2s) and the question if and how 
    % we need to correction it for normal preprocessing, because actual it 
    % lead to GM like brainstem that get lost ...
    %if debug, Yx = Yp; end  
    if ~debug, clear Ypd; end
    
    if 0
      %% display segmentation V1
      ds('l2','',vx_vol,Yp./mean(Yp(Yw(:))),(Yw*3 + Yg*2 + Yc)/3,Yp./mean(Yp(Yw(:))),(Yw*3 + Yg*2 + Yc)/3*2,360); colormap gray
      %% display segmentation V2
      ds('d2','a',vx_vol,Yd/dth,Yh,Yr./rth,Yp/pth,120); colormap gray
    end
    
        
    
    
    %% bias correction
    %  --------------------------------------------------------------------
    %  Ybiw  = initial signal intensity map (of the white matter)
    %  Ybig  = intitla signal intensity map of the gray matter
    %  Ybiws = approximated bias field map
    %          normalization in the next step
    %  --------------------------------------------------------------------
    if job.opts.bc
      useGMandCSF = [1 0]; % use also the GM (and CSF) for bias correction (CSF is problematic) 

      % approx .. 2 is fast and smooth, 1 sh
      biasstr = { % name var WMfilter CSFfitler approxres approxsmooth
        'PD' ,'Yp', 2, 3, 1.6, 2/job.opts.bc;
        'R1' ,'Yd', 3, 2, 1.6, 2/job.opts.bc;
        'R2s','Yr', 3, 2, 1.6, 1/job.opts.bc;
      };
      for bi=1:size(biasstr,1)
        stime  = cat_io_cmd(sprintf('  Bias correction (%s):',biasstr{bi,1}),'g5','',job.opts.verb-1,stime);
        eval(['Ybi = ' biasstr{bi,2} '; bith = ' biasstr{bi,2}(2) 'th;']); 
        Ybiw = Ybi .* Yw; Ybiw = cat_vol_localstat(Ybiw,Ybiw>0,1,1); 
        Ybiw = Ybi .* (Yw & Ybiw~=Ybi & abs(Ybiw-Ybi)<bith*0.2); 
        Ybiw = cat_vol_localstat(Ybiw,Ybiw>0,1,biasstr{bi,3}); 
        for ii=1:round(8/mean(vx_vol)), Ybiw = cat_vol_localstat(Ybiw,Ybiw>0,1,1); end
        if useGMandCSF(1)
          % use additonal GM
          Ybig = Ybi .* Yg; Ybig = cat_vol_localstat(Ybig,Ybig>0,1,1); 
          Ybig = Ybi .* (Yg & Ybig~=Ybi & abs(Ybig-Ybi)<bith*0.1); 
          for ii=1:round(8/mean(vx_vol)), Ybig = cat_vol_localstat(Ybig,Ybig>0,1,1); end
          Ybiw = Ybiw + (Ybig / median(Ybig(Ybig(:)>0)) * median(Ybiw(Ybiw(:)>0)));
          clear Ybig; 
        end
        if useGMandCSF(2)
          % use additonal CSF
          Ybig = Ybi .* Yc; Ybig = cat_vol_localstat(Ybig,Ybig>0,1,2); 
          Ybig = Ybi .* (Yc & Ybig~=Ybi & abs(Ybig-Ybi)<bith*0.1); 
          Ybig = cat_vol_localstat(Ybig,Ybig>0,1,biasstr{bi,4}); 
          for ii=1:round(8/mean(vx_vol)), Ybig = cat_vol_localstat(Ybig,Ybig>0,1,1); end
          Ybiw = Ybiw + (Ybig / median(Ybig(Ybig(:)>0)) * median(Ybiw(Ybiw(:)>0)));
          clear Ybig; 
        end
        %
        if job.opts.spm
          % use additonal CSF
          Ybig = Ybi .* Yh; Ybig = cat_vol_localstat(Ybig,Ybig>0,1,2); 
          Ybig = Ybi .* (Yh & Ybig~=Ybi & abs(Ybig-Ybi)<bith*0.1); 
          Ybig = cat_vol_localstat(Ybig,Ybig>0,1,biasstr{bi,4}); 
          for ii=1:round(4/mean(vx_vol)), Ybig = cat_vol_localstat(Ybig,Ybig>0,5,1); end
          %%
          Ybiw = Ybiw + (Ybig / mean(Ybiw(Ybiw(:)>0)) * mean(Ybi(Yw(:)>0)));
          if ~debug, clear Ybig; end
        end
        %%
        Ybiws = cat_vol_approx(Ybiw,'nn',vx_vol,biasstr{bi,5},struct('lfO',biasstr{bi,6}));  %#ok<NASGU>
        eval([biasstr{bi,2} 'ws = Ybiws;']);
        if debug, eval([biasstr{bi,2} 'w = Ybiw;']); end
        if ~debug, clear Ybiw Ybiws; end  
      end
      if 0
        % this is just for manual development & debugging
        %% display PD
        ds('d2','a',vx_vol,Ybiw./mean(Yp(Yw(:))),Ypws./mean(Yp(Yw(:))),Yp./mean(Yp(Yw(:))),Yp./Ypws,120); colormap gray
        %% display R1
        ds('d2','a',vx_vol,Yd./mean(Yd(Yw(:))),Ydws./mean(Yd(Yw(:))),Yd./mean(Yd(Yw(:))),Yd./Ydws,120); colormap gray
        %% display R2s
        ds('d2','a',vx_vol,Yr./mean(Yr(Yw(:))),Yrws./mean(Yr(Yw(:))),Yr./mean(Yr(Yw(:))),Yr./Yrws,120); colormap gray
      end
    end
    
    
    
    
    %% intensity normalization
    %  --------------------------------------------------------------------
    %  Ydm  = bias corrected R1 map
    %  Ypm  = bias corrected PD map
    %  Yrm  = bias corrected R2s map
    %
    %  Ydmc = intensity normalized bias corrected R1 map
    %  Ypmc = intensity normalized bias corrected PD map
    %  Yrmc = intensity normalized bias corrected R2s map
    %  --------------------------------------------------------------------
    if job.opts.in
      stime  = cat_io_cmd('  Intensity Normalization:','g5','',job.opts.verb-1,stime);

      % apply bias correction, create inverted version of the pd image
      if job.opts.bc
        Ydm = Yd ./ Ydws; 
        Ypm = Yp ./ Ypws; 
        Yrm = Yr ./ Yrws; 
      else
        Ydm = Yd ./ mean(Yd(Yw(:))); 
        Ypm = Yp ./ mean(Yp(Yw(:))); 
        Yrm = Yr ./ mean(Yr(Yw(:))); 
      end
      if ~debug, clear Yd Yp Yr; end
      Ybs = cat_vol_smooth3X(Ybm,1); 
      
      %% create synthetic T1 contrast based on the PD image
      %Ytm = 1.5 - Ypm/2; 
      Ypmi = cat_stat_histth(Ypm,99);
      Ytm  = Ybs - (Ypmi - min(Ypmi(Ybs(:)>0))) / ( max(Ypmi(Ybs(:)>0)) - min(Ypm(Ybs(:)>0)));
      Ytm  = min(cat_stat_histth( Ytm .* Ybs,99) *0.8-0.1 + 1.2*smooth3(Ybs), Ytm + (1-Ybs)); 
      if debug, clear YM; end 

      %%  the region maps are only used to estimate a the mean intensity of a tissue 
      Ywm = Yw & Yrm<mean(Yrm(Yw(:)))*1.5 & Ypg<pgth*0.5;
      Ygm = Yg & Yrm<mean(Yrm(Yw(:)))*1.5 & Ypg<pgth; 
      Ycm = Yc & Ypg<pgth;
      Yhm = ~Ybs & ~Yw & ~Yg & ~Yc & ~Ybm & Ydm>1.2 & Ypm<0.5 & Yrm>1.2;
      Ylm = ~Ybs & ~Yw & ~Yg & ~Yc & ~Ybm & Ydm>0 & Ypm>0 & Yrm>0 & Ydm<mean(Ydm(Ycm(:)))*1.5 & Ypm>mean(Ypm(Ycm(:)))/1.5 & Yrm<mean(Yrm(Ycm(:)))*1.5;
      if sum(Ylm(:))<10
        Ylm = ~Ybs & ~Yw & ~Yg & ~Yc & Ydm>min(Ydm(:)) & Ypm>min(Ypm(:)) & Yrm>min(Yrm(:)) & Ydm<mean(Ydm(Ycm(:)))*2 & Ypm>mean(Ypm(Ycm(:)))/2 & Yrm<mean(Yrm(Ycm(:)))*2;
      end  
      if ~debug, clear Yg Yc Ybs; end

      %% real thresholds [bg lm cm gm wm hm max] 
      T3th(1:4,1) = {'Ydm';'Ypm';'Yrm';'Ytm'};
      T3th{1,2} = [ 0 min([nan,min(Ydm(Ylm(:)))]) mean(Ydm(Ycm(:))) mean(Ydm(Ygm(:))) mean(Ydm(Ywm(:))) cat_stat_nanmean([nan,mean(Ydm(Yhm(:)))]) max(Ydm(:))];
      T3th{2,2} = [ 0 mean(Ypm(Ycm(:))) mean(Ypm(Ygm(:))) mean(Ypm(Ywm(:))) cat_stat_nanmean([nan,mean(Ypm(Ylm(:)))]) max(Ypm(:))];
      T3th{3,2} = [ 0 min([nan,min(Yrm(Ylm(:)))]) mean(Yrm(Ycm(:))) mean(Yrm(Ygm(:))) mean(Yrm(Ywm(:))) cat_stat_nanmean([nan,mean(Yrm(Yhm(:)))]) max(Yrm(:))];
      T3th{4,2} = [ min(Ytm(:)) min([nan,min(Ytm(Ylm(:)))]) mean(Ytm(Ycm(:))) mean(Ytm(Ygm(:))) mean(Ytm(Ywm(:))) max(Ytm(:))];
      if sum(Yhm(:)>0)<1000 && sum(Ylm(:)>0)<1000
        T3th{1,2}([2,6]) = [mean(T3th{1,2}(1:2:3)) mean(T3th{1,2}(5:2:7))];
        T3th{2,2}(5)     = mean(T3th{2,2}(4:2:6));
        T3th{3,2}([2,6]) = [mean(T3th{3,2}(1:2:3)) mean(T3th{3,2}(5:2:7))];
        T3th{4,2}(2)     = 0.02;
      end  
      Yp0 = Ycm + 2*Ygm + 3*Ywm; 
      if ~debug, clear Ylm Yhm Ycm Ygm Ywm; end
      
      %% final thresholds
      T3th{1,3} = [ 0 1/12 1/3 2/3 1   T3th{1,2}(5)+diff(T3th{1,2}(5:6)) T3th{1,2}(6)+diff(T3th{1,2}(6:7))]; 
      T3th{2,3} = [ 0 1 2/3 1/3 T3th{2,2}(5)/(1/T3th{2,2}(2)) T3th{2,2}(6)/(1/T3th{2,2}(2))]; 
      T3th{3,3} = [ 0 1/12 1/3 2/3 1   T3th{3,2}(5)+diff(T3th{3,2}(5:6)) T3th{3,2}(6)+diff(T3th{3,2}(6:7))]; 
      T3th{4,3} = [ 0 0.02 1/3 2/3 1   T3th{4,2}(6)/T3th{4,2}(5)]; 
      
      % global intensity normalization 
      for bi=1:size(T3th,1)
        %%
        [T3ths,tsi] = sort(T3th{bi,2});
        T3thx       = T3th{bi,3}(tsi);
        eval(['Ym = ' T3th{bi,1} ' + 0; if ~debug, clear ' T3th{bi,1} '; end; Ysrc = Ym;']); 
        for i=numel(T3ths):-1:2
          M = Ysrc>T3ths(i-1) & Ysrc<=T3ths(i);
          Ym(M(:)) = T3thx(i-1) + (Ysrc(M(:)) - T3ths(i-1))/diff(T3ths(i-1:i))*diff(T3thx(i-1:i));
        end
        M  = Ysrc>=T3ths(end); 
        Ym(M(:)) = numel(T3ths)/6 + (Ysrc(M(:)) - T3ths(i))/diff(T3ths(end-1:end))*diff(T3thx(i-1:i));    
        eval([T3th{bi,1} 'c = Ym;']); 
        clear Ym Ysrc M;
      end
      if 0 
        %% just for debuging and display
        ds('d2','',vx_vol,Ydmc,Ypmc,Yrmc,Ytmc,220); colormap gray
        %% 
        ds('d2','a',vx_vol,Ydm,Ydmc,Yrm,Yrmc,220); colormap gray
        %%
        ds('d2','a',vx_vol,Ypm,Ypmc*2,Ytm,Ytmc,140); colormap gray
      end
      
    else
      % no intensity normalization 
      Ydmc = Yd ./ mean(Yd(Yw(:))); 
      Ypmc = Yp ./ mean(Yp(Yw(:))); 
      Yrmc = Yr ./ mean(Yr(Yw(:))); 
      Ytmc = 1.5 - Ypmc/2; Ytmc(Ypmc<0.8 & cat_vol_morph(cat_vol_morph(Ypmc<0.5 | Yrmc>1.5,'l'),'d',1.5))=0;
      if ~debug, clear Ybs; end
    end
    
   
      
    
    
    
    %% find blood vessels
    %  --------------------------------------------------------------------
    %  Yrd, Ydd, Ytd = divergence maps
    %  --------------------------------------------------------------------
    Yb    = cat_vol_morph((Ybm | (Yp0>0.1)) & Ydmc>0 & Ydmc<1.1 & Yrmc>0 & Yrmc<1.1 & Ytmc>0 & Ytmc<1.1,'lc',max(1,min(2,1/mean(vx_vol))));   % brain mask for CSF minimum value
    Yb    = cat_vol_morph(Yb,'lo',0.5/mean(vx_vol));
    if job.opts.bvc || job.output.bv 
      stime  = cat_io_cmd('  Blood Vessel Detection:','g5','',job.opts.verb-1,stime);

      % divergence update
      Yrd = cat_vol_div(Yrmc);
      Ydd = cat_vol_div(Ydmc);
      Ytd = cat_vol_div(Ytmc);


      %% first we remove high intensity blood vessels with large area 
      %  in the save WM region only a mean filter should be used to create an
      %  corrected map, whereas in the other regions a larger and minimum
      %  based filter is applied
      Ywc   = cat_vol_morph(Yw,'lc',1) | ~cat_vol_morph(~cat_vol_morph(Yw,'o',1),'lc',max(1,min(2,1/mean(vx_vol))));       % wm mask with close microstructures
      Ywc   = Ywc | smooth3(Yrmc>0.8 & Yrmc<1.1 & Yb)>0.5 | cat_vol_morph(smooth3(Ydmc>0.7 & Ydmc<1.2 & Yb)>0.6,'o',1);
      Ywc   = cat_vol_morph(Ywc,'lc',1);
     
      %% correction maps    
      Ymin  = max(Yb.*3/12,Yb.*min(Ydmc,min(Ytmc,Yrmc)));   % minimum intensity in all images (CSF)
      Ymax  = max(Yb.*3/12,Yb.*max(Ydmc,max(Ytmc,Yrmc)));   % maximum intensity in all images (BV)
      Ymn1  = cat_vol_localstat(Ymin,Ymin>0,1,1);       % low    short  correction
      Ymin1 = cat_vol_localstat(Ymin,Ymin>0,1,2);       % storng short  correction
      Ymax1 = cat_vol_localstat(Ymax,Ymax>0,1,3);       % very close to vessel
      Ymax2 = cat_vol_localstat(Ymax1,Ymax>0,1,3);      % close to vessel
      Ymin2 = cat_vol_localstat(Ymin1,Ymin>0,1,2);      % strong medium correction
      Ymin4 = cat_vol_localstat(Ymin2,Ymin>0,2,2);      % strong long   correction


      %% specific masks 
      %  - save regions that need no corrections
      %  - save blood vessels that requires maxim correction (long distance with min)
      %  - save blood vessels that requires carfull correction (short distance with mean) 
      %  - small blood vessels that requires less agressive correction 
      YminA = Ytmc  .* (                                     Ymax<1.1 &  Ywc) + ... low  values in the WM 
              Ymn1  .* (                         Ymax>=1.1 &             Ywc) + ... high values in the WM 
              Ymin  .* (                         Ymax< 0.7 &            ~Ywc) + ... low  good values
              Ymn1  .* (                         Ymax>=0.7 & Ymax<1.0 & ~Ywc) + ... affected voxel 
              Ymn1  .* (                         Ymax>=1.0 & Ymax<1.1 & ~Ywc) + ...
              Ymin  .* (Ymin< 0.9 & Ymax2< 1.3 & Ymax>=1.1 & Ymax<1.3 & ~Ywc) + ...
              Ymin1 .* (Ymin< 0.9 & Ymax2>=1.3 & Ymax>=1.1 & Ymax<1.3 & ~Ywc) + ...
              Ymin  .* (Ymin< 0.9 & Ymax2< 1.3 & Ymax>=1.3            & ~Ywc) + ...
              Ymin2 .* (Ymin< 0.9 & Ymax2>=1.3 & Ymax>=1.3            & ~Ywc) + ...
              Ymn1  .* (Ymin>=0.9              & Ymax>=1.1 & Ymax<1.3 & ~Ywc) + ...
              Ymin4 .* (Ymin>=0.9              & Ymax>=1.3            & ~Ywc);
      

      %% find blood vessels
      Ybv1 = (Yrmc>2 | Ytmc>1.5 | Ydmc>1.5) & ~Ywc & Yb; 
      Ybv  = Yb .* max(0 , max(Yrmc.*Ydmc.*Ytmc,Ymax) - min(Yrmc.*Ydmc.*Ytmc,Ymin) - 2/3 + 2*max(Yrd.*Ydd.*Ytd,max(Yrd,max(Ydd,Ytd)))  - Ywc );
      
      if debug 
        %% display BV
        ds('d2','a',vx_vol,Ydmc*1.5,Ytmc*1.5,Ybv,Yp0/3*2,180); colormap gray
      end
    end
    %% correct blood vessels
    %  --------------------------------------------------------------------
    %  --------------------------------------------------------------------
    if job.opts.bvc
      stime  = cat_io_cmd('  Blood Vessel Correction:','g5','',job.opts.verb-1,stime);

      %% first correction 
      Ybvd = cat_vol_div(Ybv);
      Ylow = smooth3(((Ymax1+1)./(Ymin1+1))>(Ymin+1) | Ybv>2 | Ybv1)>0.2;
      %%
      %Ypc  = min(1.1.*(Ypm>0.1),max( (1.25-Ymin2) .* (Ypm>0),Ypmc - (Ylow | Ybv1) .* ~Yw .* (min(0,Ybvd) - max(0,Ybv-0.5) - max(0,(Ydmc./Ypmc) - 1.5) )));
      Ydc  = max(Ymin2,Ydmc + (Ylow | Ybv1) .* ~Yw .* (min(0,Ybvd) - max(0,Ybv-0.5) - max(0,(Ydmc./Ypmc) - 1.5) ));
      Ytc  = max(Ymin2,Ytmc + (Ylow | Ybv1) .* ~Yw .* (min(0,Ybvd) - max(0,Ybv-0.5) - max(0,(Ydmc./Ypmc) - 1.5) ));
      Yrc  = max(Ymin2,Yrmc + (Ylow | Ybv1) .* ~Yw .* (min(0,Ybvd) - max(0,Ybv-0.5) - max(0,(Ydmc./Ypmc) - 1.5) ));

      %% median fitlering
      %Ydcm = cat_vol_median3(Ydc,Ydc>0,true(size(Ydc)),0.05); 
      Ytcm = cat_vol_median3(Ytc,Ytc>0,true(size(Ydc)),0.05);
      Yrcm = cat_vol_median3(Yrc,Yrc>0,true(size(Ydc)),0.05);

      %% second correction
      % r1
      YminAm = cat_vol_median3(YminA,YminA>0,true(size(Ydc)),0.1); 
      YM   = YminAm>0.2 & ~Ywc & ((Ydmc - YminAm)>0.2 | Ydmc>1.5 | Ybv>1) & Yb; 
      Ydc2 = Ydmc; Ydc2(YM) = YminAm(YM); %Ydc2 = min(Ydc2.*Yb,Ydmc);  %.* (YminA>0.05)
      Ydc2 = cat_stat_histth(Ydc2); 
      % pd
      YM   = Ypmc>0.01 & Yb & ~Ywc & smooth3((Ypmc - (1.25-YminAm))>0.1 | Ypmc<0.3 | Ybv>1)>0.1; 
      YH   = 1.1 - 0.05*rand(size(YM)); 
      Ypc2 = Ypmc; Ypc2(YM) = min(YH(YM),1.25-YminAm(YM).*Yb(YM)); 
      Ypc2 = cat_stat_histth(Ypc2); 
      % t1 = inverse pd 
      YM   = Ytmc>0.2 & ~Ywc & ((Ytmc - Ytcm)>0.1 | Ytmc>1.1 | Ybv>0.5) & Yb; 
      Ytc2 = Ytmc; Ytc2(YM) = Ytcm(YM); Ytc2 = min(Ytc2,Ytmc);
      YM   = (Ydc2 - Ytc2)<-0.4 & Ydc2<0.5 & Yb; 
      Ytc2(YM) = min(Ytc2(YM) ,Ydc2(YM).*Yb(YM) );
      Ytc2 = cat_stat_histth(Ytc2); 
      % r2s
      YM   = Yrmc>0.2 & ~Ywc & ((Yrmc - Yrcm)>0.2 | Yrmc>1.1 | Ybv>0.5) & Yb; 
      Yrc2 = Yrmc; Yrc2(YM) = Yrcm(YM); Yrc2 = min(Yrc2,Yrmc);
      YM   = smooth3(Ydc2 - Yrc2)<-0.2 & (Ydc2<0.8 | Yrc2>1.3) & (~Ywc | (Yrd<-0.1 & Yrc2>1.3 & Ypg>0.05))>0.5 & Yb; 
      Yrc2(YM) = min(Yrc2(YM) ,Ydc2(YM).*Yb(YM) );
      Yrc2 = cat_vol_median3(Yrc2,YM | Yrc2>1.3 & Ypg>0.2); % remove small reminding WM vessels and artifacts
      Yrc2 = cat_stat_histth(Yrc2); 
      if ~debug, clear Ypg; end
      
      %% this is just for manual development & debugging
      if 0
        
        %% display BV
        ds('d2','a',vx_vol,Ybv1,Ybv,Ydcm,abs(Ydc2-Yd),200); colormap gray
        %% display PD
        ds('d2','a',vx_vol,Ypmc/2.5*1.5,Ypws./mean(Ypmc(Yw(:))),Ypmc*1.5,Ypc2*1.5,120); colormap gray
        %% display PDi
        ds('d2','a',vx_vol,2 - Ypc2*1.5,Yp/mean(Yp(Yw(:))),Ydmc*1.5,Ytc2*1.5,120); colormap gray
        %% display R1
        ds('d2','a',vx_vol,Ydws./mean(Ydws(Yw(:)))*1.5,Ydws./mean(Ydws(Yw(:))),Ydc*1.5,Ydc2*1.5,120); colormap gray
        %% display R2s
        ds('d2','a',vx_vol,Yr./mean(Yr(Yw(:))),Yrws./mean(Yr(Yw(:))),Yr./mean(Yr(Yw(:))),Yrc2,120); colormap gray
        
        
        %% create test WM surface with correct blood vessel system - Warning very slow! 
        S   = isosurface(cat_vol_smooth3X(Ytc2,0.5),0.8); 
        Sb  = isosurface(cat_vol_smooth3X(Ybv .* Yb,0.5),0.25); 
        
        % reduce surface
        SR  = reducepatch(S,500000);
        SbR = reducepatch(Sb,500000);
        
        % add some color
        SR.facevertexcdata  = zeros(size(SR.vertices,1),1); 
        SbR.facevertexcdata = ones(size(SbR.vertices,1),1); 
        
        % combine surfaces
        SRX.vertices        = [SR.vertices;SbR.vertices];
        SRX.faces           = [SR.faces;SbR.faces+size(SR.vertices,1)];
        SRX.cdata = [SR.facevertexcdata;SbR.facevertexcdata];

        % dispaly surfaces
        cat_surf_render2(SRX);
        
      end
    else
      Ypc2 = Ypmc;
      Ytc2 = Ytmc;
      Ydc2 = Ydmc;
      Yrc2 = Yrmc;
    end
    if ~debug, clear Ypmc Ytmc Ydmc Yrmc; end
    
    
     if job.opts.ss
       Ypc2 = Ypc2 .* Yb; 
       Ytc2 = Ytc2 .* Yb; 
       Ydc2 = Ydc2 .* Yb; 
       Yrc2 = Yrc2 .* Yb; 
     end
    
    
    %% write data
    %  --------------------------------------------------------------------
    stime  = cat_io_cmd('  Write Output:','g5','',job.opts.verb-1,stime);
    
    % PD 
    if job.output.pd  
      [pp,ff,ee,dd] = spm_fileparts(Vp.fname);
      Vpc = Vp;  Vpc.fname = fullfile(pp,[job.opts.prefix 'pd_' ff ee dd]); Vpc.dt(1) = 16; 
      spm_write_vol(Vpc,Ypc2); 
    end
    
    % inverted PD
    if job.output.t1 
      [pp,ff,ee,dd] = spm_fileparts(Vp.fname);
      Vtc = Vp;  Vtc.fname = fullfile(pp,[job.opts.prefix 't1_' ff ee dd]); Vtc.dt(1) = 16; 
      spm_write_vol(Vtc,Ytc2); 
    end
    
    % R1
    if job.output.r1% r1
      [pp,ff,ee,dd] = spm_fileparts(Vd.fname);
      Vdc = Vd;  Vdc.fname = fullfile(pp,[job.opts.prefix 'r1_' ff ee dd]); Vdc.dt(1) = 16; 
      spm_write_vol(Vdc,Ydc2); 
    end
    
    % R2s
    if job.output.r2s
      [pp,ff,ee,dd] = spm_fileparts(Vr.fname);
      Vrc = Vr;  Vrc.fname = fullfile(pp,[job.opts.prefix 'r2s_' ff ee dd]); Vrc.dt(1) = 16; 
      spm_write_vol(Vrc,Yrc2); 
    end
    
    % BV
    if job.output.bv 
      [pp,ff,ee,dd] = spm_fileparts(Vp.fname);
      Vbvc = Vp;  Vbvc.fname = fullfile(pp,[job.opts.prefix 'bv_' ff ee dd]); Vbvc.dt(1) = 16; 
      spm_write_vol(Vbvc,Ybv .* Yb); 
    end    
    if ~debug, clear Ypc2 Ytc2 Ydc2 Yrc2 Ybv Yb; end
    
    
    
    %% noise reduction
    if job.opts.nc
      stime  = cat_io_cmd('  Noise Correction:','g5','',job.opts.verb-1,stime); fprintf('\n')
      if job.output.pd,  cat_vol_sanlm(struct('data',cellstr(Vpc.fname),'prefix','')); end
      if job.output.t1,  cat_vol_sanlm(struct('data',cellstr(Vtc.fname),'prefix','')); end
      if job.output.r1,  cat_vol_sanlm(struct('data',cellstr(Vdc.fname),'prefix','')); end
      if job.output.r2s, cat_vol_sanlm(struct('data',cellstr(Vrc.fname),'prefix','')); end
      if job.output.bv,  cat_vol_sanlm(struct('data',cellstr(Vbvc.fname),'prefix','')); end
    end
    
    cat_io_cmd(' ','g5','',job.opts.verb-1); cat_io_cmd(' ','g5','',job.opts.verb-1,stime); % time of the last step
    cat_io_cmd(' ','','',job.opts.verb-1,stime2); fprintf('\n'); % time of the full subject
  end  
end