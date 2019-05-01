function varargout = cat_vol_isarnlm(varargin) 
% Ys = cat_vol_isarnlm(varargin)
% ______________________________________________________________________
% Iterative Spatial Adaptive mulit-Resolution Non Local Means (ISARNLM) 
% noise correction with improved filtering of parallel imaging artifacts.
%
% Filter a set of images and add the prefix 'isarnlm_'.
% Missing input will call GUI or/and use defaults. 
%
% WARNING: SPM can have problems with images with very low noise such as
%          the Brain Web Phantom with 0% noise. Although, there is more 
%          variance in real images even after noise correction please 
%          check your results.
%          Use the job.cstr parameter to modify the correction strength
%          (1=full (default), 0=none)
%          
% Input:
% job    - harvested job data structure (see matlabbatch help)
% 
% Output:
% out    - computation results, usually a struct variable.
%
% cat_vol_sanlm(job)
%   job.data   = set of images 
%   job.prefix = prefix for filtered images (default = 'isarnlm_') 
%   job.rician = noise distribution
%   job.cstr   = correction strength (1=full,0=none)
% Example:
%   cat_vol_sanlm(struct('data','','prefix','n','rician',0));
%   Ys = cat_vol_sanlm(Yo,Vm,[,verb,NCstr]); 
% 
%_______________________________________________________________________
% Christian Gaser, Robert Dahnke 
% Structural Brain Mapping Group
% University Jena
%
% $Id: cat_vol_isarnlm.m 1183 2017-09-08 16:45:08Z dahnke $
% ______________________________________________________________________
 
  if nargin == 0 
      job.data = cellstr(spm_select([1 Inf],'image','select images to filter'));
  elseif nargin == 1 && isstruct(varargin{1})
      job = varargin{1};
      if nargout>0, error('No output availailable. '); end
  elseif (nargin == 2 || nargin == 3 || nargin == 4) && isstruct(varargin{2}) && nargout==1
    stime = clock;
    if isstruct(varargin{2})
      V = varargin{2}; 
      vx_vol  = sqrt(sum(V.mat(1:3,1:3).^2));
    else
      vx_vol = varargin{2};
    end
    if nargin >= 3, verb=varargin{3};  else verb=1;  end
    if nargin == 4, NCstr=varargin{4}; else NCstr=1; end
    varargout{1} = cat_vol_sanlmX(varargin{1},'',vx_vol,struct('verb',verb,'NCstr',NCstr));
    if verb>0
      cat_io_cmd(' ','','',1); 
      fprintf('%4.0fs\n',etime(clock,stime));
    end
    return
  else
    if nargin>3 && isfield(varargin{4},'verb'), verb = varargin{4}.verb; else verb = 1; end 
    if verb, fprintf('isarnlm:\n'); stime=clock; end
    eval(sprintf('varargout{1} = cat_vol_sanlmX(varargin{1}%s);',sprintf(',varargin{%d}',2:nargin)));
    if verb, fprintf('isarnlm done in %0.0fs.\n',etime(clock,stime)); end
    return
  end
  if ~isfield(job,'data') || isempty(job.data)
     job.data = cellstr(spm_select([1 Inf],'image','select images to filter'));
     if isempty(job.data) || isempty(job.data{1}), return; end 
  else
     job.data = cellstr(job.data);
  end
  if isempty(job.data), return; end

  def.verb    = 1;           % be verbose
  def.prefix  = 'isarnlm_';  % prefix
  def.postfix = '';  
  def.NCstr   = -inf;         % 0 - no denoising, eps - light denoising, 1 - maximum denoising, inf = auto; 
  def.rician  = 0;           % use inf for GUI
  def.local   = 1;           % local weighing (only auto NCstr); 
  job = cat_io_checkinopt(job,def);
  if isinf(job.NCstr), job.NCstr = -2; end 
  %job.NCstr = max(0,min(1,job.NCstr)) + isinf(job.NCstr)*job.NCstr;          % garanty values from 0 to 1 or inf
  if isinf(job.rician), spm_input('Rician noise?',1,'yes|no',[1,0],2); end  % GUI
  
  V = spm_vol(char(job.data));

  spm_clf('Interactive'); 
  spm_progress_bar('Init',numel(job.data),'ISARNLM-Filtering','Volumes Complete');
  for i = 1:numel(job.data)
      [pth,nm,xt,vr] = spm_fileparts(deblank(V(i).fname));
      vx_vol  = sqrt(sum(V(i).mat(1:3,1:3).^2));
 
      if any(strfind(nm,job.prefix)==1), continue; end
      if job.verb==1, fprintf('ISARNLM %s:\n',spm_str_manip(deblank(V(i).fname),'a60')); end
      stime=clock;
      
      src = single(spm_read_vols(V(i)));
      % prevent NaN
      src(isnan(src)) = 0;
      src = cat_vol_sanlmX(src,'',vx_vol,job);

      V(i).fname = fullfile(pth,[job.prefix nm job.postfix '.nii' vr]);
      V(i).descrip = sprintf('%s ISARNLM filtered',V(i).descrip);

      % use at least float precision
      if V(i).dt(1)<16, V(i).dt(1) = 16; end 
      spm_write_vol(V(i), src);
      if job.verb==1, fprintf('isarnlm done in %0.0fs.\n',etime(clock,stime)); end
      spm_progress_bar('Set',i);
  end
  spm_progress_bar('Clear');
end
function [Ys,NCstr] = cat_vol_sanlmX(Y,YM,vx_vol,opt)
% Ys = cat_vol_sanlmX(Y,YM,vx_vol)
% ______________________________________________________________________
% Adaptive iterative multiresolution NLM noise correction for highres 
% images with GRAPPA or other strong (regular/non movement) artifacts.
%
% A mask can be used to filter only a specific region of the image to 
% allow faster computation. For full filtering use an empty matrix YM. 
%
%   Ys = cat_vol_sanlmX(Y,YM,vx_vol,opt)
% 
%   Ys         .. filtered image
%   Y          .. original image
%   YM         .. filter mask or empty matrix
%   vx_vol     .. voxel volume
%
%   opt.verb   .. display progess (default = 1)
%   opt.red    .. maximum number of resolution reduction (default = 2)
%   opt.iter   .. maximum number of iterations (default = 3) 
%   opt.rician .. noise type (default = 0) 
%   opt.NCstr   .. correction strength (1=full,0=none)
%   opt.SANFM  .. spatial adaptive noise filter modification (default=1)
%                 resolution and spation noise pattern depending
%                 filter strength (opt.Sth)
%   opt.Nth    .. noise threshold (default = 0.015)
%                 filter/reduce only for noise>0.015
%   opt.Sth    .. noise-signal threshold (default = 4), req. opt.SANFM
%                 lower values = less filtering of artifacts/anatomie
%
% The filter reduce high resolution images (<1.5 mm), to remove noise
% on a lower frequency level. To avoid to strong filtering of anatomical
% details, the 'filter strength' of the low resolution level also depend 
% on the resolution and the noise correction. 
% It runs multiple times because the first filtering in noisy images 
% just restores the basic anatomical pattern.
% 
% Special cases:
% * Interative filter with i iterations only on the original resolution: 
%     Ys =cat_vol_isarnlm(Yi,Yp0,vx_vol,struct('iter',i,'red',0));
% * Resolution filter without iteration: 
%     Ys = cat_vol_isarnlm(Yi,Yp0,vx_vol,struct('iter',0,'red',r));
% 
% ______________________________________________________________________
% Robert Dahnke 
% Structural Brain Mapping Group
% University Jena
%
% $Id: cat_vol_isarnlm.m 1183 2017-09-08 16:45:08Z dahnke $
% ______________________________________________________________________

  if ~exist('opt','var'), opt = struct(); end

  def.verb   = 1;     % display progess
  def.red    = 1;     % maximum number of resolution reduction
  def.iter   = 2;     % maximum number of iterations
  def.iter1  = 2;     % maximum number of iterations at full resolution
  def.rician = 0;     % noise type 
  def.cstr   = 0.5;   % correction strength 
  def.SANFM  = 1;     % spatial adaptive noise filter modification 
  def.Nth    = 0.01;  % noise threshold (filter only  
  def.fast   = 0;     % masking background?
  def.Sth    = 4;     % noise-signal threshold (lower values = less filtering of artifacts/anatomie)
  def.level  = 1;     % just for display
  def.NCstr  = -1; 
  def.local  = 1; 
  opt        = cat_io_checkinopt(opt,def);
  opt.iter   = max(1,min(10,opt.iter));  % at least one iteration (iter = 0 means no filtering)
  %opt.NCstr  = max(0,min(1,opt.NCstr)) + isinf(opt.NCstr)*opt.NCstr;  % range 0-1

  if isempty(YM), YM = true(size(Y)); end 
  YM = YM>0.5;
  
  Tth = median(Y(Y(:)>median(Y(Y>2*median(Y(:)))))); 
 
  if opt.fast 
    Y0=Y; 
    [Y,YM,BB] = cat_vol_resize({Y,YM},'reduceBrain',vx_vol,4,Y>Tth*0.2);
  end
  %Yo = Y; 
  Yi = Y .* YM;
  % just for display
  % ds('d2','',vx_vol,Y/Tth*0.95,Yi/Tth*0.95,Ys/Tth*0.95,abs(Yi-Ys)./max(Tth*0.2,Ys),90)
  
   
  iter = 0; noise = inf; Ys = Yi; noiser=1;
  while ((iter < opt.iter && opt.level>1) || (iter < opt.iter1 && opt.level==1)) && ...
      noise>opt.Nth && (opt.level<4 || noiser>1/4) && (iter==0 || mean(vx_vol)<1.5) 
    
    
    %% SANLM filtering
    if opt.level~=1
      if opt.verb, cat_io_cprintf('g5',sprintf('%3d.%d) %0.2fx%0.2fx%0.2f mm:  ',opt.level,iter+1,vx_vol)); stime = clock; end
      Ys  = Yi+0;
      YM2 = YM & Ys>Tth*0.2 & Ys<max(Ys(:))*0.98;

      cat_sanlm(Ys,3,1,opt.rician); 
    
      if opt.NCstr<0
      % adaptive local denoising 
        NCstr = opt.NCstr * 15; 

        % prepare local map
        NCs = abs(Ys - Yi) ./ max(eps,Ys); spm_smooth(NCs,NCs,2);
        NCs = max(0,min(1,NCs * abs(NCstr))); 
        opt.NCstr = -cat_stat_nanmean(NCs(:)); 

        % mix original and noise corrected image
        Ys = Yi.*(1-NCs) + Ys.*NCs; 
        clear NCs;

      elseif opt.NCstr>0
      % (adaptive) global denoising  

        NCstr = opt.NCstr * 15; 
     
        opt.NCstr = min(1,max(0,opt.NCstr));

        % mix original and noise corrected image
        Ys   = Yi*(1-opt.NCstr) + Ys*opt.NCstr; 
      end

    
      %{
      % adaptive global denoising 
      if isinf(opt.NCstr) || sign(opt.NCstr)==-1;
        YM3    = Ys>mean(Ys(:)); % object
        Tth    = mean(Ys(YM3(:)));
        NCstr  = min(1,max(0, mean( abs(Ys(YM3(:)) - Yi(YM3(:))) ./ Tth ) * 16 * min(10,max(0,abs(opt.NCstr))) ));
        NCs    = abs(Ys - Yi) ./ max(eps,Ys) * 16 * min(2,max(0,abs(opt.NCstr))); spm_smooth(NCs,NCs,2);
        NCs    = max(0,min(1,NCs));
        clear YM3; 
      else 
        NCstr  = opt.NCstr;
      end
      NCstr = max( mean(NCstr(:)) - 1*std(NCstr(:)), min ( mean(NCstr(:)) + 1*std(NCstr(:)) , NCstr)); 


      % mix original and noise corrected image
      if opt.local
        Ys = Yi.*(1-NCs) + Ys.*NCs; clear NCs; 
      else
        Ys = Yi*(1-NCstr) + Ys*NCstr; 
      end
      %}
    
      %Ys = (1-opt.NCstr) .* Yi  + opt.NCstr .* Ys; % main weighting
      %[i,txt] = feature('numCores'); i=strfind(txt,'MATLAB was assigned:');
      %fprintf(sprintf('%s',repmat('\b',1,numel('Using 8 processors '))));
      noiser = 1 - (cat_stat_nanmean(abs(Y(YM2(:))-Ys(YM2(:)))./max(Tth*0.2,Ys(YM2(:))))/sqrt(prod(vx_vol))) / noise;
      if noiser<0, noiser = noiser+1; end
      noise  = cat_stat_nanmean(abs(Y(YM2(:))-Ys(YM2(:)))./max(Tth*0.2,Ys(YM2(:))))/sqrt(prod(vx_vol));
      clear YM2;
      if opt.verb, cat_io_cprintf('g5',sprintf('  noise = %4.3f, noiser = %4.3f        %5.0fs\n',noise,noiser,etime(clock,stime))); end
    else
      noiser = 1; 
    end
 
    

    %% filtering of lower resolution level
    %  if the currect resolution is high enought
    %  important is a previous NLM on the main resolution to avoid 
    %  filtering of fine anatomical structures on lower resolutions
    %if opt.red && all(vx_vol<2.1) && sum(vx_vol<1.1)>1 && (noise>opt.Nth || iter==0) %&& def.red>0 && noiser>1/4  && iter<1
    if opt.red && all(vx_vol<2.1) && sum(vx_vol<0.8)>1 && (noise>opt.Nth || iter==0) %&& def.red>0 && noiser>1/4  && iter<1
      %%
      Yi = Ys + 0;
    
      if all(vx_vol<2)
        % first block
        [Yr,YMr,resr] = cat_vol_resize({Yi,YM},'reduceV',vx_vol,min(2.2,min(vx_vol)*2.6),32,'meanm'); YMr = YMr>0.5;
        optr       = opt;
        optr.red   = opt.red - 1; 
        %optr.iter  = opt.iter - 1;
        optr.Nth   = opt.Nth / 2; %* prod(resr.vx_vol) / prod(resr.vx_volr);
        optr.level = opt.level + 1; 
        YR  = cat_vol_resize(Yr,'dereduceV',resr,'nearest');
        Yr  = cat_vol_sanlmX(Yr,YMr,resr.vx_volr,optr);
        YRs = cat_vol_resize(Yr,'dereduceV',resr,'nearest');
        Ys  = (Yi - YR) + YRs; 
        clear YR YRs Yr YRr YRs;

        % second block
        [Yr,YMr,resr] = cat_vol_resize({Yi(2:end,2:end,2:end),YM(2:end,2:end,2:end)},...
          'reduceV',vx_vol,min(2.2,min(vx_vol)*2.6),32,'meanm'); YMr = YMr>0.5;  
        YRr = cat_vol_resize(Yr,'dereduceV',resr,'nearest');
        YR  = Yi; YR(2:end,2:end,2:end) = YRr; 
        Yr  = cat_vol_sanlmX(Yr,YMr,resr.vx_volr,optr);
        YRr = cat_vol_resize(Yr,'dereduceV',resr,'nearest');
        YRs = Yi; YRs(2:end,2:end,2:end) = YRr; 
        Ys  = Ys + (Yi - YR) + YRs; 
        clear YR YRs Yr YRr YRs; 
        Ys = Ys / 2;
        
        Yis=Ys; 
        % second block
        [Yr,YMr,resr] = cat_vol_resize({Yis(2:end,2:end,2:end),YM(2:end,2:end,2:end)},...
          'reduceV',vx_vol.*[2 2 1],min(2.2,min(vx_vol)*2.6),32,'meanm'); YMr = YMr>0.5;  
        YRr = cat_vol_resize(Yr,'dereduceV',resr,'nearest');
        YR  = Yis; YR(2:end,2:end,2:end) = YRr; 
        Yr  = cat_vol_sanlmX(Yr,YMr,resr.vx_volr,optr); %./[2 2 1]
        YRr = cat_vol_resize(Yr,'dereduceV',resr,'nearest');
        YRs = Yis; YRs(2:end,2:end,2:end) = YRr; 
        Ys  = (Yis - YR) + YRs; 
         % second block
        [Yr,YMr,resr] = cat_vol_resize({Yis,YM},...
          'reduceV',vx_vol.*[2 2 1],min(2.2,min(vx_vol)*2.6),32,'meanm'); YMr = YMr>0.5;  
        YRr = cat_vol_resize(Yr,'dereduceV',resr,'nearest');
        YR  = Yis; YR = YRr; 
        Yr  = cat_vol_sanlmX(Yr,YMr,resr.vx_volr,optr);
        YRr = cat_vol_resize(Yr,'dereduceV',resr,'nearest');
        YRs = Yis; YRs = YRr; 
        Ys  = Ys + ((Yis - YR) + YRs); 
        Ys  = Ys / 2;
        
        Yis=Ys; 
        % second block
        [Yr,YMr,resr] = cat_vol_resize({Yis(2:end,2:end,2:end),YM(2:end,2:end,2:end)},...
          'reduceV',vx_vol.*[2 1 2],min(2.2,min(vx_vol)*2.6),32,'meanm'); YMr = YMr>0.5;  
        YRr = cat_vol_resize(Yr,'dereduceV',resr,'nearest');
        YR  = Yis; YR(2:end,2:end,2:end) = YRr; 
        Yr  = cat_vol_sanlmX(Yr,YMr,resr.vx_volr,optr);
        YRr = cat_vol_resize(Yr,'dereduceV',resr,'nearest');
        YRs = Yis; YRs(2:end,2:end,2:end) = YRr; 
        Ys  = (Yis - YR) + YRs; 
        % second block
        [Yr,YMr,resr] = cat_vol_resize({Yis(2:end,2:end,2:end),YM(2:end,2:end,2:end)},...
          'reduceV',vx_vol.*[1 2 2],min(2.2,min(vx_vol)*2.6),32,'meanm'); YMr = YMr>0.5;  
        YRr = cat_vol_resize(Yr,'dereduceV',resr,'nearest');
        YR  = Yis; YR(2:end,2:end,2:end) = YRr; 
        Yr  = cat_vol_sanlmX(Yr,YMr,resr.vx_volr,optr);
        YRr = cat_vol_resize(Yr,'dereduceV',resr,'nearest');
        YRs = Yis; YRs(2:end,2:end,2:end) = YRr; 
        Ys  = Ys + ((Yis - YR) + YRs); 
        Ys = Ys / 2;
        
        Yis=Ys; 
        % second block
        [Yr,YMr,resr] = cat_vol_resize({Yis,YM},...
          'reduceV',vx_vol.*[2 1 2],min(2.2,min(vx_vol)*2.6),32,'meanm'); YMr = YMr>0.5;  
        YRr = cat_vol_resize(Yr,'dereduceV',resr,'nearest');
        YR  = Yis; YR = YRr; 
        Yr  = cat_vol_sanlmX(Yr,YMr,resr.vx_volr,optr);
        YRr = cat_vol_resize(Yr,'dereduceV',resr,'nearest');
        YRs = Yis; YRs = YRr; 
        Ys  = (Yis - YR) + YRs; 
        % second block
        [Yr,YMr,resr] = cat_vol_resize({Yis,YM},...
          'reduceV',vx_vol.*[1 2 2],min(2.2,min(vx_vol)*2.6),32,'meanm'); YMr = YMr>0.5;  
        YRr = cat_vol_resize(Yr,'dereduceV',resr,'nearest');
        YR  = Yis; YR = YRr; 
        Yr  = cat_vol_sanlmX(Yr,YMr,resr.vx_volr,optr);
        YRr = cat_vol_resize(Yr,'dereduceV',resr,'nearest');
        YRs = Yis; YRs = YRr; 
        Ys  = Ys + ((Yis - YR) + YRs); 
        clear YR YRs Yr YRr YRs; 
        
        Ys = Ys /2;
      end
    end
    
     
    

    
    %% SANFM = spatial adaptive noise filter modification 
    % - the strength filter depend on the resolution and we have to 
    %   avoid filtering on anatical frequency (lowf variable)
    % - the noise vary typically with low frequency and outlier 
    %   mostly describe anatomical structures (opt.Sth)
%     if def.SANFM 
%       if 0
%         Yh     = Ys>mean(Ys(:)); % object
%         Tth    = mean(Ys(Yh(:)));
%         YRc    = abs(Yi - Ys) / Tth * 8 * min(2,max(0,abs(opt.cstr))); spm_smooth(YRc,YRc,2);
%         YRc    = max(0,min(1,YRc));
% 
%         if sum(YRc)>0
%           YRc   = max(0.2*max(YRc(:)),min(mean(YRc(YM(:)>0.5))*3,YRc));
% 
%           [YRcr,resr] = cat_vol_resize(YRc,'reduceV',vx_vol,vx_vol*4,16,'meanm');
%           YRcr  = cat_vol_approx(YRcr,'nn',resr.vx_volr(1),1);
%           YRcr  = cat_vol_smooth3X(YRcr,2/mean(resr.vx_volr));
%           YRc   = cat_vol_resize(YRcr,'dereduceV',resr,'linear');
%           lowf  = 1; %0.5 + 0.5*min(1,max(0,mean(2 - vx_vol))); % reduce filtering on anatical frequency
%           Ys    = Yi + (Ys - Yi) .* (YM>0) .* max(0.2,min(1,YRc .* lowf));  % local filter strength weighting
%         end
%       else
%         if isinf(opt.NCstr) || sign(opt.NCstr)==-1
%           Yh     = Ys>mean(Ys(:)); % object
%           Tth    = mean(Ys(Yh(:)));
%           NCstr  = -min(1,max(0,cat_stat_nanmean(abs(Ys(Yh(:)) - Yi(Yh(:)))) * 15 * min(1,max(0,abs(opt.NCstr))) )); 
%           NC     = min(2,abs(Ys - Yi) ./ max(eps,Ys) * 15 * 2 * min(1,max(0,abs(opt.NCstr)))); 
%           NCs    = NC+ 0; spm_smooth(NCs,NCs,2); NCs = NCs .* cat_stat_nanmean(NCs(Yh(:))) / cat_stat_nanmean(NC(Yh(:)));
%           NCs  = max(0,min(1,NCs));      
%         else 
%           NCstr  = opt.NCstr;
%         end
%         NCstr = max( mean(NCstr(:)) - 2*std(NCstr(:)), NCstr); 
%         
%         %% mix original and noise corrected image
%         if opt.local
%           Ys = Yi.*(1-NCs) + Ys.*NCs; 
%         else
%           Ys = Yi*(1-NCstr) + Ys*NCstr; 
%         end
%       end
%    end



    %% prepare next iteration 
    Ys(~YM) = Y(~YM); 
    Y  = Ys;  
    Yi = Y .* (YM>0);
    iter = iter + 1; 
    
  end
  clear Y YM ; 
  
  if opt.level==1 %&& any(vx_vol<0.75)
    if opt.verb, cat_io_cprintf('g5',sprintf('%3d.%d) %0.2fx%0.2fx%0.2f mm:  ',opt.level,iter,vx_vol)); stime = clock; end
    
    cat_sanlm(Ys,3,1,opt.rician); 

    YM     = Ys>mean(Ys(:)); % object
    if opt.NCstr<0
    % adaptive local denoising 
        NCstr = opt.NCstr * 15; 

        % prepare local map
        NCs = abs(Ys - Yi) ./ max(eps,Ys); spm_smooth(NCs,NCs,2);
        NCs = max(0,min(1,NCs * abs(NCstr))); 
        opt.NCstr = -cat_stat_nanmean(NCs(:)); 

        % mix original and noise corrected image
        Ys = Yi.*(1-NCs) + Ys.*NCs; 
        clear NCs;

    elseif opt.NCstr>0
    % (adaptive) global denoising  

        NCstr = opt.NCstr * 15; 
     
        opt.NCstr = min(1,max(0,opt.NCstr));

        % mix original and noise corrected image
        Ys   = Yi*(1-opt.NCstr) + Ys*opt.NCstr; 
    end
    noiser = 1 - (cat_stat_nanmean(abs(Yi(YM(:))-Ys(YM(:)))./max(Tth*0.2,Ys(YM(:))))/sqrt(prod(vx_vol))) / noise;
    if noiser<0, noiser = noiser+1; end
    noise  = cat_stat_nanmean(abs(Yi(YM(:))-Ys(YM(:)))./max(Tth*0.2,Ys(YM(:))))/sqrt(prod(vx_vol));
      
    if opt.verb, cat_io_cprintf('g5',sprintf('  noise = %4.3f, noiser = %4.3f        %5.0fs\n',noise,noiser,etime(clock,stime))); end
    clear YM;
    
    %{
    % adaptive global denoising 
    if isinf(opt.NCstr) || sign(opt.NCstr)==-1;
      YM     = Ys>mean(Ys(:)); % object
      Tth    = mean(Ys(YM(:)));
      NCstr  = min(1,max(0, mean( abs(Ys(YM(:)) - Yi(YM(:))) ./ Tth ) * 16 * min(10,max(0,abs(opt.NCstr))) ));
      NCs    = abs(Ys - Yi) ./ max(eps,Ys) * 16 * min(2,max(0,abs(opt.NCstr))); spm_smooth(NCs,NCs,2);
      NCs    = max(0,min(1,NCs));
    else 
      NCstr  = opt.NCstr;
    end
    NCstr = max( mean(NCstr(:)) - 1*std(NCstr(:)), min ( mean(NCstr(:)) + 1*std(NCstr(:)) , NCstr)); 
    
    noiser = 1 - (cat_stat_nanmean(abs(Yi(YM(:))-Ys(YM(:)))./max(Tth*0.2,Ys(YM(:))))/sqrt(prod(vx_vol))) / noise;
    if noiser<0, noiser = noiser+1; end
    noise  = cat_stat_nanmean(abs(Yi(YM(:))-Ys(YM(:)))./max(Tth*0.2,Ys(YM(:))))/sqrt(prod(vx_vol));
      
    if opt.verb, fprintf('  noise = %4.3f, noiser = %4.3f        %5.0fs\n',noise,noiser,etime(clock,stime)); end
    clear YM;   
    
    % mix original and noise corrected image
    if opt.local
      Ys = Yi.*(1-NCs) + Ys.*NCs; 
    else
      Ys = Yi*(1-NCstr) + Ys*NCstr; 
    end
    clear Yi; 
    %}
    
  end

  % garantie positive values
  if min(Ys(:))>-0.001 && sum(Ys(:)<0)>0.01*numel(Ys(:));
    Ys = Ys - min(Ys(:)); 
  end
  
  if opt.fast
    Y0(BB.BB(1):BB.BB(2),BB.BB(3):BB.BB(4),BB.BB(5):BB.BB(6))=Ys; Ys = Y0;
  end
  
 
end
