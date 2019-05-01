function [Yb,Yl1] = cat_main_gcut_DEV(Ysrc,Yb,Ycls,Yl1,YMF,vx_vol,opt)
% This is an exclusive subfunction of cat_main.
% ______________________________________________________________________
%
% gcut+: skull-stripping using graph-cut
% ----------------------------------------------------------------------
% This routine use morphological, region-growing and graph-cut methods. 
% It starts from the WM segment and grow for lower tissue intensities.
% Atlas knowledge is used to for separate handling of the cerebellum.
% Because its high frequency structures in combination with strong
% noise or other artifacts often lead to strong underestimations.
%
% There are 4 major parameters:
%   gcutstr - strengh of skull-stripping with str=0 - more tissue, str=1 less tissue
%   vx_res - resolution limit for skull-stripping (default 1.5)
%   gcutCSF 
% Especialy the second parameter controls many subparameters for  
% different tissue thresholds, region-growing, closing and smoothing
% parameters.
% This routine have a strong relation to the previous estimated main
% partition map l1, and the blood vessel correction. Therefore, it is
% maybe useful to move it...
%
%   [Yb,Yl1] = cat_main_gcut(Ysrc,Yb,Ycls,Yl1,YMF,vx_vol,opt)
% 
%   Yb   .. updated brain mask
%   Yl1  .. updated label map
% 
%   Ysrc .. anatomical image
%   Yb   .. initial brain mask
%   Ycls .. SPM tissue classification
%   Yl1  .. CAT atlas map
%   YMF  .. subcortical/ventricular regions (for filling in surf. recon.)
%   vx_vol .. image resolutino
%   opt  .. further options
% 
% ______________________________________________________________________
%
%   Robert Dahnke (robert.dahnke@uni-jena.de)
%   Structural Brain Mapping Group (http://dbm.neuro.uni-jena.de/)
%   Department of Neurology
%   University Jena
% ______________________________________________________________________
% $Id: cat_main_gcut_DEV.m 1072 2016-11-09 15:01:16Z dahnke $

  LAB  = opt.LAB;
  NS   = @(Ys,s) Ys==s | Ys==s+1;
  voli = @(v) (v ./ (pi * 4./3)).^(1/3);           % volume > radius
  brad = double(voli(sum(Yb(:)>0).*prod(vx_vol))); % distance and volume based brain radius (brad)
  %noise   = cat_stat_nanstd(Ym(cat_vol_morph(cat_vol_morph(Ym>0.95 & Ym<1.05,'lc',1),'e')));
  Yp0  = single(Ycls{3})/255/3 + single(Ycls{1})/255*2/3 + single(Ycls{2})/255;
  rvol = [sum(round(Yp0(:)*3)==1), sum(round(Yp0(:)*3)==2), sum(round(Yp0(:)*3)==3)]/sum(round(Yp0(:)*3)>0);
  vxd  = max(1,1/mean(vx_vol)); 
  
  %% set different paremeters to modifiy the stength of the skull-stripping 
  %gc.n = max(0.05,min(0.1,noise));
  % intensity parameter
  gc.h = max(3.1,3.4 - 0.2*opt.gcutstr - 0.4*(opt.LASstr-0.5)); % 3.25, upper tissue intensity (WM vs. blood vessels)     - higher > more "tissue" (blood vessels)
  gc.l = 1.4  + 0.20*opt.gcutstr; % 1.50, lower tissue intensity (WM vs. blood vessels)     - higher > more "tissue" (blood vessels)
  gc.t = 1.3  - 0.20*opt.gcutstr; % 1.10, lower tissue intensity (WM vs. blood vessels)     - higher > more "tissue" (blood vessels)
  gc.o = 0.1  + 0.80*opt.gcutstr; % 0.50, BG tissue intensity (for high contrast CSF=BG=0!) - lower value > more "tissue"
  % distance parameter
  gc.d = brad*(8 - 6*opt.gcutstr)/mean(vx_vol);               % 3; distance  parameter for downcut - higher > more tissue
  gc.c = (0.020 - 0.04*opt.gcutstr)*mean(vx_vol);             % -0.015; growing   parameter for downcut - higher > more tissue
  gc.f = (brad/20 * opt.gcutstr * rvol(1)/0.10)/mean(vx_vol); % closing parameter             - higher > more tissue ... 8
  % smoothing parameter
  gc.s  = 0.3 + 0.40*opt.gcutstr;                             % smoothing parameter           - higher > less tissue
  gc.ss = 1.0 + 2.00*opt.gcutstr; 
  gc.lc = round(vxd*(1.9 - opt.gcutstr));
  gc.lo = round(vxd*(0.9 + opt.gcutstr));
  
  
  if opt.verb, fprintf('\n'); end
  stime = cat_io_cmd('  WM initialisation','g5','',opt.verb); dispc=1;
  %% init: go to reduces resolution 
  [Ym,Yl1,YMF,BB] = cat_vol_resize({Ysrc,Yl1,YMF},'reduceBrain',vx_vol,round(4/mean(vx_vol)),Yb);
  [Ywm,Ygm,Ycsf,Ymg,Yb] = cat_vol_resize({single(Ycls{2})/255,single(Ycls{1})/255,...
    single(Ycls{3})/255,single(Ycls{5})/255,Yb},'reduceBrain',vx_vol,round(4/mean(vx_vol)),Yb);
  Ymg = Ymg>0.05 & Ym<0.45; 
  
  clear Ycls
  Ybo=Yb;
  
  %% initial WM+ region
  Yb=Ybo; 
  YHDr = cat_vol_morph(Yl1>20 | Yl1<=0,'e',vxd*2);
  [Ybr,resT2] = cat_vol_resize(single(Yb),'reduceV',vx_vol,mean(vx_vol)*4,32); 
  Ybr = single(cat_vol_morph(Ybr>0,'e',brad/25));
  Ybr = cat_vol_resize(Ybr,'dereduceV',resT2)>0.5; 
  Yb  = Yb>0.25 & Ym>2.5/3 & Ym<gc.h/3 & Yl1<21 & Yb;  % init WM 
  Yb  = cat_vol_morph(Yb,'l'); 
  
  % if no largest object could be find it is very likeli that initial normalization failed
  if sum(Yb & mod(Yl1,2)==0)==0 || sum(Yb & mod(Yl1,2)==1)==0
    error('cat:cat_main:largestWM',['No largest WM cluster could be found: \n'...
      'Please try to set origin (AC) and run preprocessing again \n' ...
      'because it is very likeli that spatial normalization failed.']);
  end
  
  Yb  = single(Yb | (Ym>2.5/3  & Ym<gc.h/3 & Yb) | NS(Yl1,LAB.VT) | ...
    (cat_vol_morph(NS(Yl1,LAB.CB) | Ybr | NS(Yl1,LAB.HI) | ...
      NS(Yl1,LAB.HC) | NS(Yl1,LAB.BG),'e',1) & Ym>1.9/3 & Ym<1.1));     % init further WM 
  spm_smooth(Yb,Yb,gc.ss./vx_vol); Yb = Yb>(gc.s-0.2);                        % remove small dots
  Yb = single(cat_vol_morph(Yb>0,'lc'));
  
 
  
  %% region growing GM/WM (here we have to get all WM gyris!)
  Ybs = Yb;
  stime = cat_io_cmd('  GM region growing','g5','',opt.verb,stime); dispc=dispc+1;
  Yb(~Yb & (YHDr | Ym<(min(2.1,gc.l+0.5))/3 | Ym>gc.h/3 | (Ywm + Ygm)<0.5))=nan; %clear Ywm Ygm; 
  [Yb1,YD] = cat_vol_downcut(Yb,Ym,max(0,0.01+gc.c)); % this have to be not to small... 
  Yb(isnan(Yb) | YD>gc.d*vxd*2)=0; Yb(Yb1>0 & YD<gc.d*vxd*2)=1;
  for i=1:1, spm_smooth(Yb,Ybs,gc.ss./vx_vol); Yb(Ybs<(gc.s-0.1))=0; end
  Yb = single(Yb | (cat_vol_morph(Yb,'labclose',gc.lc*vxd) & Ym<gc.h/3));
   
  
   
  %% region growing CSF/GM 
  stime = cat_io_cmd('  GM-CSF region growing','g5','',opt.verb,stime); dispc=dispc+1;
  Yb(~Yb & (YHDr | Ym<gc.l/3 | Ym>gc.h/3) | Ymg)=nan;
  [Yb1,YD] = cat_vol_downcut(Yb,Ym,-0.02+gc.c*0.5);
  Yb(isnan(Yb) | YD>gc.d/2)=0; Yb(Yb1>0 & YD<gc.d)=1; 
  for i=1:2, spm_smooth(Yb,Ybs,gc.ss./vx_vol); Yb(Ybs<(gc.s-0.1))=0; end
  Yb  = single( (Yb & NS(Yl1,LAB.CB)) | cat_vol_morph(Yb,'o',max(1,min(2,4 - 0.2*gc.f* (rvol(1)/0.4) ))));
  Yb  = single(Yb | (cat_vol_morph(Yb ,'labclose',gc.lc) & Ym<gc.h/3));
  
  % Yb   = cat_vol_resize(Yb  ,'dereduceBrain',BB)>0.5; return 
  
  %% region growing - add CSF
  Yb(~Yb & (YHDr | Ym<1/3 | Ym>gc.h/3) | Ymg)=nan;
  [Yb1,YD] = cat_vol_downcut(Yb,Ym,-0.04+gc.c*0.2);
  Yb(isnan(Yb) | YD>gc.d/2)=0; Yb(Yb1>0 & YD<gc.d)=1; 
  for i=1:2, spm_smooth(Yb,Ybs,gc.ss./vx_vol); Yb(Ybs<gc.s)=0; end
  Yb  = single(Yb | (cat_vol_morph(Yb ,'labclose',gc.lc) & Ym<gc.t/3 & Ym>gc.o/3));
  Yb  = cat_vol_morph(Yb ,'labopen',gc.lo);
 
 
   
  %% region growing - add CSF regions   
  stime = cat_io_cmd('  CSF region growing','g5','',opt.verb,stime); dispc=dispc+1;
  Ygr = cat_vol_grad(Ym,vx_vol);
  Yb(~Yb & smooth3(cat_vol_morph(smooth3(Ym<0.75/3 | (Ym>1.25/3 & ~Yb) | ...
    (Ygr>0.05 & ~Yb))>0.5,'lc',vxd*2) | Ymg )>0.5)=nan; 
  [Yb1,YD] = cat_vol_downcut(Yb,Ym,-0.02+gc.c); 
  Yb(isnan(Yb) | YD>gc.d/2)=0; Yb(Yb1>0 & YD<gc.d*2 & YD>0)=1;
  for i=1:2, spm_smooth(Yb,Ybs,gc.ss./vx_vol); Yb(Ybs<(gc.s - 0.25))=0; end
  Yb = single(Yb | YMF); 
  
  % smooth / low dilated boundary 
  Ybs = single(Yb)+0; spm_smooth(Ybs,Ybs,2*gc.s./vx_vol); Yb = Yb>0.5 | (Ybs>(gc.s-0.1) & Ym<gc.h/3);
  
  %% filling of ventricles and smooth mask
  stime = cat_io_cmd('  Ventricle filling','g5','',opt.verb,stime); dispc=dispc+1; %#ok<*NASGU>
  Yb  = Yb | (cat_vol_morph(Yb ,'labclose',vxd*gc.f) & ...
    Ym>=gc.o/3 & Ym<1.25/3 & ~Ymg & Ycsf>0.75);
  Yb  = single(cat_vol_morph(Yb,'o',max(1,min(3,4 - 0.2*gc.f* (rvol(1)/0.4) ))));
  Yb  = Yb | (cat_vol_morph(Yb ,'labclose',vxd) & Ym<1.1);
  %%
  Ybs = single(Yb)+0; spm_smooth(Ybs,Ybs,3./vx_vol); Yb = Yb>0.5 | (max(Yb,Ybs)>(gc.s-0.1) & Ym<0.4); % how wide
  Ybs = single(Yb)+0; spm_smooth(Ybs,Ybs,2./vx_vol); Yb = max(Yb,Ybs)>0.4; % final smoothing
 
  %%
  Yb   = cat_vol_resize(Yb  ,'dereduceBrain',BB)>0.5;
  Yl1  = cat_vol_resize(Yl1 ,'dereduceBrain',BB);
    
  %% update Yl1 with Yb
  Yl1(~Yb)  = 0;
  [tmp0,tmp1,Yl1] = cat_vbdist(single(Yl1),Yl1==0 & Yb); clear tmp0 tmp1;

  if opt.debug
    cat_io_cmd(' ','','',opt.verb,stime); 
  else
    cat_io_cmd(' ','','',opt.verb,stime);   
%    cat_io_cmd('cleanup',dispc,'',opt.verb);
  end

end