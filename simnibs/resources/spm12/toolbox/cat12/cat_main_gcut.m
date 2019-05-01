function [Yb,Yl1] = cat_main_gcut(Ysrc,Yb,Ycls,Yl1,YMF,vx_vol,opt)
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
% $Id: cat_main_gcut.m 1152 2017-07-05 12:28:08Z dahnke $


  NS   = @(Ys,s) Ys==s | Ys==s+1;
  LAB  = opt.LAB; 
  voli = @(v) (v ./ (pi * 4./3)).^(1/3);           % volume > radius
  brad = double(voli(sum(Yb(:)>0).*prod(vx_vol))); % distance and volume based brain radius (brad)
  Yp0  = single(Ycls{3})/255/3 + single(Ycls{1})/255*2/3 + single(Ycls{2})/255;
  rvol = [sum(round(Yp0(:)*3)==1), sum(round(Yp0(:)*3)==2), sum(round(Yp0(:)*3)==3)]/sum(round(Yp0(:)*3)>0);
  %noise   = cat_stat_nanstd(Ym(cat_vol_morph(cat_vol_morph(Ym>0.95 & Ym<1.05,'lc',1),'e')));
  
  %% set different paremeters to modifiy the stength of the skull-stripping 
  %gc.n = max(0.05,min(0.1,noise));
  % intensity parameter
  gc.h = 3.5 - 0.2*opt.gcutstr + 0.2*opt.LASstr; % 3.25; upper tissue intensity (WM vs. blood vessels)     - higher > more "tissue" (blood vessels)
  gc.g = 1.9 + 0.1*opt.gcutstr; % 1.50; lower tissue intensity (WM vs. blood vessels)     - higher > more "tissue" (blood vessels)
  gc.l = 1.1 + 0.8*opt.gcutstr; % 1.50; lower tissue intensity (WM vs. blood vessels)     - higher > more "tissue" (blood vessels)
  gc.o = 0.2 + 0.8*opt.gcutstr; % 0.50; BG tissue intensity (for high contrast CSF=BG=0!) - lower value > more "tissue"
  % distance parameter
  gc.d = brad*(5 - 4*opt.gcutstr)/mean(vx_vol);               % 3.0;    distance  parameter for downcut - higher > more tissue
  gc.c = max(-0.01,(0.01 - 0.03*opt.gcutstr)*mean(vx_vol));             % -0.005; growing   parameter for downcut - higher > more tissue
  gc.f = max(1,min(4,(brad/200 / (0.7-0.4*opt.gcutstr) * rvol(1)/0.10)/mean(vx_vol))); % closing parameter   - higher > more tissue ... 8
  gc.gd = 1 + 2*opt.gcutstr;
  gc.bd = 3 + 2*opt.gcutstr; 
  % smoothing parameter
  gc.s = 0.2  + 0.30*min(0.6,opt.gcutstr);                   % 0.5;    smoothing parameter   - higher > less tissue
  
  if opt.verb, fprintf('\n'); end
  stime = cat_io_cmd('  WM initialisation','g5','',opt.verb); dispc=1;
  %% init: go to reduces resolution 
  [Ym,Yl1,YMF,BB] = cat_vol_resize({Ysrc,Yl1,YMF},'reduceBrain',vx_vol,round(4/mean(vx_vol)),Yb);
  Yp0  = single(Ycls{3})/255/3 + single(Ycls{1})/255*2/3 + single(Ycls{2})/255;
  [Yp0,Ywm,Ygm,Ycsf,Ymg,Yb] = cat_vol_resize({Yp0,single(Ycls{2})/255,single(Ycls{1})/255,...
    single(Ycls{3})/255,single(Ycls{5})/255,Yb},'reduceBrain',vx_vol,round(4/mean(vx_vol)),Yb);
  vxd  = max(1,1/mean(vx_vol)); 
  Ymg = Ymg>0.05 & Ym<0.45; 
  
  clear Ycls
  Ybo=Yb;
  %% initial WM+ region
  Yb=Ybo;
  YHDr = cat_vol_morph(Yl1>20 | Yl1<=0,'e',vxd*1);
  YGD  = cat_vbdist(max(0,1-Yp0),true(size(Yb)),vx_vol);   % something like the GWM depth/thickness
  YBD  = cat_vbdist(max(0,1-Yp0*3),true(size(Yb)),vx_vol./mean(vx_vol)); % brain depth, (simple) sulcal depth
  Yb   = Yb>0.25 & Ym>2.5/3 & Ym<gc.h/3 & Yl1<21 & Yb & YGD>gc.gd & YBD>gc.bd;  % init WM 
  Yb   = Yb | (Ym>2/3 & Ym<gc.h/3 & (YBD>10 | (YBD>2 & (NS(Yl1,LAB.CB) | NS(Yl1,LAB.HI)))));
  [Ybr,Ymr,resT2] = cat_vol_resize({single(Yb),Ym},'reduceV',1,4./vx_vol,32); 
  Ybr  = Ybr | (Ymr<0.8 & cat_vol_morph(Ybr,'lc',2)); % large ventricle closing
  Ybr  = cat_vol_resize(cat_vol_smooth3X(Ybr,2),'dereduceV',resT2)>0.9; 
  
  % if no largest object could be find it is very likeli that initial normalization failed
  if sum(Yb(:) & mod(Yl1(:),2)==0)==0 || sum(Yb(:) & mod(Yl1(:),2)==1)==0
    error('cat:cat_main:largestWM',['No largest WM cluster could be found: \n'...
      'Please try to set origin (AC) and run preprocessing again \n' ...
      'because it is very likeli that spatial normalization failed.']);
  end
  Yb   = cat_vol_morph(Yb & mod(Yl1,2)==0,'l') | ...
         cat_vol_morph(Yb & mod(Yl1,2)==1,'l') | ...
         (Ybr & Yp0>1.9 & Ym<3.5 & (NS(Yl1,LAB.CB))) | ... 
         (Ybr & Yp0<1.5 & Ym<1.5); 
  Yb  = smooth3(Yb)>gc.s;
  Yb(smooth3(single(Yb))<0.5)=0;                           % remove small dots
  Yb  = single(cat_vol_morph(Yb,'labclose',gc.f));         % one WM object to remove vbs
  
  
  %% region growing GM/WM (here we have to get all WM gyris!)
  stime = cat_io_cmd('  GM region growing','g5','',opt.verb,stime); dispc=dispc+1;
  Yb(~Yb & (YHDr | Ym<gc.g/3 | Ym>gc.h/3 | (Ywm + Ygm)<0.5 | YGD<(gc.gd-1) | YBD<(gc.bd-3)))=nan; %clear Ywm Ygm; 
  [Yb1,YD] = cat_vol_downcut(Yb,Ym,0.03+gc.c); % this have to be not to small... 
  Yb(isnan(Yb) | YD>gc.d*vxd*2)=0; Yb(Yb1>0 & YD<gc.d*vxd*2)=1;
  Yb(smooth3(single(Yb))<gc.s)=0;
  Yb = single(Yb | (cat_vol_morph(Yb,'labclose',vxd) & Ym<1.1));

  
  %% region growing CSF/GM 
  stime = cat_io_cmd('  GM-CSF region growing','g5','',opt.verb,stime); dispc=dispc+1;
  Yb(~Yb & (YHDr | Ym<gc.l/3 | Ym>gc.h/3) | Ymg)=nan; % | YBD<1
  [Yb1,YD] = cat_vol_downcut(Yb,Ym,0.01+gc.c);
  Yb(isnan(Yb) | YD>gc.d/2)=0; Yb(Yb1>0 & YD<gc.d)=1; 
  for i=1:2, Yb(smooth3(single(Yb))<gc.s)=0; end
  Yb  = single(cat_vol_morph(Yb,'o',1));
  Yb  = single(Yb | (cat_vol_morph(Yb ,'labclose',1) & Ym<1.1));
  
  %% region growing - add CSF
  Yb(~Yb & (YHDr | Ym<1/3 | Ym>gc.h/3) | Ymg)=nan; 
  [Yb1,YD] = cat_vol_downcut(Yb,Ym,-0.00+gc.c);
  Yb(isnan(Yb) | YD>gc.d/2)=0; Yb(Yb1>0 & YD<gc.d)=1; 
  for i=1:2, Yb(smooth3(single(Yb))<gc.s)=0; end
  Yb  = single(cat_vol_morph(Yb,'o',1));
  Yb  = single(Yb | (cat_vol_morph(Yb ,'labclose',1) & Ym<1.1));
  Ybox = Yb; 
  %% region growing - add CSF regions   
  stime = cat_io_cmd('  CSF region growing','g5','',opt.verb,stime); dispc=dispc+1;
  Ygr = cat_vol_grad(Ym,vx_vol); CSFth = mean(Ym(Ycsf(:)>0.8 & Ygr(:)<0.1));
  Ybb = smooth3( Ym<CSFth*0.9 | (Ym>1.5/3 & ~Yb) | (Ygr>0.15 & ~Yb))>0.5 | smooth3(Ycsf)<0.5; 
  if std(Ybb(:))>0  % no ROI in low res images 
    if sum(Ybb(:)>0.5)>0 % Ybb is maybe empty 
      Ybb = cat_vol_morph( Ybb>0.5 ,'lc',vxd);
      Yb(~Yb & smooth3( Ybb )>0.6 ) = nan;
    end 
    Yb(isnan(Yb) & cat_vol_morph(Yb>=0,'lc',2))=0;
  end
  [Yb1,YD] = cat_vol_downcut(Yb,smooth3(Ym),+0.01+gc.c); 
  Yb(isnan(Yb) | YD>gc.d/2)=0; Yb(Yb1>0 & YD<gc.d*4 & YD>0)=1;
  for i=1:2, Yb(cat_vol_smooth3X(Yb,2)<max(0.05,(gc.s - 0.2)))=0; end
  Yb  = single(cat_vol_morph(Yb,'o',2));
  Yb = Yb | YMF; 
  
  % smooth / low dilated boundary 
  Ybs = single(Yb); spm_smooth(Ybs,Ybs,4*gc.s./vx_vol); Yb   = Yb | (Ybs>(gc.s-0.2) & Ym<1.25/3);
  
  %% filling of ventricles and smooth mask
  stime = cat_io_cmd('  Ventricle filling','g5','',opt.verb,stime); dispc=dispc+1; %#ok<*NASGU>
  [Ybr,resT3] = cat_vol_resize(single(Yb),'reduceV',vx_vol,min(1,cat_stat_nanmean(vx_vol)*2),32,'meanm'); 
  Ybr = cat_vol_morph(Ybr>0.5,'labclose',gc.f/mean(resT3.vx_volr));
  Ybr = cat_vol_resize(Ybr,'dereduceV',resT3)>0.5; 
  Yb  = Yb | Ybr; clear Ybr;   % & Ym>=gc.o/3 & Ym<1.25/3 & ~Ymg & Ycsf>0.75);
  Yb  = single(cat_vol_morph(Yb,'o',max(1,min(3,4 - 0.2*gc.f* (rvol(1)/0.4) ))./cat_stat_nanmean(vx_vol)));
  Yb  = Yb | (cat_vol_morph(Yb ,'labclose',vxd*2) & Ym<1.1);
  Yb  = cat_vol_morph(Yb ,'labclose');
  Ybs = single(Yb)+0; spm_smooth(Ybs,Ybs,3./vx_vol); Yb = Yb>0.5 | (max(Yb,Ybs)>0.3 & Ym<0.4); % how wide
  Ybs = single(Yb)+0; spm_smooth(Ybs,Ybs,2./vx_vol); Yb = max(Yb,Ybs)>0.4; % final smoothing
 
  %%
  Yb   = cat_vol_resize(Yb  ,'dereduceBrain',BB)>0.5;
  Yl1  = cat_vol_resize(Yl1 ,'dereduceBrain',BB);
    
  %% update Yl1 with Yb
  Yl1(~Yb)  = 0;
  [tmp0,tmp1,Yl1] = cat_vbdist(single(Yl1),Yl1==0 & Yb); clear tmp0 tmp1;

  cat_io_cmd(' ','','',opt.verb,stime); 
%    cat_io_cmd('cleanup',dispc,'',opt.verb);

end