function [Yml,Ymg,Ycls,Ycls2,T3th] = cat_main_LAS(Ysrc,Ycls,Ym,Yb0,Yy,T3th,res,vx_vol,extopts,Tth) 
% This is an exclusive subfunction of cat_main.
% ______________________________________________________________________
%
% Local Adaptive Segmentation (LAS):
%
% This version of the local adaptive intensity correction includes a  
% bias correction that is based on a maximum filter for the WM and a mean
% filter of GM to stabilize the correction in region with less WM.
% The extension is mostly based on the assumption that the tissue next to 
% the CSF (and high divergence sulci) has to be WM (maximum, high 
% divergence) or GM. For each tissue a refined logical map is generated 
% and used to estimate the local intensity threshold.
% It is important to avoid high intensity blood vessels in the process, 
% because they will push down local WM and GM intensity - due to the CSF
% near possition of blood vessels that mostly push down GM. 
% Based on these values an intensity transformation is used. Compared to 
% the global correction this has to be done for each voxel. To save time
% only a rough linear transformation is used.
%
% Finally, a second NLM-filter is used and a refinement of WM structures
% by a divergence map 
% ______________________________________________________________________
%
%   [Yml,Ycls,Ycls2,T3th] = ...
%     cat_main_LAS(Ysrc,Ycls,Ym,Yb0,Yy,T3th,res,vx_vol,PA,template)
%
%   Yml    .. local intensity correct image
%   Ycls   .. corrected SPM tissue class map
%   Ycls2  .. ?
%   T3th   .. tissue thresholds of CSF, GM, and WM in Ysrc
%
%   Ysrc   .. (bias corrected) T1 image
%   Ym     .. intensity corrected T1 image (BG=0,CSF=1/3,GM=2/3,WM=1)
%   Ycls   .. SPM tissue class map
%   Yb0    .. brain mask
%   Yy     .. deformation map
%   res    .. SPM segmentation structure
%   vx_vol .. voxel dimensions
%   PA     .. CAT atlas map
%   template .. ?
% ______________________________________________________________________
% 
% internal maps:
%
%   Yg   .. gradient map   - edges between tissues
%   Ydiv .. divergence map - sulci, gyris pattern, and blood vessels
%   Yp0  .. label map      - tissue classes (BG=0,CSF=1,GM=2,WM=3) 
%
%   Ysw  .. save WM tissue map
%   Ybv  .. blood vessel map
%   Ycp  .. CSF / background area for distances estimation
%   Ycd  .. CSF / background distance
%   Ycm  .. CSF
%   Ygm  .. GM
%   Ywm  .. WM 
%   Yvt  .. WM next to the ventricle map 
% ______________________________________________________________________
%
%   Robert Dahnke (robert.dahnke@uni-jena.de)
%   Structural Brain Mapping Group (http://dbm.neuro.uni-jena.de/)
%   Department of Neurology
%   University Jena
% ______________________________________________________________________
% $Id: cat_main_LAS.m 1145 2017-06-20 12:18:59Z dahnke $

  % set this variable to 1 for simpler debuging without reduceBrain
  % function (that normally save half of processing time)
  verb    = extopts.verb-1;
  vxv     = 1/ mean(vx_vol);
  dsize   = size(Ysrc);
  NS      = @(Ys,s) Ys==s | Ys==s+1;                                    % function to ignore brain hemisphere coding
  LASstr  = max(eps,min(1,extopts.LASstr));  % LAS strenght (for GM/WM threshold)3 - manual correction based on R1109 (2017/02)
  LAB     = extopts.LAB;                         % atlas labels
  LABl1   = 1;                                                          % use atlas map
  cleanupstr  = min(1,max(0,extopts.gcutstr));   % required to avoid critical regions
  cleanupdist = min(3,max(1,1 + 2*cleanupstr));
    
  % set debug = 1 and do not clear temporary variables if there is a breakpoint in this file 
  dbs   = dbstatus; debug = 0; 
  for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,'cat_main_LAS'); debug = 1; break; end; end
  
  
%% ---------------------------------------------------------------------
%  First, we have to optimize the segments using further information that 
%  SPM do not use, such as the gradient, divergence and distance maps. 
%  The gradient map (average of the first derivate of the T1 map) is an 
%  edge map and independent of the image intensity. It helps to avoid PVE 
%  regions and meninges. 
%  The divergence (second derivate of the T1 map) help to identfiy sulcal
%  and gyral pattern and therefore to find WM and CSF regions for furhter 
%  corrections and to avoid meninges and blood vessels. 
%  Furhtermore, special assumption can be used. 
%  The first one is the maximum property of the WM in T1 data that allows
%  using of a maxim filter for the GM/WM region. 
%  The second is the relative stable estimation of CSF/BG that allows to 
%  estimat a distance map. Because, most regions have a thin layer of 
%  GM around the WM we can avoid overestimation of the WM by the other 
%  maps (especially the divergence). 
%  ---------------------------------------------------------------------
  fprintf('\n');
  stime = cat_io_cmd('  Prepare maps','g5','',verb); dispc=1;

  
  % brain segmentation can be restricted to the brain to save time 
  
  Yclso=Ycls; Ysrco=Ysrc;
  [Ysrc,Ym,Yb,BB] = cat_vol_resize({Ysrc,Ym,Yb0},'reduceBrain',vx_vol,round(10/mean(vx_vol)),Yb0);
  for i=1:6, Ycls{i} = cat_vol_resize(Ycls{i},'reduceBrain',vx_vol,BB.BB); end
  
  
  % helping maps (Yg = mean gradient = edge) and divergence 
  Yg    = cat_vol_grad(Ym,vx_vol);
  Ydiv  = cat_vol_div(max(0.33,Ym),vx_vol);
  Yp0   = single(Ycls{1})/255*2 + single(Ycls{2})/255*3 + single(Ycls{3})/255;
  Yb    = smooth3(Yb | cat_vol_morph(Yb,'d',2*vxv) & Ym<0.8 & Yg<0.3 & Ym>0)>0.5; % increase brain mask, for missing GM 
  
  
  %% adding of atlas information (for subcortical structures)
  %  -------------------------------------------------------------------
  if LABl1 
    stime = cat_io_cmd('  Prepare partitions','g5','',verb,stime); dispc=dispc+1;

    % map atlas to RAW space
    for i=1:5
      try
        Vl1A = spm_vol(extopts.cat12atlas{1});
        break
      catch 
        % read error in parallel processing
        pause(1)
      end
    end
    Yl1  = cat_vol_ctype(round(spm_sample_vol(Vl1A,double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),0)));
    Yl1  = reshape(Yl1,dsize);
    
    % load WM of the TPM or Dartel/Shooting Template for WMHs
    %Ywtpm = cat_vol_ctype(spm_sample_vol(res.tpm(2),double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),0)*255,'uint8');
    Vtemplate = spm_vol(extopts.templates{end}); 
    Ywtpm = cat_vol_ctype(spm_sample_vol(Vtemplate(2),double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),0)*255,'uint8');
    Ywtpm = reshape(Ywtpm,dsize); spm_smooth(Ywtpm,Ywtpm,2*vxv);
    Ywtpm = single(Ywtpm)/255;
    Yl1    = cat_vol_resize(Yl1  ,'reduceBrain',vx_vol,round(4/mean(vx_vol)),BB.BB);
    Ywtpm  = cat_vol_resize(Ywtpm,'reduceBrain',vx_vol,round(4/mean(vx_vol)),BB.BB);
    if ~debug, clear Yy; end
    
    LASmod = min(2,max(1,mean((Ym( NS(Yl1,LAB.BG) & Yg<0.1 & Ydiv>-0.05  & Ycls{1}>4)) - 2/3) * 8)); % do not reduce LASstr 
  end
  
  
  %% adaption of the LASstr depending on average basal values 
  LASstr  = min(1,max(0.01,LASstr * LASmod));                           % adaption by local BG variation
  LASfs   = 1 / max(0.01,LASstr);                                       % smoothing filter strength 
  LASi    = min(8,round(LASfs));                                        % smoothing interation (limited)
   
  
  %% helping segments
  stime = cat_io_cmd(sprintf('  Prepare segments (LASmod = %0.2f)',LASmod),'g5','',verb,stime); dispc=dispc+1;
  % Ybb = don't trust SPM to much by using Yp0 because it may miss some areas! Shood be better now with MRF.
  Ybb = cat_vol_morph((Yb & Ym>1.5/3 & Ydiv<0.05) | Yp0>1.5,'lo',vxv);
  % Ysw = save WM and blood vessels mpas
  % Ybv = possible blood vessels
  Ysw = cat_vol_morph(Ycls{2}>128 & (min(1,Ym)-Ydiv)<1.5,'lc',vxv*2) & (Ym-Ydiv)>5/6; % 1.2 
  Ybv = ((min(1,Ym) - Ydiv + Yg)>2.0 | (Ycls{5}>16 & Ym<0.6 & Ycls{1}<192)) & ...
        ~cat_vol_morph(Ysw,'d',1) & Ym>0.2;              
  % Ycp = for CSF/BG distance initialization 
  Ycp = (Ycls{3}>240 & Ydiv>0 & Yp0<1.1 & Ym<0.5) | ...                 % typcial CSF
        (Ycls{5}>8 & Ycls{2}<32 & Ym<0.6 & Ydiv>0) | ...                % venes
        ((Ym-Ydiv/4<0.4) & Ycls{3}>4 & Ycls{3}>16) | ...                  % sulcal CSF
        ...(single(Ycls{6})+single(Ycls{5})+single(Ycls{4}))>192 | ...     % save non-csf 
        (single(Ycls{6})+single(Ycls{4}))>192 | ...     % save non-csf .. class 5 with error for negative t1 values
        ~cat_vol_morph(Ybb,'lc',5) | ...                                % add background
        Ym<0.3;                                                         % but do not trust the brain mask!
  Ywd = cat_vbdist(single(Yp0>2.5),Yp0>0.5,vx_vol);                         % WM distance for skelelton
  %Ysk = cat_vol_div(min(5,Ywd),2); %clear Ywd;                      % divergence skeleton
  %Ysk = (Ym + min(0,Ysk))<0.2;                                          % binary divergence skeleton
  %Ycp = Ycp | Ysk; %clear Ysk;                                          % 
  Ycp(smooth3(Ycp)>0.4)=1;                                              % remove some meninges
  Ycd = cat_vbdist(single(Ycp),~Ycp,vx_vol);                                % real CSF distance 
  Ycd((Ym-Ydiv<2/3 | Ydiv>0.1) & Ycls{3}>4 & Ycls{3}>1) =  ... correction for sulci ... maybe a second distance estimation??=
    min(Ycd((Ym-Ydiv<2/3 | Ydiv>0.1) & Ycls{3}>4 & Ycls{3}>1),1.5);
  % we need to remove strong edge regions, because here is no GM layer between CSF and WM ???  
  Yb  = cat_vol_morph(~Ycp | (Ycls{3}>128),'lc',1);
  Ybd = cat_vbdist(single(~Yb),Yb,vx_vol);
  Yvt = (Yg+abs(Ydiv))>0.4 & smooth3(single(Ycls{1})/255)<0.5 & Ybd>20 & ...
    cat_vol_morph(Ycls{3}>8,'d',vxv) & cat_vol_morph(Ycls{2}>8,'d',vxv); 
  Yvt = smooth3(Yvt)>0.7;
  Yvt = smooth3(Yvt)>0.2;
  
  
  %% final tissue maps:  Ycm = CSF, Ygm = GM, Ywm = WM 
  Ysc = Ycp & Yb & Ycls{3}>192 & ~Ybv & Ym<0.45 & Yg<0.1;
  Ycm = Ycp & Yb & Ycls{3}>192 & ~Ybv & (Yb | Ym>1/6) & Ym<0.45 & Yg<0.25 & Ym>0; % & Ydiv>-0.05;
  %Ycm = Ycm | (Yb & (Ym-max(0,Ydiv))<0.5); 
  Ywm = (Ysw | Ycls{2}>252 | ((Ycd-Ydiv)>2 & Ydiv<0 & Ym>0.9+LASstr*0.05 & Yb) | ... % save WM 
        ((Ycd-Ydiv.*Ycd)>4 & (Ydiv<-0.01) & Yb & Ym>0.5 & Ybd<20 & Ycd>2) ) & ...
        ... ((Ycd-Ydiv*5)>3 & (Ydiv<-0.01 & (Yg + max(0,0.05-Ycd/100))<0.1) & Yb & Ym>0.4 & Ybd<20 & Ycd>2.5) ) & ... % further WM
        ~Ybv & Yb & Ybd>1 & (Ycd>1.0 | (Yvt & Yp0>2.9)) & (Yg+Ydiv<(Ybd/50) | (Ydiv-Ym)<-1); % Ybd/800 + Ycd/50
  Ygm = ~Yvt & Ybb & ~Ybv & ~Ywm & ~Ycm & Ycd>0.5 & (Ym-Ydiv-max(0,2-Ycd)/10)<0.9 & ... (Ym+Ydiv)>0.5 & ... ~Ysk & 
        (Ycls{1}>4 | (Ym>0.7 & Ycls{3}>64) | Ycd<(Ym+Ydiv)*3 ) & ...
        smooth3(Yg>(Ybd/800) & Ycls{2}<240 )>0.6; % avoid GM next to hard boundies in the middle of the brain
  Ygx = Ybb & ~Ycm & ~Ywm & Ym>1/3 & Ym<2.8/3 & Yg<0.4 & (Ym-Ydiv)>1/3 & (Ym-Ydiv)<1; Ygx(smooth3(Ygx)<0.5) = 0;
  Ygm = Ygm | Ygx; clear Ygx;
  Ygm = Ygm | (Ym>1.5/3 & Ym<2.8/3 & ~Ywm & ~Ycm & Ybb);
  Ygm(smooth3(Ygm)<0.25)=0;
  %Ygw = Ygm & smooth3(Ywm)<0.1 & smooth3(Ycm)<0.4 & Ycd>0 & Ycd<2 & Ydiv<0.4 & Ydiv>-0.3 & Yg<0.1; %& (Ydiv>-0.4 | Ycd>1.5)
  if ~debug, clear Ybv  Ycp; end %Ycd

  if debug>1
    try %#ok<TRYNC> 
      [pth,nam] = spm_fileparts(res.image0(1).fname); tpmci=0;
      tpmci=tpmci+1; tmpmat = fullfile(pth,reportfolder,sprintf('%s_%s%02d%s.mat',nam,'LAS',tpmci,'prepeaks'));
      save(tmpmat);
    end
  end
  %%

  if LABl1   
    %% ------------------------------------------------------------------
    % SPM GM segmentation can be affected by inhomogeneities and some GM
    % is missclassified as CSF/GM (Ycls{5}). But for some regions we can 
    % trust these information more
    % ------------------------------------------------------------------
    Ybd  = cat_vbdist(single(~Yb),Yb,vx_vol);
    Ycbp = cat_vbdist(single(NS(Yl1,LAB.CB)),Yb,vx_vol);                    % next to the cerebellum
    Ycbn = cat_vbdist(single(~NS(Yl1,LAB.CB)),Yb,vx_vol);                   % not to deep in the cerebellum
    Ylhp = cat_vbdist(single(mod(Yl1,2)==1 & Yb & Yl1>0),Yb,vx_vol);        % GM next to the left hemisphere 
    Yrhp = cat_vbdist(single(mod(Yl1,2)==0 & Yb & Yl1>0),Yb,vx_vol);        % GM next to the righ hemishpere
    Ybv2 = Ycls{5}>2 & Ym<0.7 & Ym>0.3 & Yb & (... 
           ((Ylhp+Ybd/2)<cleanupdist*6 & (Yrhp+Ybd/2)<cleanupdist*6) | ... % between the hemispheres next to skull                 
           ((Ycbp+Ybd/2)<cleanupdist*8 & (Ycbn+Ybd/2)<cleanupdist*8));     % between cerebrum and cerebellum next to hull
    Ybv2 = smooth3(Ybv2)>0.5;
    Ybvv = (Ym-max(0,6-abs(Ycbp-6))/50)<0.6 & Ym>0.4 & Yb & Ycbp<8 & Ycbp>1;
    
    %% subcortical map refinements
    THth = 0.8 - LASstr*0.6; %0.5; % lower more thalamus
    YTH = NS(Yl1,LAB.TH) | (cat_vol_morph(NS(Yl1,LAB.TH),'d',3) & Ym>0.5 & Ycls{1}>128);
    Ytd = cat_vbdist(single(Ym<0.45),YTH | NS(Yl1,LAB.BG),vx_vol); Ytd(Ytd>2^16)=0; % CSF distance in the TH
    Yxd = cat_vbdist(single(NS(Yl1,LAB.BG)),YTH,vx_vol); Yxd(Yxd>2^16)=0; % BG distance in the TH
    %Yyd = cat_vbdist(single(NS(Yl1,LAB.TH)),NS(Yl1,LAB.BG),vx_vol); Yyd(Yyd>2^16)=0; % TH distance in the BG
    Yss = NS(Yl1,LAB.BG) | NS(Yl1,LAB.TH); 
    Yss = Yss | (cat_vol_morph(Yss,'d',vxv*2) & Ym>2.25/3 & Ym<2.75/3 & Ydiv>-0.01); % add ihger tissue around mask
    Yss = Yss | (cat_vol_morph(Yss,'d',vxv*3) &  NS(Yl1,LAB.VT) & Yp0>1.5 & Yp0<2.3); % add lower tissue around mask
    Yss = Yss & Yp0>1.5 & (Yp0<2.75 | (Ym<(2.5+LASstr*0.45)/3 & Ydiv>-0.05)); % by intensity
    Yss = Yss | ((Yxd./max(eps,Ytd+Yxd))>THth/2 & (Yp0<2.75 | (Ym<(2.75+LASstr*0.20)/3 & Ydiv>-0.05)));  % save TH by distances - for overcorrected images
    Yss = cat_vol_morph(Yss,'o');
    Ynw = (Yxd./max(eps,Ytd+Yxd))>THth/2 | (NS(Yl1,LAB.BG) & Ydiv>-0.01);
    if ~debug, clear Ytd Yxd ; end
    % increase CSF roi
    Yvt = cat_vol_morph( (NS(Yl1,LAB.VT) | cat_vol_morph(Ycm,'o',3) ) ...
      & Ycm & ~NS(Yl1,LAB.BG) & ~NS(Yl1,LAB.TH) & Ybd>30,'d',vxv*3) & ~Yss; % ventricle roi to avoid PVE GM between WM and CSF
    Ycx = (NS(Yl1,LAB.CB) & ((Ym-Ydiv)<0.55 | Ycls{3}>128)) | (((Ym-Ydiv)<0.45 &  Ycls{3}>8)| Ycls{3}>240);
    % in the crebellum tissue can be differentated by div etc.
    Ycwm = NS(Yl1,LAB.CB) & (Ym-Ydiv*4)>5/6 & Ycd>3 & Yg>0.05;
    Yccm = NS(Yl1,LAB.CB) & Ydiv>0.02 & Ym<1/2 & Yg>0.05;
    Ybwm = (Ym-Ydiv*4)>0.9 & Ycd>3 & Yg>0.05; %Ydiv<-0.04 & Ym>0.75 & Ycd>3;
    Ybcm = Ydiv>0.04 & Ym<0.55 & Yg>0.05;
    %% correction 1 of tissue maps
    Ywmtpm = (Ywtpm.*Ym.*(1-Yg-Ydiv).*cat_vol_morph(NS(Yl1,1).*Ybd/5,'e',1))>0.6; % no WM hyperintensities in GM!
    %
    Ygm = Ygm | (Yss & ~Yvt & ~Ycx & ~Ybv2 & ~Ycwm & ~(Yccm | Ybcm));
    Ygm = Ygm & ~Ywmtpm & ~Ybvv; % no WMH area
    Ywm = (Ywm & ~Yss & ~Ybv2  & ~Ynw) | Ycwm | Ybwm; %& ~NS(Yl1,LAB.BG)
    Ywmtpm(smooth3(Ywmtpm & Ym<11/12)<0.5)=0;
    Ywm = Ywm & ~Ywmtpm & ~Ybvv & ~Yss; % no WM area
    Ycm = Ycm | ( (Ycx | Yccm | Ybcm) & Yg<0.2 & Ym>0 & Ydiv>-0.05 & Ym<0.3 & Yb ) | Ybvv;
    if ~debug, clear Ycwm Yccm Ycd; end
    % mapping of the brainstem to the WM (well there were some small GM
    % structures, but the should not effect the local segmentation to much.
    Ybs = cat_vol_morph(NS(Yl1,LAB.BS) & Ym<1.2 & Ym>0.9 & Yp0>2.5,'c',2*vxv) & Ym<1.2 & Ym>0.9 & Yp0>1.5;
    Ygm = (Ygm & ~Ybs & ~Ybv2 & ~Ywm) | Yss;
    Ywm = Ywm | (Ybs & Ym<1.1 & Ym>0.9 & Yp0>1.5) ; 
    if ~debug, clear Ycx; end
  end
  
  
  % back to original resolution for full bias field estimation
  [Ycm,Ygm,Ywm]      = cat_vol_resize({Ycm,Ygm,Ywm},'dereduceBrain',BB); % ,Ygw
  [Yvt,Yb,Yss,Ybb,Ysc,Ybs,Ybv2] = cat_vol_resize({Yvt,Yb,Yss,Ybb,Ysc,Ybs,Ybv2},'dereduceBrain',BB);
  [Ym,Yp0,Yl1] = cat_vol_resize({Ym,Yp0,Yl1},'dereduceBrain',BB);
  [Ybd,Ybvv] = cat_vol_resize({Ybd,Ybvv},'dereduceBrain',BB);
  [Yg,Ydiv] = cat_vol_resize({Yg,Ydiv},'dereduceBrain',BB);
  [Ywd] = cat_vol_resize(Ywd,'dereduceBrain',BB); % Ysk
  Ycls=Yclso; Ysrc=Ysrco; 
  clear Yclso Ysrco Ybv;
  
  if debug>1
    try %#ok<TRYNC> % windows requires this... i don't know why ... maybe the file size
      tpmci=tpmci+1; tmpmat = fullfile(pth,reportfolder,sprintf('%s_%s%02d%s.mat',nam,'LAS',tpmci,'prepeaks'));
      save(tmpmat);
    end
  end
            
  % correction for negative values 
  srcmin = min(Ysrc(:)); Ysrc = Ysrc - srcmin; T3th = T3th - srcmin; Tthc = Tth; Tthc.T3th = Tth.T3th - srcmin;
  
%% --------------------------------------------------------------------- 
%  Now, we can estimate the local peaks 
%  ---------------------------------------------------------------------
  % Estimation of the local WM threshold with "corrected" GM voxels to
  % avoid overfitting (see BWP cerebellum). 
  % CSF is problematic in high contrast or skull-stripped image should 
  % not be used here, or in GM peak estimation
  mres = 1.1; 
  stime = cat_io_cmd('  Estimate local tissue thresholds','g5','',verb,stime); dispc=dispc+1;
  Ysrcm = cat_vol_median3(Ysrc.*Ywm,Ywm,Ywm); 
  rf    = [10^5 10^4];
  T3th3 = max(1,min(10^6,rf(2) / (round(T3th(3)*rf(1))/rf(1))));
  Ysrcm = round(Ysrcm*T3th3)/T3th3;
  Ygw2 = Ycls{1}>128 & Ym>2/3-0.04 & Ym<2/3+0.04 & Ygm .*Ydiv>0.01;
  Ygw2 = Ygw2 | (Ycls{1}>128 & Yg<0.05 & abs(Ydiv)<0.05 & ~Ywm & Ym<3/4); % large stable GM areas - like the BWP cerebellum
  Ygw3 = Ycls{3}>128 & Yg<0.05 & ~Ywm & ~Ygm & Ywd<3; 
  Ygw3(smooth3(Ygw3)<0.5)=0;
  [Yi,resT2] = cat_vol_resize(Ysrcm,'reduceV',vx_vol,mres,32,'max'); % maximum reduction for the WM
  %
  if cat_stat_nanmean(Ym(Ygw3))>0.1, % not in images with to low CSF intensity (error in skull-stripped)
    Ygi = cat_vol_resize(Ysrc.*Ygw2*T3th(3)/mean(Ysrc(Ygw2(:))) + Ysrc.*Ygw3*T3th(3)/mean(Ysrc(Ygw3(:))),'reduceV',vx_vol,mres,32,'meanm');  % mean for other tissues
  else
    Ygi = cat_vol_resize(Ysrc.*Ygw2*T3th(3)/mean(Ysrc(Ygw2(:))),'reduceV',vx_vol,mres,32,'meanm'); % mean for other tissues
  end
  for xi=1:2*LASi, Ygi = cat_vol_localstat(Ygi,Ygi>0,2,1); end; Ygi(smooth3(Ygi>0)<0.3)=0;
  Yi = cat_vol_localstat(Yi,Yi>0,1,3); % one maximum for stabilization of small WM structures
  Yi(Yi==0 & Ygi>0)=Ygi(Yi==0 & Ygi>0);
  for xi=1:2*LASi, Yi = cat_vol_localstat(Yi,Yi>0,2,1); end % no maximum here!
 %
  if cat_stat_nanmean(Ym(Ygw3))>0.1 && cat_stat_nanmean(Ysrc(Ycls{6}(:)>128))>T3th(1)
    %%
     Ygw4 = Ycls{6}>240 & ~Ygw3 & ~Ygw2 & Yg<0.5 & abs(Ydiv)<0.1 & ...
        Ysrc>min(res.mn(res.lkp==6))*0.5 & Ysrc<max(res.mn(res.lkp==6))*1.5; 
     [Ygi,resTB] = cat_vol_resize(Ysrc.*Ygw4*T3th(3)/mean(Ysrc(Ygw4(:))),'reduceV',vx_vol,mres*4,16,'meanm');
     Ygi = cat_vol_approx(Ygi,'nh',resTB.vx_volr,2); Ygi = cat_vol_smooth3X(Ygi,2*LASfs);
     Ygi = Ygw4 .* max(eps,cat_vol_resize(Ygi,'dereduceV',resTB)); 
  end
  Yi(Yi==0) = Ygi(Yi==0); 
  if debug==0; clear Ygw2; end
  Yi = cat_vol_approx(Yi,'nh',resT2.vx_volr,2); Yi = cat_vol_smooth3X(Yi,2*LASfs); 
  Ylab{2} = max(eps,cat_vol_resize(Yi,'dereduceV',resT2)); 
 % Ylab{2} = Ylab{2} .* mean( [median(Ysrc(Ysw(:))./Ylab{2}(Ysw(:))),1] ); 
  if debug==0; clear Ysw; end

  %% update GM tissue map
  %Ybb = cat_vol_morph((Yb & Ysrc./Ylab{2}<(T3th(1) + 0.25*diff(T3th(1:2))) & Ydiv<0.05) | Yp0>1.5,'lo',1);
  Ygm(Ysrc./Ylab{2}>(T3th(2) + 0.90*diff(T3th(2:3)))/T3th(3))=0; % correct GM mean(T3th([2:3,3]))/T3th(3) 
  Ygm(Ysrc./Ylab{2}<(T3th(2) + 0.75*diff(T3th(2:3)))/T3th(3) & ...
      Ysrc./Ylab{2}<(T3th(2) - 0.75*diff(T3th(2:3)))/T3th(3) & ...
      Ydiv<0.3 & Ydiv>-0.3 & Ybb & ~Ywm & ~Yvt & ~Ybv2 & Ycls{1}>48)=1;
  Ywmd2 = cat_vbdist(single(Ywm),Yb);
  Ygx = Ywmd2-Ym+Ydiv>0.5 & Ym+0.5-Ydiv-Yg-Ywmd2/10>1/3 & ~Ybv2 & ... low intensity tissue
    ~(Ym-min(0.2,Yg+Ywmd2/10-Ydiv)<1/4) & Yg<Ylab{2}/T3th(3)*0.3 & Ysrc<Ylab{2}*0.9; % no real csf
  Ygx(smooth3(Ygx)<0.5)=0; 
  Ygm = ~Ywm & (Ygm | Ygx); % correct gm (intensity based)
  Ygx = (single(Ycls{1})/255 - abs(Ydiv) + min(0,Ydiv) - Yg)>0.5 & ~Ywm;
  Ygx(smooth3(Ygx)<0.5)=0; 
  Ygm = Ygm | Ygx; % correct gm (spm based)
  %%
  Ycm = ~Ygm & ~Ywm & ~Ybv2 & Yg<0.6 & (Ycm | (Yb & (Ysrc./Ylab{2})<((T3th(1)*0.5 + 0.5*T3th(2))/T3th(3)))); 
  %
  Ycp = (Ycls{2}<128 & Ydiv>0 & Yp0<2.1 & Ysrc./Ylab{2}<mean(T3th(1)/T3th(2))) | Ycm | ...      % typcial CSF
        (Ycls{5}>32 & Ycls{2}<32 & Ysrc./Ylab{2}<T3th(2)/T3th(3) & Ydiv>0) | ...                % venes
        ((Ym-Ydiv<0.4) & Ycls{3}>4 & Ycls{3}>16 & Ysrc./Ylab{2}<mean(T3th(2)/T3th(3))) | ...    % sulcal CSF
        (single(Ycls{6})+single(Ycls{5})+single(Ycls{4}))>192 | ...                             % save non-csf 
        Ysrc./Ylab{2}<T3th(1)/T3th(3);                                                          % but do not trust the brain mask!
  Ycp(smooth3(Ycp)>0.4)=1;                                                                      % remove some meninges
  Ycd = cat_vbdist(single(Ycp),~Ycp,vx_vol);  
  %%
  Ygm = Ygm & ~Ycm & ~Ywm & Ywd<5; %  & ~Ybvv  & ~Ysk
  
  Ygm = Ygm | (NS(Yl1,1) & Ybd<20 & (Ycd-Ydiv)<2 & Ycls{1}>0 & ~Ycm & Ybb & Ym>0.6 & Yg<max(0.5,1-Ybd/30)); 
  Ygm = Ygm & (Yg<0.1 | Ysrc./Ylab{2}<(T3th(2)*1/3+2/3*T3th(3))/T3th(3)); % outer high intensity GM
  %
  if LABl1
    Ygm = (Ygm | Yss) & ~Ycm & cat_vol_morph(~Ybs | ~Yvt,'e');
  end
  Ybb = cat_vol_morph(smooth3(Ygm | Ywm | Yp0>1.5 | (Ym>1.2/3 & Ym<3.1/3 & Yb))>0.6,'lo',min(1,vxv)); % clear Yp0 Yvt
  Ygm(~Ybb)=0; Ygm(smooth3(Ygm)<0.3)=0;
  Ygm(smooth3(Ygm)>0.4 & Ysrc./Ylab{2}>mean(T3th(1)/T3th(2)) & Ysrc./Ylab{2}<(T3th(2)*0.2+0.8*T3th(3)))=1;
  if debug==0; clear Ybb Ybd Yvt Ybvv Ycp Ycd Yl1 Yss; end %Ydiv Yg

  
  %% GM
%  Yi = (Ysrc./Ylab{2} .* Ygm , Ysrc./Ylab{2}.*Ywm.*Ybs.*(T3th(2)+0.8*diff(T3th(2:3)))/T3th(3));
%  Yi(Ybv2) = Ysrc(Ybv2)./Ylab{2}(Ybv2) .* T3th(2)/mean(T3th(1:2)); 
  Yi = Ysrc./Ylab{2} .* Ygm; 
  Yi = round(Yi*rf(2))/rf(2);
 % Yi(Ybv2) = Ysrc(Ybv2)./Ylab{2}(Ybv2) .* T3th(2)/mean(T3th(1:2)); % ????
  Yi(Ybs)  = Ysrc(Ybs)./Ylab{2}(Ybs)   .* T3th(2)/T3th(3); 
  Yi = cat_vol_median3(Yi,Yi>0.5,Yi>0.5); 
  Ycmx = smooth3(Ycm & Ysrc<(T3th(1)*0.8+T3th(2)*0.2))>0.9; Tcmx = mean(Ysrc(Ycmx(:))./Ylab{2}(Ycmx(:)))*T3th(3);
  Yi(Ycmx) = Ysrc(Ycmx)./Ylab{2}(Ycmx)  .* T3th(2)/Tcmx; 
  %Yii =  Ysrc./Ylab{2} .* Ycm * T3th(2) / cat_stat_nanmedian(Ysrc(Ycm(:))); 
  [Yi,Yii,Ybgx,resT2] = cat_vol_resize({Yi,Ylab{2}/T3th(3),Ycls{6}>240},'reduceV',vx_vol,1,32,'meanm');
  for xi=1:2*LASi, Yi = cat_vol_localstat(Yi,Yi>0,3,1); end
  Yi = cat_vol_approx(Yi,'nh',resT2.vx_volr,2); 
  Yi = min(Yi,Yii*(T3th(2) + 0.90*diff(T3th(2:3)))/T3th(3)); 
  Yi(Ybgx) = Yii(Ybgx)*cat_stat_nanmean(Yi(~Ybgx(:))); 
  %%
  Yi = cat_vol_smooth3X(Yi,LASfs); 
  Ylab{1} = cat_vol_resize(Yi,'dereduceV',resT2).*Ylab{2};   
  %Ylab{1}(Ygm) = Ysrc(Ygm); Ylab{1} = cat_vol_smooth3X(Ylab{1},LASfs); % can lead to overfitting
  Ycm = (single(Ycls{3})/255 - Yg*4 + abs(Ydiv)*2)>0.5 &  Ysrc<(Ylab{1}*mean(T3th([1,1:2]))/T3th(2)); %Ycm & Ysrc<(Ylab{1}*mean(T3th(1:2))/T3th(2)) & Yg<0.1;
  Ycm(smooth3(Ycm)<0.5)=0;
  Ycm(Yb & cat_vol_morph(Ysrc<mean(T3th(1:2)),'o'))=1;
  
  %% CSF & BG 
  if exist('Ygw4','var')
    Ynb = Ygw4 & cat_vol_morph(smooth3(Ycls{6})>128 | (~Yb & Yg<=0.001),'e',4*vxv); 
  else
    Ynb = cat_vol_morph(smooth3(Ycls{6})>128 | (~Yb & Yg<=0.001),'e',4*vxv); 
  end  
  [Yc,resT2] = cat_vol_resize(round(Ysrc./Ylab{2} .* (smooth3(Ycm | Ysc)>0.5) * rf(2))/rf(2),'reduceV',vx_vol,8,16,'min');% only pure CSF !!!
  [Yx,Yxm]   = cat_vol_resize({round(Ysrc./Ylab{2} .* Ynb * rf(2))/rf(2),single(Ynb)},'reduceV',vx_vol,8,16,'mean');% only pure CSF !!!
  Yx(Yxm<1) = 0; 
  %%
  Yx(Yc>0)=0; Yc(Yx>0)=0;
  if cat_stat_nanmean(Yx(Yx~=0)) < cat_stat_nanmean(Yc(Yc~=0))
    meanYx = min(median(Yc(Yc(:)>0)),median(Yx(Yx(:)>0))); 
    meanYc = max(median(Yc(Yc(:)>0)),median(Yx(Yx(:)>0))); 
  else
    meanYc = min(median(Yc(Yc(:)>0)),median(Yx(Yx(:)>0))); 
    meanYx = max(median(Yc(Yc(:)>0)),median(Yx(Yx(:)>0))); 
  end  
  stdYbc = mean([std(Yc(Yc(:)>0)),std(Yx(Yx(:)>0))]);
  %Yx = min(max(meanYx/2,Yx),min(meanYx*4,meanYc/2));
  %Yc = min(max(meanYc/2,Yx),meanYc/2);
  Yxa = cat_vol_approx(Yx ,'nh',resT2.vx_volr,16); %+(Yb>0).*stdYbc  + Yc.*meanYb/max(eps,meanYc)
  Yca = cat_vol_approx(Yc + min(max( meanYx + stdYbc , meanYc - stdYbc ),...
    Yx.*meanYc/max(eps,meanYx)),'nh',resT2.vx_volr,16); % + Yb.*meanYc/max(eps,meanYb)
  Yca = Yca*0.7 + 0.3*max(mean(Yca(:)),T3th(1)/T3th(3));
  %%
  Yxa = cat_vol_smooth3X(Yxa,LASfs*2); 
  Yca = cat_vol_smooth3X(Yca,LASfs*2); 
  Ylab{3} = cat_vol_smooth3X(cat_vol_resize(Yca,'dereduceV',resT2).*Ylab{2},LASfs*2);  
  Ylab{6} = cat_vol_smooth3X(cat_vol_resize(Yxa,'dereduceV',resT2).*Ylab{2},LASfs*2);
  if ~debug, clear Yxa Yca Yx Yc Y Ydiv; end
  
  %% local intensity modification of the original image
  % --------------------------------------------------------------------
  Yml = zeros(size(Ysrc)); 
  Yml = Yml + ( (Ysrc>=Ylab{2}                ) .* (3 + (Ysrc-Ylab{2}) ./ max(eps,Ylab{2}-Ylab{3})) );
  Yml = Yml + ( (Ysrc>=Ylab{1} & Ysrc<Ylab{2} ) .* (2 + (Ysrc-Ylab{1}) ./ max(eps,Ylab{2}-Ylab{1})) );
  Yml = Yml + ( (Ysrc>=Ylab{3} & Ysrc<Ylab{1} ) .* (1 + (Ysrc-Ylab{3}) ./ max(eps,Ylab{1}-Ylab{3})) );
  Yml = Yml + ( (Ysrc< Ylab{3}                ) .* (    (Ysrc-Ylab{6}) ./ max(eps,Ylab{3}-Ylab{6})) );
  Yml(isnan(Yml) | Yml<0)=0; Yml(Yml>10)=10;
  
  
  %% global
  if Tth.T3th(Tth.T3thx==1)==Tth.T3thx(Tth.T3thx==1) %invers
    Ymg = max(eps,(Ysrc + srcmin)./Ylab{2}); 
  else
    Ymg = (Ysrc + srcmin)./max(eps,(Ylab{2} + srcmin)); 
    Ymg = Ymg * Tth.T3th(Tth.T3thx==3)/(Tth.T3thx(Tth.T3thx==3)/3);
  end
  Ymg = cat_main_gintnorm(Ymg,Tth); 
  
  %%
  if debug>1
    try %#ok<TRYNC> % windows requires this... i don't know why
      tpmci=tpmci+1; tmpmat = fullfile(pth,reportfolder,sprintf('%s_%s%02d%s.mat',nam,'LAS',tpmci,'postpeaks'));
      save(tmpmat);
    end
  end
  
  %% fill up CSF in the case of a skull stripped image 
  if max(res.mn(res.lkp==5 & res.mg'>0.1)) < mean(res.mn(res.lkp==3 & res.mg'>0.3))
    YM   = cat_vol_morph(Yb,'d'); 
    Ymls = smooth3(max(Yml,YM*0.5));
    Yml(YM & Yml<0.5)=Ymls(YM & Yml<0.5); 
    clear Ymls YM
  end
  
  
  %% class correction and second logical class map Ycls2
  Ynwm = Ywm & ~Ygm & Yml/3>0.95 & Yml/3<1.3;
  Ynwm = Ynwm | (smooth3(Ywm)>0.6 & Yml/3>5/6); Ynwm(smooth3(Ynwm)<0.5)=0;
  Yngm = Ygm & ~Ywm & Yml/3<0.95; Yngm(smooth3(Yngm)<0.5)=0;
  Yncm = ~Ygm & ~Ywm & ((Yml/3)>1/6 | Ycls{3}>128) & (Yml/3)<0.5 & Yb;
  Ycls{2} = cat_vol_ctype(single(Ycls{2}) + (Ynwm & ~Yngm & Yp0>=1.5)*256 - (Yngm & ~Ynwm & Yp0>=2)*256,'uint8');
  Ycls{1} = cat_vol_ctype(single(Ycls{1}) - (Ynwm & ~Yngm & Yp0>=1.5)*256 + (Yngm & ~Ynwm & Yp0>=2)*256,'uint8');
  %Ycls{3} = cat_vol_ctype(single(Ycls{3}) - ((Ynwm | Yngm) & Yp0>=2)*256,'uint8');
  %Ycls{3} = cat_vol_ctype(single(Ycls{3}) + (Yb & Yml<1.1 & ~Ynwm & ~Yngm)*256,'uint8');
  Ycls{1} = cat_vol_ctype(single(Ycls{1}) - (Yb & Yml<1.1 & ~Ynwm & ~Yngm)*256,'uint8');
  Ycls{2} = cat_vol_ctype(single(Ycls{2}) - (Yb & Yml<1.1 & ~Ynwm & ~Yngm)*256,'uint8');
  Ycls2 = {Yngm,Ynwm,Yncm};
  clear Yngm Ynwm Yncm;
  
  %%
  Yml = Yml/3;
  cat_io_cmd('','','',verb,stime);
 % if debug
 %   cat_io_cmd(' ','','',verb,stime); 
 % else
 %   cat_io_cmd(' ','','',verb,stime);   
 %   cat_io_cmd('cleanup',dispc,'',verb);
 % end

end

