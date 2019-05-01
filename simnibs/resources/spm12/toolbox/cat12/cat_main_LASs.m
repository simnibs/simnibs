function [Yml,Ymg,Ycls] = cat_main_LASs(Ysrc,Ycls,Ym,Yb,Yy,Tth,res,vx_vol,extopts) 
% This is an exclusive subfunction of cat_main.
% ______________________________________________________________________
%
% Local Adaptive Segmentation (LAS):
%
% This version of the local adaptive intensity correction includes a  
% bias correction that based on a maximum filter for the WM and a mean
% filter of the GM to stabilize the correction in region with less WM.
%
% The extension based mostly on the assumption that the tissue next to 
% the CSF (and high divergence sulci) has to be WM (maximum, high 
% divergence) or GM. For each tissue a refined logical map is generated 
% and used to estimate the local intensity threshold.
%
% It is important to avoid high intensity blood vessels in the process, 
% because they will push down local WM and GM intensity. 
%
% There are further regionwise correction, e.g. , to avoid overfitting in 
% cerebellum, or adapt for age specific changes, e.g. enlarged ventricle.
%
% Based on this values a intensity transformation is used. Compared to 
% the global correciton this has to be done for each voxel. To save time
% only a rough linear transformation is used.
% ______________________________________________________________________
%
%   [Yml,Ymg,Yclsu] = ...
%     cat_main_LAS(Ysrc,Ycls,Ym,Yb0,Yy,T3th,res,vx_vol,extopts,Tth)
%
%   Yml     .. local intensity correct image
%   Ymg     .. global intensity correct image
%   Yclsu   .. corrected SPM tissue class map
%
%   Ysrc    .. (bias corrected) T1 image
%   Ycls     .. SPM tissue class map
%   Ym      .. intensity corrected T1 image (BG=0,CSF=1/3,GM=2/3,WM=1)
%   Yb      .. brain mask
%   Yy      .. deformation map for the cat atlas map
%   Tth     .. structure with tissue thresholds of CSF, GM, and WM in Ysrc
%   res     .. SPM segmentation structure
%   vx_vol  .. voxel dimensions
%   extopts .. cat options
% ______________________________________________________________________
% 
% internal maps:
%
%   Yg   .. gradient map   - edges between tissues
%   Ydiv .. divergence map - sulci, gyris pattern, and blood vessels
%   Yp0  .. label map      - tissue classes (BG=0,CSF=1,GM=2,WM=3) 
%
%   Ycm  .. CSF
%   Ygm  .. GM
%   Ywm  .. WM 
%
%   Yvt  .. WM next to the ventricle map 
%   Ygmt .. cortical thickness map
%   Ypp  .. cortical percentage position map
% ______________________________________________________________________
%
%   Robert Dahnke (robert.dahnke@uni-jena.de)
%   Structural Brain Mapping Group (http://dbm.neuro.uni-jena.de/)
%   Department of Neurology
%   University Jena
% ______________________________________________________________________
% $Id: cat_main_LASs.m 1118 2017-03-17 15:57:00Z gaser $

  % set this variable to 1 for simpler debuging without reduceBrain
  % function (that normally save half of processing time)
  verb    = extopts.verb-1;
  dsize   = size(Ysrc);
  NS      = @(Ys,s) Ys==s | Ys==s+1;          % function to ignore brain hemisphere coding
  LASstr  = max(eps,min(1,extopts.LASstr));   % LAS strenght (for GM/WM threshold)3
  LAB     = extopts.LAB;                      % atlas labels
  mvx     = mean(vx_vol);                     % mean voxel volume to correct for morphological operations  
  T3th    = Tth.T3th(3:6);                    % CSF, GM, and WM peak
  Tth.T3th(1) = min(0,Tth.T3th(1));           % correction of the background value
  
  % set debug = 1 and do not clear temporary variables if there is a breakpoint in this file 
  dbs   = dbstatus; debug = 0; 
  for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,'cat_main_LASs'); debug = 1; break; end; end
  
  
%% ---------------------------------------------------------------------
%  First, we have to optimize the segments using further information that 
%  SPM do not use, such as the gradient, divergence and distance maps. 
%  The gradient map (average of the first derivate of the T1 map) is an 
%  edge map and independent of the image intensity. It helps to avoid PVE 
%  regions and meninges. 
%  The divergence (second derivate of the T1 map) help to identfiy sulcal
%  and gyral pattern and therefore to find WM and CSF regions for further 
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
  stime = cat_io_cmd('  Prepare maps','g5','',verb); 

  % map atlas to RAW space
  % -------------------------------------------------------------------
  % avoid read error in parallel processing
  for i=1:5
    try
      Vl1A = spm_vol(extopts.cat12atlas{1}); break
    catch 
      pause(1)
    end
  end
  Yl1  = cat_vol_ctype(round(spm_sample_vol(Vl1A,double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),0)));
  Yl1  = reshape(Yl1,dsize);
  if ~debug, clear Yy; end

  
  
  % lower resolution to save time and space
  rres = 1.2; 
  [Ym,resTb] = cat_vol_resize(Ym        ,'reduceV',vx_vol,rres,32,'meanm'); 
  if any(resTb.vx_vol ~= resTb.vx_volr), Ysrco = Ysrc+0; end % save orignal image for later correction
  Ysrc       = cat_vol_resize(Ysrc      ,'reduceV',vx_vol,rres,32,'meanm'); 
  Yb         = cat_vol_resize(single(Yb),'reduceV',vx_vol,rres,32,'meanm')>0.5; 
  Yl1        = cat_vol_resize(Yl1       ,'reduceV',vx_vol,rres,32,'nearest'); 
  for i=1:6, Ycls{i} = cat_vol_ctype(cat_vol_resize(single(Ycls{i}),'reduceV',vx_vol,rres,32,'meanm')); end
  if debug, Ymo = Ym; Ybo = Yb; end 
  vx_vol = resTb.vx_volr;
  
  % help maps to detect edges (Yg) and sulci/gyris (Ydiv)
  Yg    = cat_vol_grad(Ym,vx_vol);                                                  % mean gradient map  ...  ./max(0.1,Ym)
  Ydiv  = cat_vol_div(Ym,vx_vol);                                                   % divergence map 
  Yp0   = single(Ycls{1})/255*2 + single(Ycls{2})/255*3 + single(Ycls{3})/255;      % tissue label map
  noise = std(Ym(Ycls{2}>192)); 
  
  noisef = max(0,min(0.3,noise./mvx)); 
  Ym  = cat_vol_smooth3X(Ym,noisef); 
  Yp0 = cat_vol_smooth3X(Yp0,noisef); 
  for i=1:6, Ycls{i} = cat_vol_ctype(cat_vol_smooth3X(single(Ycls{i}),noisef)); end
 
  Yb = smooth3(Yb | (cat_vol_morph(Yb,'d',2/mvx) & Ym<0.8 & Yg<0.3 & Ym>0 & Yp0>0.2))>0.5;  % increase brain mask, for missing GM 

  % correction for negative (and positive) values
  srcmin = min(Ysrc(:)); Ysrc = Ysrc - srcmin; Tth.T3th = Tth.T3th - srcmin; T3th = T3th - srcmin; 
  
  %% initial bias correction 
  %  ----------------------------------------------------------------------
  %  required especial in case of strong bias (7 Tesla) for the redefinition 
  %  of tissue segments for the LAS correction
  %  ----------------------------------------------------------------------
  if 1;; % this if is just to run the full block and keep the subblocks in debuging mode!
    stime = cat_io_cmd('  Initial bias correction','g5','',verb,stime); 
    if debug, Ym = Ymo; Yb = Ybo; end 
  
    %  background and brain distance map
    %  background may include artifacts, but no head tissue (check for 7T and motion)
    
    % ... bias correction for background detection 
    [Ymb,resT3] = cat_vol_resize(Ym .* (Ym./Ydiv<0),'reduceV',vx_vol,4,32,'max');  
    Ymb = cat_vol_approx(Ymb,'nh',resT3.vx_volr,4); Ymb = cat_vol_smooth3X(Ymb,8);   
    Ymb = Ym ./ cat_vol_resize(Ymb ,'dereduceV',resT3);  
     
    Ybg = cat_vol_smooth3X( Ycls{6}>128 & Ysrc>mean(Tth.T3th(1:2)) & ...                          % SPM background & lower boundary (real background
          (Ymb-Ydiv./Ym)<max(0.1,0.05+noise) & ...
          (Ysrc-(Ydiv./Ymb*T3th(3))<mean(T3th(1:2))) & Ydiv./Ymb<=0 & Yg./Ydiv>10000 & ...   (Ym-Ydiv./Ym)<max(0.1,0.05+noise) &   % upper boundary  
          Yg<noise & abs(Ydiv./Ymb)<0.1 & Yg./Ymb>0.05,max(0.5,min(3,1/mvx)))>0.5;              % gradients and divergence to avoid objects 
    Ybg = cat_vol_morph(Ybg,'c',3/mvx);    
    [Ybdr,Ybgr,resT3] = cat_vol_resize({single(Yb),single(Ybg)},'reduceV',vx_vol,2,32,'meanm');    % brain distance Ybd
    Ybdr  = cat_vbdist(single(Ybdr>0.1),Ybgr<=0.05,resT3.vx_volr); 
    Ybgdr = cat_vbdist(single(Ybgr>0.5),Ybdr>=0.00,resT3.vx_volr); 
    Ybd   = cat_vol_resize(Ybdr ,'dereduceV',resT3); clear Ybdr; 
    Ybgd  = cat_vol_resize(Ybgdr,'dereduceV',resT3); clear Ybgdr; 
    Ybv   = (Ym - 2*Ydiv)>1.1 & Ym>1.1 & Ycls{2}<64 & Yp0>0;                                       % blood vessels

    %% brain tissues WM, GM, CSF (Ywi, Ygi, Yci)
    stime2 = cat_io_cmd('    Brain Tissues','g5','',debug,stime);  
    Ywmh = smooth3(Yp0/3-Ym)>0.2 & Ym<2.5/3 & Ym>1.5/3 & NS(Yl1,LAB.CT)>0.5;
    Ywh = cat_vol_smooth3X(max(0,Ym-1) .* (Yp0>2.75 & Ym>1.05 & ~Ybv & ~Ywmh),4)*4;                            % biased WM
    Ym  = Ym - Ywh; clear Ywh; % hard but good working correction for uncorrected high intensity bias (WM with intensity over 1)
    Yss = (NS(Yl1,LAB.BG) | NS(Yl1,LAB.TH)) & Ym<0.99 & Ym>0.6;                                                % subcortical structures
    % WM 
    Ywi = Ysrc .* (Yp0>2.5 | ((Ym - 2*Ydiv./Ym)>1 & (Yp0 - 2*Ydiv./Yp0)>2.1 & Yp0>1.9) ) .* (~Yss & ~Ywmh & ~Ybv & Yb);  
    Ywi(smooth3(Ywi>0)<0.5)=0;
    Ywi = cat_vol_noPVE(Ywi,vx_vol,1,1);       
    % GM 
    Ygi = (Yp0 - 2*abs(Ydiv./Ym))>max(1,1.2-noise*2) & (Yp0 + 2*abs(Ydiv./Ym))<2.5 & Yb & abs(Ydiv./Ym)<0.05 & Ydiv./Yg>-0.5;                  
    Ygi = (Ygi | Yss) & ~Ywi & ~Ywmh & ~cat_vol_morph( NS(Yl1,LAB.VT) | NS(Yl1,LAB.BS) ,'d',2/mvx);            % no ventricle/brainstem
    Ygi = Ygi .* Ysrc * T3th(3)/T3th(2);                                                                       % correct for WM intensity
    Ygi(smooth3(Ygi>0)<0.5)=0;
    Ygi = cat_vol_noPVE(Ygi,vx_vol,1);                                                                                  
    % CSF
    Yci = (Ycls{3}>16 & NS(Yl1,LAB.VT) & Ym>0.5/3)  | (Ycls{3}>224 & Ym>0.8/3 );                                          
    Yci = Yci & (Ym - Ydiv./Ym)<1.5/3 & ~Ygi & ~Ywi & (Yg./Ym)<0.3 & Yg<0.1 & Ydiv./Ym>-0.02 & ~Ywmh & Ysrc>mean(Tth.T3th(1:2));  % only save CSF!
    Yci(smooth3(Yci>0)<0.5) = 0; 
    Yci = Yci .* Ysrc * T3th(3)/median(Ysrc(Yci(:)>0));                                                        % correct for WM intensity
    Yci = cat_vol_noPVE(Yci,vx_vol,1);                                                                                  
    Ygi = Ygi + Yci;                                                                                           % save memory
    if ~debug, clear Yss Yci; end

    %% head tissues 
    %  The normal head tissue below WM intensity is typically similar to the GM intensity, 
    %  even this is not true in all cases it allow a good correction even with strong inhomogeneities. 
    %  It is not clear if the high intensities head tissue can be used, becuase of the uncleare 
    %  effect of the PVE and blood vessels. But it is important do avoid high intensity tissue
    %  in the low intensity tissue mask.
    stime2 = cat_io_cmd('    Head Tissues','g5','',debug,stime2);
    Yhl = Yg./Ymb>0.2 & Yg<0.2 & Ymb<0.3 & Ydiv./Ymb>=0 & (Ycls{4}+Ycls{5}+Ycls{6})>128 & cat_vol_smooth3X(Ybg,4/mvx)<0.1;
    Yhl = smooth3(Yhl)>0.5;
    Yhl = Yhl * T3th(3)/max(eps,median([0.1;Yhl(Yhl(:)>0 & Ybd(:)<10)])); 
    Yhg = Yp0==0 & Ymb>max(0.1,2/3-Ybd/50) & Yg./Ymb<0.3 & Ymb>max(0.1,min(0.2,mean([0.2,median(Ymb(Ybg(:)))]))) & ...
      Yg<0.2 & ~Yhl & (Ymb + 2*abs(Ydiv))<max(0.8,0.6+Ybd/20) & Ybd>5 & ...
      cat_vol_smooth3X(Ybg,4/mvx)<0.1 & ~((Ym-Ydiv)>0.8 | Ym>max(0.8,1-Ybd/50));           % head GM like (mussels)
    Yhg = Yhg | (Ym>0.6 & ~Ybv & Yp0<1.5 & (Ybd ./ (Ybd + Ybgd))<0.2 & Ybgd<20); 
    Yhg(smooth3(Yhg)<0.5) = 0; 
    % remove small dots
    Yhi = Yp0==0 & ~Yhg & ~Yhg & ~Ybg & Ybd>2 & ((Ymb-Ydiv)>1.2 | Ymb>max(0.5,1.5-Ybd/50)) & Ydiv./Ymb<0.1;       % higher head tissue
    Yhi(smooth3(Yhi)<0.5) = 0;                                                                                 % remove small dots
    Yhi(cat_vol_smooth3X(Yhi,1/mvx)<0.5)=0;
    Yhi = Yhi .* Ysrc; Yhi = cat_vol_noPVE(Yhi,vx_vol,1);                                                               % correct for PVE
    Yhi = Yhi * T3th(3)/max(eps,median(Yhi(Yhi(:)>0 & Ybd(:)<10))); 
    Yhi = cat_vol_noPVE(Yhi,vx_vol,1,2);
    % correct for WM intensity
    %Yhg = Yhg | ((smooth3(Yhg | Ycls{6}>128)<0.2) & Ygi>0 & Ycls{2}<8 & (Ym - 2*Ydiv)>0.5 & Ym<0.8 & Yp0<2);    % further low intensity head tissue 
    %Yhg = cat_vol_morph(Yhg & ~Yhi & Ybd>5 & Yg<0.15 & abs(Ydiv)<0.15,'o'); 
    Yhg = Yhg .* Ysrc * T3th(3)/T3th(2); %median(Yhg(Yhg(:)>0 & Ybd(:)<8));                                   % expect GM like intensity
    Yhg = cat_vol_noPVE(Yhg,vx_vol,1,2);
    Yhl = cat_vol_noPVE(Yhl,vx_vol,1,2);
    Ygi = Ygi + Yhg + 1*Yhi + Yhl;    % ... adding of Yhi failed in some cases                                                                               % save memory
    if ~debug, clear Yhi Yhg Ybv; end  
    Ygi = cat_vol_noPVE(Ygi,vx_vol,1,2);

    %% background to avoid conflict by bad approximation 
    stime2 = cat_io_cmd('    Background','g5','',debug,stime2);
    if median(Ysrc(Ybg(:)))/T3th(3)>0.02 && median(Ysrc(Ybg(:)))/T3th(3)<T3th(1)
      Ybgi = (Ybg & Ym>0.01) .* (Ysrc * T3th(3)/max(eps,median(Ysrc(Ybg(:)))));
      Ybgi(Ysrc==0) = T3th(3);  
      [Ybgi,Yd,resT2] = cat_vol_resize({Ybgi,single(Ygi>0)},'reduceV',vx_vol,min(vx_vol*3,6),32,'meanm'); 
      Ybgi = cat_vol_noPVE(Ybgi,vx_vol,2);
      Yd   = cat_vbdist(Yd,Yd<1,resT2.vx_volr); 
      Ybgi = cat_vol_approx(Ybgi,'nh',resT2.vx_volr,4); Ybgi = cat_vol_smooth3X(Ybgi,4); 
      Ybgi = max(eps,cat_vol_resize(Ybgi,'dereduceV',resT2));
      Yd   = max(eps,cat_vol_resize(Yd,'dereduceV',resT2));
      Ygi = Ygi + Ybgi .* (Yd>20 & Ygi==0); clear Ybgi Yd;
    else
      Ygi = Ygi + Ybg * T3th(3); 
    end

    %% estimate correction on lower resolution
    stime2 = cat_io_cmd('    Bias Field','g5','',debug,stime2);
    [Ywi,resT2] = cat_vol_resize(Ywi,'reduceV',vx_vol,min(vx_vol*2,2),32,'max');          % maximum reduction for the WM
    [Ygi,Ybdr]   = cat_vol_resize({Ygi,Yp0},'reduceV',vx_vol,min(vx_vol*2,2),32,'meanm');  % mean for the rest
    Ywi = cat_vol_localstat(Ywi,Ywi>0,1,3);                                               % maximum for stabilization of small WM structures
    Ywi(Ywi==0 & Ygi>0)=Ygi(Ywi==0 & Ygi>0); clear Ygir;                                  % mixing of both maps

    % Ywi = cat_vol_noPVE(Ywi,vx_vol,1);   
    Yis = cat_vol_approx(Ywi,'nh',resT2.vx_volr,4); Yis = cat_vol_smooth3X(Yis,4);        % first bias field map 
    Ywi(Ywi>Yis*1.1 | (Ywi>Yis*1.01 & ~Ybdr) | (Ywi<Yis*0.5 & Ywi>0)) = 0; 
    Yd  = cat_vbdist(single(Ybdr>1),Ybdr<=1,resT2.vx_volr);                                 % distance map the brain tissue 
    Ywi(Yd>8) = Yis(Yd>8); clear Yis Yd Ybr;                                              % use the first bias field far away from the brain
    Ywi = cat_vol_approx(Ywi,'nh',resT2.vx_volr,2); Ywi = cat_vol_smooth3X(Ywi,2);        % ... to use less filtering here
    Ywi = max(eps,cat_vol_resize(Ywi,'dereduceV',resT2));
    Ybf = Ywi; 
    Ym  = ( Ysrc./max(eps,Ywi) ) / mean( (Ysrc(Yp0(:)>2.9 & ~Ywmh(:)) ./ max(eps,Ywi(Yp0(:)>2.9 & ~Ywmh(:))))) * T3th(3); 
    Ym  = cat_main_gintnorm(Ym,Tth); 
    %%{
    if ~debug, clear Ywi; end
    % hard but good working correction for uncorrected high intensity bias (WM with intensity over 1)
    Ybv = (Ym - 2*Ydiv)>1.05 & Ym>1.05 & Ycls{2}<192;              
    Ywh = cat_vol_smooth3X(max(0,Ym-1) .* (Yp0>2.75 & Ym>1.05 & ~Ybv),4)*4;
    Ym = Ym - Ywh; 
    clear Ybv Ywh;
    %}
    cat_io_cmd(' ','','',debug,stime2);   
  end


  %% Optimization of the tissue segments
  %  ----------------------------------------------------------------------
  %  The correction of the SPM tissue classification is required especially
  %  in cases of failed bias correction, that strongly based on the T1
  %  intensities that were corrected before. So this is main part of LAS.
  %  ----------------------------------------------------------------------
  if 1;; % also this if is just for the debuging mode and contain the second block the correction of the segmentation 
    %% brain segmentation can be restricted to the brain to save time 
    stime = cat_io_cmd('  Prepare partitions','g5','',verb,stime); 
    [Ymr,Yb,BB] = cat_vol_resize({Ym,Yb},'reduceBrain',vx_vol,round(10/mvx),Yb);
    [Ygr,Ydivr,Yp0]   = cat_vol_resize({Yg,Ydiv,Yp0},'reduceBrain',vx_vol,BB.BB);
    Yl1               = cat_vol_resize(Yl1          ,'reduceBrain',vx_vol,round(4/mvx),BB.BB);
    Yclsr = cell(size(Ycls)); for i=1:6, Yclsr{i} = cat_vol_resize(Ycls{i},'reduceBrain',vx_vol,BB.BB); end
    %Ywtpm = cat_vol_resize(Ywtpm,'reduceBrain',vx_vol,round(4/mvx),BB.BB);


    % adaption of the LASstr depending on average basal values 
    LASmod  = min(2,max(0,mean((Ymr( NS(Yl1,LAB.BG) & Ygr<0.1 & Ydivr>-0.05  & Yclsr{1}>4)) - 2/3) * 8));
    LASstr  = min(1,max(0.05,LASstr * LASmod)); clear LASmod                 % adaption by local BG variation
    LASfs   = 1 / max(0.05,LASstr);                                          % smoothing filter strength 
    LASi    = min(8,round(LASfs));                                           % smoothing interation (limited)


    %  GM thickness (Ygmt) and percentage possition map (Ypp) estimation
    %  The Ypp and Ygmt maps are used to refine the GM especially to correct
    %  highly myelinated GM regions. Using allow to avoid overcorrections.
    fastppi = 1; % 1 mm: 15 vs. 40 seconds ... it thing fast is ok here 20161014 
    if fastppi
      [Ygmt,Ypp,Ywmd,Ycsfd] = cat_vol_pbt( (Yp0 + (Ymr*3 .* (Yp0>0)))/2 , ...
        struct('resV',vx_vol,'verb',0,'dmethod','vbdist','method','pbt2x') );
    else
      [Ygmt,Ypp,Ywmd,Ycsfd] = cat_vol_pbt( (Yp0 + (Ymr*3 .* (Yp0>0)))/2 , ...
        struct('resV',vx_vol,'verb',0) ); %#ok<*UNRCH>
    end
    Ygmtroi = Ygmt>0 & Ygmt<6 & NS(Yl1,LAB.CT); 
    GMTstat = [cat_stat_nanmedian(Ygmt(Ygmtroi(:))) cat_stat_nanstd(Ygmt(Ygmtroi(:)))]; 
    [D,I]   = cat_vbdist(single(Ygmt>0.1),Yp0>0,vx_vol); Ygmt = Ygmt(I); clear D I;  % full GMT map for the whole brain
    if ~debug, clear Ywmd Ygmtroi; end



    %% As far as SPM segmentation is not optimal we need some refinements.
    %  ----------------------------------------------------------------------
    stime = cat_io_cmd('  Prepare segments','g5','',verb,stime);
    if 1; 
      % use SPM segment 
      Ycm = Yp0>0.5 & Yp0<=1.5;
      Ygm = Yp0>1.5 & Yp0<=2.5;
      Ywm = Yp0>2.5;

      % corrections for original intensies 
      Ywm = Ywm & Ymr>2/3 & Ymr<1.3; 
      Ygm = Ygm & Ymr>1/3 & Ymr<1.0; 
      Ycm = Ycm & Ymr<1.5/3; % no lower limit here

      % add subcortical structures
      Yss = (NS(Yl1,LAB.BG) | NS(Yl1,LAB.TH)) & Ymr>2/3 & Ymr<2.75/3 & Ygr<0.1 & Ydivr>-0.05; 
      Yss = smooth3(Yss)>0.5;
      Ygm = Ygm | Yss;
      Ywm = Ywm & ~Yss & ~(NS(Yl1,LAB.BG) | NS(Yl1,LAB.TH));

      % blood vessel correction
      Ybv = (Ymr-Ydivr)>1.2 & Ymr>0.9 & Yclsr{2}<192;
      Ywm = Ywm & ~Ybv; 
      Ygm = Ygm & ~Ybv; 
      Ycm = Ycm & ~Ybv; 

      % correction by divergence to avoid gyri in the GM segment 
      Ymlwm = Ygm & (Ymr - 2*Ydivr + (Ygmt-3)/10 + 0.1*NS(Yl1,LAB.CB))>0.9 & ...
              ~Yss & Ymr>(2.25 + min(0.25,max(0,mvx)))/3  & ...
              cat_vol_morph((Ymr - 2*Ydivr + 0.2*NS(Yl1,LAB.CB))>1.75/3,'e',1/mvx);  
      Ymlwm = Ymlwm | (smooth3(Yp0>1.5 & Ymr>2.5/3 & Ymr<3.2/3 & (Ydivr<0 | Ygr<0.1) & ~Ybv)>0.3 & Ymr>2.5/3);
      Ymlwm = (Ymlwm & ~Yss & (Ymr - Ydivr)>0.9) | (cat_vol_morph(Ymlwm,'e') & (Ymr - Ydivr)>0.8 & ~Yss); 
      Ygm   = Ygm & ~Ymlwm; 
      Ywm   = Ywm | Ymlwm; 
      if ~debug, clear Ymlwm; end
      
      % correction for WMHs
      Yvt  = cat_vol_morph( NS(Yl1,LAB.CB) | NS(Yl1,LAB.BS) | NS(Yl1,LAB.PH) | ...
                            NS(Yl1,LAB.BG) | NS(Yl1,LAB.TH) | NS(Yl1,LAB.HC)  ,'d',5/mvx); 
      Ywmh = smooth3( (Yp0/3-Ymr)>0.1 & Ymr<2.5/3 & Ymr>1.5/3 & ~Yss & ~Yvt)>0.5;
      Ywmh = cat_vol_morph( NS(Yl1,LAB.VT) | Ywmh ,'c',4) & ~NS(Yl1,LAB.VT); 
      Ywm  = Ywm & ~Ywmh; 
      Ygm  = Ygm & ~Ywmh; 
      if ~debug, clear Yss; end


      % added undetected GM around WM 
      Ymlgm = ~Ywm & Yb & cat_vbdist(single(Ywm),Yb,vx_vol)<3 & ... 
        Ymr>(1.25 + min(0.25,max(0,mvx)))/3 & Ymr<2.5/3 &  ...
        abs(Ydivr)<0.1 & cat_vol_morph(NS(Yl1,LAB.CT),'d',2/mvx) & ...
        ~cat_vol_morph(NS(Yl1,LAB.VT) | Ywmh,'d',2/mvx); 
      Ygm   = Ygm | Ymlgm; 
      Ycm   = Ycm & ~Ymlgm;
      if ~debug, clear Ymlgm; end

      % correct for uncorrected WM bias
      Ywh  = cat_vol_smooth3X(Ywm & Ymr>1.1 & ~Ybv,3)>0.05;
      Ygmc = Ywm & (Ywh & Ymr<2.75/3); 
      Ycmc = Ygm & (Ywh & Ymr<0.75); 
      if ~debug, clear Ywh; end
      Ygm  = (Ygm & ~Ycmc) | Ygmc; 
      Ywm  = Ywm & ~Ygmc; 
      Ycm  = Ycm | Ycmc;
      if ~debug, clear Ygmc Ycmc; end

      %% correction by divergence to avoid sulci in the GM segment
      Ypvec = Ygm & (Ymr - 2*Ydivr - (Ygmt)/10 - Ypp/10)<0.5 & Ymr<1.25/3;  
      Ygm   = Ygm & ~Ypvec; 
      Ycm   = Ycm | Ypvec; 
      if ~debug, clear Ypvec; end

      
      % myelinated GM  -  this is the major goal of LAS !!!
      Ymgm = ~(cat_vbdist(single(Ypp<0.1),Yb,vx_vol)>(min(2.5,max(1,GMTstat(1)+GMTstat(2))))) & ~Yvt & ...
                cat_vol_smooth3X(Ygmt,4)<=min(2,max(1,GMTstat(1)-GMTstat(2)/2)); 
      Ymgm = Ymgm | (Yclsr{1}>0 & Ycsfd<min(2,max(1,GMTstat(1)-GMTstat(2)/2)) & ~Ycm & ~Yvt);
      if ~debug, clear Ycsfd; end
      Ymgm = Ymgm & Yp0>2.0 & Ymr<0.95 & Ymr>0.8 & Ymr<1 & Ywm & ~Ybv & ~Ywmh & Ymr>2/3 & ~Yvt & ~Ycm & NS(Yl1,LAB.CT); 
      Yhd  = cat_vbdist(single(~Yb),Yb,vx_vol);
      Ymgm = Ymgm | (Yp0>2.0 & Ymr>0.5 & Yhd<min(2,max(1,GMTstat(1)-GMTstat(2)/2)) & ~Yvt & ~Ybv);
      if ~debug, clear Yvt Yhd Ywmh Ybv; end
      Ygm  = Ygm | Ymgm; 
      Ywm  = Ywm & ~Ymgm;
      Ycm  = Ycm & ~Ymgm;
      if ~debug, clear Ymgm; end

      % parahippocampale gyrus and hippocampus
      Yhcg  = (NS(Yl1,LAB.PH) | NS(Yl1,LAB.HC)) & (Ymr - Ydivr)>2.5/3;
      Ywm  = Ywm | Yhcg; 
      Ygm  = Ygm & ~Yhcg;
      if ~debug, clear Yhcg; end

      % Ycsf
      Yngm = Ygm & Ypp<0.1 & Ymr<0.5 & Ygmt<=GMTstat(1); 
      Ygm  = Ygm & ~Yngm;
      if ~debug, clear Yngm Ygmt Ypp; end

      %% remove sinus (rectus) and brainstem
      Ysr = cat_vbdist(single(NS(Yl1,LAB.CB)),NS(Yl1,LAB.CT) | NS(Yl1,LAB.BV),vx_vol)<30 & ...
            cat_vol_morph(NS(Yl1,LAB.BV),'d',10/mvx) & ...
            (NS(Yl1,LAB.CT)) & Ymr>0.5/3 & Ymr<1.8/3 & ~cat_vol_morph(NS(Yl1,LAB.CT) & Ymr>2/3,'d',1);
      Ysr(smooth3(Ysr)<0.5)=0; 
      Ygm = Ygm & ~Ysr & ~cat_vol_morph(NS(Yl1,LAB.BS),'d',2/mvx); 
      if ~debug, clear Ysr; end

      %% non CSF
      Ycm = Ycm & (Ymr-Ydivr)<.35 & Yclsr{3}>240 & abs(Ydivr)<0.1 & Ygr<0.2; 
      Ycm(smooth3(Ycm)<0.5)=0;  
      Ygm = Ygm | (~Ycm & (Ymr-Ydivr)>.4 & Ygr<0.2 & Ymr<0.5 & Yb); 
      Ygm(smooth3(Ygm | Ywm)<0.5)=0;  
      if ~debug, clear Ymr Ydivr Ygr Yclsr; end 
    end
    

    %% back to original resolution for full bias field estimation
    [Ycm,Ygm,Ywm] = cat_vol_resize({Ycm,Ygm,Ywm},'dereduceBrain',BB); 
    [Yp0,Yl1]     = cat_vol_resize({Yp0,Yl1},'dereduceBrain',BB);
    [Yb]          = cat_vol_resize({Yb},'dereduceBrain',BB);
  
  
    %% new head tissues
    Ynb  = smooth3(Ycls{6})>128 & Ysrc<mean(Tth.T3th(2:3)) & Ym<1/6; % ... zero background
    Yhdh = Yp0==0 & (Ym - Ydiv./Ym)>0.5 & (Ym - Yg)>0.8 & (Ym>1.2 | smooth3(Yg)>0.1) & ~Ybg & Ybd>2 & (Ym>1.2 | Ym>max(0.4,(1-Ybd/300)));
    Yhdm = Yp0==0 & Ym<1.0 & Ym>0.5; 
    Yhdm = Yhdm | (cat_vol_morph(Yp0>1.5,'d',3) & Yp0<1.5 & cat_vol_smooth3X(Yhdh>0.5,8)>0.1 & ~Yhdh & Ym>0.5 & Ym<0.9);
    Yhdm = Yhdm | ((smooth3(Yhdh | Ycls{6}>128)<0.2) & ~Ygm & ~Ywm & ~Ycm & (Ym - 2*Ydiv)>0.5 & Ym>0.6);   
    Yhdm = Yhdm & (Yg./Ym)<0.5 & abs(Ydiv./Ym)<0.1 & ~Ybg & ~Yhdh;
    Yhdm = Yhdm | (Ym>0.6 & Yp0<1.5 & (Ybd ./ (Ybd + Ybgd))<0.2 & Ybgd<20);
    Yhdm((smooth3(Yhdm) - smooth3(Yhdh)/1.5)<0.6) = 0;
    if ~debug, clear Ybg Yg; end 
  end
  
  
  
  
  %% Estimation of local peaks and creation of normalized T1 maps
  %  --------------------------------------------------------------------- 
  %  Estimation of the local WM threshold with "corrected" GM voxels to
  %  avoid overfitting (see BWP cerebellum). 
  %  CSF is problematic in high contrast or skull-stripped image should 
  %  not be used here, or in GM peak estimation
  %  ---------------------------------------------------------------------
  if 1;; % the last if block for debuging mode 
    mres  = 1.1; 
    stime = cat_io_cmd('  Estimate local tissue thresholds','g5','',verb,stime);  if debug, fprintf('\n'); end
    stime2 = cat_io_cmd('    WM intensity','g5','',debug);  
    %%
    Ysrcm = Ywm & ~(NS(Yl1,LAB.BG) | NS(Yl1,LAB.TH));
    Ysrcm(smooth3(Ysrcm>0)<0.5)=0;
    Ysrcm = Ysrc.*Ysrcm;
    Ysrcm = cat_vol_localstat(Ysrcm,Ysrcm>0,1,3);
    Ysrcm2 = Ysrcm; 
    YM = cat_vol_morph(Ysrcm>0,'d') & Ysrcm==0; Ysrcm2(YM)=Ysrc(YM);
    Ysrcm2 = cat_vol_localstat(Ysrcm2,Ysrcm2>0,2,3); Ysrcm(YM) = Ysrcm2(YM);
    Ysrcm = cat_vol_noPVE(Ysrcm,vx_vol,2) * T3th(3)/mean(Ysrcm(Ysrcm(:)>0)); 
    
    % major correction outside the brain 
    Ygmw = Ygm & Ym<0.9 & Ym>0.6 & abs(Ydiv)<0.1;
    Ygi = cat_vol_noPVE(Ysrc .* (Ygmw & ~Ysrcm)  * T3th(3)/mean(Ysrc(Ygmw(:))),vx_vol,2) + ...
          cat_vol_noPVE(Ysrc .* (Ycm  & ~Ysrcm)  * T3th(3)/mean(Ysrc(Ycm(:))),vx_vol,2) + ...
          ... cat_vol_noPVE(Ysrc .* Yhdh * T3th(3)/median(Ysrc(Yhdh(:)>0 & Ybd(:)<15)))
          cat_vol_noPVE(Ysrc .* (Yhdm & ~Ysrcm)  * T3th(3)/mean(T3th(2)),vx_vol,2);
    Ygi = cat_vol_noPVE(Ygi,vx_vol,2);  Ygi(smooth3(Ygi>0)<0.5)=0;     
    Ysrcm = Ysrcm + Ygi.*(Ysrcm==0); Ysrcm = cat_vol_noPVE(Ysrcm,vx_vol,2); 
    Ysrcm(Ybd>20) = Ybf(Ybd>20); 
    % fine correction
%    Ygi         = cat_vol_resize(Ygi,'reduceV',vx_vol,mres,32,'meanm'); % maximum reduction for the WM
    [Yi ,resT2] = cat_vol_resize(Ysrcm,'reduceV',vx_vol,mres,32,'meanm'); % maximum reduction for the WM
    if ~debug, clear Ysrcm Ydiv; end
    Yi(smooth3(Yi>0)<0.5)=0; Yi(smooth3(Yi>0)<0.5)=0;
    
 %   Yi(Yi==0 & Ygi>0)=Ygi(Yi==0 & Ygi>0); if ~debug, clear Ygi; end
    Yi = cat_vol_approx(Yi,'nh',resT2.vx_volr,2); Yi = cat_vol_smooth3X(Yi,LASfs*2); 
    Ylab{2} = max(eps,cat_vol_resize(Yi,'dereduceV',resT2)); 
    if ~debug, clear Yi; end


    %% GM
    stime2 = cat_io_cmd('    GM intensity','g5','',debug,stime2);  
    Yi = Ysrc ./ max(eps,Ylab{2}) .* (Ygm | Yhdm);
    Ybs = NS(Yl1,LAB.BS) & Ym<1.1 & Ym>0.9 & Yp0>2.5 & Ym>0.9;
    Yi(Ybs)  = Ysrc(Ybs) ./ max(eps,Ylab{2}(Ybs))   .* T3th(2)/T3th(3); clear Ybs;
    Yi = cat_vol_noPVE(Yi,vx_vol,1);
    Yi(Ybd>20 & Yi==0) = 2/3; 
    if ~debug, clear Yhdm; end
    [Yir,Ygmr,resT2] = cat_vol_resize({Yi,Ygm},'reduceV',vx_vol,mres,32,'meanm'); if ~debug, clear Yi; end
    Yir = cat_vol_noPVE(Yir,vx_vol,1);
    Yir = cat_vol_approx(Yir,'nh',resT2.vx_volr,2); 
    Yir = cat_vol_smooth3X(Yir,LASfs); 
    Ylab{1} = cat_vol_resize(Yir,'dereduceV',resT2).*Ylab{2};   
    if ~debug, clear Yir Ybd Yl1; end


    %% CSF & BG 
    stime2 = cat_io_cmd('    CSF and BG intensity','g5','',debug,stime2);  
    Ynb = Ynb | smooth3( (Ycls{4}>128 | Ycls{6}>128) & Ym<median(Ym(Ycls{4} & Ym<0.4)))>0.5; if ~debug, clear Ym; end
    Ynb = cat_vol_morph(Ynb,'e',2/mvx);
    Ynb = Ynb & Ysrc~=0; 
    [Yc,resT2] = cat_vol_resize(Ysrc ./ max(eps,Ylab{2}) .* (smooth3(Ycm)>0.5),...
       'reduceV',vx_vol,8,16,'min');% only pure CSF !!!
    Ynbr = cat_vol_resize(Ysrc ./ max(eps,Ylab{2}) .* Ynb,...
       'reduceV',vx_vol,8,16,'meanm');
    Ynbr(Yc>0)=0; Yc(Ynbr>0)=0;
    for xi=1:2*LASi, Ynbr = cat_vol_localstat(Ynbr,Ynbr>0,2,1); end
    for xi=1:2*LASi, Yc  = cat_vol_localstat(Yc,Yc>0,2,1); end
    Ynba = cat_vol_approx(Ynbr,'nh',resT2.vx_volr,2); 
    Yca  = cat_vol_approx(Yc  ,'nh',resT2.vx_volr,2); % + min(max( meanYnb + stdYbc , meanYc - stdYbc ),...
    if ~debug; clear Yc Ynbr Ycls; end
    % Yca  = max(Yca,Ynba*1.5); 
    Yca   = max(Yca,mean(Tth.T3th(1:2))/T3th(3)); 
    Ynba  = min(Yca/1.5,Ynba); 
    Yca   = max(Yca,Ynba + 0.01*1.5); 

    Ynba = cat_vol_smooth3X(Ynba,LASfs*2); 
    Yca  = cat_vol_smooth3X(Yca,LASfs*2); % * meanYc/mean(Yca(:)); 
    Ylab{3} = cat_vol_smooth3X(cat_vol_resize(Yca,'dereduceV',resT2).*Ylab{2},LASfs*2);  
    Ylab{6} = cat_vol_smooth3X(cat_vol_resize(Ynba,'dereduceV',resT2).*Ylab{2},LASfs*2);
    clear Ynba Yca; 
    
    
    %% back to original resolution
    if any(resTb.vx_vol ~= resTb.vx_volr), 
      Ysrc = Ysrco; clear Ysrco; 
      for i=1:6, if ~isempty(Ylab{i}), Ylab{i} = cat_vol_resize(Ylab{i},'dereduceV',resTb); end; end
      Yb  = cat_vol_resize(Yb,'dereduceV',resTb);
      Yp0 = cat_vol_resize(Yp0,'dereduceV',resTb);
      Ywm = cat_vol_resize(Ywm,'dereduceV',resTb);
    end 
    
    
    %% local intensity modification of the original image
    % --------------------------------------------------------------------
    stime2 = cat_io_cmd('    Intensity mapping','g5','',debug,stime2);  
    Yml = zeros(size(Ysrc));  
    Yml = Yml + ( (Ysrc>=Ylab{2}                ) .* (3.0 + (Ysrc-Ylab{2}) ./ max(eps,Ylab{2}-Ylab{1})    ) );
    Yml = Yml + ( (Ysrc>=Ylab{1} & Ysrc<Ylab{2} ) .* (2.0 + (Ysrc-Ylab{1}) ./ max(eps,Ylab{2}-Ylab{1})    ) );
    Yml = Yml + ( (Ysrc>=Ylab{3} & Ysrc<Ylab{1} ) .* (1.0 + (Ysrc-Ylab{3}) ./ max(eps,Ylab{1}-Ylab{3})    ) );
    Yml = Yml + ( (Ysrc>=Ylab{6} & Ysrc<Ylab{3} ) .* (1/3 + (Ysrc-Ylab{6}) ./ max(eps,Ylab{3}-Ylab{6})*2/3) );
    Yml = Yml + ( (Ysrc< Ylab{6}                ) .* (0.0 + (Ysrc        ) ./ max(eps,Ylab{6}*3          )) );
    Yml(isnan(Yml) | Yml<0)=0; Yml(Yml>10)=10;
    Yml = Yml/3;

    %% global
    Ymg = Ysrc ./ max(eps,Ylab{2}) * Tth.T3th(5);
    Ymg = cat_main_gintnorm(Ymg,Tth); 
    if ~debug, clear Ylab Ysrc; end

    % fill up CSF in the case of a skull stripped image 
    if max(res.mn(res.lkp==5 & res.mg'>0.1)) < mean(res.mn(res.lkp==3 & res.mg'>0.3))
      YM   = cat_vol_morph(Yb,'d'); 
      Ymls = smooth3(max(Yml,YM*0.5));
      Yml(YM & Yml<0.5)=Ymls(YM & Yml<0.5); 
      clear Ymls YM
    end
    cat_io_cmd(' ','','',debug,stime2);  
  end
  
  
  
  
  %% final corrections
  stime = cat_io_cmd('  Final corrections','g5','',verb,stime);  
  % prepare class correction 
  Yp0n = cat_vol_smooth3X(Yml,noisef) .* Yb;
  Yp0n = min(Yp0n,2-(smooth3(Ywm)>0.3));
  Yp0n(Yp0n>7/6)=0;
  Yp0n(Yp0<1.1 & Yp0n<1.5/3 & Yp0>1/3) = 1/3; % set CSF
  Yp0n(smooth3(Yp0n<1.5/3)>0.5 & Yp0n>1/3) = 1/3;
  Yp0n(cat_vol_morph(Yp0n>1.25/3,'labopen',1/mvx)==0 & Yp0n>1/3)=1/3;   % remove small tissue dots 
  Yp0n(smooth3((Yml<1/6 | Yml>0.45) & Yp0n<0.34 & Yp0<1.5)>0.5)=0;      % remove meninges
  Yp0n(cat_vol_morph(Yp0n>0.5/3,'labopen',1/mvx)==0 & Yp0n<=0.34)=0;
  Yp0n(cat_vol_morph(Yp0n>0.5/3,'labclose',2/mvx) & Yp0n<=0.34)=1/3;
  Yp0n = Yp0n*3;
  Yp0n = cat_vol_median3(Yp0n,Yp0n>0,Yp0n>0);
  if ~debug, clear Yp0; end
  
  %% update classes
  Yp0toC  = @(Yp0,c) 1-min(1,abs(Yp0-c));
  Ycls{1} = cat_vol_ctype(Yp0toC(Yp0n,2)*255,'uint8');
  Ycls{2} = cat_vol_ctype(Yp0toC(Yp0n,3)*255,'uint8');
  Ycls{3} = cat_vol_ctype(Yp0toC(Yp0n,1)*255,'uint8');
  Ycls{6} = cat_vol_ctype(cat_vol_morph(smooth3(1-max(0,Ymg*3) - Yp0n)>0.75,'lc')*255,'uint8'); 
  Ycls{5} = cat_vol_ctype((Ymg<4/6 & ~Yp0n & ~Ycls{6})*255,'uint8'); 
  Ycls{4} = cat_vol_ctype((Ymg>5/6 & ~Yp0n & ~Ycls{6})*255,'uint8'); 
  cat_io_cmd('','','',verb,stime);


end
function Ygi = cat_vol_noPVE(Ygi,vx_vol,mres,dist,iter)
  if ~exist('dist','var'), dist=2; end
  if ~exist('iter','var'), iter=1; end
  if ~exist('mres','var'), mres=1; end
  
  if any(vx_vol<0.5)
    YM = Ygi>0;
    [Ygi,resT2] = cat_vol_resize(Ygi,'reduceV',vx_vol,mres,32,'meanm');
  end
  for i=1:iter
    Ygi      = cat_vol_median3(Ygi,Ygi>0,Ygi>0);
    [Ygistd,Ygimean] = cat_vol_localstat(Ygi,Ygi>0,dist,4); 
    Ygx      = Ygi<(Ygimean-Ygistd/4) | Ygi>(Ygimean+Ygistd/2); 
    Ygi(Ygx) = Ygimean(Ygx);
    clear Ygx;
  end
  if any(vx_vol<0.5)
    Ygi = cat_vol_approx(Ygi,'nh',resT2.vx_volr,4); 
    Ygi = cat_vol_resize(Ygi,'dereduceV',resT2);
    Ygi(~YM)=0;
  end
end