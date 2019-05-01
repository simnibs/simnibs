function [Ya1,Ycls,YBG,YMF] = cat_vol_partvol(Ym,Ycls,Yb,Yy,vx_vol,extopts,Vtpm,noise)
% ______________________________________________________________________
% Use a segment map Ycls, the global intensity normalized T1 map Ym and 
% the atlas label map YA to create a individual label map Ya1. 
% The atlas contain main regions like cerebrum, brainstem, midbrain,
% cerebellum, ventricle, and regions with blood vessels. 
%
% This function try to solve the following problems:
%  1) Finding of the cerebrum, the cerebellum, the head, blood vessels, 
%     brain skin and other mayor structures based on atlas (YA) and 
%     tissue class information (Yp0). 
%     To do this it is important to use data from the T1-map (Ym) that
%     use the same intensity scaling as the segment map Yp0, but have 
%     more informations about partial volume regions.
%  2) Set Partions:
%     2.1) Find biggest WM part of each region.
%     2.2) Align the nearest region class for other voxel
%     2.3) Finding and Filling of the ventricle and the Basalganglia
%     2.4) Find blood vessels
%     2.5) Brain extraction
%     2.6) Side alignment
% ______________________________________________________________________
%
% Structure:
%
%   [vol,Ya1,Yb,YMF] = cat_vol_partvol(YA,Yp0,Ym,Yl0,opt)
%
%   INPUT:  YA = 3D-volume with brain regions (altas map)
%           Yp0  = 3D-volume with tissue propability map (CSF=1,GM=2;WM=3)
%           Ym   = intensity normalized T1 image (BG=0,CSF=1/3,GM=2/3,WM=1)
%           Yl0  = spm-classes 4-6 (intracranial=1,skull=2,background=3)
%           opt
%            .res    = resolution for mapping
%            .vx_vol = voxelsize
%            .LAB    = label of Ya1 map (see LAB definition below)
%            
%
%   OUTPUT: vol = structure with volumes
%           Ya1 = individual label map 
%           Yb  = brain mask
%           YMF = filling mask for ventricle and subcortical structures
%
% ______________________________________________________________________
% Structural Brain Mapping Group, University Jena, Germany
% Robert Dahnke
% 2013/05
%
% $Id: cat_vol_partvol.m 1149 2017-06-23 10:31:47Z dahnke $

% ______________________________________________________________________
%
% Development comments:
%
%   Was ist neu im Vergleich zu anderen?
%   - Zuweisung durch Dartel mit hoher Genauigkeit m??glich 
%   - Erweiterung von SPM/VBM durch MainROIs (Seiten, Lappen, ...)
%   - Verbesserung der SPM/VBM durch bessere Enfernung von unerw??nschtem
%     Gewebe (ON, Blutgef????e ...)
%   - Blutgef????e k??nnnen als erweitere Masken f??r fMRI genutzt werden um
%     Seiteneffekte besser ausblenden zu k??nnen.
%  [- Beliebige Atlanten k??nnen genutzt werden.]
%
%  Todo:
%   - Besserer Atlas
%   - BV vs. HD - gl??tten in dilated HD region
%   - F??llen von CSF l??cken bei LAB~=BV und Ym<1.2 und LAB==NV?
%
% ______________________________________________________________________



% ----------------------------------------------------------------------
% fast partitioning for B3C[, and LAS]
% ----------------------------------------------------------------------
% VBM atlas atlas map to find important structures for the LAS and the
% skull-stripping, which are the subcortical GM regions and the cerebellum.
% Maybe also WM hyperintensity have to be labeled here as a region without
% local correction - actual clear WMHs are handeled as GM.
% ----------------------------------------------------------------------
  
  % definition of ROIs 
  
%   LAB.CT =  1; % cortex
%   LAB.MB = 13; % MidBrain
%   LAB.BS = 13; % BrainStem
%   LAB.CB =  3; % Cerebellum
%   LAB.ON = 11; % Optical Nerv
%   LAB.BG =  5; % BasalGanglia 
%   LAB.TH =  9; % Hypothalamus 
%   LAB.HC = 19; % Hippocampus 
%   LAB.VT = 15; % Ventricle
%   LAB.NV = 17; % no Ventricle
%   LAB.BV =  7; % Blood Vessels
%   LAB.NB =  0; % no brain 
%   LAB.HD = 21; % head
%   LAB.HI = 23; % WM hyperintensities
%   LAB.PH = 25; % Gyrus parahippocampalis

  LAB     = extopts.LAB;
  BVCstr  = extopts.BVCstr; 
  WMHCstr = extopts.WMHCstr; 
  verb    = extopts.verb-1;
  debug   = extopts.verb>2;
  PA      = extopts.cat12atlas;
  vx_res  = mean([max(vx_vol) min(vx_vol)]); % cat_get_defaults('extopts.vx_res'); 
  

  %% map atlas to RAW space
  if verb, fprintf('\n'); end
  stime = cat_io_cmd('  Atlas -> subject space','g5','',verb); dispc=1;
  VA = spm_vol(PA{1});
  YA = cat_vol_ctype(spm_sample_vol(VA,double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),0));
  YA = reshape(YA,size(Ym));
 % clear Yy;
  
  Yp0  = (single(Ycls{1})*2/255 + single(Ycls{2})*3/255 + single(Ycls{3})/255) .* Yb; 
  Yp0A = single(spm_sample_vol(Vtpm(1),double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),1))*2 + ...
         single(spm_sample_vol(Vtpm(2),double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),1))*3 + ...
         single(spm_sample_vol(Vtpm(3),double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),1))*1;
  Yp0A = reshape(Yp0A,size(Ym));   
  
  % work on average resolution
  Ym0 = Ym; 
  [Ym,YA,Yp0,Yb,Yp0A,BB] = cat_vol_resize({Ym,YA,Yp0,Yb,Yp0A},'reduceBrain',vx_vol,2,Yb);
  [Ym,Yp0,Yb,Yp0A,resTr] = cat_vol_resize({Ym,Yp0,Yb,Yp0A},'reduceV',vx_vol,vx_res,64);
  [YA]              = cat_vol_resize(YA ,'reduceV',vx_vol,vx_res,64,'nearest'); 
  
  spm_smooth(Ym,Ym,0.9./vx_vol);
  
  vx_vol = resTr.vx_volr; 
  vxd    = 1/mean(vx_vol); 
  
  % prepare maps
  [tmp0,tmp1,YS] = cat_vbdist(single(mod(YA,2)) + single(YA>0)); YS=~mod(YS,2); clear tmp0 tmp1;  % side map
  YA(mod(YA,2)==0 & YA>0)=YA(mod(YA,2)==0 & YA>0)-1;                    % ROI map without side
  YA   = cat_vol_ctype(cat_vol_median3c(single(YA),Yp0>0));
  Yg   = cat_vol_grad(Ym,vx_vol);
  Ydiv = cat_vol_div(Ym,vx_vol); Ymo=Ym; 
  Ym   = Ym*3 .* (Yb);
  Yb   = Yb>0.5;

  
  
  %% Create individual mapping:
  stime = cat_io_cmd('  Major structures','g5','',verb,stime); dispc=dispc+1;
 % noise = double(max(0.02,min(0.1,mean(Yg(cat_vol_morph(Yp0>2.8,'lc')))/3)));
  noise = double(noise);
  
  % Major structure mapping:
  % Major structure mapping with downcut to have a better alginment for 
  % the CB and CT. Simple setting of BG and TH as GM structures. 
  Ya1  = zeros(size(Ym),'single');
  
  % Basal Ganglia
  Ybg  = zeros(size(Ym),'single');
  Ybgd = cat_vbdist(single(YA==LAB.BG),Yb,vx_vol); 
  Yosd = cat_vbdist(single(YA==LAB.TH | YA==LAB.VT | YA==LAB.HC  | YA==LAB.BS | (YA==LAB.CT & Ym>2.9)),Yb,vx_vol); 
  Ybg(smooth3(Yosd>3 & Ybgd<5  & Ym>1.9 & Ym<2.85 & Yg<4*noise & ((Ybgd<1 & Ydiv>-0.01) | (Ydiv>-0.01+Ybgd/100)))>0.7)=1;
  Ybg(smooth3((Ybg==0 & Yp0>2.8 & Ym>2.8 & YA==LAB.CT) | Ym>2.9 | YA==LAB.TH | YA==LAB.HC | Yosd<2 | ...
    (Ybg==0 & Yp0<1.25) | (Ybg==0 & Ybgd>8) | (Ybg==0 & Ydiv<-0.01+Ybgd/200))>0.3)=2;
  Ybg(Ybg==0 & Ybgd>0 & Ybgd<10)=1.5; 
  Ybg = cat_vol_laplace3R(Ybg,Ybg==1.5,0.005)<1.5 & Ym<2.9 & Ym>1.8 & Ydiv>-0.02;
  Ya1(Ybg)=LAB.BG;                                                      % basal ganglia
  Ya1(YA==LAB.TH & Ym>1.9 & Ym<2.85 & Ydiv>-0.1)=LAB.TH;                % thalamus
  Ya1(YA==LAB.HC & Ym>1.9 & Ym<2.85 & Ydiv>-0.1)=LAB.HC;                % hippocampus
  Ya1(((Yp0>2.5 & Ym>2.5 & (YA==LAB.CT | YA==LAB.BG)) | (Yp0>1.5 & Ym<3.5 & Ybgd>1 & Ybgd<8 & (Ybgd>4 | Ydiv<-0.02+Ybgd/200))) & Ya1==0)=LAB.CT;  % cerebrum
  Ya1((Yp0>2.0 & Ym>2.0) & YA==LAB.CB)=LAB.CB;                          % cerebellum
  Ya1((Yp0>2.0 & Ym>2.0) & YA==LAB.BS)=LAB.BS;                          % brainstem
  Ya1((Yp0>2.0 & Ym>2.0) & YA==LAB.ON)=LAB.ON;                          % optical nerv
  Ya1((Yp0>2.0 & Ym>2.0) & YA==LAB.MB)=LAB.MB;                          % midbrain
  clear Ybg Ybgd; 
  %% region-growing
  Ya1(Ya1==0 & Yp0<1.9)=nan; 
  Ya1 = cat_vol_downcut(Ya1,Ym,4*noise); Ya1(isinf(Ya1))=0; 
  Ya1 = cat_vol_median3c(Ya1,Yb);                                       % smoothing
  Ya1((Yp0>1.75 & Ym>1.75 & Yp0<2.5 & Ym<2.5) & Ya1==LAB.MB)=0;         % midbrain correction
  
  
  
  %% Ventricle:
  % Ventricle estimation with a previous definition of non ventricle CSF
  % to have a second ROI in the region-growin. Using only the ventricle
  % ROI can lead to overgrowing. Using of non ventrilce ROI doesn't work
  % because dartle failed for large ventricle. 
  % It is important to use labopen for each side!
  stime = cat_io_cmd('  Ventricle detection','g5','',verb,stime); dispc=dispc+1;
  Ynv = cat_vol_morph(cat_vol_morph(~Yb,'d',4*vxd) | (YA==LAB.CB | YA==LAB.BS),'d',2) | cat_vol_morph(YA==LAB.NV,'e',1);
  Ynv = single(Ynv & Ym<2 & ~cat_vol_morph(Yp0<2 & (YA==LAB.VT) & Yg<0.2,'d',4*vxd));
  Ynv = smooth3(round(Ynv))>0.5; 
  % between thamlamus
  Ynv = Ynv | (cat_vol_morph(Ya1==LAB.TH,'c',10) & Yp0<2) | YA==LAB.CB | YA==LAB.BS;
  Ynv = smooth3(Ynv)>0.8;
  Yvt = single(smooth3(Yp0<1.5 & (YA==LAB.VT) & Yg<0.25 & ~Ynv)>0.7); 
  Yvt(Yvt==0 & Ynv)=2; Yvt(Yvt==0 & Ym>1.8)=nan; Yvt(Yvt==0)=1.5;
  Yvt2 = cat_vol_laplace3R(Yvt,Yvt==1.5,0.005);
 
  warning('off','MATLAB:cat_vol_morph:NoObject');
  Yvts1 = cat_vol_morph(Yvt2<1.5 & YS==0,'lo',1);
  if sum(Yvts1)==0, Yvts1 = cat_vol_morph(smooth3(Yvt2<1.5 & YS==0)>0.6,'l'); end
  if sum(Yvts1)==0, Yvts1 = smooth3(Yvt2<1.5 & YS==0)>0.5; end
  Yvts2 = cat_vol_morph(Yvt2<1.5 & YS==1,'lo',1); 
  if sum(Yvts2)==0, Yvts2 = cat_vol_morph(smooth3(Yvt2<1.5 & YS==1)>0.6,'l'); end
  if sum(Yvts2)==0, Yvts2 = smooth3(Yvt2<1.5 & YS==1)>0.5; end
  warning('on','MATLAB:cat_vol_morph:NoObject');
  
  Yvt = smooth3((Yvts1 | Yvts2 | (YA==LAB.VT & Ym<1.7)) & Yp0<1.5 & Ym<1.5)>0.5; 
  Ya1(Yvt)=LAB.VT; 
  clear Yvts1 Yvts2;

  
  
  %% Blood vessels
  % For this we may require the best resolution!
  % first a hard regions growing have to find the real WM-WM/GM region
  if BVCstr
    stime = cat_io_cmd('  Blood vessel detection','g5','',verb,stime); dispc=dispc+1;
    Ywm = Yp0>2.5 & Ym>2.5 & Yp0<3.1 & Ym<4;                              % init WM 
    Ywm = Ywm | (cat_vol_morph(Ywm,'d') & Ym<3.5); 
    %%
    Ywm = single(cat_vol_morph(Ywm,'lc',2));                              % closing WM               
    Ywm(smooth3(single(Ywm))<0.5)=0;                                      % remove small dots
    Ywm(~Ywm & (Yp0<0.5 | Ym<1.2 | Ym>4))=nan;                            % set regions growing are
    [Ywm1,YDr] = cat_vol_downcut(Ywm,Ym,2*noise*(1-BVCstr/2));            % region growing
    Ywm(Ywm==-inf | YDr>20)=0; Ywm(Ywm1>0)=1; clear Ywm1                  % set regions growing
    % smoothing
    Ywms = smooth3(single(Ywm)); Yms=smooth3(Ym);                         
    Ywm(Ywms<0.5)=0; Ywm(Ywms>0.5 & Yb & (Ym-Yms)<0.5)=1;                 
    Ywm(Ywms<0.5 & Yb & (Ym-Yms)>0.5)=0; clear Ywms                       
    %% set blood vessels
    Ybv=cat_vol_morph( (Ym>3.75-(0.5*BVCstr) & Yp0<2+(0.5*BVCstr)) | ... % high intensity, but not classified as WM (SPM)
      (Yms>2.5 & (Ym-Yms)>0.6) | ...                                     % regions that strongly change by smoothing
      (Ym>2.5-(0.5*BVCstr) & Ywm==0) | ...                               % high intensity, but not classified as WM (SPM)
      (Ym>2.5-(0.5*BVCstr) & Yp0<2+(0.5*BVCstr) & Ya1==0 & YA==LAB.CT),'c',1) & ...
      cat_vol_morph(Ya1==LAB.CT,'d',2) & ~cat_vol_morph(Ya1==LAB.HC,'d',2) & ...
      cat_vol_morph((Ya1==0 | Ya1==LAB.CT | Ya1==LAB.BV | Ym>1.5) & Ya1~=LAB.VT & Yp0<2.5,'e',1) & ... avoid subcortical regions
      ~Ywm;  clear Ywm 
    Ybb = cat_vol_morph(Yp0>0.5,'lc',1); 
    Ybv = (Ybv & Ybb) | smooth3(Yp0<0.5 & Ybb)>0.4; clear Ybb; 
    %% smoothing
    Ybvs = smooth3(Ybv);
    Ybv(Ybvs>0.3 & Ym>2.5 & Yp0<2.5)=1; Ybv(Ybvs>0.3 & Ym>3.5 & Yp0<2.9)=1;
    Ybv(Ybvs<0.2 & Ym<4-2*BVCstr)=0; clear Yvbs;
    Ya1(Ybv)=LAB.BV; clear Ybv 
  end

  

  %% WMH (White Matter Hyperintensities):
  % WMHs can be found as GM next to the ventricle (A) that do not belong 
  % to a subcortical structure (A) or there must be a big difference 
  % between the tissue SPM expect and the real intensity 'Yp0 - Ym' (C).
  % Furthermore no other Sulic (=near other CSF) should be labeld (D).
  % ####################################################################
  % There can also be deep GM Hyperintensities! 
  % ####################################################################
  % ds('l2','',vx_vol,Ym,Ywmh,Ym/3,Ym/3,90)
  % ####################################################################
  % UPDATE: 
  % 1) add fast Shooting (6:1.5:1.5 mm resolution) for better initialization
  % 2) update values ... seperate detection of ventricular lession and 
  %    other 
  % ####################################################################
  try 
    Yp0e = Yp0.*cat_vol_morph(Yb,'e',2); 
    vols = mean([sum(round(Yp0e(:))==1) sum(round(Yp0e(:))==1 & Yvt(:))] / sum(round(Yp0e(:))>0.5));

    % only if there is a lot of CSF and not to much noise
    volth = 0.05; 
    if vols(1)>volth %&& noise<0.10 
      stime = cat_io_cmd(sprintf('  WMH detection (WMHCstr=%0.02f)',WMHCstr),'g5','',verb,stime); dispc=dispc+1;
      
      %WMHCstr = WMHCstr * (1 + (vols(1)-0.10)*2); 
      
      YBG2 = cat_vol_morph(Ya1==LAB.BG,'d',1); 
      
      % cortex 
      Yvto = cat_vol_morph(Yvt,'o',3/mean(vx_vol)); 
      Ywmh = single(smooth3(cat_vol_morph(Yvto,'d',2/mean(vx_vol)) & Ym<2.25 & cat_vol_morph(YA==LAB.CT,'e',2) &...
        ~(cat_vol_morph(YA==LAB.HC & Ym>1.5,'d',4*vxd) & Ym>1.5))>0.5); clear Yvto;  % ventricle
      Ywmh(smooth3((Ym.*Yp0A - Ym.*Ym)>2-WMHCstr/5+0.05-vols+noise & Ym<2.8 & Ym>2.8)>0.5 & cat_vol_morph(YA==LAB.CT,'e',2))=1; % WMH
      Ywmh(smooth3(~cat_vol_morph(~(Yp0>2.2 | Yvt),'lo',0) & Yp0<2.8 & Ym<2.8 & Yp0>1.5 & Ym>1.5)>0.6)=1; % WMH wholes
      Ywmh((Ywmh==0 & Ym>2.8) | YBG2 | YA==LAB.BG | Ya1==LAB.TH | YA==LAB.TH)=-inf;
      Ywmh(cat_vol_morph(Ya1==LAB.HC,'d',3) & Ywmh==0)=2; 
      Ywmh(Yvt2>1.75 & Yvt2<2.2 | (Ywmh==0 & Ym<1.5 & ~cat_vol_morph(Yvt,'d',2)))=2;

      % == dieser abschnitt ist noch in der entwicklung ==
     % Ywmh( smooth3((smooth3(~cat_vol_morph(Yp0<2.5,'lc',1))>0.5 | smooth3(Yp0 + Yp0A - Ym - Yg*2)>=(2.5+0.05-vols+noise)) & ... 
     %   Ym>2 & Ym<max(2.5,2.9-noise) & cat_vol_morph(YA==LAB.CT,'e',2))>0.5)=1;
      if mean(vx_vol)<1.5
        %%
        Ynwmh = cat_vol_morph(YBG2 | YA==LAB.BG | Ya1==LAB.BG | Ya1==LAB.TH | YA==LAB.TH | Ya1==LAB.VT | YA==LAB.HC,'c',8);
        Ynwmh = Ynwmh | cat_vol_smooth3X(Ynwmh,4)>0.4;  
        Ynwmh = cat_vol_morph(Ynwmh,'d',1) & cat_vol_morph(Ynwmh | Yp0>2.5,'c',1);  
        Yoc  = cat_vol_morph(((Ym<=2 & Yp0<=2 & ~cat_vol_morph(Yp0>2.5,'lc')) | Yp0<1) & ~cat_vol_morph(Ynwmh,'d',2),'l',0.5); 
        Ygmd = cat_vbdist(single(Yoc),Yp0>=1,vx_vol); % abstand zum CSF/GM bereich
        Ywmm = cat_vol_localstat(Ym,cat_vol_morph(Yp0>2.2,'lc'),2,3); % lokaler wm threshold 
        Ywmm = cat_vol_localstat(Ywmm,Ywmm>0,1,1); 
        %%
        Ywmhsm = cat_vol_smooth3X( ((Ywmm - Ym)>max(0.15,min(0.5,noise/2))) &  Ygmd>3  & cat_vol_morph(YA==LAB.CT,'e',2) & ...
          (cat_vol_morph(~cat_vol_morph(Yp0<2.5,'l'),'e') | cat_vol_morph(cat_vol_morph(~cat_vol_morph(Yp0<2.9,'e'),'l'),'e',2)) & ...
          cat_vol_morph(YA==LAB.CT,'e',2),0.6)>0.5;  % kleine wmhs
        Ywmhsm = Ywmhsm | ((Yp0A.*Ygmd/4) > Yp0 & ~Ynwmh & Yp0<2.5 & Yp0>1.5); 
        %%
       
        %%
        [YSr,red1]  = cat_vol_resize(YS,'reduceV',vx_vol,vx_vol*4,16,'meanm'); 
        YSdr   = cat_vbdist(single(YSr>0),true(size(YSr)),red1.vx_volr) + cat_vbdist(single(YSr<1),true(size(YSr)),red1.vx_volr); 
        YSD    = cat_vol_resize(YSdr,'dereduceV',red1);
         Ywmhsm = cat_vol_morph(Ywmhsm & YSD>15 & ~Ynwmh,'l',[inf -20/mean(vx_vol)]); 
        %%
        %Ywmh(Ywmhsm>0)=1;
        Ywmh(Ynwmh & Ywmh==1)=0; 
        Ywmh(Ya1==LAB.VT) = 2; 
      end
      %  
     % Ywmh(Ywmh==1 & smooth3(Ywmh==1)<0.65 - vols(1) - noise)=0; % more WMHs if lot of CSF and less noise 
      %Ywmh((Ywmh==0 & Ym>2.5) | YBG2 | Ya1==LAB.TH)=-inf;
      %Ywmh = cat_vol_downcut(Ywmh,(3-Ym)/3,noise*WMHCstr/2,vx_vol); 
      %Ywmh = Ywmh .* (cat_vol_morph(Ywmh,'l',[inf -30/mean(vx_vol)])>0); 
      %%
      Ywmh(Ywmh==2 & smooth3(Ywmh==2)<0.1+WMHCstr/2)=0;
      Ywmh(Ywmh==1 & smooth3(Ywmh==1)<0.1+WMHCstr/2)=0;
      Ywmh(Ywmh<0 & (YA==LAB.CT | Ya1==LAB.CT) & Ym<2.95)=0;
      Ywmh = cat_vol_downcut(Ywmh,(3-Ym)/3,noise*WMHCstr/2,vx_vol); 

      %%
      %{
      Ywmh(isinf(Ywmh))=0;
      Ywmh((Ywmh==0 & Ym>2.75) | YBG2 | Ya1==LAB.TH)=nan; Ywmh(Ywmh==0)=1.5;
      Ywmh = cat_vol_laplace3R(Ywmh,Ywmh==1.5,0.005);
      Ywmh(cat_vol_morph(YS==1,'d',3*vxd) & cat_vol_morph(YS==0,'d',3*vxd))=2; % not for the CC
      Ynwmh = ~smooth3(cat_vol_morph(Ya1==LAB.VT | Ya1==LAB.TH | YBG2,'c',4))>0.5;
      Ywmh = smooth3(Ywmh<1.1+WMHCstr/4 & Ya1~=LAB.VT & Ynwmh)>0.75-WMHCstr/2;
      %%
      Ywmh = smooth3(cat_vol_morph(Ywmh,'c',1) & Ym<2.95)>(0.75-WMHCstr/2);
      %Ywmh(Ywmhsm>0)=1; 
      %}
      Ywmh = Ywmh==1; 
      %% entfernen zu kleiner wmhs
      [Ywmhl,num] = spm_bwlabel(double(Ywmh>0));
      lhst = hist(Ywmhl(:),1:num); lhstind = 1:num; lhstind(lhst>9)=[];
      for lhsti=1:numel(lhstind), Ywmh(Ywmhl==lhstind(lhsti))=0; end
      %%
      Ya1(Ywmh)=LAB.HI;
    elseif vols(1)<volth
      stime = cat_io_cmd(sprintf('  NO WMH detection (CSF ~%0.0f%%%%)',vols(1)*100),'g5','',verb,stime); dispc=dispc+1;
    elseif noise>0.10 
      stime = cat_io_cmd(sprintf('  NO WMH detection (too noisy ~%0.2f)',noise),'g5','',verb,stime); dispc=dispc+1;
    else
      stime = cat_io_cmd(sprintf('  NO WMH detection (CSF ~%0.0f%%%% and too noisy ~%0.2f)',...
        vols(1)*100,noise),'g5','',verb,stime); dispc=dispc+1;
    end
  end
  %{
   Yvt2(Yvt2>1.45 & Yvt2<1.55)=inf; Yvt2=round(Yvt2);
  Yvt2(Yvt2>3 & Ym<2.5 & Ym>1.5 & YA~=LAB.CB & YA~=LAB.BS & Ya1~=LAB.BG & Ya1~=LAB.TH)=1.5;
  Yvt2(Yvt2==1.5 & smooth3((Yp0 - Ym)>0.6 & Ym<2.75 & Ym>1.75)>0.7)=1; 
  Yvt2(Yvt2>3)=1.5;
  Yvt2 = cat_vol_laplace3R(Yvt2,Yvt2==1.5 & Ym<2.75 & YA~=LAB.BS & Ya1~=LAB.BG & Ya1~=LAB.TH,0.005); %Yvt2 = round(Yvt2);
  %}
  if ~debug, clear Ywmh Yvt Ynwmh Yvt2; end

  
  
  %% Closing of gaps between diffent structures:
  stime = cat_io_cmd('  Closing of deep structures','g5','',verb,stime); dispc=dispc+1;
  Yvtd2 = cat_vol_morph(Ya1==LAB.VT,'d',2*vxd) & Ya1~=LAB.VT;
  % CT and VT
  Yt = cat_vol_morph(Ya1==LAB.VT,'d',2*vxd) & ...
       cat_vol_morph(Ya1==LAB.CT,'d',2*vxd) & Ya1==0 ;
  Ya1(Yt & Yp0<=1.5 & ~Ynv)=LAB.VT; Ya1(Yt & Yp0>1.5)=LAB.CT; 
  % WMH and VT
  Yt = cat_vol_morph(Ya1==LAB.HI,'d',1*vxd) & Yvtd2 & ~Ynv & Ya1==0;
  Ya1(Yt &  Ym<=1.25)=LAB.VT; Ya1(Yt & Ym>1.25 & Ym<2.5)=LAB.HI; 
  % TH and VT
  Yt = cat_vol_morph(Ya1==LAB.TH,'d',1*vxd) & Yvtd2;
  Ya1(Yt & Ym<=1.5)=LAB.VT; Ya1(Yt & Ym>1.5 & Ym<2.85)=LAB.TH; 
  % BG and VT
  Yt = cat_vol_morph(Ya1==LAB.BG,'d',1*vxd) & Yvtd2;
  Ya1(Yt & Ym<=1.5)=LAB.VT; Ya1(Yt & Ym>1.5 & Ym<2.85)=LAB.BG;
  % no bloodvessels next to the ventricle, because for strong atrophy
  % brains the WM structures can be very thin and may still include 
  % strong bias
  Ya1(Ya1==LAB.BV & cat_vol_morph(Ya1==LAB.VT,'d',3*vxd))=0;
  clear Yt Yh Yvtd2 Yw
 
  
  
  %% complete map
  [tmp0,tmp1,Ya1] = cat_vbdist(Ya1,Yb); clear tmp0 tmp1;
  
  % consider gyrus parahippocampalis
  Ya1(YA==LAB.PH) = LAB.PH;
  
  %% side aligment using laplace to correct for missalignments due to the normalization
  stime = cat_io_cmd('  Side alignment','g5','',verb,stime); dispc=dispc+1;
  YBG  = Ya1==LAB.BG | Ya1==LAB.TH;
  YMF  = Ya1==LAB.VT | Ya1==LAB.BG | Ya1==LAB.TH | Ya1==LAB.HI; 
  YMF2 = cat_vol_morph(YMF,'d',2*vxd) | Ya1==LAB.CB | Ya1==LAB.BS | Ya1==LAB.MB;
  Ymf  = max(Ym,smooth3(single(YMF2*3))); 
  Yt = cat_vol_smooth3X(YS==0,6)<0.9 & cat_vol_smooth3X(YS==1,6)<0.9 & ~YMF2 & Yp0>0 & Ym<3.1 & (Yp0<2.5 | Ya1==LAB.BV);
  Ys = (2-single(YS)) .* single(smooth3(Yt)<0.4);
  Ys(Ys==0 & (Ym<1 | Ym>3.1))=nan; Ys = cat_vol_downcut(Ys,Ymf,0.1,vx_vol); 
  [tmp0,tmp1,Ys] = cat_vbdist(Ys,Ys==0);
  clear YMF2 Yt YS tmp0 tmp1;
  
  %% YMF for FreeSurfer fsaverage
  Ysm  = cat_vol_morph(Ys==2,'d',1.75*vxd) & cat_vol_morph(Ys==1,'d',1.75*vxd);
  YMF  = cat_vol_morph(Ya1==LAB.VT | Ya1==LAB.BG | Ya1==LAB.HI | (Ya1==LAB.TH & smooth3(Yp0)>2),'c',3) & ~Ysm; 
  %YMF  = YMF | (cat_vol_morph(YA==LAB.CT & YBG,'c',6) & ~Ysm); 
  YMF  = Ym<=2.5  & cat_vol_morph(YMF | Ym>2.3,'c',1) & cat_vol_morph(YMF,'d',2);
  YMF  = smooth3(YMF)>0.5;
  clear Ysm; 
  
  
  %% back to original size
  stime = cat_io_cmd('  Final corrections','g5','',verb,stime); dispc=dispc+1;
  Ya1 = cat_vol_resize(Ya1,'dereduceV',resTr,'nearest'); Ya1 = cat_vol_median3c(Ya1,Ya1>1 & Ya1~=LAB.BV);
  Ys  = cat_vol_resize(Ys ,'dereduceV',resTr,'nearest'); Ys  = 1 + single(smooth3(Ys)>1.5);
  YMF = cat_vol_resize(YMF,'dereduceV',resTr);
  YBG = cat_vol_resize(YBG,'dereduceV',resTr);
  
  Ya1 = cat_vol_resize(Ya1,'dereduceBrain',BB); Ya1 = cat_vol_ctype(Ya1);
  Ys  = cat_vol_resize(Ys ,'dereduceBrain',BB); [tmp0,tmp1,Ys] = cat_vbdist(Ys,Ya1>0); clear tmp0 tmp1;
  YMF = cat_vol_resize(YMF,'dereduceBrain',BB); 
  YBG = cat_vol_resize(YBG,'dereduceBrain',BB); 
  Ym  = Ym0; clear Ym0;

  % final side alignment
  Ya1(Ya1>0)=Ya1(Ya1>0)+(Ys(Ya1>0)-1);
 
  
  % class correction
  % YBG is smoothed a little bit and (B) reset all values that are related
  % with GM/WM intensity (Ym<2.9/3) (A)
  Yclssum = single(Ycls{1})+single(Ycls{2})+single(Ycls{3});
  YBGs    = min( max(0,min(255, 255 - cat_vol_smooth3X(Ya1==1 & Ycls{2}>round(2.9/3),0.8) .* single(Ycls{2}) )), ... (A)
                 max(0,min(255, 255 * cat_vol_smooth3X(YBG .* (Ym<=2.9/3 & Ym>2/3) ,0.5) )) ); % (B)
  Ycls{1} = cat_vol_ctype(single(Ycls{1}) + YBGs .* (single(Ycls{2})./max(eps,Yclssum)));
  Ycls{2} = cat_vol_ctype(single(Ycls{2}) - YBGs .* (single(Ycls{2})./max(eps,Yclssum)));
  clear YBGs Yclssum; 
 
  if debug
    cat_io_cmd(' ','','',verb,stime); 
  else
    cat_io_cmd(' ','','',verb,stime); 
    %cat_io_cmd('cleanup',dispc,'',verb); 
  end
  
end
%=======================================================================
function Yg = cat_vol_grad(Ym,vx_vol)
% ----------------------------------------------------------------------
% gradient map for edge description
% ----------------------------------------------------------------------
  [gx,gy,gz] = cat_vol_gradient3(Ym); 
  Yg = abs(gx./vx_vol(1))+abs(gy./vx_vol(2))+abs(gz./vx_vol(3)); 
  Yg = Yg ./ max(eps,Ym);
end
%=======================================================================

%=======================================================================
function Ydiv = cat_vol_div(Ym,vx_vol)
% ----------------------------------------------------------------------
% Diverence helps to identify all gyri that should not be in the GM, but 
% helps to improve the WM. Divergence estimation is very memory intensive 
% so it is better to limit the resolution.
% ----------------------------------------------------------------------
  [Ymr,resT2] = cat_vol_resize(Ym,'reduceV',vx_vol,1.5,32);
  [gx,gy,gz]  = cat_vol_gradient3(max(2/3,Ymr)); 
  Ydivr = smooth3(divergence(gy./vx_vol(1),gx./vx_vol(1),gz./vx_vol(3))); clear gx gy gz Ymr;
  Ydiv  = cat_vol_resize(Ydivr,'dereduceV',resT2); 
end
%=======================================================================