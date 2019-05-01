function [Ycls,Yp0b] = cat_main_cleanup(Ycls,prob,Yl1b,Ymb,extopts,inv_weighting,vx_vol,indx,indy,indz)
%  -----------------------------------------------------------------
%  final cleanup 2.0
%  
%  First we need to describe our region of interest. Blood vessels and 
%  menignes occure in the sulci and next to the skull. Therefore we 
%  use the brainmask and the label map to identify special regions.
%  -----------------------------------------------------------------
  dbs   = dbstatus; debug = 0; 
  for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,'cat_main_cleanup'); debug = 1; break; end; end

  
  LAB  = extopts.LAB;
  vxv  = 1/max(vx_vol);           % use original voxel size!!!
  NS   = @(Ys,s) Ys==s | Ys==s+1; % remove side alignment from atlas maps

  cleanupstr  = min(1,max(0,extopts.cleanupstr * 1/(inv_weighting+1) / max(1,mean(vx_vol)) ));
  cleanupdist = min(2,max(0,1 + 2*extopts.cleanupstr));

  stimec = cat_io_cmd(sprintf('Final cleanup (gcutstr=%0.2f)',cleanupstr));
  fprintf('\n');
  stime = cat_io_cmd('  Level 1 cleanup (ROI estimation)','g5','',extopts.verb); %dispc=1;
  
  %% estimate the ROI
  % ------------------------------------------------------------------
  % This part removes menignes next to the skull and between large 
  % structes.
  % ------------------------------------------------------------------
  Yvt  = cat_vol_morph(NS(Yl1b,LAB.VT) | NS(Yl1b,LAB.BG),'d',vxv*3);  % ventricle ... no cleanup here
  Yp0  = single(prob(:,:,:,1))/255*2 + single(prob(:,:,:,2))/255*3 + single(prob(:,:,:,3))/255;
  Ybd  = cat_vbdist(single(~cat_vol_morph(Yp0>0,'lc',vxv)),true(size(Yp0)),vx_vol);
  Ybd  = cat_vbdist(single(~cat_vol_morph(Yp0>1.5 | Ybd>8,'lc',vxv)),true(size(Yp0)),vx_vol);
  Ycbp = cat_vol_morph(NS(Yl1b,LAB.CB),'d',cleanupdist*vxv);          % next to the cerebellum
  Ycbn = cat_vol_morph(NS(Yl1b,LAB.CB),'e',0.2*cleanupdist*vxv);      % not to deep in the cerebellum
  Ylhp = cat_vol_morph(Yl1b==1 & Yp0<2.1,'d',cleanupdist*vxv*2);      % GM next to the left hemisphere 
  Yrhp = cat_vol_morph(Yl1b==2 & Yp0<2.1,'d',cleanupdist*vxv*2);      % GM next to the righ hemishpere
  Yroi = Ybd<cleanupdist*2 | ...                                      % next to the brain mask
         (~Ycbn & Ycbp & (Ylhp | Yrhp)) | ...                         % between the cortex and the cerebellum                       
         (Ylhp & Yrhp) | ...                                          % between left and right hemisphere
         NS(Yl1b,LAB.VT) | ...                                        % in the ventricle 
         (NS(Yl1b,LAB.BS) & Ycbp);                                    % between brainstem and crebellum
  Yrbv = Yp0>0 & Ybd<6 & cat_vol_morph( (Ylhp & Yrhp) | (~Ycbn & Ycbp & (Ylhp | Yrhp)),'d',4);
  Yroi = (Yroi | Yrbv) & ~NS(Yl1b,LAB.BS) & ~Ycbn; 
  % bv
%     Ycd  = cat_vbdist(single(Yp0>2.5 & Yl1b==LAB.CB),Ylhp & Yrhp,vx_vol);
%     Ylhd = cat_vbdist(single(Yp0>2.5 & Yl1b==1),Ylhp & Yrhp,vx_vol);
%     Yrhd = cat_vbdist(single(Yp0>2.5 & Yl1b==1),Ylhp & Yrhp,vx_vol);
%     Ybvx = single(min(cat(4,Ylhd,Yrhd,Ycd),[],4)<8 & (min(cat(4,Ylhd,Yrhd,Ycd),[],4))>3 & Ybd<5 & Ymb>0.3); 
%     Ywm  = single(smooth3(Yp0>2)>0.5); Ybvx(Ymb<0.67 | Yp0==0)=nan;
%     Ywm  = cat_vol_downcut(Ywm,Ymb,0.1);
%     Ybvx(smooth3(Ybvx)<0.7)=0; Ybvx(smooth3(Ybvx)<0.5)=0; Ybvx(Ywm & Ybvx==0)=2; Ybvx(Yp0==0)=nan;
%     Ybvx = cat_vol_downcut(Ybvx,Ymb,0.05);

  if ~debug, clear Ycbp Ycbn Ylhp; end

  %% roi to change GM or WM to CSF or background
  stime = cat_io_cmd('  Level 1 cleanup (brain masking)','g5','',extopts.verb,stime); %dispc=dispc+1;
  Yrw = Yp0>0 & Yroi & Ymb>1.1+Ybd/20 & ~NS(Yl1b,LAB.CB);             % basic region with cerebellum
  Yrw = Yrw | smooth3(Yrw)>0.4-0.3*cleanupstr;                        % dilate region
  Ygw = cat_vol_morph(Yp0>=1.9 & ~Yrw,'lo',0); % even one is to much in adrophic brains :/ 
  Yrw = Yrw | (Yp0>1 & Yroi & ~Ygw);                                  % further dilation
  Yrw = Yrw & ~Yvt & ~cat_vol_morph(Ygw,'d',1); 
  Yrw(smooth3(Yrw)<0.5+0.2*cleanupstr)=0; 
  Yrw(smooth3(Yrw)<0.5-0.2*cleanupstr)=0;                             % only larger objects
  if ~debug, clear Ygw Yroi; end

  %% update brain masks and class maps
  Ybb = cat_vol_morph((Yp0>0 & ~Yrw) | Ybd>2,'lo',2/vxv);          
  Ybb(cat_vol_smooth3X(Ybb,2)>0.4 & ~Yrw)=1;
  Ybb = cat_vol_morph(Ybb | Ybd>3,'lc',1/vxv); 
  Ybb = single(Ybb); spm_smooth(Ybb,Ybb,0.6./vx_vol); Ybb = Ybb>1/3;

  %% correct to background
  for i=1:3, prob(:,:,:,i)=min(prob(:,:,:,i),uint8(Ybb*255)); end
  % correct to CSF
  prob(:,:,:,1)=min(prob(:,:,:,1),uint8(~(Ybb & Yrw)*255));
  prob(:,:,:,2)=min(prob(:,:,:,2),uint8(~(Ybb & Yrw)*255));
  prob(:,:,:,3)=max(prob(:,:,:,3),uint8( (Ybb & Yrw)*255));    
  if ~debug, clear Yrw; end



  %% cleanup of meninges
  % ------------------------------------------------------------------
  % This removes meninges next to the brain... works quite well.
  clear Yrg Yrw Yroi
  stime = cat_io_cmd('  Level 2 cleanup (CSF correction)','g5','',extopts.verb,stime); %dispc=dispc+1;
  Yp0 = single(prob(:,:,:,1))/255*2 + single(prob(:,:,:,2))/255*3 + single(prob(:,:,:,3))/255;
  YM  = single(cat_vol_morph((prob(:,:,:,1) + prob(:,:,:,2))>(160 + 32*cleanupstr) & ...
         ~cat_vol_morph(Yp0>1 & Yp0<1.5+cleanupstr/2,'o',vxv)  ,'l')); 
  YM2 = cat_vol_morph(YM,'o',min(1,0.7/max(vx_vol)));
  YM(NS(Yl1b,1) & YM2==0)=0;
  spm_smooth(YM,YM,0.6./vx_vol); % anisotropic smoothing!
  YM  = ( (YM<0.1*cleanupstr) ) & Ybb & ~Yvt & Ymb>0.25;
  prob(:,:,:,1)=min(prob(:,:,:,1),uint8(~YM*255));
  prob(:,:,:,2)=min(prob(:,:,:,2),uint8(~YM*255));
  prob(:,:,:,3)=max(prob(:,:,:,3),uint8( (YM | (Ybb & Yp0==0))*255));
  Yp0  = single(prob(:,:,:,1))/255*2 + single(prob(:,:,:,2))/255*3 + single(prob(:,:,:,3))/255;
  


  %% cleanup WM 
  % ------------------------------------------------------------------
  % the idea was to close WMH ... but its not stable enough yet
  %{
  Ywmh = false(size(Yp0)); 
  for p0thi=2.1:0.2:2.9
    Ywmh = Ywmh | ~cat_vol_morph(Yp0<p0thi,'l') & (Yp0<p0thi); 
  end
  Ywmh = smooth3(Ywmh)>0.1 & NS(Yl1b,1) & Yp0>=2 & Yp0<3; 
  Yl1b(Ywmh) = LAB.HI + ~mod(Yl1b(Ywmh),2); 
  clear Ywmh;
  %}
  % correction later depending on WMHC



%%
  % ------------------------------------------------------------------
  % cleanup in regions with PVE between WM and CSF without GM
  % ------------------------------------------------------------------
  stime = cat_io_cmd('  Level 3 cleanup (CSF/WM PVE)','g5','',extopts.verb,stime); %dispc=dispc+1;
  Yp0  = single(prob(:,:,:,1))/255*2 + single(prob(:,:,:,2))/255*3 + single(prob(:,:,:,3))/255;
  Ybs  = NS(Yl1b,LAB.BS) & Ymb>2/3;
  YpveVB = cat_vol_morph(NS(Yl1b,LAB.VT) | Ybs,'d',2);                % ventricle and brainstem
  YpveCC = cat_vol_morph(Yl1b==1,'d',3*vxv) & cat_vol_morph(Yl1b==2,'d',3*vxv) & ...
           cat_vol_morph(NS(Yl1b,LAB.VT),'d',2);                      % corpus callosum
  Ynpve  = smooth3(NS(Yl1b,LAB.BG) | NS(Yl1b,LAB.TH))>0.3;            % no subcortical structure 
  Yroi = (YpveVB | YpveCC) & ~Ynpve & ...
         cat_vol_morph(Yp0==3,'d',2) & cat_vol_morph(Yp0==1,'d',2) & ...
         Yp0<3 & Yp0>1 & ...
         smooth3((Yp0<3 & Yp0>1) & ~cat_vol_morph(Yp0<3 & Yp0>1,'o',1))>0.1;
  clear YpveVB YpveCC Ybs Ynpve;         
  Yncm = (3-Yp0)/2.*Yroi; 

  for i=1:3, Ycls{i}=zeros(size(Ycls{i}),'uint8'); end
  Ycls{1}(indx,indy,indz) = min(prob(:,:,:,1),uint8(~Yroi*255));
  Ycls{2}(indx,indy,indz) = cat_vol_ctype(single(prob(:,:,:,2)).*~Yroi + (Yroi - Yncm)*255,'uint8');
  Ycls{3}(indx,indy,indz) = cat_vol_ctype(single(prob(:,:,:,3)).*~Yroi + Yncm*255,'uint8');

  Yp0b = cat_vol_ctype(single(Ycls{1})*2/3 + single(Ycls{2}) + single(Ycls{3})*1/3,'uint8');
  Yp0b = Yp0b(indx,indy,indz); 

  cat_io_cmd(' ','','',extopts.verb,stime);   
  fprintf('%5.0fs\n',etime(clock,stimec));
