function [Ym,Yt,Ybg,WMth,bias] = cat_run_job_APP_init(Ysrco,vx_vol,verb)
%  _____________________________________________________________________
%  The rough bias correction is a subfunction of cat_run_rob.
% 
%  All tissues (low gradient areas) should have a similar intensity.
%  A strong smoothing of this approximation is essential to 
%  avoid anatomical filtering between WM and GM that can first 
%  be seen in overfitting of the subcortical structures.
%  However, this filtering will overcorrect head tissue with
%  a typical intensity around GM.
%  _____________________________________________________________________
%  Robert Dahnke
%  $Id: cat_run_job_APP_init1070.m 1147 2017-06-21 09:31:53Z gaser $


%    ds('l2','',0.5,Yo/WMth,Yg<0.2,Yo/WMth,Ym,80)

  rf = 10^9; 
  bfsmoothness = 3; 
  if verb, fprintf('\n'); end
  
  stime = cat_io_cmd('  Initialize','g5','',verb);
  msize = 222; %round(222 ./ max(size(Ysrco).*vx_vol) .* min(size(Ysrco).*vx_vol));  

  [Ysrc,resT3] = cat_vol_resize(Ysrco,'reduceV',vx_vol,min(1.2,cat_stat_nanmean(vx_vol)*2),msize,'meanm'); 

  % correction for negative backgrounds (MT weighting)
  WMth = roundx(single(cat_stat_nanmedian(Ysrc(Ysrc(:)>cat_stat_nanmean( ...
          Ysrc(Ysrc(:)>cat_stat_nanmean(Ysrc(:))))))),rf); 
  BGth = max( min(Ysrc(:))*0.7 + 0.3*WMth ,...
    cat_stat_nanmean(Ysrc(Ysrc(:)<cat_stat_nanmean(Ysrc(:))))); BGth = roundx(BGth,rf); 

  Ysrc = Ysrc - BGth; Ysrco = Ysrco - BGth; BGth2 = BGth; 
  Yg   = cat_vol_grad(Ysrc,resT3.vx_volr) ./ max(eps,Ysrc); 
  
  WMth = roundx(single(cat_stat_nanmedian(Ysrc(Yg(:)<0.2 & Ysrc(:)>cat_stat_nanmean( ...
           Ysrc(Yg(:)<0.2 & Ysrc(:)>cat_stat_nanmean(Ysrc(:))))))),rf); 
  BGth = max( min(Ysrc(:))*0.7 + 0.3*cat_stat_nanmean(Ysrc(:)) ,...
    cat_stat_nanmean(Ysrc(Ysrc(:)<cat_stat_nanmean(Ysrc(:))))); BGth = roundx(BGth,rf); 
  Ym   = (Ysrc - BGth) ./ (WMth - BGth);
  
  Ydiv = cat_vol_div(Ym,resT3.vx_volr/2) ./ (Ym+eps); % lower resolution is 8 times faster 
  
  
  %% background
  stime = cat_io_cmd('  Estimate background','g5','',verb,stime);
  Ybg = ((Yg.*Ym)<cat_vol_smooth3X(Ym,2)*1.2) & Ym>0.1; 
  Ybg([1,end],:,:)=0; Ybg(:,[1,end],:)=0; Ybg(:,:,[1,end])=0; Ybg(smooth3(Ybg)<0.5)=0;
  [Ybg,resT2] = cat_vol_resize(single(Ybg),'reduceV',resT3.vx_volr,2,32,'meanm'); 
  Ybg([1,end],:,:)=0; Ybg(:,[1,end],:)=0; Ybg(:,:,[1,end])=0; Ybg = Ybg>0.5;
  Ybg  = cat_vol_morph(Ybg,'lc',8);
  Ybg  = cat_vol_smooth3X(Ybg,2); 
  Ybg  = cat_vol_resize(Ybg,'dereduceV',resT2)<0.5;    
  BGth = roundx(mean(Ysrc(Ybg(:))),rf);
  Ym   = (Ysrc - BGth) ./ (WMth - BGth);
  
  %% first WM inhomogeneity with low tissue boundary (may include CSF > strong filtering for IXI175)
  stime = cat_io_cmd('  Initial correction','g5','',verb,stime);
  Yms  = cat_vol_smooth3X( min(2 .* ~Ybg,Ym .* (Ydiv>-0.2) .* ~Ybg .* (Ym>0.1)),16*mean(vx_vol));     % this map is to avoid CSF in the mask!
  Yms  = (Yms ./ mean(Yms(~Ybg(:)))) * WMth;
  Yms  = cat_vol_smooth3X( min(Yms*1.5 .* ~Ybg,Ysrc .* ~Ybg),16*mean(vx_vol));
  Yms  = (Yms ./ mean(Yms(~Ybg(:)))) * WMth;
  Yt   = Ysrc>max(BGth,Yms*0.3) & Ysrc<Yms*2 & Ysrc<WMth*(1+Yms/WMth*2) & Yg<0.9 & Ydiv<0.2 & ...
         Ydiv>-0.6 & smooth3(Ysrc./(Yms+eps).*Yg.*Ydiv<-0.2)<0.3 & ~Ybg; Yt(smooth3(Yt)<0.5)=0;
  Ywi  = (Ysrc .* Yt) ./ max(eps,Yt);  
  [Ywi,resT2] = cat_vol_resize(Ywi,'reduceV',resT3.vx_volr,cat_stat_nanmean(resT3.vx_volr)*2,32,'max'); 
  for i=1:1, Ywi = cat_vol_localstat(Ywi,Ywi>0,2,3); end % only one iteration!
  for i=1:4, Ywi = cat_vol_localstat(Ywi,Ywi>0,2,1); end
  Ywi  = cat_vol_approx(Ywi,'nn',resT2.vx_volr,4);
  Ywi  = cat_vol_smooth3X(Ywi,bfsmoothness.*mean(vx_vol)); % highres data have may stronger inhomogeneities 
  Ywi  = cat_vol_resize(Ywi,'dereduceV',resT2);    
  Ybc  = Ysrc./Ywi;
  WMt2 = roundx(cat_stat_nanmedian(Ybc(Yg(:)<0.2 & Ybc(:)>0.9)),rf); 
  Ywi  = Ywi * WMt2;
  
  %% background update
  stime = cat_io_cmd('  Refine background','g5','',verb,stime);
  Ybg = ((Yg.*(Ysrc./Ywi))<cat_vol_smooth3X(Ysrc./Ywi,2)*1.2) & Ysrc./Ywi>0.2; 
  Ybg([1,end],:,:)=0; Ybg(:,[1,end],:)=0; Ybg(:,:,[1,end])=0; Ybg(smooth3(Ybg)<0.5)=0;
  [Ybg,resT2] = cat_vol_resize(single(Ybg),'reduceV',resT3.vx_volr,2,32,'meanm'); 
  Ybg([1,end],:,:)=0; Ybg(:,[1,end],:)=0; Ybg(:,:,[1,end])=0; Ybg = Ybg>0.5;
  Ybg  = cat_vol_morph(Ybg,'lc',8);
  Ybg  = cat_vol_smooth3X(Ybg,2); 
  Ybg  = cat_vol_resize(Ybg,'dereduceV',resT2)<0.5 & Ysrc<min(WMth*0.2,BGth*0.8+0.2*WMth);
  Ybg  = cat_vol_morph(Ybg,'lo');

  %% second WM inhomogeneity with improved Yt with higher lower threshold (avoid CSF and less filtering)
  stime = cat_io_cmd('  Final correction','g5','',verb,stime);
  Yt   = Ysrc>max(BGth,Yms*0.3)  & Ysrc./(Ywi+eps)>0.2 & Ysrc./(Ywi+eps)<1.2 & Ysrc./(Ywi+eps)<Yms/WMth*2 & Yg<0.9 & Ydiv<0.2 & Ydiv>-0.6 & ...
         smooth3(Ysrc./(Yms+eps).*Yg.*Ydiv<-0.1)<0.1 & ~Ybg; Yt(smooth3(Yt)<0.5)=0;
  Yt   = Yt | (~Ybg & Ysrc>BGth/2 & Ysrc>Yms*0.5 & Ysrc<Yms*1.2 & Ydiv./(Yg+eps)<0.5 & ((Ysrc./(Ywi+eps)>0.3 & Yg>0.1 & Ydiv<0) | (~Ybg & Ysrc./(Ywi+eps)>0.6)) & Ysrc./(Ywi+eps)<1.2); 
  Yt(smooth3(Yt)<0.7)=0;
  Ywi  = (Ysrc .* Yt) ./ max(eps,Yt);  
  [Ywi,resT2] = cat_vol_resize(Ywi,'reduceV',resT3.vx_volr,cat_stat_nanmean(resT3.vx_volr)*2,32,'max'); 
  for i=1:1, Ywi = cat_vol_localstat(Ywi,Ywi>0,2,3); end % only one iteration!
  for i=1:4, Ywi = cat_vol_localstat(Ywi,Ywi>0,2,1); end
  Ywi  = cat_vol_approx(Ywi,'nn',resT2.vx_volr,4);
  Ywi  = cat_vol_smooth3X(Ywi,bfsmoothness.*mean(vx_vol)); %.*mean(vx_vol)); % highres data have may stronger inhomogeneities 
  Ywi  = cat_vol_resize(Ywi,'dereduceV',resT2);    
  Ybc  = Ysrc./Ywi;
  WMt2 = roundx(cat_stat_nanmedian(Ybc(Yg(:)<0.2 & Ybc(:)>0.9)),rf); 
  Ywi  = Ywi * WMt2;
  bias = std(Ywi(:))/mean(Ywi(:)); 
  
  %% BG inhomogeneity (important for normalization of the background noise)
  %[Ybc,Ygr,resT2] = cat_vol_resize({Ysrc./Ywi,Yg},'reduceV',resT3.vx_volr,cat_stat_nanmean(resT3.vx_volr)*4,16,'meanm'); 
  %Ybc  = cat_vol_morph(Ybc<BGth/WMth*2 & Ygr<0.05,'lc',2);
  %Ybc  = cat_vol_resize(smooth3(Ybc),'dereduceV',resT2)>0.5; 
  stime = cat_io_cmd('  Background correction','g5','',verb,stime);
  [Ybc,resT2] = cat_vol_resize(single(Ysrc .* Ybg),'reduceV',resT3.vx_volr,max(8,min(16,cat_stat_nanmean(resT3.vx_volr)*4)),16,'min'); 
  Ybc  = cat_vol_localstat(Ybc,Ybc>0,2,2);
  Ybc  = cat_vol_localstat(Ybc,Ybc>0,2,1);
  %Ybc  = cat_vol_approx(Ybc,'nn',resT2.vx_volr,4); % no aproximation to correct only in the background! 
  Ybc  = cat_vol_smooth3X(Ybc,4);
  Ybc  = cat_vol_resize(Ybc,'dereduceV',resT2); 


  %% back to original size
  stime = cat_io_cmd('  Final scaling','g5','',verb,stime);
  [Ywi,Ybc] = cat_vol_resize({Ywi,Ybc},'dereduceV',resT3); 
  Yg        = cat_vol_resize(Yg,'dereduceV',resT3); 
  [Yt,Ybg]  = cat_vol_resize({single(Yt),single(Ybg)},'dereduceV',resT3); Yt = Yt>0.5; Ybg = Ybg>0.5;
  Ysrc      = Ysrco; clear Ysrco;

  %% intensity normalization (Ybc is the average background noise)
  % in data with strong inhomogeneities (7T) the signal can trop below the noise level 
  Ym   = (Ysrc - min(BGth,min(Ybc/2,Ywi/20))) ./ (Ywi - min(BGth,min(Ybc/2,Ywi/20))); 
  Wth  = single(cat_stat_nanmedian(Ym(Yg(:)<0.2 & Ym(:)>cat_stat_nanmean( Ym(Yg(:)<0.2 & Ym(:)>cat_stat_nanmean(Ym(:))))))); 
  [WIth,WMv] = hist(Ym(Yg(:)<0.2 & Ym(:)>Wth*0.5 & Ym(:)<Wth*1.5),0:0.01:2);
  WIth = find(cumsum(WIth)/sum(WIth)>0.8,1,'first'); WIth = roundx(WMv(WIth),rf); 
  Ym   = Ym ./ WIth; 
  % update WMth
  Ysrc = Ysrc + BGth2;
  [WIth,WMv] = hist(Ysrc(Yg(:)<0.2 & Ym(:)>Wth*0.5 & Ym(:)<Wth*1.5),1000);
  WMth = find(cumsum(WIth)/sum(WIth)>0.7,1,'first'); WMth = roundx(WMv(WMth),rf); 
  
  cat_io_cmd(' ','','',verb,stime); 
end
%=======================================================================
function r = roundx(r,rf)
  r(:) = round(r(:) * rf) / rf;
end
%=======================================================================
