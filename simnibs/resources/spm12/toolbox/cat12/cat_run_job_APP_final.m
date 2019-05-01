function  [Ym,Yp0,Yb] = cat_run_job_APP_final(Ysrco,Ym,Yb,Ybg,vx_vol,gcutstr,verb)
%  _____________________________________________________________________
%  The final bias correction is a subfunction of cat_run_job.
% 
%  The affine registration, especially spm_preproc8 requires a very good 
%  masking! Because this is also required for the Unified Segmentation
%  a wider mask with a complete brain is important.
%  _____________________________________________________________________
%  Robert Dahnke
%  $Id: cat_run_job_APP_final.m 1184 2017-09-08 17:02:40Z dahnke $


  % if there is a breakpoint in this file set debug=1 and do not clear temporary variables 
  dbs   = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,'cat_run_job_APP_final'); debug = 1; break; end; end
  zeroBG = (cat_stat_nanmean(Ysrco(Ybg(:)>0))/cat_stat_nanmean(Ysrco(Yb(:)>0 & Ym(:)>0.9)))<0.4;

% ds('l2','m',0.5,Ym*0.7+0.3,Yb,Ysrc/WMth,Ym,80)

  if verb, fprintf('\n'); end
  stime = cat_io_cmd('  Initialize','g5','',verb);
  msize = 222; %round(222 ./ max(size(Ysrco).*vx_vol) .* min(size(Ysrco).*vx_vol));  
  
  [Ysrc,Ym,resT3] = cat_vol_resize({Ysrco,Ym},'reduceV',vx_vol,min(1.5,min(vx_vol)*2),msize,'meanm'); 
  [Yb,Ybg]        = cat_vol_resize({single(Yb),single(Ybg)},'reduceV',vx_vol,min(1.5,min(vx_vol)*2),msize,'meanm'); 
  if debug, Ybo = Yb; end %#ok<NASGU>
  Ybg = Ybg>0.5;
  
  Yg   = cat_vol_grad(Ym,resT3.vx_volr) ./ max(eps,Ym); 
  Ydiv = cat_vol_div(Ym,resT3.vx_volr) ./ Ym;
  Ygs  = smooth3(Yg);
  
  % greater mask - distance based brain radius (brad)
  [dilmsk,resT2] = cat_vol_resize(Yb,'reduceV',resT3.vx_volr,mean(resT3.vx_volr)*2,32); 
  dilmsk  = cat_vbdist(dilmsk,true(size(dilmsk)))*mean(resT2.vx_volr); %resT2.vx_volr);
  dilmsk  = dilmsk - cat_vbdist(single(dilmsk>0),true(size(dilmsk)))*mean(resT2.vx_volr); %,resT2.vx_volr);
  dilmsk  = cat_vol_resize(smooth3(dilmsk),'dereduceV',resT2); 
  brad    = -min(dilmsk(:));
  dilmsk  = dilmsk / brad; 

  voli  = @(v) (v ./ (pi * 4./3)).^(1/3);                        % volume > radius
  brad  = double(mean([brad,voli(sum(Yb(:)>0).*prod(vx_vol))])); % distance and volume based brain radius (brad)
  
  
  % thresholds
  rf   = 6; 
  Hth  = round2(cat_stat_nanmean(Ym(Ym(:)>0.4 & Ym(:)<1.2  & Ygs(:)<0.2 & ~Yb(:) & Ydiv(:)<0.05 & Ydiv(:)>-0.5 & dilmsk(:)>0 & dilmsk(:)<10)),rf); % average intensity of major head tissues
  if isnan(Hth), Hth = 0.8; end
  GMth = round2(cat_stat_nanmean(Ym(Ym(:)>0.2  & Ym(:)<0.9      & Ygs(:)<0.2 & Yb(:) & Ydiv(:)<0.1 & Ydiv(:)>-0.1)),rf);  % first guess of the GM intensity
  CMth = round2(cat_stat_nanmean(Ym(Ym(:)>0.05 & Ym(:)<GMth*0.5 & Ygs(:)<0.2 & Yb(:) & Ydiv(:)>-0.10)),rf);  % first guess of the CSF intensity
  %WMth = cat_stat_nanmean(Ym(Ym(:)>0.8 & Ym(:)<1.2 & Ygs(:)<0.2 & ~Yb(:) & Ydiv(:)>-0.05)); 
  BGth = round2(cat_stat_nanmean(Ym(Ybg(:))),rf); 
  if isnan(CMth), CMth=mean([BGth,GMth]); end
  
  
  %% Skull-Stripping
  % intensity parameter
  gc.h = 1.3  - 0.20*gcutstr; % upper tissue intensity (WM vs. blood vessels)     - higher > more "tissue" (blood vessels)
  gc.l = 0.6  + 0.20*gcutstr; % lower tissue intensity (WM vs. blood vessels)     - higher > more "tissue" (blood vessels)
  gc.o = 0.1  + 0.10*gcutstr; % BG tissue intensity (for high contrast CSF=BG=0!) - lower value > more "tissue"
  % distance parameter
  gc.d = (brad*4 - brad*2*gcutstr)/mean(vx_vol);   % distance  parameter for downcut   - higher > more tissue
  gc.c = (0.01 - 0.02*gcutstr)*mean(vx_vol);         % growing  parameter for downcut    - higher > more tissue
  gc.f = (brad/5 - brad/10*gcutstr)/mean(vx_vol);    % closing   parameter               - higher > more tissue (higher as in gcut2, because of unkown ventricle)
  % smoothing parameter
  gc.s = -0.1 + 0.20*gcutstr;         % smoothing parameter                   - higher > less tissue

  stime = cat_io_cmd('  Skull-Stripping','g5','',verb,stime);
  Yb = (dilmsk<0.5 & Ym<1.5) & (Ym>(GMth*gc.l + (1-gc.l))) & Yg<0.3 & Ydiv<0.2 & Ym+Ydiv<gc.h & ~Ybg & Yb;
  Yb = cat_vol_morph(Yb & ~Ybg,'l',1); Yb(smooth3(Yb)<0.5)=0; 
  
  %% the hull is required to remove large scull tissues such as muscle in apes\monkeys 
  [hull,resT2] = cat_vol_resize(single(Yb),'reduceV',resT3.vx_volr,mean(resT3.vx_volr)*4,32); 
  hull  = cat_vol_morph(hull,'labclose',3); 
  hull  = cat_vbdist(single(hull<=0),true(size(hull)))*mean(resT2.vx_volr);
  hull  = cat_vol_resize(cat_vol_smooth3X(hull,2),'dereduceV',resT2); 
  %%
  Yb2 = single(smooth3(cat_vol_morph(cat_vol_morph(Yb & hull/max(hull(:))>0.2 & ~Ybg,'l')>0,'d',max(hull(:))*0.3) )>0.5 + gc.s); 
  
  Yb  = single(smooth3( cat_vol_morph(Yb & Yb2 & ~Ybg,'lo',2) |  cat_vol_morph(Yb & Yb2 & ~Ybg & hull/max(hull(:))>0.2,'lo') | ...
     cat_vol_morph(Yb & Yb2 & ~Ybg & hull/max(hull(:))>0.4,'l') )>0.5 + gc.s); 
  %if ~debug, clear hull; end 
  [dilmsk2,resT2] = cat_vol_resize(single(Yb),'reduceV',resT3.vx_volr,mean(resT3.vx_volr)*2,32); 
  dilmsk2  = cat_vbdist(dilmsk2,true(size(dilmsk2)))*mean(resT2.vx_volr); %resT2.vx_volr);
  dilmsk2  = dilmsk2 - cat_vbdist(single(dilmsk2>0),true(size(dilmsk2)))*mean(resT2.vx_volr); %,resT2.vx_volr);
  dilmsk2  = cat_vol_resize(smooth3(dilmsk2),'dereduceV',resT2); 
  dilmsk2  = dilmsk2 / brad; 
  %% WM growing
  Yb = cat_vol_morph(Yb,'d');
  Yb(Yb<0.5 & (dilmsk2>0.1 | Ym>1.05 | Ym<mean([1,GMth]) | (Yg.*Ym)>0.5))=nan;
  [Yb1,YD] = cat_vol_downcut(Yb,Ym,0.01); 
  Yb(isnan(Yb))=0; Yb((YD)<gc.d/2)=1; Yb(isnan(Yb))=0;
  Yb = smooth3(Yb)>0.2 + gc.s;
  %% the hull is required to remove large scull tissues such as muscle in apes\monkeys 
  [hull,resT2] = cat_vol_resize(single(Yb),'reduceV',resT3.vx_volr,mean(resT3.vx_volr)*4,32); 
  hull  = cat_vol_morph(hull,'labclose',3); 
  hull  = cat_vbdist(single(hull<=0),true(size(hull)))*mean(resT2.vx_volr);
  hull  = cat_vol_resize(cat_vol_smooth3X(hull,2),'dereduceV',resT2); 
  %% ventricle closing
  [Ybr,resT2] = cat_vol_resize(single(Yb),'reduceV',resT3.vx_volr,mean(resT3.vx_volr)*4,32); 
  Ybr = cat_vol_morph(Ybr>0,'lc',2*gc.f/mean(resT2.vx_volr));
  Ybr = cat_vol_resize(smooth3(Ybr),'dereduceV',resT2)>0.5 + gc.s; 
  Yb = single(Yb | (Ym<1.2 & Ybr & hull/max(hull(:))>0.3));
  %% GWM growing
  Yb(Yb<0.5 & (dilmsk2>0.4  | Ym>0.95 | Ym<GMth | (Yg.*Ym)>0.8))=nan;
  [Yb1,YD] = cat_vol_downcut(Yb,Ym,0.02 + gc.c); 
  Yb(isnan(Yb))=0; Yb((YD)<gc.d)=1; Yb(isnan(Yb))=0;
  Yb = smooth3(Yb)>0.5 + gc.s; 
  Yb = single(Yb | (Ym>0.2 & (Ym<0.9 | hull/max(hull(:))<0.2) & cat_vol_morph(Yb,'lc',max(1,min(3,0.2*gc.f)))));
  %% GM growing
  Yb(Yb<0.5 & (dilmsk2>0.5  | Ym>0.9 | Ym<CMth | (Yg.*Ym)>0.9))=nan;
  [Yb1,YD] = cat_vol_downcut(Yb,Ym,0.01 + gc.c);
  Yb(isnan(Yb))=0; Yb((YD)<gc.d)=1; Yb(isnan(Yb))=0; clear Yb1 YD; 
  Yb(smooth3(Yb)<0.5 + gc.s)=0;
  Yb = single(Yb | (Ym>0.1 & (Ym<0.9 | hull/max(hull(:))<0.2) & cat_vol_morph(Yb,'lc',max(1,min(3,0.1*gc.f)))));
  %% CSF growing (add some tissue around the brain)
  Yb(Yb<0.5 & (dilmsk2>0.6  | Ym< gc.o | Ym>1.1 | (Yg.*Ym)>0.9))=nan;
  [Yb1,YD] = cat_vol_downcut(Yb,Ym,-0.01 + gc.c); Yb(isnan(Yb))=0; 
  Yb(isnan(Yb))=0; Yb(YD<gc.d)=1; Yb(isnan(Yb))=0; clear Yb1 YD; 
  Yb(smooth3(Yb)<0.7 + gc.s)=0;
  Yb  = single(cat_vol_morph(Yb,'o',max(1,min(3,4 - 0.2*gc.f))));
  Yb  = Yb | (Ym>0.2 & Ym<0.5 & cat_vol_morph(Yb,'lc',max(1,min(3,0.2*gc.f))));
  Yb  = cat_vol_smooth3X(Yb,2)>0.3 + gc.s;
  Ymo = Ym; 

  %% ---------------------------------------------------------
  %  improve bias correction:
  %  Also here it is important that the bias field is very smooth
  %  to avoid overcorrections. In contrast to the first
  %  correction we will try to separate between tissues...
  %  ---------------------------------------------------------
  stime = cat_io_cmd('  Tissue classification','g5','',verb,stime);
  Ym   = Ymo;
  Yg   = cat_vol_grad(Ym,vx_vol);%./max(eps,Ym)
  Ygs  = smooth3(Yg);
  Ydiv = cat_vol_div(Ym,vx_vol);

      % WM Skeleton:  Ydiv./Yg<-1    
% nicht WM: Ydiv.*Yg<-0.02

  %% tissue classes WM, GM, subcortical GM (Ybm), CSF, head tissue (Yhm) 
  Ywm  = ((Ym-max(0,(Ydiv+0.01))*10)>(GMth*0.1+0.9) | Ym>(GMth*0.1+0.9) | ...
          (Ydiv./Yg<0.5 & ((Ym>0.9 & Yg>0.1 & Ydiv<0) | (Ym>0.9)) & Ym<1.2) ) & ...
          Ym<1.3 & Yg<0.6 & Ygs<0.9 & Yb & Ydiv<0.1 & Ydiv>-0.5; 
  CSFD = cat_vbdist(single(Ym<(CMth*0.5+0.5*GMth)),Yb,vx_vol);
  Ywm  = Ywm & CSFD>2;
  Ywm  = smooth3(Ywm)>0.5;      
  % subcotical GM 
  Ybm  = ((Ym-max(0,(Ydiv+0.01))*10)<0.98) & Ym<0.98 & dilmsk<-brad*0.3 & ... 
         Ym>(GMth*0.6+0.4) & Yb & Ygs<0.2 & Yg<0.2 & ... Ydiv<0.1 & Ydiv>-0.02 & ... Yg<0.1 & 
         ~(Ydiv./Yg<0.5 & ((Ym>0.9 & Yg>0.1 & Ydiv<0) | (Ym>0.95)) & Ym<1.2);
       %& ~Ywm;  
  Ybm  = smooth3(Ybm)>0.5;
  Ybm  = cat_vol_morph(Ybm,'o',1);
  % cortical GM 
  Ygm  = Ym<(GMth*0.3+0.7) & Ym>(CMth*0.6+0.4*GMth) & Yg<0.4 & Yb & Ydiv<0.4 & Ydiv>-0.3 & ~Ywm & ~Ybm; % & (Ym-Ydiv*2)<GMth;  
  Ygm(smooth3(Ygm)<0.3 | ~cat_vol_morph(Ywm,'d',3/mean(vx_vol)))=0;
  Ygm(CSFD<3 & Ym>(CMth*0.5+0.5*GMth) & ~Ywm & Ym<(CMth*0.5+0.5*GMth))=0; 
  % CSF
  Ycm  = Ym<(CMth*0.5+0.5*GMth) & Yg<0.1 & Yb & ~Ygm & dilmsk<-brad*0.3; 
  Ycm  = smooth3(Ycm)>0.5;
  % head tissue
  Yhm  = Ym>max(mean([CMth,GMth]),Hth*0.2) & Ym<1.2 & Yg<0.8 & cat_vol_smooth3X(Yb,2)<0.1 & Ydiv<0.6 & Ydiv>-0.6;
  Yhm  = smooth3(Yhm)>0.5;

  %% refine
  %Ygm  = Ygm | (cat_vol_morph(Ywm,'d',3) & ~Ybm & ~Ywm & (Ym-Ydiv*2)<GMth & ...
  %   ~Ycm & smooth3((Ym + Yg)<(CMth*0.8+0.2*GMth))<0.5) & Ym>(CMth*0.9+0.1*GMth) & Ym<(GMth*0.2+0.8);;
  
  %% masking of the original values and local filtering
  stime = cat_io_cmd('  Filtering','g5','',verb,stime);
  fi   = round2(max(3,min(resT3.vx_volr)*3)/3); 
  Ywm  = Ysrc .* Ywm; Ywm  = cat_vol_localstat(Ywm,Ywm>0,1,3); % PVE
  Ycm  = Ysrc .* Ycm; Ycm  = cat_vol_localstat(Ycm,Ycm>0,1,2);
  for i=1:fi-1, Ywm = cat_vol_localstat(Ywm,Ywm>0,2,1); end
  for i=1:fi-1, Ycm = cat_vol_localstat(Ycm,Ycm>0,2,1); end
  Ybm  = Ysrc .* Ybm; for i=1:fi, Ybm = cat_vol_localstat(Ybm,Ybm>0,2,1); end
  Ygm  = Ysrc .* Ygm; for i=1:fi, Ygm = cat_vol_localstat(Ygm,Ygm>0,2,1); end
  Yhm  = Ysrc .* Yhm; for i=1:fi, Yhm = cat_vol_localstat(Yhm,Yhm>0,2,1); end

  % estimate intensity difference bettween the tissues
  Ywmr = cat_vol_resize(Ywm,'reduceV',resT3.vx_volr,8,16,'meanm'); 
  Ybmr = cat_vol_resize(Ybm,'reduceV',resT3.vx_volr,8,16,'meanm'); 
  Ygmr = cat_vol_resize(Ygm,'reduceV',resT3.vx_volr,8,16,'meanm'); 
  Ycmr = cat_vol_resize(Ycm,'reduceV',resT3.vx_volr,8,16,'meanm'); 
  bmth = mean(Ybmr(Ybmr(:)>0 & Ywmr(:)>0))/mean(Ywmr(Ybmr(:)>0 & Ywmr(:)>0));
  gmth = mean(Ygmr(Ygmr(:)>0 & Ywmr(:)>0))/mean(Ywmr(Ygmr(:)>0 & Ywmr(:)>0));
  cmth = mean(Ycmr(Ycmr(:)>0 & Ywmr(:)>0))/mean(Ywmr(Ycmr(:)>0 & Ywmr(:)>0));
  Ywmr = cat_vol_resize(Ywm,'reduceV',resT3.vx_volr,16,8,'meanm'); 
  Yhmr = cat_vol_resize(Yhm,'reduceV',resT3.vx_volr,16,8,'meanm'); 
  hmth = mean(Yhmr(Yhmr(:)>0 & Ywmr(:)>0))/mean(Ywmr(Ywmr(:)>0 & Ywmr(:)>0));
  hmth = min(max(hmth,mean([cmth,gmth])),mean([gmth,1]));
  % if something failed use global thresholds
  if isnan(bmth), bmth = GMth; end
  if isnan(gmth), gmth = GMth; end
  if isnan(cmth), cmth = CMth; end
  if isnan(hmth), hmth = GMth; end
  clear Ywmr Ybmr Ygmr Yhmr Ycmr; 

  Yhm = Yhm .* (dilmsk>20);  % to avoid near skull tissue
  
  %% estimate bias fields
  stime = cat_io_cmd('  Bias correction','g5','',verb,stime);
  Ywi = sum( cat(4,Ywm,Ygm/gmth,Ybm/bmth,Ycm/cmth,Yhm/hmth),4) ./ sum( cat(4,Ywm>0,Ybm>0,Ygm>0,Ycm>0,Yhm>0),4 );
  if ~zeroBG
    Ybg2 = Ybg(:) & Yg(:)<(cat_stat_nanmean(Yg(Ybg(:))) + 2*(cat_stat_nanstd(Yg(Ybg(:))))); 
    Ywi(Ybg2) = Ysrc(Ybg2); clear Ybg2;
  end
  %%
  [Ywi,resT2]  = cat_vol_resize(Ywi,'reduceV',resT3.vx_volr,min(4,min(resT3.vx_volr)*2),32,'meanm'); 
  for i=1:4, Ywi=cat_vol_localstat(Ywi,Ywi>0,2,1); end
  Ywi   = cat_vol_approx(Ywi,'nn',resT2.vx_volr,2);
  Ywi   = cat_vol_smooth3X(Ywi,4.*mean(vx_vol)); 
  Ywi   = cat_vol_resize(Ywi,'dereduceV',resT2);

  %% background noise
  if zeroBG
    stime = cat_io_cmd('  Background correction','g5','',verb,stime);
    %Ybc  = cat_vol_morph(smooth3(Ym<mean([BGth,CMth]) & Ym<CMth & Ygs<0.05 & ~Yb & dilmsk2>8)>0.5,'lo',3); 
    [Ybc,resT2] = cat_vol_resize(Ysrc .* Ybg,'reduceV',resT2.vx_volr,max(8,max(16,cat_stat_nanmean(resT2.vx_volr)*4)),16,'min'); 
    Ybc  = cat_vol_localstat(Ybc,Ybc>0,2,2);
    for i=1:1, Ybc  = cat_vol_localstat(Ybc,Ybc>0,2,1); end
    %Ybc2 = cat_vol_approx(Ybc,'nn',resT2.vx_volr,4); % no aproximation to correct only in the backgound!  
    %Ybc2 = cat_vol_smooth3X(Ybc2,4);
    Ybc  = cat_vol_smooth3X(Ybc,2);
    Ybc  = cat_vol_resize(Ybc,'dereduceV',resT2); 
    %Ybc2 = cat_vol_resize(Ybc2,'dereduceV',resT2); 
  else
    % correction for negative backgrounds (MT weighting)
    [x,y]=hist(Ysrc(:),200); cx = cumsum(x)/sum(x);
    Ybc = max([min(Ysrc(:)),y(find(cx>0.01,1,'first'))]); 
  end

  %% back to original size
  stime = cat_io_cmd('  Final scaling','g5','',verb,stime);
  Ywi   = cat_vol_resize(Ywi,'dereduceV',resT3); 
  %if zeroBG, Ybc = cat_vol_resize(Ybc,'dereduceV',resT2); end
  %Ybc2 = cat_vol_resize({Ybc2},'dereduceV',resT3); 
  [Yg,Ygs]  = cat_vol_resize({Yg,Ygs},'dereduceV',resT3); 
  Yb   = cat_vol_resize(Yb,'dereduceV',resT3)>0.5; 
  Yp0 = cat_vol_resize(((Ywm>0)*3 + (Ygm>0)*2 + (Ybm>0)*2.3 + (Ycm>0)*1 + (Yhm>0)*2.7 + (Ybg>0)*0.5)/3,'dereduceV',resT3);
  Ysrc = Ysrco; clear Ysrco;

  %%  Final intensity scaling
  Ym   = (Ysrc - Ybc) ./ (Ywi - Ybc); % correct for noise only in background
 % Ym   = (Ysrc - Ybc) ./ (Ywi - Ybc2 + Ybc); % correct for noise only in background
  Wth  = single(cat_stat_nanmedian(Ym(Ygs(:)<0.2 & Yb(:) & Ym(:)>0.95))); 
  [WIth,WMv] = hist(Ym(Ygs(:)<0.2 &  Yb(:) & Ym(:)>mean([GMth,Wth]) & Ym(:)<Wth*1.1),0:0.01:2);
  WIth = find(cumsum(WIth)/sum(WIth)>0.8,1,'first'); WIth =  round2(min([mean(Ym(Ym(:)>0.5)),WMv(WIth)]),rf);  
  %[BIth,BMv] = hist(Ym(Ym(:)<mean([BGth,CMth]) & Yg(:)<0.2),-1:0.01:2);
  %BIth = find(cumsum(BIth)/sum(BIth)>0.02,1,'first'); BIth = round2(BMv(BIth),rf);  
  Ym   = Ym ./ WIth; 
  
  cat_io_cmd(' ','','',verb,stime); 
end
function X = round2(X,N)
  if ~exist('N','var'), N = 0; end
  X = round(X*10^N)/10^N; 
end