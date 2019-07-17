function  [Ymc,Ymi] = cat_run_job_APP_SPM(Po,vout,vx_vol,verb,fstr)
%  _____________________________________________________________________
%  The final bias correction is a subfunction of cat_run_job.
%  
%  [Ymc,Ymi] = cat_run_job_APP_SPM(Po,vout,vx_vol,verb,fwhm)
%
%  _____________________________________________________________________
%  Robert Dahnke
%  $Id: cat_run_job_APP_SPM.m 1212 2017-11-10 12:39:37Z dahnke $


  % if there is a breakpoint in this file set debug=1 and do not clear temporary variables 
  dbs   = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end
  
  %resMth  = 0.5; % general resolution limitation 
  resTth  = 1.0; % tissue intensity estimation resolution limitation 
  
  if ~exist('verb','var'), verb = 1; end
  if ~exist('fstr','var'), fstr = 0.5; end
  fwhm    = max(1.5,min(4.5,4.5 - 3*fstr)); 
  resMth  = 1 - 0.75*fstr; 

  if verb>0, fprintf('\n'); end
  stime = cat_io_cmd('    Initialize','g5','',verb); 
  
  Vo   = spm_vol(Po);
  Yo   = single(spm_read_vols(Vo));
  Ym   = vout.Ym; 
  Ymo  = Ym; 
  Ycls = vout.Ycls; 
  res  = vout.res; 
  if ~debug; clear vout; end 

  % general resolution limitation 
  [Ym,resTM] = cat_vol_resize(Ym,'reduceV',vx_vol,resMth,32,'meanm');
  Yo         = cat_vol_resize(Yo,'reduceV',vx_vol,resMth,32,'meanm');
  for i=1:6, Ycls{i} = cat_vol_resize(Ycls{i},'reduceV',vx_vol,resMth,32); end
  vx_vol = resTM.vx_volr;

 
  %% prepare data
  Yp0   = single(Ycls{1})/255*2 + single(Ycls{2})/255*3 + single(Ycls{3})/255;
  Yb    = cat_vol_morph(cat_vol_morph(Yp0>0.6,'lo'),'lc',1); 

  % global intensity normalization 
  if any( min(vx_vol*2,resTth)./vx_vol >= 2 )
    Ymr = cat_vol_resize(Ym,'reduceV',vx_vol,min(vx_vol*2,resTth),32,'meanm');
    Ybr   = cat_vol_resize(single(Yb),'reduceV',vx_vol,min(vx_vol*2,resTth),32,'meanm')>0.5;
    Yclsr = cell(size(Ycls)); for i=1:6, Yclsr{i} = cat_vol_resize(Ycls{i},'reduceV',vx_vol,min(vx_vol*2,resTth),32); end
    [Ymr,Ybr,T3th,Txth,inv_weighting,noise] = cat_main_gintnorm(Ymr,Yclsr,Ybr,vx_vol,res);
    clear Ymr Ybr Yclsr; 
    Ymi = cat_main_gintnorm(Ym,Txth); 
  else
    [Ymi,Yb,T3th,Txth,inv_weighting,noise] = cat_main_gintnorm(Ym,Ycls,Yb,vx_vol,res);
  end

  
  % gradient & divergence maps
  stime = cat_io_cmd('    Prepare measures','g5','',verb,stime);
  Yg    = cat_vol_grad(Ymi,vx_vol*2); gth = mean(Yg(Yb(:)));
  Ydiv  = cat_vol_div(Ymi,vx_vol/2);
  Ycsfd = cat_vbdist(single( ((~Yb | Ymi<0.5) & (Yp0<1.5 | Ycls{3}>128) &  Ymi<0.6 )| ~Yb),true(size(Yp0)),vx_vol); % csf distance

  %% brain and head distances
  [Yp0r,Ymir,Ybr,resT1] = cat_vol_resize({Yp0,Ymi.*Yb,Yb},'reduceV',vx_vol,2,32,'meanm');
  Ybd   = cat_vbdist(min(Yp0r,1) ,true(size(Yp0r)),resT1.vx_volr); 
  Yhd   = cat_vbdist(single(~Ybr),true(size(Yp0r)),resT1.vx_volr); Yhd  = Yhd/max(Yhd(Yhd(:)<10^6));
  Yd    = cat_vbdist(single(Ymir<0.55)); Yd  = Yd/max(Yd(:));
  Ydi   = cat_vbdist(single(Yd>0.8),Ymir>0.5); Ydi = Ydi/max(Ydi(Ydi(:)<10^4));
  Ybd   = cat_vol_resize(Ybd,'dereduceV',resT1); % brain distance
  Yhd   = cat_vol_resize(Yhd,'dereduceV',resT1); % head distance
  Ydi   = cat_vol_resize(Ydi,'dereduceV',resT1); % head distance
  if ~debug, clear Yd; end
 
  %% == tissues ==
  
  % subcortical structures
  Yss  = Ymi>0.6 & (Ymi + max(0,Ydi*30-3) < 0.98 ) & Ymi<0.98 & Yhd./(Yhd + Ydi)>0.7 & Ydi<0.2; 
  Ybgi = cat_vol_localstat(Yss.*Ymi,Yss,2,1);
  Yss(abs(((Yss.*Ymi)./Ybgi)-Yss)>0.1 | abs(((Yss.*Ymi)./Ybgi)-Yss)==0)=0;
  if ~debug, clear Ybgi; end 
  Yss(smooth3(Yss)<0.6)=0; 
  Yss  = cat_vol_morph(Yss,'l',[6 0.2])>0;
  [Yssc,resT1] = cat_vol_resize(Yss,'reduceV',vx_vol,4,32,'meanm');
  Yssc  = cat_vol_morph(Yssc,'c',8); 
  Yssc  = cat_vol_resize(Yssc,'dereduceV',resT1);
  Yss   = Yss | (Yssc & Ymi<2.8/3 & Ymi>0.5 & cat_vol_morph(Yss,'d',4));
  Yss(smooth3(Yss)<0.6)=0; 
  if ~debug, clear Yssc; end 
  
  %% cortex close to head
  Yct = Yb & Ycsfd>0 & Ycsfd<2.5 & Yhd<0.2 & cat_vol_smooth3X(Ycls{6}>128,16)>0.01 & Ymi>0.5 & (Yg./Ymi)<4*gth & Ycls{2}<192 & Ycls{3}<192;
  Yct(smooth3(Yct)<0.5) = 0;
  Ycw = Yb & Ycsfd>0 & Ycsfd<4 & Yhd<0.2 & cat_vol_smooth3X(Ycls{6}>128,16)>0.01 & Ymi>0.9 & Ycls{1}<128 & Ycls{3}<128 & ~Yct & Yb;
  Ycw(smooth3(Ycw)<0.5) = 0;
  
  %% WM
  Ywm = cat_vol_morph(cat_vol_morph(Ycls{2}>8 & Ycls{5}<192,'l'),'lc'); 
  Ywm = Yb & (Ywm | cat_vol_morph(Ymi>.95 & Ymi<1.2 & Yp0>2 & (Yp0<1.8 | Yp0>2.2) ,'lc')) &  ...
        Ycls{1}<240 & Ycls{3}<32 & Yg<0.5 & abs(Ydiv)<0.5 & Ycsfd>2 & ~Yss & ~Yct;
  Ywm = Ywm | Ycw | (Yb & Ymi-Ydiv+((Ycsfd-3)/10)>0.9 & Ydiv<-0.05 & Ycls{1}<240 & Ycls{3}<32 & ~Yss & ~Yct); 
  Ywm = cat_vol_morph(Ywm,'l',[3 0.1])>0; Ywm(smooth3(Ywm)<0.5)=0;
  Ywme1 = cat_vol_morph(Ywm,'e'); 
  Ywme2 = cat_vol_morph(Ywme1,'e'); 
  Ymiw = cat_vol_localstat(Ymi,Ywm & ~Ywme2,2,3);
  Ygmb = (Ywm & ~Ywme2 & Ymi./max(eps,Ymiw)<0.9 & Ycsfd<3); 
 
  
  %% GM
  Ygm  = Yb & ~Ywme1 & Ycls{1}>64 & Ycls{2}<192 & Ycls{3}<128 & Yg./Ymi<gth*3 & abs(Ydiv)<gth*3 & (Ycsfd<5);  
  Ygm  = Ygm | Yct | Ygmb | (Ycls{1}>192 & Ycls{2}<8 & Yg./Ymi<gth*2 & ~Ywme1) | ...
         (Yss & abs(Ydiv)<gth & Yg./Ymi<gth*1 & Ymi<0.95 & Ycls{1}>4 & Ycls{2}<252 ); %  | ...
         %(Ycls{1}>64 & Ycls{2}<240 & Ycls{3}<128 & Ycsfd<3 & Yhd<3 & ~Ywm & abs(Ydiv)<gth*3); 
  Ygmi = cat_vol_localstat(Ygm.*Ymi,Ygm,2,1);
  Ygm(abs(((Ygm.*Ymi)./Ygmi)-Ygm)>gth*mean(vx_vol)*3 | abs(((Ygm.*Ymi)./Ygmi)-Ygm)==0)=0;
  Ygm(smooth3(Ygm)<0.4)=0; 
  Ywm  = Ywm & ~Ygm; 
  
  %% CM (this did not work yet) 
  %  
  if 0
    Ycm = Yb & ~Ywm & ~Ygm & Yg<gth*2 & Ymi>0.01 & (Ymi-Ydiv)<0.5 & Ycls{3}>128 & cat_vol_smooth3X(Yb,4)>0.95; 
    Ycm(smooth3(Ycm)<0.6)=0;
  end
  
  
  %% HM 
  %  there is now way to use high intensity information from the scull, but
  %  it is possbile to use the large areas of mussels 
  Yhm = Ybd>2 & smooth3( Yg>gth*2 | Ymi>1.2 | Ydiv<-0.1)>0.5; 
  Yhm = Yg./Ymi<gth & ~Yb & abs(Ydiv)<0.2 & Ycls{6}<128 & Ym<T3th(3)*1.5 & (Ybd>10 | Ycls{5}>128) & Ycls{4}<128 & ...
    Ydiv>-0.2 & Ydiv<0.2 & ~Yhm & cat_vol_smooth3X(Ycls{6}>128,4)<0.5 & Ymi>gth;
  Yhmi = cat_vol_localstat(Yhm.*Ymi,Yhm,1,2);
  Yhm(abs(((Yhm.*Ymi)./Yhmi)-Yhm)>0.5 | abs(((Yhm.*Ymi)./Yhmi)-Yhm)==0)=0;
  Yhm(smooth3(Yhm)<0.6)=0; 
  
  
  %% HBG
  %  In MR protocols with high intensity background (e.g., MP2Rage and R1)
  %  we can use the intensity information for our bias field estiamtion
  if mean(Yo(Ycls{6}(:)>192)) > T3th(2)
    Yhbg = Ycls{6}>96 & Yo > mean(T3th(1)) & ~isnan(Yo) & ~isinf(Yo) & Ydiv>-0.1 & Yg<gth*2;
    Yhbg = cat_vol_morph(Yhbg,'c',2); 
    Yhm(Yhbg)=0; Ywm(Yhbg)=0; Ygm(Yhbg)=0; 
    if exist('Ycm','var'), Ycm(Yhbg)=0; end
    Yhbgi = cat_vol_approx(Yhbg .* Yo,'nn',vx_vol,4,struct('lfO',fwhm));
    
    %%
    Ya = Ycls{2}<240 & Yp0>1.5 & Yo./(Yhbgi*T3th(3)/cat_stat_nanmean(Yhbgi(Ywm(:))))>1.05 & Yhd<0.2; Ya(smooth3(Ya)<0.5)=0;
    Ya = cat_vol_smooth3X(Ya,2)>0.1 & (Yo-Yg)/T3th(3)>1.1 ;
    Ywm(Ya)=0; Ygm(Ya)=1;
  end
  
  %%
  stime = cat_io_cmd('    Smooth values','g5','',verb,stime);
  Ywmi = cat_vol_median3(Ywm.*Yo,Ywm,Ywm,0.2);
  if inv_weighting % PVE filting 
    Ywmi = cat_vol_localstat(Ywmi,Ywm,1,2);
  else
    Ywmi = cat_vol_localstat(Ywmi,Ywm,1,3);
  end
  Ygmi = cat_vol_median3(Ygm.*Yo,Ygm,Ygm); Ygmi = cat_vol_localstat(Ygmi,Ygm,1,1); 
  Yhmi = cat_vol_median3(Yhm.*Yo,Yhm,Yhm); Yhmi = cat_vol_localstat(Yhmi,Yhm,1,1);
  if exist('Ycm','var')
    %Ycmi = cat_vol_median3(Ycm.*Yo,Ycm,Ycm); Ycmi = cat_vol_localstat(Ycmi,Ycm,1,1);
    Ycmi = cat_vol_approx(Ycm.* Yo,'nn',vx_vol,4,struct('lfO',fwhm));
  end
  
  % threshold of head tissues 
  hdth = res.mn(res.lkp==5); hdth(hdth<T3th(2) | hdth>T3th(3)) = []; 
  if isempty(hdth), hdth=T3th(2); else hdth = cat_stat_nanmedian(hdth); end

  % add tissues intensity maps
  Ytmi = Ywmi + Ygmi*(T3th(3)/T3th(2)) + Yhmi*(T3th(3)/hdth); 
  if exist('Ycm','var')
    Ytmi = Ytmi + Ycm.*Ycmi*(T3th(3)/cat_stat_nanmean(Ym(Ycm(:) & Yg(:)<gth/2 & Ycls{3}(:)>192)));
  end
  % final smoothing of brain/head tissues
  for i=1:4, Ytmi = cat_vol_localstat(Ytmi,Ytmi>0,2,1); end
  if exist('Yhbg','var')
    Ytmi = Ytmi + Yhbg.*Yhbgi*(cat_stat_nanmean(Ym(Ywm(:)))/cat_stat_nanmean(Ym(Yhbg(:) & Ybd(:)<50 & Yg(:)<0.5)));
  end
  
  %% bias field estimation 
  %    ds('l2','',vx_vol,Ymi,Ywm,Ymi,Ymc,80)
  stime = cat_io_cmd('    Estimate bias field','g5','',verb,stime);
  Ywi = cat_vol_approx(Ytmi,'nn',vx_vol,3,struct('lfO',fwhm));
  if 0
    %%
    [Ywir,Ybr,resTr] = cat_vol_resize({Ywi,single(Yb)},'reduceV',vx_vol,3,16,'meanm');
    meanY   = cat_stat_nanmedian(Ywir(Ybr(:)>0.5));
    Ywir    = Ywir / meanY; 
    Ygr     = cat_vol_grad(Ywir); 
    lfO     = min( 0.49 , max( 0.0001 , min(  mean(resTr.vx_volr)/10 , median(Ygr(Ywir(:)>0)) / 4*fwhm ))); 
    YM      = cat_vol_smooth3X(Ybr,24/mean(resTr.vx_volr(:))); YM = YM/max(YM(:));
    Ywir    = cat_vol_laplace3R(Ywir,YM>0.8,double(lfO)); 
    Ywir    = Ywir * meanY;
    [Ywir,YM] = cat_vol_resize({Ywir,max(0,YM-0.8)/0.2},'dereduceV',resTr);
    Ywir    = Ywir.*YM  +  Ywi.*(1-YM); 
  end
  %if exist('Ya','var'); Ywi(Ya) = Yo(Ya)*(T3th(3)/T3th(2)); Ywi = cat_vol_laplace3R(Ywi,cat_vol_smooth3X(Ya,2)>0.1,0.2); end
  if debug
    Ymc = Yo ./ Ywi; Ymc = cat_main_gintnorm(Ymc * (mean(Ym(Ywm(:)))/mean(Ymc(Ywm(:)))) ,Txth); 
  end

  %%
  Ywi = cat_vol_resize(Ywi,'dereduceV',resTM); 
  Ywm = cat_vol_resize(Ywm,'dereduceV',resTM)>0.9; 
  Yo  = single(spm_read_vols(Vo));
  Ymc = Yo ./ Ywi; Ymc = Ymc * (mean(Ymo(Ywm(:)))/mean(Ymc(Ywm(:))));
  cat_io_cmd(' ','','',verb,stime); 

  if nargout>1
    Ymi = cat_main_gintnorm(Ymc * (mean(Ym(Ywm(:)))/mean(Ymc(Ywm(:)))) ,Txth); 
  end
end