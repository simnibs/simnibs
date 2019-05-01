function [Ym,Yb,T3th3,Tth,inv_weighting,noise,cat_warnings] = cat_main_gintnorm(Ysrc,Ycls,Yb,vx_vol,res,Yy,extopts)
% This is an exclusive subfunction of cat_main.
% ______________________________________________________________________
% Global intensity normalization based on tissue thresholds estimated as 
% median intensity in the SPM tissue maps refined by edge (gradient) 
% information. Class propability should be higher than 50% (=128) to 
% avoid problems by the PVE or bias regions like basal ganglia or the CSF.
% Especialy correct CSF estimation can be problematic, because it is
% strongly influenced by the PVE and other tissues like blood vessels 
% and meninges. This structures with GM like intensity will cause a to 
% high global CSF value.
% For CSF, and WM we can use low gradient thesholds to avoid the PVE, but
% for GM this can lead to strong problems because to low thresholds will
% only give large GM areas like the basal ganlia, that have often a to high
% intensity. 
%
%   [Ym,Yb,T3th3,Tth,inv_weighting,noise,cat_warnings] =
%     cat_main_gintnorm(Ysrc,Ycls,Yb,vx_vol,res)
%
%   Ym      .. intensity normalized image
%   Yb      .. brain mask
%   T3th3   .. [CSF,GM,WM] peak intensity
%   Tth     .. structure for inverse function cat_main_gintnormi
%   inv_weighting .. true in T2/PD images
%   noise   .. first guess of the noise level
%   cat_waring .. structure with warings
%
%   Ysrc    .. the original (noise/bias corrected) image
%   Ycls    .. SPM classification [GM,WM,CSF,HD1,HD2,BG]
%              (6 cells with uint8 classes images)
%   Yb      .. brain mask
%   vx_vol  .. voxel resolution of the images
%   res     .. SPM segmentation structure
% ______________________________________________________________________
%
%   Robert Dahnke (robert.dahnke@uni-jena.de)
%   Structural Brain Mapping Group (http://dbm.neuro.uni-jena.de/)
%   Department of Neurology
%   University Jena
% ______________________________________________________________________
% $Id: cat_main_gintnorm.m 1159 2017-08-04 14:08:14Z dahnke $
  dbs   = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end

 
  if isstruct(Ycls)
    
    %% final peaks and intesity scaling
    %  -----------------------------------------------------------------
    T3th  = Ycls.T3th;
    T3thx = Ycls.T3thx;

    % intensity scaling
    Ym    = Ysrc; 
    
    if all(T3th==T3thx), return; end

    isc   = 1;
    %T3th  = interp1(T3th,1:1/isc:numel(T3th)*isc,'spline');  %pchip');
    %T3thx = interp1(T3thx,1:1/isc:numel(T3th)*isc,'spline'); %pchip');

    for i=2:numel(T3th)
      M = Ysrc>T3th(i-1) & Ysrc<=T3th(i);
      Ym(M(:)) = T3thx(i-1) + (Ysrc(M(:)) - T3th(i-1))/diff(T3th(i-1:i))*diff(T3thx(i-1:i));
    end
    M  = Ysrc>=T3th(end); 
    Ym(M(:)) = numel(T3th)/isc/6 + (Ysrc(M(:)) - T3th(i))/diff(T3th(end-1:end))*diff(T3thx(i-1:i));    
    Ym = Ym / 3; 

    return
  end
    
  clsint  = @(x) round( sum(res.mn(res.lkp==x) .* res.mg(res.lkp==x)') * 10^5)/10^5;
  clsints = @(x,y) [round( res.mn(res.lkp==x) * 10^5)/10^5; res.mg(res.lkp==x-((y==0)*8))']; 
 
  inv_weighting = 0;
  if nargout==7
    cat_warnings = struct('identifier',{},'message',{});
  end
  vxv    = 1/cat_stat_nanmean(vx_vol);
  res.mn = round(res.mn*10^5)/10^5; 
  
  if cat_get_defaults('extopts.subfolders')
    reportfolder  = 'report';
  else
    reportfolder  = '';
  end
  
  %% initial thresholds and intensity scaling
  T3th3 = [clsint(3) clsint(1) clsint(2)];
  BGth  = min(mean(Ysrc(Ycls{6}(:)>192)),clsint(6));
  
  %% -------------------------------------------------------------------
  %  intensity checks and noise contrast ratio (contrast part 1)
  %  -------------------------------------------------------------------
  % relation between the GM/WM and CSF/GM and CSF/WM contrast has to be
  % greater that 3 times of the maximum contrast (max-min).
  clear Yn
  checkcontrast = @(T3th,minContrast) ...
    abs(diff(T3th([1,3]))) < (max(T3th(:))-min(T3th(:)))*minContrast || ...
    abs(diff(T3th(1:2)))   < (max(T3th(:))-min(T3th(:)))*minContrast || ...
    abs(diff(T3th(2:3)))   < (max(T3th(:))-min(T3th(:)))*minContrast;
   
  if checkcontrast(T3th3,1/9) && exist('cat_warnings','var') % contrast relation
    cat_warnings = cat_io_addwarning(cat_warnings,...
      'CAT:cat_main:LowContrast',...
      sprintf(['The contrast between different tissues is relative low! \n' ...
           '  (BG=%0.2f, CSF=%0.2f, GM=%0.2f, WM=%0.2f)\n'],BGth,T3th3),numel(cat_warnings)==0);
  end
  
  if T3th3(1)>T3th3(3) && T3th3(2)>T3th3(3) && T3th3(1)>T3th3(2) % invers (T2 / PD)
    cat_warnings = cat_io_addwarning(cat_warnings,...
      'CAT:cat_main:InverseContrast',...
      sprintf(['Inverse tissue contrast! \n' ...
           '(BG=%0.2f, CSF=%0.2f, GM=%0.2f, WM=%0.2f)\n'],BGth,T3th3(1:3)),numel(cat_warnings)==0);
    T3th3(1) = max( max(clsints(3,0)) , mean(Ysrc(Ycls{3}(:)>240)));     
    
    % first initial scaling for gradients and divergence
    if abs(diff( abs(diff( T3th3/diff(T3th3([3,1])) )) ))>0.4 % %T2
      T3th  = [min(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:)))) ...
               min( T3th3(3)*0.8+0.2*min(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:)))) , BGth ) ...
               T3th3 ...
               ([T3th3(3)*0.5+0.5*min(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:)))) ... WM
                cat_stat_nanmean([T3th3(3), BGth ]) ... WM
                T3th3(2)*0.5 + 0.5*T3th3(1)... % CSF/GM
                max(T3th3) + abs(diff(T3th3([1,3])/2)) ... % CSF / BG
                 ]) ];
      T3thx = [0,0.05, 1,2,3.2, 1.1, 1.0, 1.75, 0.8]; 

      [T3th,si] = sort(T3th);
      T3thx     = T3thx(si);
    else 
      T3th  = [min(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:)))) ...
               min( T3th3(3)*0.2+0.8*min(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:)))) , BGth ) ...
               T3th3 ...
               ([T3th3(3)*0.5+0.5*min(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:)))) ... WM
                cat_stat_nanmean([T3th3(3), BGth ]) ... WM
                T3th3(2)*0.5 + 0.5*T3th3(1)... % CSF/GM
                max(T3th3) + abs(diff(T3th3([1,3])/2)) ... % CSF / BG
                max(T3th3(end) + abs(diff(T3th3([1,numel(T3th3)])/2)) , ...
                 max(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:)))) ) ]) ];
      T3thx = [0,0.05, 1,2,3, 2.0, 1.0, 1.75, 0.8, 0.2];

      [T3th,si] = sort(T3th);
      T3thx     = T3thx(si);
    end
    
    Ym = Ysrc+0; 
    for i=2:numel(T3th)
      M = Ysrc>T3th(i-1) & Ysrc<=T3th(i);
      Ym(M(:)) = T3thx(i-1) + (Ysrc(M(:)) - T3th(i-1))/diff(T3th(i-1:i))*diff(T3thx(i-1:i));
    end
    M  = Ysrc>=T3th(end); 
    Ym(M(:)) = numel(T3th)/6 + (Ysrc(M(:)) - T3th(i))/diff(T3th(end-1:end))*diff(T3thx(i-1:i));    
    Ym = Ym / 3; 
    %%
    Yg    = cat_vol_grad(Ym,vx_vol);
    Ydiv  = cat_vol_div(Ym,vx_vol);
    
    %% tissues for bias correction
    Ycm   = (Ym + Yg - Ydiv + single(Ycls{2})/255)<2/3 | ...
            (Ycls{3} + Ycls{6} + Ycls{4} + Ycls{5})>128; 
    Ycd   = cat_vbdist(single(Ycm));
    Ybd   = cat_vbdist(cat_vol_morph(single((Ycls{6} + Ycls{4} + Ycls{5})>128),'lo',1));
    
    Ywm  = (single(Ycls{2})/255 - Yg - Ydiv - max(0,3-Ycd-Ybd/40)/2)>0.7 | ... 
           (Ym-Yg-Ydiv-max(0,3-Ycd-Ybd/40)/2)>0.8 & Ycls{1}+Ycls{2}>240;
    Ywm(smooth3(Ywm)<0.3)=0;
    Ygm  = (single(Ycls{1})/255 - abs(Ydiv)*8)>0.5 | ...
           (single(Ycls{2}+Ycls{1})>240 & max(0,2-Ycd - max(0,Ybd/2-10))>0 & abs(Ydiv)<0.1);

    %% bias correction
    [Yi,resT2] = cat_vol_resize(Ysrc.*Ywm,'reduceV',vx_vol,1,16,'min');
    Yig        = cat_vol_resize(Ysrc.*Ygm./median(Ysrc(Ygm(:)))*median(Ysrc(Ywm(:))),'reduceV',vx_vol,1,16,'meanm');
    Yi = max(Yi,Yig); Yi(Yig>0) = min(Yig(Yig>0),Yi(Yig>0));
    Yi = cat_vol_localstat(Yi,Yi>0,1,2);
    for xi=1:2, Yi = cat_vol_localstat(Yi,Yi>0,2,1); end
    Yi = cat_vol_approx(Yi,'nh',resT2.vx_volr,2); Yi = cat_vol_smooth3X(Yi,2); 
    Yi = cat_vol_resize(Yi,'dereduceV',resT2)./median(Yi(Ycls{2}>192));  
    Ysrcr = round(Ysrc ./ Yi * 10^5)/10^5; % * T3th3(3) * 1.05
    if debug==0, clear Yg Ydiv Yn Yi; end
    
    %% final thresholds
    if 0  % old 
      T3th  = [min(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:)))) ...
               min( T3th3(3)*0.2+0.8*min(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:)))) , BGth ) ...
               T3th3 ...
               ([T3th3(3)*0.3+0.7*min(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:)))) ...
                cat_stat_nanmean([T3th3(3),BGth]) ...
                clsint(2)*0.8 ...
                max(T3th3) + abs(diff(T3th3([1,numel(T3th3)])/2)) ...
                max(T3th3(end) + abs(diff(T3th3([1,numel(T3th3)])/2)) , ...
                max(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:)))))]) ];
      T3thx = [0,0.05, 1,2,3, 2.9, 2.5, 2.0, 1.0, 0.7];
    end
    if abs(diff( abs(diff( T3th3/diff(T3th3([3,1])) )) ))>0.4 % %T2
      T3th  = [min(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:)))) ...
               min( T3th3(3)*0.8+0.2*min(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:)))) , BGth ) ...
               T3th3 ...
               ([T3th3(3)*0.5+0.5*min(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:)))) ... WM
                cat_stat_nanmean([T3th3(3), BGth ]) ... WM
                T3th3(2)*0.5 + 0.5*T3th3(1)... % CSF/GM
                max(T3th3) + abs(diff(T3th3([1,3])/2)) ... % CSF / BG
                 ]) ];
      T3thx = [0,0.05, 1,2,3.2, 1.1, 1.0, 1.75, 0.8]; 
    else
      T3th  = [min(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:)))) ...
             min( T3th3(3)*0.2+0.8*min(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:)))) , BGth ) ...
             T3th3 ...
             ([T3th3(3)*0.5+0.5*min(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:)))) ... WM
              cat_stat_nanmean([T3th3(3), BGth ]) ... WM
              T3th3(2)*0.5 + 0.5*T3th3(1)... % CSF/GM
              max(T3th3) + abs(diff(T3th3([1,3])/2)) ... % CSF / BG
              max(T3th3(end) + abs(diff(T3th3([1,numel(T3th3)])/2)) , ...
               max(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:)))) ) ]) ];
      T3thx = [0,0.05, 1,2,3, 2.0, 1.0, 1.75, 0.8, 0.2];
    end

    
    [T3th,si] = sort(T3th);
    T3thx     = T3thx(si);
    
    
    inv_weighting = 1;
    
  elseif T3th3(1)<T3th3(2) && T3th3(2)<T3th3(3) % T1
    %%
    BGmin = min(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:)))); 
    T3th3(1) = min( min(clsints(3,0)) , mean(Ysrc(Ycls{3}(:)>240))); 
    BGcon = max([BGmin*1.1,T3th3(1) - cat_stat_nanmean(diff(T3th3)),median(Ysrc(Ycls{6}(:)>128))]);
    %T3th3 = [max( min(res.mn(res.lkp==3 & res.mg'>0.3/sum(res.lkp==3)))*.05 + .95*max(res.mn(res.lkp==2 & res.mg'>0.3/sum(res.lkp==2))) , ...
    %              min(res.mn(res.lkp==3 & res.mg'>0.3/sum(res.lkp==3)))) ...
    %         max(res.mn(res.lkp==1 & res.mg'>0.1)) ...
    %         max(res.mn(res.lkp==2 & res.mg'>0.1))];
    T3th  = [BGmin ... minimum
             BGcon ... cat_stat_nanmean background (MT contrast with strong background noise)
             T3th3 ... csf gm wm 
             max(T3th3) + abs(diff(T3th3([1,numel(T3th3)])/2)) ... higher
             max(T3th3(end) + abs(diff(T3th3([1,numel(T3th3)])/2)) , ... maximum
              max(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:))))) ];
    T3thx = [0,0.05,1:5];
    Ysrcr = round(Ysrc*10^5)/10^5; 

    
    
  elseif T3th3(1)>T3th3(3) && T3th3(2)<T3th3(3) && T3th3(1)>T3th3(2)
  % This is a very special case of T2 weighting (of neonates) with
  % BG < GM < WM < CSF that need special correction for the sulcal 
  % CSF that mostly have WM like intensity due to the PVE. Hence, 
  % there are a lot of miss-labeled voxel that need further correction 
  % that follow in the second part of this file.
  %
  % Because, I observed a lot of problems by the SPM bias correction, I
  % increased the biasfwhm and the biasreg. Although, this did not work, 
  % the indroduce bias was smoother and could be corrected here very well. 
  
    cat_warnings = cat_io_addwarning(cat_warnings,...
      'CAT:cat_main:InverseContrast',...
      sprintf(['Inverse tissue contrast that require strong modifications! \n' ...
           'In case of "BG<GM<WM<CSF", the CSF in sulci got WM-like intensities \n' ...
           'due to the PVE and require severe correction that may fail! \n' ...
           '  (BG=%0.2f, CSF=%0.2f, GM=%0.2f, WM=%0.2f)\n'],BGth,T3th3(1:3)),numel(cat_warnings)==0);
    
    Sth   = clsint(2);      
    Ym    = (Ysrc - BGth) / ( Sth - BGth); 
    Yp0   = single(Ycls{3})/255/3 + single(Ycls{1})/255*2/3 + single(Ycls{2})/255;
    Yg    = cat_vol_grad(Ym,vx_vol);
    Ydiv  = cat_vol_div(Ym,vx_vol);
    [Yp0r,resTh] = cat_vol_resize(Yp0,'reduceV',vx_vol,max(vx_vol)*3,16,'meanm'); 
    Yhd   = cat_vbdist(min(max(0,1-Yp0r*3),1)); 
    Yhd   = cat_vol_resize(Yhd,'dereduceV',resTh);   
    Yhd   = Yhd ./ max(Yhd(Yhd(:)>0)); 
    Ysrco = Ysrc+0;

    T3th2 = [clsint(3) ...
             min( [clsint(1) , ...
                   cat_stat_nanmean(Ysrco( Ysrco(:)>(BGth*0.8+0.2*Sth) & Ysrco(:)<cat_stat_nanmean(Sth) & ...
                    Yp0(:)>1.9/3 & Yp0(:)<2.2/3 & Ydiv(:)<0.05 & Yg(:)>0.05 & Yg(:)<0.15))  ]) ... 
             max( [clsint(2) median( Ysrco(Ycls{2}(:)>192) ) ]) ...
            ];
          
     
    %%      
    if ~exist('Yy','var') && isfield(res,'Twarp'), Yy = res.Twarp; end
    LAB = extopts.LAB;
    if ~exist('Yy','var')
      PA  = extopts.cat12atlas;
      VA  = spm_vol(PA{1});
      YA  = cat_vol_ctype(spm_sample_vol(VA,double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),0));
      YA  = reshape(YA,size(Ym));
      YA(mod(YA,2)==0 & YA>0)=YA(mod(YA,2)==0 & YA>0)-1;  
    else
      YA = ones(size(Yg)); 
    end
    if ~debug, clear Yy; end
          
          
    %% bias correction
    Ym   = (Ysrco - BGth) / ( Sth - BGth); 
    Ycbg = Yg<0.1 & (Ym<1 | Ycls{1}>Ycls{2} & Ydiv>=0) & YA==LAB.CB; 
    Ycbw = Yg<0.1 & Ycls{2}>Ycls{1} & Ydiv<0 & YA==LAB.CB & ~Ycbg; 
    Ybg  = cat_vol_morph(cat_vol_morph(cat_vol_smooth3X(Ym<1.02 & Yp0>2/3 & Yp0<1 & Yg<cat_stat_nanmean(Yg(:))/2,1.2)>0.8 & ...
      Yhd>0.5,'o',1),'d',4) & Yp0>1.5/3 & Ym<1.02 & Yg<cat_stat_nanmean(Yg(:)) & cat_vol_morph(YA==LAB.BG | YA==LAB.TH,'d',2); % subcortical structures
    Ygw  = (Yp0>1.9/3 | (Yp0>0 & Ym>0.7 & Ym<1.4)) & cat_vol_morph(cat_vol_morph(smooth3(Ycls{1}>240 & Yg<0.1 & abs(Ydiv)<0.01)<0.3,'c',1),'e',3) & ...
      Ym<(clsint(3)/clsint(2)*0.8 + 0.2) & Ym>(clsint(1)/clsint(2)*0.6) & Yg<0.4 & ~Ybg & ~Ycbg; 
    Ygw = cat_vol_morph(Ygw,'l',[inf 0.05])>0;
    Ygw  = Ygw | (Yg<0.1 & Ycls{2}>Ycls{1} & Ycls{2}>Ycls{3} & (Yp0>1.5/3 | Ym>0.95 & Ym<1.5) & ~Ybg & ~Ycbg) | YA==LAB.BS | Ycbw;
    Ygw = cat_vol_morph(Ygw,'l',[inf 0.1])>0;
    Ygw  = Ym .* Ygw;
    Ygw  = cat_vol_median3(Ygw,Yp0>0);
    Ygm  = Ym .* (cat_vol_morph(smooth3(~Ygw & Yg<0.10 & abs(Ydiv)<0.05 & Ycls{1}>Ycls{2} & Ycls{1}>Ycls{3} & Yhd<0.4)>0.5,'o',2) | Ycbg )/ ...
            (T3th2(2) / T3th2(3)); 
    %Ycm  = Ym .* cat_vol_morph(smooth3(~Ygw & ~Ygm & Yg<0.1 & (Ycls{3}>Ycls{2} & Ycls{3}>Ycls{1})>0.5) | (Ym>1.2 & Yb)); 
    %Ycm  = cat_vol_localstat(Ycm,Ycm>0,1,3);
    %Ycm  = Ycm / mean(Ycm(Ycm(:)>0)); %  * (T3th2(1) / T3th2(3)); 
    Ygw  = Ygw + (Ygm .* (Ygw==0)) + ((Ym .* Ybg .* (Ygw==0)) / (cat_stat_nanmean(T3th2(2)) / T3th2(3)) ); % + (Ycm .* (Ygm==0)); 
    Ygw2 = Ym .* (Yg<max(0.1,min(0.2,cat_stat_nanmean(Yg(:)))) & ( ((Yp0>2.9/3 | Ym>0.9) & Ym<1.5 & Ygw>0) | Ycbw | ...
                  (abs(Ydiv)<0.1 & Ycls{2}/2>Ycls{1} & Ycls{2}/2>Ycls{3}) ) & ~Ybg & ~Ycbg);
    %% field approximation
    [Ywi,Ywi2,resT2] = cat_vol_resize({Ygw,Ygw2},'reduceV',vx_vol,max(vx_vol)*2,16,'max'); 
    for i=1:1, Ywi2 = cat_vol_localstat(Ywi2,Ywi2>0,1,3); end
    for i=1:1, Ywi = cat_vol_localstat(Ywi,Ywi>0,2,3); end % only one iteration!
    for i=1:1, Ywi = cat_vol_localstat(Ywi,Ywi>0,1,1); end
    Ywi  = cat_vol_approx(Ywi,'nn',resT2.vx_volr,2);
    Ywi  = cat_vol_smooth3X(Ywi,1); Ywi(Ywi2>0)=Ywi2(Ywi2>0);
    Ywi  = cat_vol_smooth3X(Ywi,1); % highres data have may stronger inhomogeneities 
    Ywi  = cat_vol_resize(Ywi,'dereduceV',resT2);    
    Ywi  = Ywi / cat_stat_nanmean(Ywi(Ygw2>0)); 
    Ysrc  = Ysrco ./ Ywi * ( cat_stat_nanmean(Ysrco(Ygw2>0)/T3th2(3)) / cat_stat_nanmean(Ywi(Ygw2>0)) );
    
    
    
    %% first initial scaling for gradients and divergence
    %T3th3 = [cat_stat_nanmean( res.mn(res.lkp==3 & res.mg'>0.1)) ...
    %         min( [res.mn(res.lkp==1 & res.mg'>0.1) , ...
    %               cat_stat_nanmean(Ysrc( Ysrc(:)>(BGth*0.8+0.2*Sth) & Ysrc(:)<cat_stat_nanmean(Sth) & ...
    %                Yp0(:)>1.9/3 & Yp0(:)<2.2/3 & Ydiv(:)<0.05 & Yg(:)>0.05 & Yg(:)<0.15))  ]) ... 
    %         max( [sum(res.mn(res.lkp==2) .* res.mg(res.lkp==2)') median( Ysrc(Ycls{2}(:)>192) ) ]) ...
    %        ];
    clear T3th T3thx; 
    T3th = [ min(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:)))) ...
             min( T3th3(2)*0.5+0.5*min(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:)))) , BGth ) ...
             ...
             min(Ysrc( Ysrc(:)>(BGth*0.8+0.2*Sth) & Ysrc(:)<cat_stat_nanmean(Sth) & ...
                    Yp0(:)>1.9/3 & Yp0(:)<2.2/3 & Ydiv(:)<0.05 & Yg(:)>0.05 & Yg(:)<0.15)) ... head
             T3th3 ...
             ...
             (T3th3(3)*0.8 + 0.2*max(res.mn(res.lkp==3 & res.mn>max(res.mn(res.lkp==2)))) ) ...
             (T3th3(3)*0.5 + 0.5*max(res.mn(res.lkp==3 & res.mn>max(res.mn(res.lkp==2)))) ) ...
             (T3th3(3)*0.2 + 0.8*max(res.mn(res.lkp==3 & res.mn>max(res.mn(res.lkp==2)))) ) ...
             ...
             max(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:))))]; ...
           
    T3thx = [0,0.05,1.75, 1,2,3,  3.1,2.0,1.1,  1.0]; 

    [T3th,si] = sort(T3th);
    T3thx     = T3thx(si);
    
    Ysrcr = round(Ysrc*10^5)/10^5; 
    
    inv_weighting = 1;
  
    if debug
        Ym = Ysrcr+0;
        for i=2:numel(T3th)
          M = Ysrc>T3th(i-1) & Ysrc<=T3th(i);
          Ym(M(:)) = T3thx(i-1) + (Ysrc(M(:)) - T3th(i-1))/diff(T3th(i-1:i))*diff(T3thx(i-1:i));
        end
        M  = Ysrc>=T3th(end); 
        Ym(M(:)) = numel(T3th)/6 + (Ysrc(M(:)) - T3th(i))/diff(T3th(end-1:end))*diff(T3thx(i-1:i));    
        Ym = Ym / 3; 
        ds('l2','',vx_vol,Ysrc/T3th3(3), round(Ym*3),Ysrc/T3th3(3),Ym,60)
    end
  

  else    
    error('CAT:cat_main:badTissueContrast',...
      sprintf('Bad tissue contrast (C=%0.2f, G=%0.2f, W=%0.2f)\n',...
        T3th3(1),T3th3(2),T3th3(3)),numel(cat_warnings)==0); %#ok<SPERR>
  end

  
  
  
  %% intensity scaling for gradient estimation
  Tth.T3th  = T3th;
  Tth.T3thx = T3thx;
    
  Ym = Ysrcr+0; 
  for i=2:numel(T3th)
    M = Ysrcr>T3th(i-1) & Ysrcr<=T3th(i);
    Ym(M(:)) = T3thx(i-1) + (Ysrcr(M(:)) - T3th(i-1))/diff(T3th(i-1:i))*diff(T3thx(i-1:i));
  end
  M  = Ysrcr>=T3th(end); 
  Ym(M(:)) = numel(T3th)/6 + (Ysrcr(M(:)) - T3th(i))/diff(T3th(end-1:end))*diff(T3thx(i-1:i));    
  Ym = Ym / 3; 
 
  
 
  %% new initial segment threshold
  if ~exist('Yg','var'), Yg  = cat_vol_grad(Ym,vx_vol)./max(eps,Ym); end
  T3th  = [median(Ysrcr(Ycls{3}(:)>192 & Yg(:)<0.20 & Ym(:)<0.45)) ...
           median(Ysrcr(Ycls{1}(:)>192 & Yg(:)<0.20)) ...
           median(Ysrcr(Ycls{2}(:)>192 & Yg(:)<0.10))];
  Ynw   = cat_vol_localstat(Ysrc,Ycls{3}>192,2,4);
  Ync   = cat_vol_localstat(Ysrc,Ycls{2}>192,2,4); 
  noise = round(min(cat_stat_nanmean(Ynw(Ynw(:)>0)),cat_stat_nanmean(Ync(Ync(:)>0))) / min(abs(diff(T3th(1:3)))) * 10^6)/10^6; 
  clear Ynw Ync;
 
  if debug==2
    [pth,nam] = spm_fileparts(res.image0(1).fname);
    tmpmat = fullfile(pth,reportfolder,sprintf('%s_%s%02d%s.mat',nam,'write',1,'gintnorm00'));
    save(tmpmat,'Ysrc','Ycls','Yb','vx_vol','res','T3th','T3thx','Yg','Ym','noise');
  end
  
  
  
  
  

  %% -------------------------------------------------------------------
  %  check modality (contrast part 2)
  %  -------------------------------------------------------------------
  %  It is possible to invert T2 and PD images based on the SPM class 
  %  information, but actual there is no time to develope and proof this 
  %  function in detail, due to the most other functions ...
  %  -------------------------------------------------------------------
  if T3th(1)<T3th(2) && T3th(2)<T3th(3)
  %  -------------------------------------------------------------------
  %  standard T1 contrast
  %  -------------------------------------------------------------------
  %  For T1 data SPM mean tissue values were not always correct. 
  %  Especially, high and low contrast images or images with incomplete
  %  inhomogeneity correction can have bad peaks (ADHD200/..NYC..14). 
  %  So it is better to use the SPM segments and add some further 
  %  knowledge (gradient & divergence) to refine these segments and 
  %  estimate the median value of the segment that is typcialy more 
  %  stable than the mean value. 
  %  -------------------------------------------------------------------
   
    % check SPM segmentation
    if exist('cat_warnings','var')
      Ymx = single(Ycls{1})/255*2/3 + single(Ycls{2})/255+ single(Ycls{3})/255*1/3;  
      Ygw = Yb & ((Ycls{1}+Ycls{2})>128);
      Ymp0diff = sqrt(cat_stat_nanmean(Ym(Ygw(:)) - Ymx(Ygw(:)))^2); 
      if Ymp0diff>0.10 && debug
        cat_warnings = cat_io_addwarning(cat_warnings,...
          'CAT:cat_main:badSPMsegment',sprintf(...
          ['SPM segmentation does not fit to the image (RMS(Ym,Yp0)=%0.2f).\n'...
           'This can be an alignment problem (check origin), ' ...
           'untypical subjects (neonates, non-human),\n'...
           'bad image contrast (C=%0.2f,G=%0.2f,W=%0.2f), \n'...
           'low image quality (NCR~%0.2f), or something else ...'],Ymp0diff,T3th,noise),numel(cat_warnings)==0); 
      end
      clear Ymx;
    end
    
    
    %% skull-stripping warning
    skulltest = (median(Ysrc(Ycls{5}(:)>192 & Ysrc(:)>T3th(2))) < ... 
       median(Ysrc(Ycls{3}(:)>192 & Ysrc(:)>0))); 
    if exist('cat_warnings','var') &&  (isnan(skulltest) || skulltest)
      
      % Skull-Stripped images can of course lead to problems with to strong
      % brain masks, but the bigger problem here is that the CSF intensity 
      % threshold were maybe affected. 
     
      % If a skull-stripping was used, we will use this as initial mask 
      % that we close and dilate a little bit. 
      % Now, the original image can be corrected in the stripped area, 
      % because some images have missing points (slicewise). Becuase of 
      % the gaussian functions a hard boundary is better.
      if Ymp0diff<0.05 && numel(Ysrc>0)/numel(Ysrc)<0.8
        Yb    = smooth3(cat_vol_morph(cat_vol_morph(Ysrc>0,'lc',3),'d'))>0.5;
        CSFth = min([nanmedian(Ysrc(Ycls{3}(:)>240 & Ysrc(:)>0)), ... 
                     nanmedian(Ysrc(Ycls{3}(:)>192 & Ysrc(:)>0)), ... 
                     nanmedian(Ysrc(Ycls{3}(:)>128 & Ysrc(:)>0)), ...
                     cat_stat_nanmean(Ysrc(Ysrc>0))*0.5])*0.9; % 
        Ysrc  = cat_vol_laplace3R(max(CSFth,Ysrc),Yb & Ysrc==0,0.2) .* Yb;
         cat_warnings = cat_io_addwarning(cat_warnings,...
           'CAT:cat_main:SkullStripped',...
           'Skull-stripped input image detected! Try boundary cleanup.',numel(cat_warnings)==0);  
      else
         cat_warnings = cat_io_addwarning(cat_warnings,...
           'CAT:cat_main:SkullStripped',...
           'Skull-stripped input image?',numel(cat_warnings)==0); 
      end
    end
    

    %% segment refinement and median peak estimation 
    %  -----------------------------------------------------------------
    Yg    = cat_vol_grad(Ym,vx_vol);
    Ydiv  = cat_vol_div(Ym,vx_vol);
    %noise = estimateNoiseLevel(Ym,Ycls{2}>192); 
    
    Yb2   = cat_vol_morph(Yb & Ym>0.5,'e',2*vxv); 
    gth   = max(0.06,min(0.3,noise*6));
    %Ybm   = cat_vol_morph(Ycls{6}>240 & Ysrc<min(T3th),'lc'); 
    BGmin = min(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:)))); 
    BGcon = max([BGmin*1.1,T3th3(1) - cat_stat_nanmean(diff(T3th3)),median(Ysrc(Ycls{6}(:)>128))]);
    BMth  = max(BGmin,min(BGcon,T3th(1) - diff(T3th(1:2)))); %max(0.01,cat_stat_nanmedian(Ysrc(Ybm(:))));
    Ywm   = (Ycls{2}>128  & Yg<gth) | ((Ym-Ydiv*2)>(1-0.05*cat_stat_nanmean(vx_vol)) & Yb2); % intensity | structure (neonate contast problem)
    Ycm   = smooth3((Ycls{3}>240 | Ym<0.4) & Yg<gth*3 & Yb & ~Ywm & Ycls{1}<8 & Ysrc>BMth & Ym<0.5)>0.5; % important to avoid PVE!

    % If SPM get totaly wrong maps due to bad image orientations our 
    % segment were incorrect too (or empty) and peak estimation fail.
    % I try to use the kmeans, but in WM it is affected by WMHs, in 
    % CSF by blood vessels and meninges and in GM noise and subcortical
    % structures were problematic. In ADHD/..NYC..14 the basal structes 
    % get the average peak and the cortex was detected as CSF. There 
    % were much more images with smaller problems ...
    Ysrcr  = round( Ysrc.*10^5 ) / 10^5;
    WMth   = cat_stat_nanmedian(Ysrcr(Ywm(:))); % kmeans3D(Ysrc(Ycls{2}(:)>192 & Yg(:)<gth),1); % GM/WM WM  
    CSFth  = cat_stat_nanmedian(Ysrcr(Ycm(:))); % kmeans3D(Ysrc(Ycls{3}(:)>64 & Yg(:)>gth & Yb(:)),2); % CSF CSF/GM
      %  0.05 <<<<< BMth + 4*cat_stat_nanstd(Ysrc(Ybm(:)))
    Ybg    = cat_vol_morph(Yg<0.10 & Yb & Ysrc<WMth*(1-0.03*cat_stat_nanmean(vx_vol)) & Ysrc>CSFth*1.5 & Ycls{3}<64,'o',2);
    Ygm    = ~Ybg & Yg<0.4 & Ysrc<min(clsint(2)*0.8+clsint(1)*0.2,WMth+0.5*diff([CSFth,WMth])) & Yg<gth*2 & Ycls{1}>32 & ~Ywm & Ycls{2}<64 & ...
      Ysrc>(CSFth+0.1*diff([CSFth,WMth])) & ~Ywm & ~Ycm & Yb & abs(Ydiv)<0.2; 
    %Ygm   = Ygm | (Ycls{1}>64 & Ybg & ~Ywm);
    GMth   = cat_stat_nanmedian(Ysrcr(Ygm(:))); %kmeans3D(Ysrc(Ygm(:)),3); % CSF/GM GM GM/WM
    T3th_cls  = round([CSFth(1) GMth(1) WMth(1)]*10^4)/10^4;
    %clear Ybg
   %
    if any(isnan(T3th_cls)) 
      fprintf('\n');
      error('CAT:cat_main:cat_pre_gintnorm:nobrain',...
        'Bad SPM-Segmentation. Check image orientation!');
    end
    % median tissue peaks
    
   
    if debug==2
      tmpmat = fullfile(pth,reportfolder,sprintf('%s_%s%02d%s.mat',nam,'write',1,'gintnorm01'));
      save(tmpmat,'Ysrc','Ycls','Yb','vx_vol','res','T3th','Yg','Ydiv','Ym',...
       'Yb2','gth','Ybm','BMth','Ywm','Ygm','Ycm','Ybg','T3th_cls','T3th','noise');
    end
    
    
    %% final peaks and intesity scaling
    %  -----------------------------------------------------------------
    T3th3 = T3th_cls;
    T3th  = [min(Ysrcr(~isnan(Ysrcr(:)) & ~isinf(Ysrcr(:)))) BMth min(BGth,mean([BMth,T3th3(1)])) T3th3 ...
              T3th3(end) + diff(T3th3([1,numel(T3th3)])/2) ... WM+
              max(T3th3(end)+diff(T3th3([1,numel(T3th3)])/2) , ... max
              max(Ysrcr(~isnan(Ysrcr(:)) & ~isinf(Ysrcr(:))))) ];
    T3thx = [0,0.02,0.05,1:5];


    % intensity scaling
    Ym    = Ysrc; 
    isc   = 1;
    %T3th  = interp1(T3th,1:1/isc:numel(T3th)*isc,'spline');  %pchip');
    %T3thx = interp1(T3thx,1:1/isc:numel(T3th)*isc,'spline'); %pchip');

    for i=2:numel(T3th)
      M = Ysrc>T3th(i-1) & Ysrc<=T3th(i);
      Ym(M(:)) = T3thx(i-1) + (Ysrc(M(:)) - T3th(i-1))/diff(T3th(i-1:i))*diff(T3thx(i-1:i));
    end
    M  = Ysrc>=T3th(end); 
    Ym(M(:)) = numel(T3th)/isc/6 + (Ysrc(M(:)) - T3th(i))/diff(T3th(end-1:end))*diff(T3thx(i-1:i));    
    Ym = Ym / 3; 
    
    Tth.T3th  = T3th;
    Tth.T3thx = T3thx;
  elseif T3th3(1)>T3th3(3) && T3th3(2)>T3th3(3)
    %% reestimation of brain mask
    Yb  = Ym>0.8 & Ym<1.2 & (Ycls{5}<64); Yb  = single(cat_vol_morph(Yb,'lo',1));
    [Ybr,Ymr,Ycls5,resT2] = cat_vol_resize({single(Yb),Ym,single(Ycls{5})/255},'reduceV',vx_vol,2,32); 
    Ybr(~Ybr & (Ymr<2.5/3 | Ymr>3.2/3 | Ycls5>0.5))=nan; 
    [Ybr1,Ydr] = cat_vol_downcut(Ybr,Ymr,0.03); Ybr(Ydr<100)=1; Ybr(isnan(Ybr))=0;
    Ybr(~Ybr & (Ymr<1.9/3 | Ymr>3.2/3 | Ycls5>0.5))=nan; 
    [Ybr1,Ydr] = cat_vol_downcut(Ybr,Ymr,0.01); Ybr(Ydr<100)=1; Ybr(isnan(Ybr))=0;
    Ybr(~Ybr & (Ymr<1/3 | Ymr>2.5/3 | Ycls5>0.5))=nan; 
    [Ybr1,Ydr] = cat_vol_downcut(Ybr,Ymr,-0.01); Ybr(Ydr<100)=1; Ybr(isnan(Ybr))=0;
    Ybr = Ybr>0 | (Ymr<0.8 & cat_vol_morph(Ybr,'lc',6) & Ycls5<0.02); % large ventricle closing
    Ybr = cat_vol_morph(Ybr,'lc',2);                 % standard closing
    Yb  = cat_vol_resize(cat_vol_smooth3X(Ybr,2),'dereduceV',resT2)>0.4; 
    clear Ybr Ymr;
    %% filtering
    %YM  = Ysrc<Tth.T3th(5)/1.2; 
    %Ym(YM(:)) = Ysrc(YM(:)) / (Tth.T3th(5)/1.2);    
    YM  = (smooth3(Ysrc<Tth.T3th(5)/1.2) & smooth3(Ysrc>Tth.T3th(4))) | Ym>2; 
    Ym = cat_vol_median3(Ym,YM,Ym<1.5,0.1); 
    cat_sanlm(Ym,1,3)
  elseif T3th3(1)>T3th3(3) && T3th3(2)<T3th3(3) && T3th3(1)>T3th3(2)
     %% filtering
    Ybb  = cat_vol_morph(cat_vol_morph(Yp0>0.5/3,'c',4),'d',2);
    if debug, Ymo = Ym; end 
    if 1
      %% new approach
      Yswm = Ycls{2}>Ycls{3} & Ycls{2}/32>Ycls{1} & Yg<0.2 & abs(Ydiv)<0.1; 
      % identify CSF areas next to high intensity regions
      Ycm = smooth3(Ym<0.5 & Ycls{3}>Ycls{2} & Ycls{3}>Ycls{1} & Ysrc>(T3th3(1)*0.5+0.5*(T3th3(3))))>0.2 & (Ysrc>(T3th3(1)*0.5+0.5*(T3th3(3))));
      Ycm = cat_vol_morph(Ycm | (Yp0<1.2 & (Ym<0.5 | Ym>1.2) & cat_vol_morph(Ybb,'d',1)),'c',1);
      Ycm = single(Ycm); Ycm(Ycm==0 & (Ysrc>(BGth*0.2+0.8*T3th3(2)) & Ysrc<(T3th3(2)*0.5+0.5*T3th3(3)) | ~Ybb | Yswm)) = nan;
      [Ycm,YD]= cat_vol_downcut(Ycm,Ysrc./T3th3(3),-0.2); Ycm(YD>50)=0;  Ycm(smooth3(Ycm)<0.4)=0; 
      
      %% identify the WM 
      Ywm = ~Ycm & Ym>0.85 & Ym<1.05 & (Ycls{2}*2>Ycls{3} |  Ycls{1}*2>Ycls{3} | Yp0>1.2); Ywm(smooth3(Ywm)<0.5)=0;  
      Ywm = cat_vol_morph(Ywm,'l',[inf 0.001])>0;
      Ywm = single(Ywm); Ywm((Ywm==0 & Ym<0.8) | Ycm | ~Ybb | Ysrc>(T3th3(1)*0.5+0.5*(T3th3(3)))) = nan;
      [Ywm,YD] = cat_vol_downcut(Ywm,Ysrc./T3th3(3),0.02);  Ywm(YD>400)=0; Ywm(smooth3(Ywm)<0.55)=0; clear YD; 
      Ywm = single(Ywm); Ywm((Ywm==0 & Ym<0.8) | ~Ybb | Ysrc>(T3th3(1)*0.9+0.1*(T3th3(3)))) = nan;
      [Ywm,YD] = cat_vol_downcut(Ywm,Ysrc./T3th3(3),0.01); Ywm(YD>800)=0; Ywm(smooth3(Ywm)<0.55)=0; clear YD; 
      Ywm(Ysrc<(T3th3(2)*0.2+0.8*(T3th3(3))))=0;
      Ywm = cat_vol_morph(Ywm,'l',[inf 0.001])>0;
      %
      Ywm = single(Ywm); Ywm(Ywm==0 & (Ym<0.7 | Ycm | Ysrc<(T3th3(2)*0.2+0.8*(T3th3(3)))) | ~Ybb) = nan;
      [Ywm,YD] = cat_vol_downcut(Ywm,Ysrc./T3th3(3),0.01); Ywm(YD>800)=0; clear YD; 
      Ywm = cat_vol_morph(Ywm,'l',[inf 0.001])>0;
      %
      Ywm = single(Ywm); Ywm(Ywm==0 & (Ym<0.7 | Ycm | Ysrc<(T3th3(2)*0.5+0.5*(T3th3(3)))) | ~Ybb) = nan;
      [Ywm,YD] = cat_vol_downcut(Ywm,Ysrc./T3th3(3),-0.001); Ywm(YD>800)=0; clear YD; 
      Ywm(Ysrc<(T3th3(2)*0.5+0.5*(T3th3(3))))=0;
      Ywm = cat_vol_morph(Ywm,'l',[inf 0.001])>0;
      Ywmd = cat_vbdist(single(Ywm),Ybb);
      
%%
      Ycmx = (Ywmd .* Ym)>1 | ~Ybb | Ycm; 
      Ycmx = Ycmx & (Yg>0.01 & abs(Ydiv)>0.0001) & cat_vol_smooth3X(Yp0)<0.995 & ~Ywm ; 
      Ycmx(smooth3(Ycmx)<0.5)=0;
      Ycmx = cat_vol_morph(Ycmx,'l',[inf 0.01])>0;
      
      
      %% improve the CSF classification 
      Ycm = Ycmx | smooth3(Ym<0.5 & Ycls{3}>Ycls{2} & Ycls{3}>Ycls{1} & Ysrc>(T3th3(1)*0.5+0.5*(T3th3(3))))>0.2 & (Ysrc>(T3th3(1)*0.5+0.5*(T3th3(3))));
      Ycm = single(Ycm); Ycm(Ycm==0 & (~Ybb | Ywm | Ysrc<( T3th3(2)*0.8+0.2*(T3th3(3)))) ) = nan;
      [Ycm,YD]= cat_vol_downcut(Ycm,Ysrc./T3th3(3),0.01); Ycm(YD>50)=0; Ycm(smooth3(Ycm)<0.4)=0; clear YD; 
      Ycm = cat_vol_morph(Ycm,'l',[inf 0.1])>0;
      %%
      %Ycm = single(Ycm); Ycm((Ycm==0 & Ysrc>(T3th3(1)*0.5+0.5*(T3th3(3))) ) | ~Yb | smooth3(Ywm)>0 | Ysrc<(T3th3(2)*0.9+0.1*(T3th3(3)))) = nan;
      %[Ycm,YD]= cat_vol_downcut(Ycm,Ysrc./T3th3(3),0.005); Ycm(YD>800)=0; Ycm(smooth3(Ycm)<0.4)=0; clear YD; 
      %Ycm = Ycm | (Yp0<1.5 & Ym<0.5 & Yb); 

      %% simultan region growing of CSF and WM tissue
      Ypx = single(Ywm*3 + Ycm); Ypx(~Ybb | (Ym>0.6 & Ym<0.7)) = nan; Ypx = cat_vol_median3c(Ypx,Ybb); 
      [Ypx,YD] = cat_vol_downcut(single(Ypx),Ysrc/T3th3(3),-0.01); Ypx(YD>50)=0; Ypx = cat_vol_median3c(Ypx,Ybb); Ypx(~Ybb) = nan;
      [Ypx,YD] = cat_vol_downcut(single(Ypx),Ysrc/T3th3(3),-0.01); Ypx(YD>50)=0; clear YD; 

      %% cortical thickness based modification 
      Ys = Ypx==3 & Ym>2.5/3;
      [Ysr,Ybr,resTx] = cat_vol_resize({single(Ys),Ybb},'reduceV',vx_vol,1,16,'meanm'); 
      [Ygmt,Ypp] = cat_vol_pbt( single(1 + (Ybb.* (Ym>0.5)) + Ys) , ...
        struct('verb',0,'dmethod','eidist','method','pbt2x') );

      %% final correction
      YM = cat_vol_smooth3X(Ypp,0.5)<0.3 & Ygmt>1 & Ym>2/3 & Ypx<2; 
      YM = YM | (Ypx==1 & Ym>0.7 & Ym<0.9) | (Ycm & Ym>0.7);
      YM = YM | (cat_vol_smooth3X(cat_vol_morph(~Ywm & Ygmt>0.1 & Ygmt<1.2 & Ypx<2,'c',1),1)>0.5 & Ywm); 
      Ym(YM) = 0.5 + max(0,1 - Ym(YM))/2; 
      %%
      %Ym = cat_vol_median3(Ym,Yb & Ypp<0.5,Ypp<0.5,0.05); 
      %Ym = cat_vol_median3(Ym,Yb & Ypp>0.5,Ypp>0.5,0.05);
      %Ym = cat_vol_median3(Ym,Yb & Ym>0.8 & Ym>1.1,Yb,0.1);
      cat_sanlm(Ym,1,3)
    else
      % old approach ...
      %YM  = Ysrc<Tth.T3th(5)/1.2; 
      %Ym(YM(:)) = Ysrc(YM(:)) / (Tth.T3th(5)/1.2); 
      YM = (smooth3(Ysrc<Tth.T3th(5)/1.2) & smooth3(Ysrc>Tth.T3th(4))) | Ym>2; 
      Ym = cat_vol_median3(Ym,YM,Ym<1.5,0.1);
    end
  end
    
  
  
  
  
  %% if there was a warning we need a new line 
  if nargout==7 && numel(cat_warnings)>1, fprintf('\n'); cat_io_cmd(' ','','',1); end

end