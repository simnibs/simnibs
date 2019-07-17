function [Ygmt,Ypp,Ywmd,Ycsfdc] = cat_vol_pbt(Ymf,opt)
% ______________________________________________________________________
%
% Cortical thickness and surface position estimation. 
% 
%   [Ygmt,Ypp,Ywmd,Ycsfd] = cat_vol_pbt(Ymf,opt)
%  
%   Ygmt:      GM thickness map 
%   Ypp:       percentage position map
%   Ywmd:      WM distance map
%   Ycsfd:     CSF distance map
% 
%   Ymf:       tissue segment image or better the noise, bias, and 
%              intensity corrected 
%
%   opt.resV   voxel resolution (only isotropic)
%   opt.method choose of method {'pbt2x','pbt2'} with default=pbt2x as 
%              the method that is described in the paper.
% ______________________________________________________________________
%
%   Dahnke, R; Yotter R; Gaser C.
%   Cortical thickness and central surface estimation.
%   NeuroImage 65 (2013) 226-248.
% ______________________________________________________________________
%
%   Robert Dahnke (robert.dahnke@uni-jena.de)
%   Structural Brain Mapping Group (http://dbm.neuro.uni-jena.de/)
%   Department of Neurology
%   University Jena
%
% ______________________________________________________________________
% $Id: cat_vol_pbt.m 1218 2017-11-17 16:27:43Z dahnke $ 


% default variables and check/set function  
  if ~exist('opt','var'), opt=struct(); end

  def.resV      = 1;
  def.dmethod   = 'eidist';
  def.method    = 'pbt2x';  % pbt is worse ... just for tests!
  def.debug     = cat_get_defaults('extopts.verb')>2;
  def.verb      = cat_get_defaults('extopts.verb')-1;
  opt           = cat_io_checkinopt(opt,def);
  opt.resV      = mean(opt.resV);
  
  minfdist = 0; 
  


  %% Distance maps
  if (sum(round(Ymf(:))==Ymf(:)) / numel(Ymf))>0.9, bin=1; else bin=0; end
  
  %  WM distance 
  %  Estimate WM distance Ywmd and the outer CSF distance Ycsfdc to correct
  %  the values in CSF area are to limit the Ywmd to the maximum value that 
  %  is possible within the cortex.
  %  The increasement of this area allow a more accurate and robust projection. 
  %  cat_vol_eidist used speed map to align voxel to the closer gyrus
  %  that is not required for the correction map.
  if opt.verb, fprintf('\n'); end
  stime = cat_io_cmd('    WM distance: ','g5','',opt.verb); stime2=stime;
  YMM = cat_vol_morph(Ymf<1.5,'e',1) | isnan(Ymf);
  switch opt.dmethod
    case 'eidist' 
      % [D,I] = vbm_vol_eidist(B,L,[vx_vol,euclid,csf,setnan,verb])
      YM  = max(0,min(1,(Ymf-2))); YM(YMM) = nan; Ywmd   = cat_vol_eidist(YM,max(0.5,min(1,Ymf/2)),[1 1 1],1,1,0,opt.debug); 
      YM  = max(0,min(1,(Ymf-1))); YM(YMM) = nan; Ycsfdc = cat_vol_eidist(YM,max(1.0,min(1,Ymf/2)),[1 1 1],1,1,0,opt.debug); 
    case 'vbdist'
      YM  = max(0,min(1,(Ymf-2))); Ywmd   = max(0,cat_vbdist(single(YM>0.5),~YMM)-0.5); 
      YM  = max(0,min(1,(Ymf-1))); Ycsfdc = max(0,cat_vbdist(single(YM>0.5),~YMM)-0.5); 
  end
  clear YMM; 
  if ~bin
    % limit the distance values outside the GM/CSF boudary to the distance possible in the GM
    YM  = Ywmd>minfdist & Ymf<=1.5; Ywmd(YM) = Ywmd(YM) - Ycsfdc(YM); Ywmd(isinf(Ywmd)) = 0; clear Ycsfdc;
    % smoothing of distance values inside the GM
    YM  = Ywmd>minfdist & Ymf> 1.5; YwmdM = Ywmd; YwmdM = cat_vol_localstat(YwmdM,YM,1,1); Ywmd(YM) = YwmdM(YM);
    % smoothing of distance values outside the GM
    YM  = Ywmd>minfdist & Ymf<=1.5; YwmdM = Ywmd; for i=1:2, YwmdM = cat_vol_localstat(YwmdM,YM,1,1); end; Ywmd(YM) = YwmdM(YM);
    % reducing outliers in the GM/CSF area
    YM  = Ywmd>minfdist & Ymf< 2.0; YwmdM = Ywmd; YwmdM = cat_vol_median3(YwmdM,YM,YM); Ywmd(YM) = YwmdM(YM); clear YwmdM YM;
  end
  
  minfdist = 1; 
  %  CSF distance
  %  Similar to the WM distance, but keep in mind that this map is
  %  incorrect in blurred sulci that is handled by PBT
  stime = cat_io_cmd('    CSF distance: ','g5','',opt.verb,stime);
  YMM = cat_vol_morph(Ymf<1.5,'e',1) | cat_vol_morph(Ymf>2.5,'e',1) | isnan(Ymf); % this was dilate???
  switch opt.dmethod
    case 'eidist'
      YM  = max(0,min(1,(2-Ymf)));   YM(YMM) = nan; Ycsfd = cat_vol_eidist(YM,max(0.5,min(1,(4-Ymf)/2)),[1 1 1],1,1,0,opt.debug); 
      YM  = max(0,min(1,(3-Ymf)));   YM(YMM) = nan; Ywmdc = cat_vol_eidist(YM,max(1.0,min(1,(4-Ymf)/2)),[1 1 1],1,1,0,opt.debug); 
      YM  = max(0,min(1,(2.7-Ymf))); YM(YMM) = nan; Ywmdx = cat_vol_eidist(YM,max(1.0,min(1,(4-Ymf)/2)),[1 1 1],1,1,0,opt.debug)+0.3;
    case 'vbdist'
      YM  = max(0,min(1,(2-Ymf)));   Ycsfd = max(0,cat_vbdist(single(YM>0.5),~YMM)-0.5); 
      YM  = max(0,min(1,(3-Ymf)));   Ywmdc = max(0,cat_vbdist(single(YM>0.5),~YMM)-0.5); 
      YM  = max(0,min(1,(2.7-Ymf))); Ywmdx = max(0,cat_vbdist(single(YM>0.5),~YMM)-0.2); 
  end
  Ywmdc = min(Ywmdc,Ywmdx);
  clear YMM;
  if ~bin
    YM = Ycsfd>minfdist & Ymf>=2.5; Ycsfd(YM) = Ycsfd(YM) - Ywmdc(YM); Ycsfd(isinf(-Ycsfd)) = 0; clear Ywmdc;
    YM = Ycsfd>minfdist & Ymf< 2.5; YcsfdM = Ycsfd; YcsfdM = cat_vol_localstat(YcsfdM,YM,1,1); Ycsfd(YM) = YcsfdM(YM);
    YM = Ycsfd>minfdist & Ymf>=2.5; YcsfdM = Ycsfd; for i=1:2, YcsfdM = cat_vol_localstat(YcsfdM,YM,1,1); end; Ycsfd(YM) = YcsfdM(YM);
    YM = Ycsfd>minfdist & Ymf> 2.0; YcsfdM = Ycsfd;  YcsfdM = cat_vol_median3(YcsfdM,YM,YM); Ycsfd(YM) = YcsfdM(YM); clear YcsfdM YM;
  end  


  %% PBT thickness mapping 
  %  PBT is the default thickness estimation, but PBT2x is the optimized
  %  version that use both sulci and gyris refinements, because not only 
  %  thin sulci can blurred. PBT2x is furthermore the methode that is
  %  described in the paper.
  iter = 1/mean(opt.resV);
  if strcmp(opt.method,'pbt2x')  
    % Estimation of the cortical thickness with sulcus (Ygmt1) and gyri 
    % correction (Ygmt2) to create the final thickness as the minimum map
    % of both.
    stime = cat_io_cmd('    PBT2x thickness: ','g5','',opt.verb,stime);
    
    % estimate thickness with PBT approach
    Ygmt1 = cat_vol_pbtp(Ymf,Ywmd,Ycsfd);  
    Ygmt2 = cat_vol_pbtp(4-Ymf,Ycsfd,Ywmd); 
   
    % avoid meninges !
    Ygmt1 = min(Ygmt1,Ycsfd+Ywmd);
    Ygmt2 = min(Ygmt2,Ycsfd+Ywmd);  
    
    % median filter to remove outliers
    Ygmt1 = cat_vol_median3(Ygmt1,Ygmt1>0,Ygmt1>0);
    Ygmt2 = cat_vol_median3(Ygmt2,Ygmt2>0,Ygmt2>0); 
    
    % estimation of Ypp for further GM filtering without sulcul blurring
    Ygmt  = min(cat(4,Ygmt1,Ygmt2),[],4);
    YM    = Ymf>=1.5 & Ymf<2.5; Ypp = zeros(size(Ymf),'single'); Ypp(Ymf>=2.5)=1;
    Ypp(YM) = min(Ycsfd(YM),Ygmt(YM) - Ywmd(YM)) ./ (Ygmt(YM) + eps); Ypp(Ypp>2) = 0;
    YM  = (Ygmt<=opt.resV & Ywmd<=opt.resV & Ygmt>0); Ypp(YM) = (Ymf(YM)-1)/2; 
    Ygmts = Ygmt; for i=1:iter, Ygmts = cat_vol_localstat(Ygmts,Ygmt1>0,1,1); end; Ygmt(Ygmts>0) = Ygmts(Ygmts>0);
    
    % filter result
    Ygmts = Ygmt1; for i=1:iter, Ygmts = cat_vol_localstat(Ygmts,(Ygmt>1 | Ypp>0.1) & Ygmt>0 & (Ygmt>1 | Ymf>1.8),1,1); end; Ygmt1(Ygmts>0) = Ygmts(Ygmts>0); 
    Ygmts = Ygmt2; for i=1:iter, Ygmts = cat_vol_localstat(Ygmts,(Ygmt>1 | Ypp>0.1) & Ygmt>0 & (Ygmt>1 | Ymf>1.8),1,1); end; Ygmt2(Ygmts>0) = Ygmts(Ygmts>0);

    % mix result 
    % only minimum possible, because Ygmt2 is incorrect in blurred sulci 
    Ygmt  = min(cat(4,Ygmt1,Ygmt2),[],4); %clear Ygmt1 Ygmt2; 
    Ygmts = Ygmt; for i=1:iter, Ygmts = cat_vol_localstat(Ygmts,(Ygmt>1 | Ypp>0.1) & Ygmts>0 & (Ygmt>1 | Ymf>1.8),1,1); end; Ygmt(Ygmts>0) = Ygmts(Ygmts>0);
  else
    % Estimation of thickness map Ygmt and percentual position map Ypp.
    stime = cat_io_cmd('    PBT2 thickness: ','g5','',opt.verb,stime);
    
    % estimate thickness with PBT approach
    Ygmt  = cat_vol_pbtp(Ymf,Ywmd,Ycsfd);   
    
    % filter result
    Ygmts = Ygmt; for i=1:iter, Ygmts = cat_vol_localstat(Ygmts,Ygmt>0 & Ymf>1.5,1,1); end; Ygmt(Ygmts>0) = Ygmts(Ygmts>0);
  end
  
  
  %% Estimation of a mixed percentual possion map Ypp.
  stime = cat_io_cmd('    Final Corrections: ','g5','',opt.verb,stime);
  YM  = Ymf>=1.5 & Ymf<2.5 & Ygmt>eps;
  Ycsfdc = Ycsfd; Ycsfdc(YM) = min(Ycsfd(YM),Ygmt(YM) - Ywmd(YM)); 
  Ypp = zeros(size(Ymf),'single'); Ypp(Ymf>=2.5)=1;
  Ypp(YM) = Ycsfdc(YM) ./ (Ygmt(YM) + eps); 
  Ypp(Ypp>2) = 0;
  YM  = (Ygmt<=opt.resV & Ywmd<=opt.resV & Ygmt>0); Ypp(YM) = (Ymf(YM)-1)/2 - 0.2; % correction of voxel with thickness below voxel resolution
  Ypp(isnan(Ypp)) = 0; 
  Ypp(Ypp<0) = 0; 
  
  %% Final corrections for position map with removing of non brain objects.
  % ds('d2','',1,Ymf/3,Ywmd/3,Ygmt/5,Ypp,70)
  
  % Final corrections for thickness map with thickness limit of 10 mm. 
  % Resolution correction of the thickness map after all other operations, 
  % because PBT actually works only with the voxel-distance (isotropic 1 mm)
  Ygmt = Ygmt*opt.resV; 
  Ygmt(Ygmt>10) = 10; 
  if exist('Ycsfdc','var'), Ycsfdc = Ycsfdc*opt.resV; end
  if exist('Ywmd','var'), Ywmd = Ywmd*opt.resV; end
  
  cat_io_cmd(' ','g5','',opt.verb,stime);
  if opt.debug, cat_io_cmd(' ','','',opt.debug,stime2); end
end