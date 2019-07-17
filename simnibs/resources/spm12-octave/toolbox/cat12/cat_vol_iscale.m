function [TI,varargout] = cat_vol_iscale(T,action,vx_vol,varargin)
% CAT Prprocessing Intensity Scaling Functions
% ______________________________________________________________________
% Set of functions for intensity scaling of an image T. 
%
% actions:
%   - findhead:   [TI,H]  = cat_vol_iscale(T,'findhead' ,vx_vol);
%   - findbrain:  [TI,B]  = cat_vol_iscale(T,'findbrain',vx_vol);
%   - gCGW:       [TI,tp] = cat_vol_iscale(T,'gCGW',vx_vol,...);
% 
%   [TI,varargout] = cat_vol_iscale(T,action,vx_vol,varargin)
% 
%   T            = original 3d-volume
%   action       = used method {'findhead'|'findbrain'|'gCGW'}
%                = test method {'test_findhead_findbrain'};
%   vx_vol       = 1x3 matrix with the voxel size (default [1 1 1])
%  
%   TI           = intensity scaled volumen with BG=0 and WM=1.
%   varargout{1} = head mask (for action 'findhead')
%   varargout{1} = brain mask (for action 'findbrain')
%   varargout{1} = tissue peaks (BG,CSF,GM,WM,high) 
% 
% ______________________________________________________________________
% Robert Dahnke 2012_10
% Structural Brain Mapping Group
% University Jena
%
% $Id: cat_vol_iscale.m 1103 2017-01-13 10:01:32Z gaser $
% ______________________________________________________________________

  
  if ~exist('T','var')       
    error('MATLAB:cat_vol_iscale:input','ERROR: Need input image!\n'); 
  else 
    T = single(T);
  end
  if ~exist('action','var')
    error('MATLAB:cat_vol_iscale:input','ERROR: Unkown action ''%s''!\n',action);
  end
  if ~exist('vx_vol','var')
    vx_vol=ones(1,3);
  end
  
  
  switch action
    case 'findhead'
    % UPDATE durch findbrain nötig
    % __________________________________________________________________
    
      %if nargin>3, redres = varargin{1}; end
      
      % noise estimation is only possible in edge-free regions
      [gx,gy,gz] = cat_vol_gradient3(T); G=abs(gx)+abs(gy)+abs(gz); 
      G=G./max(eps,T); G(isinf(G) | isnan(G) | G<0)=0; clear gx gy gz;
      [Tr,Gr,resT] = cat_vol_resize({T,G},'reduceV',vx_vol,12,24);

      Mr  = Gr<mean(Gr(:)) & Tr>mean(Tr(:));  
      Gth = cat_stat_nanstat1d(Gr(Mr(:)),'median');
      Mr  = cat_vol_morph(Tr>mean(Tr(Gr(:)<Gth)) & Gr<Gth & ~(Tr<mean(Tr(:)) & Gr<Gth),'o')>0; 
      Tth = cat_stat_nanstat1d(Tr(Mr(:)),'median');
      clear GMr TMr; 

      % background estimation
      Hr = cat_vol_smooth3X(cat_vol_morph(Tr>Tth*0.3,'ldc',8),1)>0.5;  %0.1 %2 % 2
      
      TI = T./Tth*0.7;
      varargout{1} = cat_vol_resize(Hr,'dereduceV',resT);
    
      
      
    case 'findbrain'
    % Finds the brain and scale the intensity of a T1 MR image based on 
    % Eikonal distance and object proberties. 
    % __________________________________________________________________

      if nargin>3, redres = varargin{1}; else redres = 4; end
      
   tic 
      % save original image 
      TO = T+0;
      
      %% reduce to ultra highres images to save time
      [T,resTO] = cat_vol_resize(T,'reduceV',vx_vol,1.2); vx_vol=resTO.vx_volr;
      
      % initial noise correction
      TSDs = [round(size(TO) * 2/5);round(size(TO) * 3/5)]; 
      TSD  = single(TO(TSDs(1):TSDs(2),TSDs(3):TSDs(4),TSDs(5):TSDs(6)));
      TSD  = TSD/cat_stat_nanstat1d(TSD,'median');
      [gx,gy,gz] = cat_vol_gradient3(TSD); GTSD=abs(gx)+abs(gy)+abs(gz); 
      TSD  = cat_vol_localstat(TSD,GTSD<0.2,1,4);
      n1   = cat_stat_nanstat1d(TSD(TSD>0),'mean');
      nc   = min(1.0,max(1/2,n1*10));
      TS   = cat_vol_smooth3X(single(T),nc);
      clear TSD TSDs nc GTSD gx gy gz; 
      % gradient map
      [gx,gy,gz] = cat_vol_gradient3(TS); G=abs(gx)+abs(gy)+abs(gz); 
      G=G./max(0,TS-mean(TS(:)/3)); G(isinf(G) | isnan(G) | G<0)=10; clear gx gy gz;

      
      % Intensity scaling for bad contrast pictures, but only for positiv 
      % minima, because heavy background noise ("sMT01-0003-00001-000176-01_MT")
      % can cause problems in background detection.
      [Tr,Gr] = cat_vol_resize({max(0,TS-(G.*TS)),G},'reduceV',vx_vol,4);
      minT = max(0,min(Tr(:))); T=T-minT; TS=TS-minT; TO=TO-minT; Tr=Tr-minT; 
      Mr   = Gr<mean(Gr(:)) & Tr>mean(Tr(:));                                                  
      Gth  = cat_stat_nanstat1d(Gr(Mr(:)),'median');
      Mr   = cat_vol_morph(Tr>mean(Tr(Gr(:)<Gth)) & Gr<Gth/2 & ~(Tr<mean(Tr(:)) & Gr<Gth/2),'o')>0;
%      Tth  = cat_stat_nanstat1d(Tr(Mr(:)),'median');
      Tth  = peak(Tr(Mr(:))/max(Tr(Mr(:))))*max(Tr(Mr(:)));
      clear GMr TMr; 

      
      
      % Addaptive noise correction:
      % We filter in every case, because although high field strength
      % produce lower noise they have higher bias that lead to stronger
      % noise after bias correction for regions with lower intensity.
      % For noise estimation a part of the WM region should be used, by
      % the labopen function that may remove one hemissphere, but allows
      % a better WM peak (Tth) estimation.
      M2  = smooth3(TS>Tth*2/3 & TS<Tth*1.6 & G<Gth/4)>0.5; M2max=max(TS(M2(:)>0));
      Tth = peak(TS(M2>0)/M2max)*M2max*0.5 + 0.5*cat_stat_nanstat1d(TS(M2(:)),'median');
      M2  = cat_vol_morph(smooth3(TS>Tth*2/3 & TS<Tth*1.6 & G<Gth/3)>0.5,'l');
      TSD = cat_vol_localstat(T./Tth,M2,1,4); 
      n2  = cat_stat_nanstat1d(TSD(TSD>0),'mean');
      nc  = min(1.0,max(1/2,n2*10));
      TS  = cat_vol_smooth3X(single(T),nc); 
      clear TSD;
      % gradient map update after smoothing
      [gx,gy,gz] = cat_vol_gradient3(TS); G=abs(gx)+abs(gy)+abs(gz); 
      G=G./max(0,TS-mean(TS(:)/3)); G(isinf(G) | isnan(G) | G<0)=10; clear gx gy gz;

      % test functions
      % fprintf(' %0.2f %0.2f | %0.3f %0.3f',n1,n2,Gth,Tth); return
      % ds('l2','',vx_vol,TS/Tth,M2,G>Gth,TS/Tth,round(size(TS,3)*3/5)); pause(2); return
    
      %% background and distance estimation to find parts of brain tissue
    
      % background estimation
      [Hr,resTr] = cat_vol_resize(TS./G>0.5 & TS/Tth>0.2,'reduceV',vx_vol,8,8,'nearest');
      Hr = cat_vol_morph(Hr,'dc',3);
      H = cat_vol_resize(Hr,'dereduceV',resTr)>0.5; clear Hrr
      
      % distance estimation to find the brain 
      [Tr,Hr,resT] = cat_vol_resize({min(1.5,TS/Tth),H},'reduceV',vx_vol,redres,32,'mean');
      [Gr]         = cat_vol_resize(G,'reduceV',vx_vol,redres,32,'max');
      HDr = cat_vbdist(single(~Hr),true(size(Hr)),resT.vx_volr);
      [Mr,Dr] = cat_vol_downcut(single(~cat_vol_morph(Hr,'e')), ...
        1-Gr,0.1,resT.vx_volr,[double(eps('single')) 1]);
      Dr = Dr .* HDr/cat_stat_nanstat1d(HDr(:),'max'); %Dr = max(0,Dr - cat_stat_nanstat1d(Dr(:),'max')/2);
      Dr(isnan(Dr) | isinf(Dr))=0; clear Mr; 
      D  = cat_vol_resize(cat_vol_smooth3X(Dr),'dereduceV',resT);
  %  toc  
   
    
    
      %% BIAS Correction:
    tic
      TSO=TS; it=1; itmax=8; CV=0; CVO=inf; CVP=1; CVth=0.001; bs=3;
      while it<=itmax && CV<CVO-CVth && CV<CVP-CVth
       
        
      % Initial correction:
      % ----------------------------------------------------------------
      %tic
      if it>=1
        % GM/WM boundary estimation to classify the WM and GM segemnts:
        res1 = 2.5; %1.6;  
        [Gr,Dr,resT1] = cat_vol_resize({G,D},'reduceV',vx_vol,res1,32,'meanm');
        TSr  = cat_vol_resize(TS,'reduceV',vx_vol,res1,32,'max');
        M2r  = cat_vol_resize(M2,'reduceV',vx_vol,res1,32,'meanm')>0.5;
        TSXr = TSr .* cat_vol_morph(Gr<Gth*2 & TSr/Tth<1.5 & TSr/Tth>0.1,'lo'); TSXr=cat_vol_median3(TSXr,TSXr>0);
        
        % very rough maximum based correction for better GWM estimation
        [WIrr,resT2] = cat_vol_resize(TSXr,'reduceV',vx_vol,4,32,'max'); TSXrr=WIrr>0;
        WIrr = cat_vol_approx(WIrr,'nh',resT1.vx_volr,8); WIrr = cat_vol_smooth3X(WIrr,max(2,min(4,Tth/std(WIrr(TSXrr(:)))))); 
        WIr  = cat_vol_resize(WIrr,'dereduceV',resT2); clear WIrr; 
        WIr  = WIr * median(TSr(M2r(:) & Gr(:)<Gth & TSr(:)>WIr(:)*0.95) ./ WIr(M2r(:) & Gr(:)<Gth & TSr(:)>WIr(:)*0.95));
        
        %% estimating of the GM-WM threshold GWM
        GWMr = TSr .* (smooth3(Gr>Gth/2 & TSr./WIr<1.2 &  Gr<Gth*2 & TSXr>0)>0.5); 
        GWMr = TSr .* (GWMr>WIr*(0.6-2*(1 - mean(GWMr(GWMr(:)>0)./WIr(GWMr(:)>0))))) .* (GWMr<WIr*0.95);
        [GWMrr,resT2] = cat_vol_resize(GWMr,'reduceV',resT1.vx_volr,8,32,'meanm'); TSXrr=GWMrr>0;
        for i=1:4, GWMrr = cat_vol_localstat(GWMrr,GWMrr>0,2,1); end
        GWMrr = cat_vol_approx(GWMrr,'nh',resT2.vx_volr,8); GWMrr = cat_vol_smooth3X(GWMrr,max(2,min(4,Tth/std(GWMrr(TSXrr(:))))));
        GWMr = cat_vol_resize(GWMrr,'dereduceV',resT2); clear GWMrr;
         
        
        % rought maximum based bias correction 
        TSXOr = cat_vol_morph(Gr<Gth*2 & TSr./WIr<1.2 &  TSr>GWMr*0.95 & TSr<GWMr*2,'l') .* TSr; TSXOr=cat_vol_median3(TSXOr,TSXOr>0);
        [WIrr,resT2] = cat_vol_resize(TSXOr,'reduceV',resT1.vx_volr,4,32,'max');  TSXrr=WIrr>0;
        WIrr = cat_vol_localstat(WIrr,WIrr>0,1,3);
        WIrr = cat_vol_approx(WIrr,'nh',resT2.vx_volr,8); WIrr = cat_vol_smooth3X(WIrr,max(2,min(4,Tth/std(WIrr(TSXrr(:))))));
        WIr  = cat_vol_resize(WIrr,'dereduceV',resT2); clear WIrr;
        WIr  = WIr * median(TSr(M2r(:) & Gr(:)<Gth & TSr(:)>WIr(:)*0.95) ./ WIr(M2r(:) & Gr(:)<Gth & TSr(:)>WIr(:)*0.95));
        
        %% update of the GWM threshold
        GWMr = TSr .* (smooth3(Gr>Gth/4 & Gr<Gth*2 & TSr>WIr*0.4 & TSr<WIr*0.95)>0.5);
        [GWMrr,resT2] = cat_vol_resize(GWMr,'reduceV',resT1.vx_volr,8,32,'meanm'); TSXrr=GWMrr>0;
        for i=1:2, GWMrr = cat_vol_localstat(GWMrr,GWMrr>0,2,1); end
        GWMrr = cat_vol_approx(GWMrr,'nh',resT2.vx_volr,8); GWMrr = cat_vol_smooth3X(GWMrr,max(2,min(4,Tth/std(GWMrr(:)))));
        GWMr = cat_vol_resize(GWMrr,'dereduceV',resT2); clear GWMrr;
  
        
        % final inital maximum based bias correction WI
        TSXOr = cat_vol_morph(Gr<Gth*2 & TSr./WIr<1.2 & TSr>min(WIr*0.9,(GWMr*0.5+0.5*WIr)) & TSr./WIr<GWMr*2,'l') .* TSr; TSXOr=cat_vol_median3(TSXOr,TSXOr>0);
        [WIrr,resT2] = cat_vol_resize(TSXOr,'reduceV',resT1.vx_volr,4,32,'max'); TSXrr=WIrr>0;
        WIrr = cat_vol_localstat(WIrr,WIrr>0,1,3);
        WIrr = cat_vol_approx(WIrr,'nh',resT2.vx_volr,8); WIrr = cat_vol_smooth3X(WIrr,max(2,min(4,Tth/std(WIrr(TSXrr(:))))));
        WIr  = cat_vol_resize(WIrr,'dereduceV',resT2); clear WIrr;
        WIr  = WIr * median(TSr(M2r(:) & Gr(:)<Gth & TSr(:)>WIr(:)*0.95) ./ WIr(M2r(:) & Gr(:)<Gth & TSr(:)>WIr(:)*0.95));

        %% rough skull-stipping
        % use the head distance map D to find a central tissue (low gradient) with WM intensity
        [Grr,Trr,GWMrr,Drr,M2rr,resT4] = cat_vol_resize({Gr,TSr./WIr,GWMr./WIr,Dr,M2r},'reduceV',vx_vol,2,32,'meanm');
        Brr  = cat_vol_morph(Trr<1.2 & Trr>max(0.8,GWMrr) & Grr<Gth & M2rr); 
        [L,num] = spm_bwlabel(double(Brr),6);
        HST = hist(L(L(:)>0),1:num); [Hs,Hi] = sort(HST,'descend'); 
        Hi(min(numel(Hi), max(2,find(Hs<Hs(1)/4,1,'first'))):end) = [];  
        maxD=0; maxL=0; 
        for l=1:numel(Hi)
          sumD=sum(Drr(L(:)==Hi(l)))/sum(L(:)==Hi(l));
          if sumD>maxD; maxD=sumD; maxL=Hi(l); end; 
        end
        Brr  = L==maxL; 
        Brr  = cat_vol_morph(Brr | (cat_vol_morph(Brr>0.5,'d',round(4 / mean(resT4.vx_volr))) & Trr>1/2 & Trr<1.2),'dc',10); 
        Br   = cat_vol_resize(cat_vol_smooth3X(Brr),'dereduceV',resT4)>0.5; clear Trr Brr;
        
        % dereduce maps
        %GWM  = cat_vol_resize(GWMr./WIr,'dereduceV',resT1);
        WI   = cat_vol_resize(WIr,'dereduceV',resT1);
        B    = cat_vol_resize(cat_vol_smooth3X(Br),'dereduceV',resT1)>0.5;  
      
        
      
        TI=TS./WI;
%        [TI,T3th,T3] = cat_vol_CGWscale(TS./WI,G,B,Gth,resT1.vx_vol);
        Tr  = cat_vol_resize(TI,'reduceV',vx_vol,res1,32,'max');
        Br  = Br>0.5 & Tr>4/5 & Tr<8/6 & Gr<Gth*1.5; 
        Br  = single(cat_vol_morph(Br,'l')); 
        Br(~Br & (Tr<5/6 | Tr>8/6))=-inf; Br = single(cat_vol_smooth3X(cat_vol_downcut(Br,Tr, 0.02/mean(resTr.vx_volr))>0,1)>0.5);
        Br(~Br & (Tr<2/3 | Tr>8/6))=-inf; Br = single(cat_vol_smooth3X(cat_vol_downcut(Br,Tr, 0.01/mean(resTr.vx_volr))>0,1)>0.5);
        Br(~Br & (Tr<1/3 | Tr>2/3))=-inf; Br = single(cat_vol_smooth3X(cat_vol_downcut(Br,Tr,-0.01*mean(resTr.vx_volr))>0,1)>0.5);
        Br(~Br & (Tr<1/6 | Tr>1/2))=-inf; Br = cat_vol_smooth3X(cat_vol_downcut(Br,Tr,-0.02*mean(resTr.vx_volr))>0,1)>0.5;
        [Trr,Brr,resTBr] = cat_vol_resize({Tr,Br},'reduceV',vx_vol,4,32); Brr=Brr>0.5;
        Brr = cat_vol_morph(Brr | (cat_vol_morph(Brr,'ldc',4) & Trr<8/6),'lo',2);
        Br  = (Br.*Tr)>0.5 | (cat_vol_resize(cat_vol_smooth3X(Brr),'dereduceV',resTBr)>0.5 & Tr<1.05);
        B    = cat_vol_resize(cat_vol_smooth3X(Br),'dereduceV',resT1)>0.5;   
      end
      %toc,tic   
      %}  
        
      %% Fine correction
      % ----------------------------------------------------------------
        % tissue segmentation and intensity scaling
        %res1 = 1.6;
        
        %TI = cat_vol_iscale(TS./WI,'gCGW',vx_vol,T3th); 
  %{      
        % final maximum based bias correction WI3
        %GM  = B & cat_vol_resize(cat_vol_smooth3X(T3r{2},0.5),'dereduceV',resT1)>0.5;  
        %GMC = B & cat_vol_resize(cat_vol_smooth3X(cat_vol_morph(T3r{2},'c',4),0.5),'dereduceV',resT1)>0.5;  
        GM  = B & T3{2};  
        [T3r,resT3] = cat_vol_resize(T3{2},'reduceV',resT1.vx_volr,4,32,'max');
        GMC = B & cat_vol_resize(cat_vol_smooth3X(cat_vol_morph(T3r,'c',4),0.5),'dereduceV',resT3)>0.5;  
        WM  = cat_vol_morph(B & TI>0.9  & TI<1.5 & G<Gth*1.5,'l'); %,'d');& TS./WI>GWM
        WM  = WM | (~GM& GMC & TI<1.5 & G<Gth*1.5 & cat_vol_smooth3X(WM,2)>0.1) & ~T3{2}; % & TS./WI>GWM 
        [WI3r,resT6] = cat_vol_resize(TSO.*WM,'reduceV',resT1.vx_volr,1.6,32,'max'); WMrr=WI3r>0;
        WI3r = cat_vol_localstat(WI3r,WI3r>0,round(3/mean(resT6.vx_volr)),3);
        WI3r = cat_vol_approx(WI3r,'linear',resT6.vx_volr,8); WI3r = cat_vol_smooth3X(WI3r,max(2,min(6,Tth/std(WI3r(WMrr(:))))));
        
        WI3 = cat_vol_resize(WI3r,'dereduceV',resT6); %clear WI3r;
        [TI,T3th,T3] = cat_vol_CGWscale(TSO./WI3,G,B,Gth,resT1.vx_vol);
        WM  = cat_vol_morph(B & TI>11/12 & TI<8/6 & G<Gth*1.5,'l');
        WM  = WM | cat_vol_smooth3X( (B & ~GM & GMC & TI>5/6 & TI<1.2 & G<Gth),0.5)>0.5;     
 %}       
        %{
        [Tr,Br,Gr,BOr,resTr] = cat_vol_resize({TI,single(WM.*B),G,B},'reduceV',vx_vol,res1,32);
        Br  = BOr>0.5 & Tr>9/12 & Tr<8/6 & Gr<Gth*1.5; 
        Br  = single(cat_vol_morph(Br,'l')); 
        Br(~Br & (Tr<5/6 | Tr>8/6))=-inf; Br = single(cat_vol_smooth3X(cat_vol_downcut(Br,Tr, 0.020/mean(resTr.vx_volr))>0,1)>0.5);
        Br(~Br & (Tr<5/6 | Tr>8/6))=-inf; Br = single(cat_vol_smooth3X(cat_vol_downcut(Br,Tr, 0.010/mean(resTr.vx_volr))>0,1)>0.5);
        Br(~Br & (Tr<2/3 | Tr>8/6))=-inf; Br = single(cat_vol_smooth3X(cat_vol_downcut(Br,Tr, 0.005/mean(resTr.vx_volr))>0,1)>0.5);
        Br(~Br & (Tr<1/3 | Tr>2/3))=-inf; Br = cat_vol_smooth3X(cat_vol_downcut(Br,Tr,-0.01*mean(resTr.vx_volr))>0,1)>0.5;
        [Trr,Brr,resTBr] = cat_vol_resize({Tr,Br},'reduceV',vx_vol,4,32); Brr=Brr>0.5;
        Brr = cat_vol_morph(Brr | (cat_vol_morph(Brr,'ldc',4) & Trr<7/6),'lo',2);
        Br  = (Br.*Tr)>0.5 | (cat_vol_resize(cat_vol_smooth3X(Brr),'dereduceV',resTBr)>0.5 & Tr<1.05);
        B   = cat_vol_resize(cat_vol_smooth3X(Br),'dereduceV',resTr)>0.5;
        %}
        
        if it==1, WI3O=ones(size(TS)); else WI3O=WI; end
     
        WM  = cat_vol_morph(B & TI>11/12 & TI<8/6 & G<Gth*1.5,'l');
        
        %% bias measurements and updates
        T = TS; TS = TSO./WI; TS = TS ./ median(TS(WM(:))); TS=TS.*Tth; it=it+1;
       
        CV   = std( TS(WM(:)>0))/mean( TS(WM(:)>0));% + std( TS(B(:)>0 & GM(:)>0))/mean( TS(B(:)>0 & GM(:)>0));
        CVP  = std(  T(WM(:)>0))/mean(  T(WM(:)>0));% + std(  T(B(:)>0 & GM(:)>0))/mean(  T(B(:)>0 & GM(:)>0));
        CVO  = std(TSO(WM(:)>0))/mean(TSO(WM(:)>0));% + std(TSO(B(:)>0 & GM(:)>0))/mean(TSO(B(:)>0 & GM(:)>0));
     %  toc
        if CVP<CV || CVO<CVP, TS=T; WI3=WI3O; break; end
        fprintf('%0.3f >> %0.3f > %0.3f\n',CVO,CVP,CV);
      end
      TS = TSO; %T = T/median(T(WM(:)>0)); 
      %fprintf('Tissue Peaks: %0.2f %0.2f %0.2f\n',T3th);
      
  %  toc
     
      %% scull-stipping on low res
    tic; clear TSO WI GWM;
      res1=1.5;
      [Tr,Br,Gr,BOr,resTr] = cat_vol_resize({TI,single(WM.*B),G,B},'reduceV',vx_vol,res1,32);
     % Tr  = cat_vol_CGWscale(Tr,Gr,BOr,Gth,vx_vol);
      Br  = BOr>0.5 & Tr>5/6 & Tr<8/6 & Gr<Gth*1.5; 
      Br  = single(cat_vol_morph(Br,'l')); 
      Br(~Br & (Tr<5/6 | Tr>8/6))=-inf; Br = single(cat_vol_smooth3X(cat_vol_downcut(Br,Tr, 0.010/mean(resTr.vx_volr))>0,1)>0.5);
      Br(~Br & (Tr<2/3 | Tr>5/6))=-inf; Br = single(cat_vol_smooth3X(cat_vol_downcut(Br,Tr, 0.001/mean(resTr.vx_volr))>0,1)>0.5);
      Br(~Br & (Tr<1/3 | Tr>2/3))=-inf; Br = cat_vol_smooth3X(cat_vol_downcut(Br,Tr,-0.01*mean(resTr.vx_volr))>0,1)>0.5;
      [Trr,Brr,resTBr] = cat_vol_resize({Tr,Br},'reduceV',vx_vol,4,32); Brr=Brr>0.5;
      Brr = cat_vol_morph(Brr | (cat_vol_morph(Brr,'lc',1) & Trr<7/6),'lo',2);
      Br  = (Br.*Tr)>0.5 | (cat_vol_resize(cat_vol_smooth3X(Brr),'dereduceV',resTBr)>0.5 & Tr<1.05);
      B   = cat_vol_resize(cat_vol_smooth3X(Br),'dereduceV',resTr)>0.5;
  %  toc,tic 

    %%  
      %T   = cat_vol_resize(T  ,'dereduceV',resTO);
      G   = cat_vol_resize(G  ,'dereduceV',resTO);
      B   = cat_vol_resize(B  ,'dereduceV',resTO)>0.5;
      WI3 = cat_vol_resize(WI3,'dereduceV',resTO);
      H   = cat_vol_resize(H  ,'dereduceV',resTO)>0.5;
      WM  = cat_vol_resize(WM ,'dereduceV',resTO)>0.5;

      % segment image
      %{
      %   ds('l2','',vx_vol,T,p0T,TO./WI2,T/Tth,90)
      WM   = cat_vol_morph(WM,'e');
      TI   = cat_vol_CGWscale(TO./WI3,G,B,Gth,vx_vol);
      TSD = cat_vol_localstat(TI,WM,1,4); 
      n3  = cat_stat_nanstat1d(TSD(TSD>0),'mean');
      nc  = min(0.9,max(0.3,n3*10));
      BV   = B & ((TI>7/6 & ~cat_vol_morph(TI>5/6,'l')) | TI>9/6); BV = 2*cat_vol_morph(BV,'d') - BV;
      p0T  = max(B,min(3,round(cat_vol_smooth3X(TI,nc)*3).*B));
      p0T(BV & p0T) = min(BV(BV & p0T),p0T(BV & p0T)); 
      p0T(cat_vol_morph(p0T==3,'l'))=3; % close small WM wholes
     %}
    %toc
      varargout{1} = B;
      varargout{2} = WM;
%      varargout{3} = H;
%      varargout{4} = cat_stat_nanstat1d(TO(WM(:)),'mean'); 
%      varargout{5} = min(1.2,max(0,(TI*3)-1)/2); %(TO./WI3);
%      varargout{5} = varargout{5}  / cat_stat_nanstat1d(varargout{5}(WM(:)),'mean');
 %     varargout{6} = p0T;
    case 'test_findhead_findbrain' 
      opt.fnamelenght = 40;
      
      if isempty(T), T = cat_io_checkfilelist(spm_select(Inf,'image','select raw images')); end
      V = spm_vol(T); 
      
      fprintf('FindBrainTest: \n');
      for subj = 1:numel(V)
        % display subject name and home directory
        [pp,ff] = spm_fileparts(V(subj).fname); [pp,hh]=spm_fileparts(pp); 
        fn = [hh filesep ff]; clear pp ff;
        fprintf(1,'%s',fliplr(sprintf(sprintf('%% %ds',opt.fnamelenght),...
          fliplr(fn(1:min(numel(fn),opt.fnamelenght)))))); clear space fn; 
      
        % test
        vx_vol  = sqrt(sum(V(subj).mat(1:3,1:3).^2)); 
        stime1  = clock; [TH,MH] = cat_vol_iscale(single(spm_read_vols(V(subj))), ...
          'findhead' ,vx_vol,4); 
        stime2  = clock; [TB,MB] = cat_vol_iscale(single(spm_read_vols(V(subj))), ...
          'findbrain',vx_vol,4); 
        ds('cat_vol_iscale','',vx_vol,TH,MH,TB,MB,round(size(TH,3)/9*5)); 
        pause(1); fprintf('\n');        
      end
      
    case 'gCGW'
    % __________________________________________________________________
    
%  function tp=tissue_peaks(T,GWM)
%  [gx,gy,gz] = cat_vol_gradient3(T); G=abs(gx)+abs(gy)+abs(gz); G=G./T; clear gx gy gz;
%  tp(3) = cat_stat_nanmedian(T(GWM(:) & T(:)>0.6 & T(:)<1.5 & G(:)<mean(G(GWM(:)))));
%  tp(2) = cat_stat_nanmedian(T(GWM(:) & T(:)>0.4 & T(:)<tp(3)));
%  tp(1) = max(0, min( tp(2)/2 , tp(2)-diff(tp(2:3)) ));
%
      if ndims(varargin{1})>2
        if max(varargin{1})==3
          for i=1:3, tp2=cat_stat_nanstat1d(T(varargin{1}(:)==i),'median'); end
          tp2(4) = min(4,tp2(3)+diff(tp2(2:3)));
        else
          B=varargin{1};

          [gx,gy,gz] = cat_vol_gradient3(T); G=abs(gx)+abs(gy)+abs(gz);
          G=G./T; G(isinf(G) | isnan(G) | G<0)=0; clear gx gy gz;
          tpa = cat_stat_nanstat1d(T(B & G<0.2 & T>mean(T(B(:)>0))),'median'); T = T / tpa; 

          [Tr,Gr,Br] = cat_vol_resize({T,G,B},'reduceV',vx_vol,2);

          
          % inital tissues values
          WMr    = Br & Gr<0.2 & Tr>0.9 & Tr<1.5;                   
          tp0(3) = cat_stat_nanstat1d(Tr(WMr(:)),'median');
          
          GMr    = Br & Gr<0.3 & Tr<0.9 & Tr>0.4;                   
          tp0(2) = cat_stat_nanstat1d(Tr(GMr(:)),'median');
          
          CMr    = Br & Gr<0.4 & Tr<max(0.1,tp0(2)-diff(tp0(2:3))); 
          tp0(1) = cat_stat_nanstat1d(Tr(CMr(:)),'median');
          
          tp0(1) = max(0.05,min(tp0(1),tp0(2)-diff(tp0(2:3))));
          tp0(4) = min(1.5,tp0(3)+diff(tp0(2:3)));

          
          % median tissue values
          WMr    = Br & Gr<0.2 & Tr>mean(tp0(2:3)) & Tr<tp0(4);              
          tp1(3) = cat_stat_nanstat1d(Tr(WMr(:)),'median');
          
          GMr    = Br & Gr<0.2 & Tr>mean(tp0(1:2)) & Tr<mean(tp0(2:3));      
          tp1(2) = cat_stat_nanstat1d(Tr(GMr(:)),'median');
          
          CMr    = Br & Gr<0.4 & Tr<mean(tp0(1:2));                          
          tp1(1) = cat_stat_nanstat1d(Tr(CMr(:)),'median');
          
          tp1(1) = max(0.05,min(tp1(1),tp1(2)-diff(tp1(2:3))));
          tp1(4) = min(1.5,tp0(3)+diff(tp0(2:3)));

          
          % peak tissue values
          %  tp2=kmeans3D(T(B(:)>0),3,100,tp1(1:3));
          WMr    = Br & Tr>mean(tp1(2:3)) & Tr<tp1(4);              
          tp2(3) = peak(Tr(WMr(:)));
          
          GMr    = Br & Tr>mean(tp1(1:2)) & Tr<mean(tp1(2:3));      
          tp2(2) = peak(Tr(GMr(:)));
          
          CMr    = Br & Tr<mean(tp1(2));                            
          tp2(1) = peak(Tr(CMr(:)));
          tp2(1) = min(tp2(1),tp2(2)-diff(tp2(2:3)));
          tp2(4) = min(1.5,tp2(3)+diff(tp2(2:3)));

          tp2=tp1;
          tp2=tp2*tpa;
        end
      else
        tp2 = varargin{1};
      end
      
      % Intensity Normalized Image
      if size(tp2,2)<4,  tp2(4) = min(tp2(3).*1.5,tp2(3)+diff(tp2(2:3))); end % add WM+
      BG=true(size(T)); BG(3:end-2,3:end-2,3:end-2)=0; 
      tp2=[cat_stat_nanmean(T(BG(:))),tp2];                               % add background
      T = T - tp2(1); tp2=tp2-tp2(1);
      %{
      try
        imax = double(intmax('uint16')); isc=1000; 
        lut  = [-1/3 0:1/3:10]; lut2=[-1/3 tp2(1:4) tp2(5)+(0:1/3:10)];
        lut2 = lut2(1:numel(lut));
        TI   = single(intlut(uint16(isc*T),uint16(isc * ...
                 interp1(lut2,lut,0:1/isc:imax/isc,'pchip'))))/isc;
        clear lut lut2 imax isc;  
      catch %#ok<CTCH>
      % if inlut is not available... only temporarly!!!
      % ################################################################
     
      % ################################################################
      end
      end
      %}
        TI = T; 
        isc = 2;
        tp2 = interp1(tp2,1:1/isc:5,'pchip');
        for i=2:numel(tp2)
          M = T>tp2(i-1) & T<=tp2(i);
          TI(M(:)) = (i-2)/isc/3 + (T(M(:)) - tp2(i-1))/diff(tp2(i-1:i))/isc/3;
        end
        M  = T>=tp2(end); 
        TI(M(:)) = numel(tp2)/isc/3 + (T(M(:)) - tp2(i))/diff(tp2(end-1:end))/isc/3;
      
      varargout{1} = tp2;
    otherwise
  end
end


function p=peak(T,ss)
  if ~exist('ss','var'), ss=0.01; end
  H=hist(T(:),0:ss:2.00); H=smooth(H,20); %'rloess',20);
  [v,p]=max(H(:)); p=p*ss;
end

function HS = smooth(H,span)
  window = ones(span,1)/span; 
  HS = convn(H,window,'same');
end

%{ 
% working version
function [TI,T3th,T3r] = cat_vol_CGWscale(T,G,B,Gth,vx_vol)
  %% T=cat_vol_smooth3X(T,0.5);

  % Los gehts mit dem WM, das wir ja nun schon durch die Biaskorrektur
  % kennen und das sich um 1.00 bewegt.
  T7{6}   = cat_vol_morph(B & G<Gth/2 & T>0.9 & T<=1.3,'l');        
  T7th(6) = min(1.3,max(0.95,cat_stat_nanstat1d(T(T7{6}),'median'))); 

  % Nun können wir uns die Kannte zum GM anschauen, die sich in D<2 zum 
  % knapp unterhalb des WM. Da diese beide Seiten gleichartig betrifft, 
  % sollte man so recht gut den GWM peak bestimmen können. Das da teil-
  % weise bissel CSF mit dabei sollst schnuppe sein.
  T7{5}   = (G>Gth & G<0.5 & B & cat_vol_morph(T>T7th(6),'d',1) & T<T7th(6));
  T7th(5) = max(0.75*T7th(6),min(0.95*T7th(6),cat_stat_nanstat1d(T(T7{5}),'median'))); 
  T7th(7) = 2*T7th(6) - T7th(5);

  % Nun zum GM das wir durch 
  T7{4}   = G>Gth/2 & G<0.5 & B & T>0.3 & T<T7th(5);
  T7th(4) = max(0.5*T7th(6),min(0.90*T7th(6),cat_stat_nanstat1d(T(T7{4}),'median'))); 
  T7{4}   = G<Gth*2 & B & T>T7th(4)-std(T(T7{4}(:))) & T<(T7th(4)*0.2 + T7th(6)*0.8);
  T7th(4) = max(0.5*T7th(6),min(0.90*T7th(6),cat_stat_nanstat1d(T(T7{4}),'median'))); 
   
  T7{3}   = (G>Gth*2 & G<1 & B & cat_vol_morph(T<T7th(4)*0.95 & B,'d',1) & T<T7th(4)*0.95);
  T7th(3) = max(0.25*T7th(6),min(0.70*T7th(6),cat_stat_nanstat1d(T(T7{3}),'median'))); 
  
  % CSF
  T7{2} = (G<Gth/2 | (G.*T<Gth*2 & B & T<T7th(2))) & B & T<(T7th(4)/2);
  if sum(T7{2}(:))>10
    T7{2}   = cat_vol_morph(T7{2},'l');
    T7th(2) = max(0.1*T7th(6),min(0.60*T7th(6),cat_stat_nanstat1d(T(T7{2}),'median'))); 
  else
    T7th(2) = T7th(3) - T7th(4)-diff(T7th([4,6]));
  end
  T7th(1) = T7th(2)/2;
  
  % Als dem WM glaub ich, beim GM und CSF ist es schwieriger.
  % Ist viel CSF da (Ventrikel), ist es meist eine brauchere Information
  % als das GM da dieses halt durch den PVE und Krankheiten in
  % Mitleidenschaft gezogen sein kann. Im gesunden fall haben wir aber
  % meist zu wenig CSF als das man damit was machen könne, weshalb die
  % Abschätzung des CSFs über den GM Wert dann schon wieder besser ist.
  T3th(3) = T7th(6);
  T3th(2) = T7th(4);
  T3th(1) = max(0,min(T7th(2),T7th(4)-diff(T7th([4,6]))));
  
  TI   = cat_vol_iscale(T,'gCGW',vx_vol,T3th); 
  p0T  = max(B,min(3,round(TI*3).*B));
  for i=2:2, T3th(i) = median(T(p0T==i)); end
  TI   = cat_vol_iscale(T,'gCGW',vx_vol,T3th); 
  T3r  = T7(2:2:end);
end
function p=peak(T,ss)
  if ~exist('ss','var'), ss=0.01; end
  H=hist(T(:),0:ss:2.00); H=smooth(H,20); %'rloess',20);
  [v,p]=max(H(:)); p=p*ss;
end

function HS = smooth(H,span)
  window = ones(span,1)/span; 
  HS = convn(H,window,'same');
end
%}
