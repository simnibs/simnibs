function cat_vol_correct_slice_scaling(varargin)
% ______________________________________________________________________
% Correction of slice artifacts for a set of images.
%
% WARNING: This kind of artifacts are untypical for 98% of the data!
%          Apply this filter only for images with slice artifacts, 
%          because this filter can introduce inhomogeneities!
%
% To control the filter direction and number of iterations we test the 
% if one gradient is untpyical high. If this is the case the filter is
% applyed for this direction as long as all gradient get more similar.
% In most cases only 1-3 interations are necessary.
% The filter use information from the foreground (object/tissue) to 
% estimate the correction filed. Background information are used to 
% stabilize the estimation of the WM threshhold. 
%
% This is only a slice-bias correction and the image will include   
% inhomogeneity after correction that requires further correction by 
% standard approaches like N3 or the SPM and CAT proprocessings!
%
%   cat_vol_correct_slice_scaling(job)
% 
%   job.data
%   job.prefix      .. file prefix (default = 'slicecorr_')
%   job.s           .. filtersize (default = 12)
%   job.iter        .. maximum number of iterations (default = 5);
%   job.verb        .. display progress in the command window
%   job.lb          .. lower intensity boundary
%   job.ub          .. upper intensity boundary
%
% ______________________________________________________________________
% Robert Dahnke 
% Structural Brain Mapping Group
% University Jena
%
% $Id: cat_vol_correct_slice_scaling.m 1118 2017-03-17 15:57:00Z gaser $
% ______________________________________________________________________

% ______________________________________________________________________
% Further comments: 
% - There is a private test script 'cat_tst_BWPsliceartifact.m' that use
%   BWP (brain web phantom) or real data to simulate slice artifacts for
%   each direction. 
% - It is also possible to use the background for correction, but due to
%   the low values the division by zero is critical. It is possible to 
%   smooth i and i+1 separate, but this will increase calculatin time by
%   factor 3 and the improvements are only a little. 
% - In images with a high intensity low gradient tissues higher than the 
%   WM like FAT were overcorrected, although the intensity boundies were
%   correct. I expect that this depend on the slicewise correction.
% - The most problematic cases is given by images, were the slice 
%   artifact is in the anatomic x-direction, because the slices around 
%   the corpus callosum have untpyical high GM tissue. This slices trend 
%   to be brighter than they should. Using the background as control 
%   parameter for the WM threshold helps to reduce this problem. 
% - The boundaries have to be wide to allow also corrections of block 
%   artifacts like in IBSR1:5_8 were slice 1-x are normal and x+1 to end 
%   are 50% brighter.
% 
% - Other methods:
%   - Slice correction also by 'caret_command -volume-bias-correction'?
%     - not full automatic, no direction / interation criteria
%   - SLED log correction?
%   - 

% ______________________________________________________________________



% check / set intput
  spm_clf('Interactive'); 
 
  
  % data
  if nargin == 0 
    job.data = cellstr(spm_select([1 Inf],'image','select images to filter'));
  else
    job = varargin{1};
  end
  if ~isfield(job,'data') || isempty(job.data)
    job.data = cellstr(spm_select([1 Inf],'image','select images to filter'));
  else
    job.data = cellstr(job.data);
  end
  if isempty(job.data), return; end

 
  % filter size
  if ~isfield(job,'s') 
    %job.s = min(32,max(8,spm_input('filter size? (8-32)','+1','i',12,1)));
    job.s = 12;
  end  
  
  if ~isfield(job,'method')
    job.method = 1;
  end
  
  
  % maximum iteration number
  if ~isfield(job,'iter')
    %job.iter = min(10,max(1,spm_input('maximum iterations? (1-10)','+1','i',3,1)));
    job.iter = 10;
  end
  
  
  % tissue limits for WM==1
  if ~isfield(job,'lb')
    lb  = 0.5;
  else
    lb  = job.lb;
  end
  if ~isfield(job,'ub')
    ub  = 2.0;
  else
    ub  = job.ub;
  end
  
  
  % force at least one iteration in one direction processing
  if ~isfield(job,'force')
    force = 3;
  else
    force = job.force;
  end
  
  
  % file prefix
  if ~isfield(job,'prefix') 
    job.prefix = spm_input('file prefix','+1','s',...
      sprintf('slicecorrS%02.0f_',job.s),1);
  end
  if strcmp(job.prefix,'auto')
    job.prefix = sprintf('slicecorrF%dS%02.0f_',force,job.s);
  end
  
  
  % processing information
  if ~isfield(job,'verb')
    job.verb = 1;
  end
  
  
% start processing ...  
  V = spm_vol(char(job.data));
 
  spm_progress_bar('Init',numel(job.data),'Slice-Filtering','Volumes Complete');
  if job.verb, fprintf('Correct Slice Scaling:\n'); stime0=clock; end
  for j = 1:length(V)
    if job.verb, fprintf('  %s:',V(j).fname); stime=clock; end
    
    Yo = single(spm_read_vols(V(j))); Y=Yo;
    [gx,gy,gz] = cat_vol_gradient3(Y); 
          
    % estimate slice error and create a gradient map Yg (tissue) and Ygb
    % (background) without the critical direction 
    calc1 = cat_stat_nanmean([ abs(gx(:)),abs(gy(:)),abs(gz(:))]) / ...
            cat_stat_nanmean([ abs(gx(:));abs(gy(:));abs(gz(:))]);
    calc  = (calc1 > 1.01) & (calc1==max(calc1));
    
    % intern noise estimation and correction
    signal = cat_stat_nanmedian(Y(Y(:)>cat_stat_nanmedian(Y(:))));
    signal = cat_stat_nanmedian(Y(Y>cat_stat_nanmedian(Y(Y(:)>signal & Y(:)<2*signal))));  
    noise  = cat_vol_localstat(Y,Y<signal*0.2,1,4); noise = median(noise(noise(:)>0));
    NSR    = double(noise / signal); sx = repmat(NSR * 5,1,3); sx(calc==1) = 0;
    Y = double(Y); spm_smooth(Y,Y,sx); Y = single(Y);
    
    %%
    [gx,gy,gz] = cat_vol_gradient3(Y);
    Yg   = ((calc1(1)<1.05) .* abs(gx) + (calc1(2)<1.05) .* abs(gy) + (calc1(3)<1.05) .* abs(gz)) .* ...
            3/(sum(calc1<1.05)) ./ max(eps,Y);       
    Ygb  = ((calc1(1)<1.05) .* abs(gx) + (calc1(2)<1.05) .* abs(gy) + (calc1(3)<1.05) .* abs(gz)) .* ...
            3/(sum(calc1<1.05)) ./ max(signal*0.3,Y);    
          
    % estimate thresholds for Yg, Ygb, and Y (WM threshold)      
    gth  = cat_stat_nanmean(Yg(Yg(:)<cat_stat_nanmedian(Yg(:))));
    Yth  = cat_stat_nanmean(Y(Y>cat_stat_nanmean(Y(Yg(:)>0 & Yg(:)<gth))));                       
    Ythw = cat_stat_nanmean(Y(Y>cat_stat_nanmean(Y(Yg(:)>0 & Yg(:)<gth & Y(:)>Yth*0.5))));  

    % test segments
    Ywm  = smooth3(Yg<gth & Y>Ythw*0.2)>0.5;  
    Ybg  = smooth3(Ygb./Yg<1.0 & Y<Ythw*NSR*4 & Yo>0)>0.5;  %smooth3(Y<Ythw*0.3*NSR2 & Ygb<Ythw*0.3*NSR2)>0.5; 

    Ythb = cat_stat_nanmean(Y(Ybg(:))); 
    %Ythw = cat_stat_nanmedian(Y(Ywm(:))); 
     
    % estimate slice error to find the affected slice  
    calc1 = cat_stat_nanmean(abs([ gx(Ywm(:)),gy(Ywm(:)),gz(Ywm(:)) ])) / ...
            cat_stat_nanmean(abs([ gx(Ywm(:));gy(Ywm(:));gz(Ywm(:)) ]));
    calc3 = calc1; calc5 = calc1; 
    calc  = (calc1 > 1.01) & (calc1==max(calc1));
    clear gx gy gz;

    if force
      calcdir = force;
      calc    = [0 0 1];
      calc1   = ones(1,3); calc1(force) = inf;
      calc3   = calc1;
      calc5   = calc1; 
    else
      calcdir = find(calc==1);
    end
    
    % adapt boundaries for artifact strenth
    lb  = max(0.3,lb  ./ max(calc1));
    ub  = max(2.0,ub  .* max(calc1));

    % ds('l2','',[1 1 1],Y,Ywm + (1-Ybg),Yg,Ywm,90)
    
    %%
    if any(isnan(calc)) || all(calc==0)
      fprintf('\n  No slice artifact detected! No correction! \n'); 
      desc = 'not-slicecorreted';
    else
      if job.verb==2;  xyz='xyz'; fprintf('\n\t%s: %0.3f',xyz(calc==1),mean(abs(calc3-1))); end
      Y1 = Y;
      Y2 = Yo;
      it = 0;



      while it < job.iter
        it  = it + 1;
        Yf  = Y1;
        Yof = Y2;
        switch calcdir
          case 1 % x-direction
            for i = 2:V(j).dim(1) % forward
              [Yf(i,:,:),Yof(i,:,:)] = correctslice(Yof(i,:,:),Yf(i,:,:),Yf(i-1,:,:),Ybg(i,:,:),Ybg(i-1,:,:),lb,ub,Ythw,Ythb,job.s,job.method,1);
            end   
            % global intensity correction
            Yof = Yof ./ (mean(Yf(Ywm(:)))./mean(Y1(Ywm(:))));
            Yf  = Yf  ./ (mean(Yf(Ywm(:)))./mean(Y1(Ywm(:))));
            for i = V(j).dim(1)-1:-1:1 % backward
              [Yf(i,:,:),Yof(i,:,:)] = correctslice(Yof(i,:,:),Yf(i,:,:),Yf(i+1,:,:),Ybg(i,:,:),Ybg(i+1,:,:),lb,ub,Ythw,Ythb,job.s,job.method,1);
            end   
          case 2 % y-direction 
            for i = 2:V(j).dim(2)
              [Yf(:,i,:),Yof(:,i,:)] = correctslice(Yof(:,i,:),Yf(:,i,:),Yf(:,i-1,:),Ybg(:,i,:),Ybg(:,i-1,:),lb,ub,Ythw,Ythb,job.s,job.method,2);
            end   
            Yof = Yof ./ (mean(Yf(Ywm(:)))./mean(Y1(Ywm(:))));
            Yf  = Yf  ./ (mean(Yf(Ywm(:)))./mean(Y1(Ywm(:))));
            for i = V(j).dim(2)-1:-1:1
              [Yf(:,i,:),Yof(:,i,:)] = correctslice(Yof(:,i,:),Yf(:,i,:),Yf(:,i+1,:),Ybg(:,i,:),Ybg(:,i+1,:),lb,ub,Ythw,Ythb,job.s,job.method,2);
            end   
          case 3 % z-direction 
            for i = 2:V(j).dim(3)
              [Yf(:,:,i),Yof(:,:,i)] = correctslice(Yof(:,:,i),Yf(:,:,i),Yf(:,:,i-1),Ybg(:,:,i),Ybg(:,:,i-1),lb,ub,Ythw,Ythb,job.s,job.method,3);
            end 
            Yof = Yof ./ (mean(Yf(Ywm(:)))./mean(Y1(Ywm(:))));
            Yf  = Yf  ./ (mean(Yf(Ywm(:)))./mean(Y1(Ywm(:))));
            for i = V(j).dim(3)-1:-1:1
              [Yf(:,:,i),Yof(:,:,i)] = correctslice(Yof(:,:,i),Yf(:,:,i),Yf(:,:,i+1),Ybg(:,:,i),Ybg(:,:,i+1),lb,ub,Ythw,Ythb,job.s,job.method,3);
            end
        end
        % global intensity correction
        Yof = Yof ./ (mean(Yf(Ywm(:)))./mean(Y1(Ywm(:))));
        Yf  = Yf  ./ (mean(Yf(Ywm(:)))./mean(Y1(Ywm(:))));

        % estimate artifact size to handly iteration
        [calc2,calc4] = estimateError(Yf,Ywm,Ythw);

        if it==1, imp=0.9; else imp=1; end
        if mean(abs(calc2-1))/mean(abs(calc3-1))<imp && mean(abs(calc4-1))/mean(abs(calc5-1))<imp;
          if job.verb==2, fprintf(' > %0.3f',mean(abs(calc2-1))); end
          Y1 = Yf; Y2 = Yof;
        else
          if job.verb==2, fprintf(' > %0.3fx',mean(abs(calc2-1))); end
          if it==1, fprintf('\n  No slice artifact detected! No correction! \n'); end
          break;
        end;
        if (mean(abs(calc2-1))/mean(abs(calc3-1))>0.95) || ...
           (mean(abs(calc4-1))/mean(abs(calc5-1))>0.95) || calc4(calcdir)<1.01
          break;
        end
        calc3 = calc2; calc5 = calc4;
      end
      desc = 'slicecorreted';
    end
   
    % write result
    % ------------------------------------------------------------------
    [pth, nam, ext, num] = spm_fileparts(V(j).fname);
    Vc = V; Vc(j).fname = fullfile(pth, [job.prefix nam ext num]);
    Vc(j).descrip = sprintf('%s<%s',desc,Vc(j).descrip); 
    spm_write_vol(Vc(j),Y2); 
    
    spm_progress_bar('Set',j);
    if job.verb, fprintf('\t%6.2fs\n',etime(clock,stime)); end
  end
  spm_progress_bar('Clear');
  if job.verb, fprintf('done (%3.2f Minute(s)).\n',etime(clock,stime0)/60); end
end
function [Ythi,Ythwi] = estimateWMth(Y,Yth,lb,ub)
  Ythi  = cat_stat_nanmean(Y(Y>cat_stat_nanmean(Y(Y(:)>0))));             % object threshold (lower boundary)
  Ythwi = cat_stat_nanmean(Y(Y>cat_stat_nanmean(Y(Y(:)>Ythi))));        % WM threshold for displaying
  Ythwi = max(Ythi*lb,min(Yth*ub,Ythwi)); 
  Ythi  = Ythwi .* lb; 
end
function [cimgt,cimgo] = correctslice(imgo,imgi,imgj,imgbi,imgbj,lb,ub,Ythw,Ythb,s,method,dim)

% ds('l2','',[1 1 1],imgi,imgi,imgj/Ythwj,imgi/Ythwi,1)

  imgt  = imgi;

  imgb = mean(cat(dim,imgbj,imgbi),dim)==1; % average background
  
  % background intensity for special cases
  Ythbi = mean(imgi(imgb(:))); 
  Ythbj = mean(imgj(imgb(:))); 
  
  [Ythi,Ythwi] = estimateWMth(imgi,Ythw,lb,ub);
  [Ythj,Ythwj] = estimateWMth(imgj,Ythw,lb,ub);
  
  lb1 = 0.1;
  if Ythwi<Ythw*lb1 || Ythwj<Ythw*lb1 || Ythbi<Ythb*lb || Ythbj<Ythb*lb  % no correction 
    cimgt = imgt;
    cimgo = imgo;  
  else % default correction
    Ythwj = min(Ythbj .* (Ythw/Ythb) .* ub,max( Ythwj,Ythbj .* (Ythw/Ythb) .* lb*2)); 
    
    imgi  = max(Ythi,imgi); 
    imgj  = max(Ythj,imgj);

    img   = imgi ./ max(eps,imgj);

    cimg  = smoothslice(img,s,method,dim);
    cimg  = cimg .* (Ythwj ./ max(eps,Ythw));

    cimgt = imgt ./ max(eps,cimg);
    cimgo = imgo ./ max(eps,cimg);
  end
end
function [cimgt,cimgo] = correctslice2(imgo,imgi,imgj,imgbi,imgbj,lb,ub,Ythw,Ythb,s,method,dim)
  imgt  = imgi;
  %%
  img  = mean(cat(dim,imgj,imgi),dim);      
  imgb = mean(cat(dim,imgbj,imgbi),dim)==1; 
  
  imgbth = median(img(imgb(:))); 
  imgbth = max(Ythb*lb,min(Ythb*ub,imgbth));
  
  %imgw = img>Ythw*0.3 & img<Ythw*4;
  imgb(img(:)>imgbth*1.5)=0; imgb(img(:)<imgbth*0.2)=0;
  
  sx = repmat(s,1,3); sx(dim) = 0;
  
  [Ythi,Ythwi] = estimateWMth(imgi,Ythw,lb,ub,imgbi,Ythw/Ythb); 
  [Ythj,Ythwj] = estimateWMth(imgj,Ythw,lb,ub,imgbj,Ythw/Ythb);

  %Ythwi = mean(imgi(imgw(:)));
  %Ythwj = mean(imgj(imgw(:)));
 
  Ythbi = mean(imgi(imgb(:))); 
  Ythbj = mean(imgj(imgb(:))); 
  
  %Ythwi = max(Ythbi.*(Ythw/Ythb)*0.5,min(Ythbi.*(Ythw/Ythb)*8,Ythwi));
  %Ythwj = max(Ythbj.*(Ythw/Ythb)*0.5,min(Ythwj.*(Ythw/Ythb)*8,Ythwi));  
  
  imgix=imgi; imgix(imgb)=imgix(imgb).*(Ythwi/Ythbi); imgix = max(Ythwi*0.9,imgix);
  imgjx=imgj; imgjx(imgb)=imgjx(imgb).*(Ythwj/Ythbj); imgjx = max(Ythwj*0.9,imgjx);
  
  % spm-smoothing
  imgix=imgix-Ythwi; spm_smooth(imgix,imgix,sx); imgix=imgix+Ythwi;
  imgjx=imgjx-Ythwj; spm_smooth(imgjx,imgjx,sx); imgjx=imgjx+Ythwj;
  
  %imgi  = max(Ythi,imgi); 
  %imgj  = max(Ythj,imgj);

  %img   = imgi ./ max(eps,imgj);
  %cimg  = smoothslice(img,s,method,dim);
  %cimg  = cimg .* (Ythwj ./ max(eps,Ythw));
  
  cimg  = imgix ./ max(eps,imgjx);
 % cimg  = smoothslice(cimg,s,method,dim);
  cimg  = cimg .* (Ythwj ./ max(eps,Ythw));
  
  cimgt = imgt ./ max(eps,cimg);
  cimgo = imgo ./ max(eps,cimg);
end
function [calc2,calc4] = estimateError(Yf,Ywm,Ythw)
  [gx,gy,gz] = cat_vol_gradient3(Yf); 
  Yg    = (abs(gx) + abs(gy) + abs(gz)) ./ max(eps,Yf); 
  gth   = cat_stat_nanmedian(Yg(Yg(:)<cat_stat_nanmedian(Yg(:))));
  Ywm2  = cat_vol_morph(smooth3(Yf>Ythw*0.9 & Yf<Ythw*1.2 & Yg<gth & Yg>0)>0.5,'close',2);
  calc2 = cat_stat_nanmean([ abs(gx(Ywm(:))),abs(gy(Ywm(:))),abs(gz(Ywm(:)))]) / max(eps,...
          cat_stat_nanmean([ abs(gx(Ywm(:)));abs(gy(Ywm(:)));abs(gz(Ywm(:)))]));
  calc4 = cat_stat_nanmean([ abs(gx(Ywm2(:))),abs(gy(Ywm2(:))),abs(gz(Ywm2(:)))]) / max(eps,...
          cat_stat_nanmean([ abs(gx(Ywm2(:)));abs(gy(Ywm2(:)));abs(gz(Ywm2(:)))]));
end
function cimg = smoothslice(img,s,method,dim)
  % filtering by image reduction, smoothing and reinterpolation
  % 3d-data required 
  % not nice, but better than conv2
 if method==1
    switch dim
      case 1
        cimg = repmat(img,[3,1,1]);
        [cimg,IR1] = cat_vol_resize({cimg},'reduceV',1,s,2,'mean');
        cimg = smooth3(cimg); 
        cimg = cat_vol_resize({cimg},'dereduceV',IR1); 
        cimg = smooth3(cimg); 
        cimg = cimg(2,:,:); 
      case 2
        cimg = repmat(img,[1,3,1]);
        [cimg,IR1] = cat_vol_resize({cimg},'reduceV',1,s,2,'mean');
        cimg = smooth3(cimg); 
        cimg = cat_vol_resize({cimg},'dereduceV',IR1); 
        cimg = smooth3(cimg); 
        cimg = cimg(:,2,:); 
      case 3
        cimg = repmat(img,[1,1,3]); 
        [cimg,IR1] = cat_vol_resize({cimg},'reduceV',1,s,2,'mean');
        cimg = smooth3(cimg); 
        cimg = cat_vol_resize({cimg},'dereduceV',IR1); 
        cimg = smooth3(cimg); 
        cimg = cimg(:,:,2); 
    end
 elseif method==2 % spm-smoothing approach - bad boundary properies, even if I correct for the mean intensity
    sx = repmat(s,1,3); sx(dim) = 0; ofs = mean(img(:));
    cimg = double(img-ofs); spm_smooth(cimg,cimg,sx); cimg = single(cimg+ofs);
 else % christian all smoothing approach - not realy smooth
    x = [-s:s];
    x = exp(-(x).^2/(2*(s).^2));
    x = x/sum(x);
    cimg = conv2(img,x'*x,'same');
  end
end