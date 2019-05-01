function varargout=cat_tst_calc_kappa(P,Pref,opt)
% ______________________________________________________________________
% Estimates Kappa for a set of input images P to one or an equal number of
% reference image(s) Pref. Use realgignment if image properties does not 
% match.
%
%??[txt,val] = cat_tst_calc_kappa(P,Pref,opt)
% 
% P             .. list of images
% Pref          .. ground truth segmentation
% opt 
%  .methodname  .. just to display            (default = datasubpath)
%  .verb        .. verbose level              (default = 1)
%  .realign     .. force realignment          (default = 0)
%  .realignres  .. resolution of realignment  (default = 1.5)
%  .diffimg     .. write difference image     (default = 0)
%  .testcase    .. evalution of specific label maps (default = 'auto')
%  .finishsound .. bong                       (default = 0)
%                   ...
%  .spaces      .. length of filename field   (default = 50)               
% ______________________________________________________________________
% based on cg_calc_kappa by Christian Gaser
%
% Robert Dahnke 
% Structural Brain Mapping Group
% University Jena
%
% $Id: cat_tst_calc_kappa.m 1136 2017-06-11 14:04:58Z dahnke $
% ______________________________________________________________________
%#ok<*AGROW>
%#ok<*ASGLU>

% ______________________________________________________________________
% ToDo:
%  - opt.fname .. csv-result file % not yet
%  - single segment comparison with input of the first image c1 or p1
%  - csv file export
%  - xml file update / export
%  - update for 4 class p0 case of CAT12
%  - interpolate to gt resolution
%  
%  - slicewise and tissue class based evaluation
% ______________________________________________________________________


% set defaults, get files:
%  spm('Defaults','FMRI');
  
  
% initialize output
  txt = '';
  val = struct('fname','','path','','name','', ...
              'BE',struct('kappa',[],'accuracy',[],'FP','','FN','', ...
                   'sensit_all',[],'sensit',[],'specif',[],'dice',[],'jaccard',[]),... 
              'SEG',struct('kappa',[],'rms',[],'kappaGW',[],'rmsGW',[]));
  if nargout>0, varargout{1}=''; end
  if nargout>1, varargout{2}={''}; end
  if nargout>2, varargout{3}=val; end 
 
  
% first input - test data
  if ~exist('P','var')
    P = spm_select(Inf,'image','Select images to compare'); 
  else
    if isa(P,'cell'), if size(P,1)<size(P,2), P=P'; end; P=char(P); end
  end
  if isempty(P), return; end
    
  
% second input - ground truth
  V = spm_vol(P);
  n = numel(V);   % number of test cases
  if ~exist('Pref','var')
    Pref = spm_select([1 n],'image','Select reference mask');
    Vref = spm_vol(Pref); 
  else
    Pref = cellstr(Pref);
    if size(Pref,1)<size(Pref,2), Pref=Pref'; end; 
    Vref = spm_vol(char(Pref));
  end
  if isempty(V) || isempty(Vref), return; end 
  
  
% default parameter settings  
  if ~exist('opt','var'), opt=struct(); end
  def.methodname  = ['(' spm_str_manip(spm_str_manip(P(1,:),'h'),'l20') ')'];
  def.verb        = 2;
  def.realign     = 0;
  def.realignres  = 1.5; 
  def.diffimg     = 0;
  def.testcase    = 'auto';
  def.spaces      = 70; 
  def.finishsound = 0; 
  def.allkappa    = 1; 
  opt = cat_io_checkinopt(opt,def);
  
  
% check how we can compare the images:
  vol  = single(spm_read_vols(Vref(1))); 
  switch opt.testcase
    case 'slices'
      h=hist(vol(:),0:255);
      ncls = min(3,sum(h>0)-1);
    case 'binary'
      ncls = 1; 
    case 'IBSR'
      ncls = 3;
    case 'p03'
      ncls = 3; 
    case 'p04',
      ncls = 4;
    otherwise
      ncls = max(round(vol(:))); 
      if     ncls==255, ncls=1; 
      elseif ncls==254, ncls=3; % IBSR
      elseif ncls==4,   ncls=3; % default 3-class label images with CSF,GM and WM (and
      end
  end
  clear vol;
  
% rating system and color output 
  MarkColor   = cat_io_colormaps('marks+',40); 
  setnan      = [0 nan];
  evallinearb = @(x,best,worst,marks) min(marks,max( 1,(abs(best-x) ./ ...
    abs(diff([worst,best]))*(marks-1)+1))) + setnan(isnan(x)+1); 
  estr = sprintf('%s\n%s\n\n',spm_str_manip(P(1,:),'h'),Vref(1).fname);

  Yp0toC  = @(Yp0,c) 1-min(1,abs(Yp0-c));
% loop  
  for nc=1:(ncls>1 && nargout>2)+1  
    
    
  % create header of output table 
    if opt.verb>1
      fprintf('cat_tst_calc_kappa with %d classes.\n',ncls);
    end
    switch ncls
      case 0, txt{1}='Error ground truth empty!'; continue
      case 1, tab = {['File ' sprintf(sprintf('%%%ds',opt.spaces-4),opt.methodname)],...
                    'kappa','jaacard','dice','sens.','spec.','FP(F)','FN(N)','N/(P+N)','RMS'};
              txt{1} = sprintf(sprintf('\\n%%%ds%%6s%%8s%%8s%%8s%%8s%%8s%%8s%%8s%%8s%%8s\\n',opt.spaces),...
                estr,tab{1},tab{2},tab{3},tab{4},tab{5},tab{6},tab{7},tab{8},tab{9},tab{10});
              k = zeros(n,9);
      case 3, tab = {['File ' sprintf(sprintf('%%%ds',opt.spaces-4),opt.methodname)],...
                    'K(C)','K(G)','K(W)','K(CGW)','K(B)','RMS(C)','RMS(G)','RMS(W)','RMS(CGW)','RMS(B)'};
              txt{1} = sprintf(sprintf('\\n%%%ds%%s%%8s%%8s%%8s%%8s%%8s |%%8s%%8s%%8s%%8s%%8s\\n',opt.spaces),...
                estr,tab{1},tab{2},tab{3},tab{4},tab{5},tab{6},tab{7},tab{8},tab{9},tab{10},tab{11}); 
              k = zeros(n,10);
    end
    txt{2} = ''; 
    if opt.verb && ~isempty(txt{1}) && opt.verb>1, fprintf(txt{1}); end

    
  % data evaluation
    for i=1:n 
      %% for all test cases
      [pth, name] = fileparts(V(i).fname); 
      val(i).fname = V(i).fname;
      val(i).path  = pth;
      val(i).name  = name;
      fnamestr     = [spm_str_manip(pth,sprintf('k%d',max(0,min(floor(opt.spaces/3),opt.spaces-numel(name)-1)-1))),'/',...
                      spm_str_manip(name,sprintf('k%d',opt.spaces - floor(opt.spaces/3)))];
      
      % if only one ground-truth image is give use this, otherwise their
      % should be a gound-truth for each image
      if numel(Vref)==numel(V), Vrefi=i; else Vrefi=1; end                    
      if numel(Vref)==numel(V) || i==1
        vol1 = single(spm_read_vols(Vref(Vrefi)));
      end  
      vx_vol  = sqrt(sum(Vref(Vrefi).mat(1:3,1:3).^2)); 
      
      % realginment
      if any(V(i).dim ~= Vref(Vrefi).dim) %|| any(V(i).mat(:) ~= Vref(Vrefi).mat(:)) 
        [pp,ff,ee] = spm_fileparts(V(i).fname);
        if opt.realign
          Vir = [tempname ee];
          try
            copyfile(V(i).fname,Vir);
            spm_realign(char([{Vref(Vrefi).fname};{Vir}]),...
              struct('sep',opt.realignres,'rtm',0,'interp',4,'graphics',0,'fwhm',opt.realignres));
            fprintf(repmat('\b',1,73*3+1))
          catch
            delete(Vir);
          end  
        else
          Vir = V(i).fname; 
        end
        [V(i),vol2] = cat_vol_imcalc(Vir,Vref(Vrefi),'i1',struct('interp',6,'verb',0));
      else
        vol2 = single(spm_read_vols(V(i)));
      end

      
      %% ds('l2','',1,vol2/ncls,vol1/ncls,vol2/ncls,vol1/ncls,126)
      switch opt.testcase
        case 'slices'
          %%
          vol1 = round(vol1); 
          %vol2 = round(vol2); 
          
          xsum = shiftdim(sum(sum(vol1,2),3),1); xslices = find(xsum>max(xsum)*.1);
          ysum = sum(sum(vol1,1),3);             yslices = find(ysum>max(ysum)*.1);
          zsum = shiftdim(sum(sum(vol1,1),2),1); zslices = find(zsum>max(zsum)*.1);
          mask = false(size(vol1));
          for xi=1:numel(xslices), mask(xslices(xi),:,:) = true; end
          for yi=1:numel(yslices), mask(:,yslices(yi),:) = true; end
          for zi=1:numel(zslices), mask(:,:,zslices(zi)) = true; end
          %%
          switch ncls
            case 1
              [kappa_all, kappa, accuracy_all, accuracy, sensit_all, sensit, specif, confusion, dice, jaccard] = ...
                cg_confusion_matrix( uint8(vol1(mask(:))>0) + 1, uint8( round(vol2(mask(:))) == max(vol1(:))-1) + 1, 2);
              
              % rms for GM class
              rms    = sqrt( cat_stat_nanmean( ( (vol1(mask(:))>0) - Yp0toC(vol2(mask(:)),2) ).^2 ) );

              FP     = confusion(1,2); FN = confusion(2,1);
              k(i,:) = [kappa_all,jaccard(1),dice(1),sensit(1),sensit(2),FP,FN,FN/(FN+FP),rms];
              txti   = sprintf(sprintf('%%%ds:%%8.4f%%8.4f%%8.4f%%8.4f%%8.4f%%8.0f%%8.0f%%8.4f%%8.4f\\n',opt.spaces),...
                fnamestr,k(i,:)); 
              
              val(i).BE  = struct('kappa',kappa_all,'accuracy',accuracy_all, ...
                           'FP',FP,'FN',FN, ...
                           'sensit_all',sensit_all,'sensit',sensit(1),'specif',specif(1),'dice',dice(1),'jaccard',jaccard(1),'rms',rms); 
              colori = mean(kappa_all);
            case 3
              vol1o = vol1; vol1 = (vol1o==1) + (vol1o==3)*2 + (vol1o==6)*3 + (vol1o==5)*4;
%%
              if opt.allkappa
                [kappa_all,kappa] = cg_confusion_matrix( uint8(round(vol1(mask(:))+1)) ,uint8(round(vol2(mask(:))+1)), 4); 
                kappa_all = [kappa(2:4)' kappa_all kappa(1)]; 
              else
                for c=1:2, kappa_all(1,c) = cg_confusion_matrix(uint8((round(vol1(mask(:)))==c)+1),uint8((round(vol2(mask(:)))==c)+1), 2); end
                c=3;       kappa_all(1,c) = cg_confusion_matrix(uint8((round(vol1(mask(:)))==c)+1),uint8((round(vol2(mask(:)))>=c)+1), 2); 
                bth=0.5;   kappa_all(1,5) = cg_confusion_matrix(uint8((vol1(mask(:))>=bth)+1     ),uint8((vol2(mask(:))>=bth)+1     ), 2); 
                kappa_all(1,4) = mean(kappa_all(1,1:3)); 
              end
              % rms 
              rms = calcRMS(vol1(mask(:)),vol2(mask(:))); rms = [rms(1:3) mean(rms(1:3)) rms(4)];
              k(i,:) = [kappa_all,rms];
              txti   = sprintf(sprintf('%%%ds:%%8.4f%%8.4f%%8.4f%%8.4f%%8.4f |%%8.4f%%8.4f%%8.4f%%8.4f%%8.4f\\n', ...
                opt.spaces),fnamestr,k(i,:)); 

              val(i).SEG = struct('kappa',kappa_all(1:3),'rms',rms(1:3),'kappaGW',kappa_all(4),'rmsGW',rms(4));
             
              %val(i).BE  = struct('kappa',kappa_all,'accuracy',accuracy_all, ...
              %             'FP',FP,'FN',FN, ...
              %             'sensit_all',sensit_all,'sensit',sensit(1),'specif',specif(1),'dice',dice(1),'jaccard',jaccard(1),'rms',rms); 
              
              colori = mean(kappa_all(1,2:3)); %  colori = kappa_all(4);
          end
        otherwise
          %%



          if opt.diffimg
            [pp,ff,ee] = spm_fileparts(V(i).fname);
            [pr,fr]    = spm_fileparts(Vref(Vrefi).fname);

            Vd = Vref(Vrefi);
            %Vd.fname = fullfile(pp,['diffimg.' strrep(ff,'-','') '-' strrep(fr,'-','') '.' genvarname(strrep(opt.methodname,'/','-')) ee]);
            Vd.fname = fullfile(pr,['diffimg.' ff '.' fr '.' ...
              strrep(genvarname(strrep(['XNT',opt.methodname],'/','-')),'XNT','') ee]);
            spm_write_vol(Vd,vol2-vol1);
          end




          %% class-based evaluation
          switch ncls
            case 1
              %if length(Vref)==n,  vol1 = spm_read_vols(Vref(i))/255+1;
              %else                 vol1 = spm_read_vols(Vref(i));
              %end

              maxv=max((vol1(:))); if maxv==255, vol1=vol1/maxv; else vol1=vol1/maxv; end

              [kappa_all, kappa, accuracy_all, accuracy, sensit_all, sensit, specif, confusion, dice, jaccard] = ...
                cg_confusion_matrix(uint8((round(vol1(:))>0)+1), uint8((round(vol2(:))>0)+1), 2);

              rms    = sqrt( cat_stat_nanmean( ( ( vol1(:) - vol2(:) ).^2 ) ));
              
              FP     = confusion(1,2); FN = confusion(2,1);
              k(i,:) = [kappa_all,jaccard(1),dice(1),sensit(1),sensit(2),FP,FN,FN/(FN+FP),rms];
              txti   = sprintf(sprintf('%%%ds:%%8.4f%%8.4f%%8.4f%%8.4f%%8.4f%%8.0f%%8.0f%%8.4f%%8.4f\\n',opt.spaces),...
                fnamestr,k(i,:)); 

              val(i).BE  = struct('kappa',kappa_all,'accuracy',accuracy_all, ...
                           'FP',FP,'FN',FN, ...
                           'sensit_all',sensit_all,'sensit',sensit(1),'specif',specif(1),'dice',dice(1),'jaccard',jaccard(1)); 
              colori = mean(kappa_all);
            case 3
              maxv=max((vol1(:))); 
               switch opt.testcase
                 case 'IBSR'
                   vol1=round(vol1);
                   if maxv>3
                     vol1=round(((vol1-64)/(maxv-64)) * 3);
                   end
                 otherwise
                   if maxv>4, vol1=round(vol1); vol1=vol1/maxv*3; end
               end
               
              if 0
                % temporare
                % bei dem BWP test bei fsl gibts einen ungekl??rten versatz
                if ~isempty(strfind(upper(V(i).fname),'FSL')) && ~isempty(strfind(upper(V(i).fname),'BWP'))
                  vol2(2:end,:,:)=vol2(1:end-1,:,:);
                end
              end

              if opt.allkappa
                [kappa_all,kappa] = cg_confusion_matrix( uint8(round(vol1(:)+1)) ,uint8(round(vol2(:)+1)), 4); 
                kappa_all = [kappa(2:4)' kappa_all kappa(1)]; 
              else
                for c=1:2, kappa_all(1,c) = cg_confusion_matrix(uint8((round(vol1(:))==c)+1),uint8((round(vol2(:))==c)+1), 2); end
                c=3;       kappa_all(1,c) = cg_confusion_matrix(uint8((round(vol1(:))==c)+1),uint8((round(vol2(:))>=c)+1), 2); 
                bth=0.5;   kappa_all(1,5) = cg_confusion_matrix(uint8((vol1(:)>=bth)+1     ),uint8((vol2(:)>=bth)+1     ), 2); 
                kappa_all(1,4) = mean(kappa_all(1,1:3)); 
              end
              
              rms = calcRMS(vol1,vol2); rms = [rms(1:3) mean(rms(1:3)) rms(4)];
              k(i,:) = [kappa_all,rms];
              txti   = sprintf(sprintf('%%%ds:%%8.4f%%8.4f%%8.4f%%8.4f%%8.4f |%%8.4f%%8.4f%%8.4f%%8.4f%%8.4f\\n', ...
                opt.spaces),fnamestr,k(i,:)); 

              val(i).SEG = struct('kappa',kappa_all(1:3),'rms',rms(1:3),'kappaGW',kappa_all(4),'rmsGW',rms(4));
              switch opt.testcase
                case 'IBSR'
                  colori = mean(kappa_all(2:3));
                otherwise
                  colori = kappa_all(4);
              end
            otherwise
              if numel(Vref)==numel(V), Vrefi=i; else Vrefi=1; end
              vol1 = single(spm_read_vols(Vref(Vrefi))); 
              vol2 = single(spm_read_vols(V(i)));

              for c=1:ncls, kappa_all(i,c) = cg_confusion_matrix(uint8((round(vol1(:))==c)+1),uint8((round(vol2(:))==c)+1), 2); end
              colori = mean(kappa_all);
          end
          if exist('Vo','var') && exist(Vo.fname,'file')
            txti = [txti(1:end-1) 'i' txti(end)]; 
            delete(Vo.fname); 
            clear Vo;
          end
      end
      
      %%    
      if opt.verb
        if ncls==1 && ~strcmp(opt.testcase,'slices')
          cat_io_cprintf(MarkColor(round(min(40,max(1,evallinearb(colori,1.00,0.80,6)/10*40))),:),txti); 
        else
          cat_io_cprintf(MarkColor(round(min(40,max(1,evallinearb(colori,0.95,0.65,6)/10*40))),:),txti); 
        end
      end; 
      txt{2}=[txt{2} txti]; tab=[tab;[{name},num2cell(k(i,:))]]; 
    end
   
    
    %% conclustion
    if numel(n)
      switch ncls
        case 1, txt{3} = sprintf(sprintf( ...
                           ['\\n%%%ds:%%8.4f%%8.4f%%8.4f%%8.4f%%8.4f%%8.0f%%8.0f%%8.4f%%8.4f\\n', ...
                               '%%%ds:%%8.4f%%8.4f%%8.4f%%8.4f%%8.4f%%8.0f%%8.0f%%8.4f%%8.4f\\n\\n'], ...
                               opt.spaces,opt.spaces),'mean',mean(k,1),'std',std(k,1,1));
        case 3, txt{3} = sprintf(sprintf( ...
                           ['\\n%%%ds:%%8.4f%%8.4f%%8.4f%%8.4f%%8.4f |%%8.4f%%8.4f%%8.4f%%8.4f%%8.4f\\n' ...
                               '%%%ds:%%8.4f%%8.4f%%8.4f%%8.4f%%8.4f |%%8.4f%%8.4f%%8.4f%%8.4f%%8.4f\\n\\n'], ...
                               opt.spaces,opt.spaces),'mean',mean(k,1),'std',std(k,1,1));     
      end
      if opt.verb>1 && n>1, fprintf(txt{3}); end; 
      tab = [tab;[{'mean'},num2cell(mean(k,1));'std',num2cell(std(k,1,1))]];                   
    end
    
    % export
    if nc==1
      if nargout>0, varargout{1}=txt'; end
      if nargout>1, varargout{2}=tab; end
      if nargout>2, varargout{3}=val; end
    else
      if nargout>0, varargout{1}=[varargout{1};txt']; end
      if nargout>1, varargout{2}{nc}=tab; end
      if nargout>2, varargout{3}{nc}=val; end
    end
    ncls=1;
  end
  
  if opt.finishsound 
    load gong.mat; 
    soundsc(y(5000:25000),Fs)
  end
end
function rms=calcRMS(v1,v2)
% boundary box????
  v1(v1>3)=3;
  v2(v2>3)=3;
 
  for ci=1:3
    c1 = (v1-(ci-1)).* (v1>(ci-1) & v1<ci) + ((ci+1)-v1).*(v1>=ci & v1<(ci+1));
    c2 = (v2-(ci-1)).* (v2>(ci-1) & v2<ci) + ((ci+1)-v2).*(v2>=ci & v2<(ci+1));
    rms(1,ci) = sqrt(cat_stat_nanmean((c1(:)-c2(:)).^2));
  end
  
  rms(1,4) = sqrt(cat_stat_nanmean((v2(:)-v1(:)).^2));
end
function varargout = cg_confusion_matrix(reference, classified, n_class)
% compute statistic from confusion matrix
% [kappa_all, kappa, accuracy_all, accuracy, sensit_all, sensit, specif, confusion] = cg_confusion_matrix(reference, classified, n_class)

  % get sure that image is integer
  
  if nargin < 3
    n_class = max(classified);
  end

  % build confusion matrix
  confusion = zeros(n_class,n_class);
  for i = 1:n_class
    for j = 1:n_class
      confusion(i,j) =  length(find(round(reference)==i & round(classified)==j));
    end
  end

  N = sum(confusion(:));
  kappa    = zeros(size(confusion,1),1,'single');
  sensit   = zeros(size(confusion,1),1,'single');
  specif   = zeros(size(confusion,1),1,'single');
  accuracy = zeros(size(confusion,1),1,'single');

  sum_col  = sum(confusion,1);
  sum_row  = sum(confusion,2);

  Pc = 0;
  for i = 1:n_class
    sum_row_x_col = sum_row(i)*sum_col(i);

    % calculate a..d of confusion matrix
    a = confusion(i,i);
    b = sum_col(i) - a;
    c = sum_row(i) - a;
    d = N - (a + b + c);

    specif(i) = d/(b+d);
    sensit(i) = a/(a+c);
    accuracy(i) = 1-(b+c)/N;
    dice(i)     = d/(0.5*(d+d+b+c)); % Shattuck 2008, Online resource for validation of brain segmentation methods
    jaccard(i)  = d/(d+b+c);         % Shattuck 2008, Online resource for validation of brain segmentation methods

    kappa(i) = (N*confusion(i,i) - sum_row_x_col)/(N*sum_row(i) - sum_row_x_col + eps);
    Pc = Pc + sum_row_x_col/N^2;
  end

  P0 = sum(diag(confusion))/N;

  kappa_all = (P0-Pc)/(1-Pc);
  sensit_all = P0;
  accuracy_all = P0;

  varargout{1} = kappa_all;
  varargout{2} = kappa;
  varargout{3} = accuracy_all;
  varargout{4} = accuracy;
  varargout{5} = sensit_all;
  varargout{6} = sensit;
  varargout{7} = specif;
  varargout{8} = confusion;
  varargout{9} = dice;
  varargout{10} = jaccard;
end