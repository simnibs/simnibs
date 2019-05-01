function varargout=cat_vol_resize(T,operation,varargin)
% ______________________________________________________________________
% 
%   varargout=cat_vol_resize(T,operation,varargin)
%
% Examples:
%
% - Resizing to the half resolution 
%
%
% - Resizing of the image resolution to an lower and more isotropic resolution:
%   [TIr,TIGr,Br,Gr,resTr] = cat_vol_resize({TI,TIG,single(B),G./TI},'reduceV',vx_vol,2,64); 
%   TV = cat_vol_resize(TV,'dereduceV',resT);
%
% - Removing/readding of background:
%   [Tr,BB] = cat_vol_resize(T ,'reduceBrain'  ,vx_vol,10);      
%    T      = cat_vol_resize(Tr,'dereduceBrain',BB);
%
% - Interpolation for low slice resolution to obtain an more isotropic resolution
%   [Tr,aniso] = cat_vol_resize(T ,'aniso2iso',vx_vol,method);         
%   TV         = cat_vol_resize(Tr,'iso2aniso',aniso,method);  
% ______________________________________________________________________
% Robert Dahnke
% Structural Brain Mapping Group
% University Jena
% ______________________________________________________________________
% $Id: cat_vol_resize.m 1212 2017-11-10 12:39:37Z dahnke $ 
  
  if nargin==0, help cat_vol_resize; return; end
  if isempty(T), varargout{1} = T; return; end
  if ndims(T)>2, TI=T; clear T; T{1}=TI; end %else varargout{1}=T; end 
  if nargin<2, error('ERROR: cat_vol_resolution: not enought input!\n'); end

  
  switch lower(operation)
    % REDUCE & DEREDUCE
    % __________________________________________________________________
    case 'reduce'
      if numel(varargin)<1, method='linear'; else method=varargin{2}; end
      for i=1:numel(T)
        if mod(size(T{i},1),2)==1, T{i}(end+1,:,:)=T{i}(end,:,:); end
        if mod(size(T{i},2),2)==1, T{i}(:,end+1,:)=T{i}(:,end,:); end
        if mod(size(T{i},3),2)==1, T{i}(:,:,end+1)=T{i}(:,:,end); end
        [Rx,Ry,Rz] = meshgrid(single(1.5:2:size(T{i},2)),single(1.5:2:size(T{i},1)),single(1.5:2:size(T{i},3)));
        varargout{i} = cat_vol_interp3f(single(T{i}),Rx,Ry,Rz,method);
      end
        
      
    case 'dereduce'
      if numel(varargin)<2, method='linear'; else method=varargin{2}; end
      sD = varargin{1}/2+0.25;
      [Rx,Ry,Rz] = meshgrid(single(0.75:0.5:sD(2)),single(0.75:0.5:sD(1)),single(0.75:0.5:sD(3)));
      for i=1:numel(T)
        if islogical(T{i}), varargout{i} = cat_vol_interp3f(smooth3(single(imageExpand(T{i}))),Rx,Ry,Rz,method)>0.5;
        else                varargout{i} = cat_vol_interp3f(single(imageExpand(T{i})),Rx,Ry,Rz,method);
        end
      end
    
          
      
      
    % REDUCEV & DEREDUCEV
    % __________________________________________________________________
    case 'reducev'  
      sizeT=size(T{1});
      if numel(varargin)<1, vx_vol  = [1 1 1];  else vx_vol  = round(varargin{1}*100)/100; end
      if numel(varargin)<2, vx_volr = 2;        else vx_volr = round(varargin{2}*100)/100; end; 
      if numel(varargin)<3, minSize = 32;       else minSize = varargin{3}; end; 
      if numel(varargin)<4, method  = 'linear'; else method  = varargin{4}; end
      
      if numel(vx_vol)==1,  vx_vol =repmat(vx_vol ,1,3); end; 
      if numel(vx_volr)==1, vx_volr=repmat(vx_volr,1,3); end; vx_volr = max(vx_volr,vx_vol); 
      if numel(minSize)==1, minSize=min(sizeT,repmat(minSize,1,3)); end
      
      ss = floor(vx_volr ./ vx_vol); sizeTr = floor(sizeT ./ ss);
      ss = floor(sizeT ./ max(sizeTr,minSize));
      vx_volr = vx_vol .* ss; vx_red = ones(1,3)./ss;
%      vx_red  = vx_vol./min(max(vx_vol,minSize),vx_volr); vx_red  = ceil(vx_red*4)/4;
%      vx_red  = vx_red .* (1+(sizeT.*vx_red < minSize));  vx_red  = ceil(vx_red*4)/4;
%      vx_volr = vx_vol./vx_red;
      %minvoxcount = min(8,mean(ss)/2);
      minvoxcount = mean(ss);
      for i=1:numel(T)
        if any(vx_red<=0.5) % 0.65
          if mod(size(T{i},1),2)==1 && vx_red(1)<=0.75, T{i}(end+1,:,:)=T{i}(end,:,:); end
          if mod(size(T{i},2),2)==1 && vx_red(2)<=0.75, T{i}(:,end+1,:)=T{i}(:,end,:); end
          if mod(size(T{i},3),2)==1 && vx_red(3)<=0.75, T{i}(:,:,end+1)=T{i}(:,:,end); end
     
          %ss=floor([1 1 1]./vx_red); % <=0.5); % 0.65
          if strcmp(method,'max')
            varargout{i} = zeros(floor(size(T{i})./ss),'single');
            counter = varargout{i};
            nsize = floor(size(T{i})./ss).*ss;
            for ii=1:ss(1)
              for jj=1:ss(2)
                for kk=1:ss(3)
                  Tadd = T{i}(ii:ss(1):nsize(1),jj:ss(2):nsize(2),kk:ss(3):nsize(3));
                  Tadd(isnan(Tadd(:))) = 0;
                  counter = counter + (Tadd>0);
                  varargout{i} = max(varargout{i}, ...
                    (T{i}(ii:ss(1):nsize(1),jj:ss(2):nsize(2),kk:ss(3):nsize(3))>0) .* ...
                     T{i}(ii:ss(1):nsize(1),jj:ss(2):nsize(2),kk:ss(3):nsize(3)));
                    
                end
              end
            end
            varargout{i}(counter(:)<minvoxcount) = 0;
          elseif strcmp(method,'min') 
            varargout{i} = nan(floor(size(T{i})./ss),'single');
            counter = varargout{i};
            nsize = floor(size(T{i})./ss).*ss; T{i}(T{i}<eps)=nan; 
            for ii=1:ss(1)
              for jj=1:ss(2)
                for kk=1:ss(3)
                  Tadd = T{i}(ii:ss(1):nsize(1),jj:ss(2):nsize(2),kk:ss(3):nsize(3));
                  Tadd(isnan(Tadd(:))) = 0;
                  counter = counter + (Tadd>0);
                  varargout{i} = min(varargout{i}, ...
                    T{i}(ii:ss(1):nsize(1),jj:ss(2):nsize(2),kk:ss(3):nsize(3)));
                    
                end
              end
            end
            varargout{i}(counter(:)<minvoxcount) = 0;

          elseif strcmp(method,'meanm')
            varargout{i} = zeros(floor(size(T{i})./ss),'single');
            counter = varargout{i};
            nsize = floor(size(T{i})./ss).*ss;
            for ii=1:ss(1)
              for jj=1:ss(2)
                for kk=1:ss(3)
                  Tadd = T{i}(ii:ss(1):nsize(1),jj:ss(2):nsize(2),kk:ss(3):nsize(3));
                  Tadd(isnan(Tadd(:))) = 0;
                  varargout{i} = varargout{i} + Tadd;
                  counter = counter + (Tadd~=0);
                  clear Tadd;
                end
              end
            end
            varargout{i}(counter(:)<minvoxcount) = 0;
            varargout{i}(counter(:)>0) = varargout{i}(counter(:)>0) ./ counter(counter(:)>0);   
            varargout{i}(isnan(varargout{i})) = 0;
         elseif strcmp(method,'stdm')
            % mean
            meanx = zeros(floor(size(T{i})./ss),'single');
            counter = meanx;
            nsize = floor(size(T{i})./ss).*ss;
            for ii=1:ss(1)
              for jj=1:ss(2)
                for kk=1:ss(3)
                  Tadd = T{i}(ii:ss(1):nsize(1),jj:ss(2):nsize(2),kk:ss(3):nsize(3));
                  Tadd(isnan(Tadd(:))) = 0;
                  meanx = meanx + Tadd;
                  counter = counter + (Tadd~=0);
                  clear Tadd;
                end
              end
            end
            meanx(counter(:)<minvoxcount) = 0;
            meanx(counter(:)>0) = meanx(counter(:)>0) ./ counter(counter(:)>0); 
            meanx(isnan(meanx)) = 0;
            % std
            varargout{i} = zeros(floor(size(T{i})./ss),'single');
            for ii=1:ss(1)
              for jj=1:ss(2)
                for kk=1:ss(3)
                  Tadd = T{i}(ii:ss(1):nsize(1),jj:ss(2):nsize(2),kk:ss(3):nsize(3));
                  Tadd(isnan(Tadd(:))) = 0;
                  varargout{i} = varargout{i} + (Tadd - meanx).^2;
                  counter = counter + (Tadd~=0);
                  clear Tadd;
                end
              end
            end
            varargout{i}(counter(:)<minvoxcount) = 0;
            varargout{i}(counter(:)>0) = sqrt(varargout{i}(counter(:)>0) ./ counter(counter(:)>0));   
            varargout{i}(isnan(varargout{i})) = 0;
         elseif strcmp(method,'median')
            varargout{i} = zeros([floor(size(T{i})./ss),prod(ss)],'single'); 
            %zeros([floor(size(T{i})./ss),prod(size(T{i}) ./ floor(size(T{i})./ss))],'single');
            medi=1; nsize = floor(size(T{i})./ss).*ss;
            for ii=1:ss(1)
              for jj=1:ss(2)
                for kk=1:ss(3)
                  Tadd = T{i}(ii:ss(1):nsize(1),jj:ss(2):nsize(2),kk:ss(3):nsize(3));
                  Tadd(isnan(Tadd(:))) = 0;
                  varargout{i}(:,:,:,medi) = Tadd; medi=medi+1;
                  clear Tadd;
                end
              end
            end
            varargout{i} = median(varargout{i},4);
         elseif strcmp(method,'cat_stat_nanmean') || strcmp(method,'meannan')
            varargout{i} = zeros(floor(size(T{i})./ss),'single');
            counter = varargout{i};
            nsize = floor(size(T{i})./ss).*ss;
            for ii=1:ss(1)
              for jj=1:ss(2)
                for kk=1:ss(3)
                  Tadd = T{i}(ii:ss(1):nsize(1),jj:ss(2):nsize(2),kk:ss(3):nsize(3));
                  Tadd(isnan(Tadd(:))) = 0;
                  varargout{i} = varargout{i} + Tadd;
                  counter = counter + ~isnan(T{i}(ii:ss(1):nsize(1),jj:ss(2):nsize(2),kk:ss(3):nsize(3)));
                  clear Tadd;
                end
              end
            end
            varargout{i} = varargout{i} ./ counter;       
          elseif strcmp(method,'mean')
            varargout{i} = zeros(floor(size(T{i})./ss),'single');
            nsize = floor(size(T{i})./ss).*ss;
            for ii=1:ss(1)
              for jj=1:ss(2)
                for kk=1:ss(3)
                  varargout{i} = varargout{i} + T{i}(ii:ss(1):nsize(1),jj:ss(2):nsize(2),kk:ss(3):nsize(3));
                end
              end
            end
            varargout{i} = varargout{i} ./ prod(ss);
          else
            [Rx,Ry,Rz] = meshgrid(single(1+0.5*(ss(2)-1):ss(2):size(T{i},2)),...
                                  single(1+0.5*(ss(1)-1):ss(1):size(T{i},1)),...
                                  single(1+0.5*(ss(3)-1):ss(3):size(T{i},3)));
            varargout{i} = cat_vol_interp3f(single(T{i}),Rx,Ry,Rz,method);
          end
          
          if islogical(T{i}), varargout{i} = varargout{i}>0.5; end
%           if any(vx_red<=0.5)
%             [varargout{i},resT] = cat_vol_resize(varargout{i},'reduceV',vx_vol.*ss,vx_volr,minSize,method);
%             resT.vx_red = resT.vx_red .* ss;
%           else
            resT.vx_red = ss; resT.vx_volr=vx_volr;
 %         end
        else 
          varargout{i} = T{i}; 
          resT.vx_red=[1 1 1]; resT.vx_volr=vx_vol;
        end

      end
      varargout{i+1}.vx_red  = resT.vx_red;
      varargout{i+1}.vx_vol  = vx_vol;
      varargout{i+1}.vx_volr = resT.vx_volr;
      varargout{i+1}.sizeT   = sizeT;
      varargout{i+1}.sizeTr  = size(varargout{1});
     
      
    case 'dereducev'
      vx_red = varargin{1}.vx_red;
      sD     = varargin{1}.sizeT./vx_red + 0.5;
      if numel(varargin)<2, method='linear'; else method=varargin{2}; end;

      [Rx,Ry,Rz]   = meshgrid(single(0.5+0.5/vx_red(2):1/vx_red(2):sD(2)),...
                              single(0.5+0.5/vx_red(1):1/vx_red(1):sD(1)),...
                              single(0.5+0.5/vx_red(3):1/vx_red(3):sD(3)));

      for i=1:numel(T)  
        T{i}(isnan(T{i})) = 0; 
        if islogical(T{i}) && any(vx_red>1), varargout{i} = cat_vol_smooth3X(cat_vol_interp3f(single(T{i}),Rx,Ry,Rz,method),mean(vx_red))>0.5;
        else                                 varargout{i} = cat_vol_interp3f(single(T{i}),Rx,Ry,Rz,method);
        end
      end
      
      
    % ANISO2ISO & ISO2ANISO
    % __________________________________________________________________
    case 'aniso2iso'
      vx_vola = varargin{1};
      if numel(varargin)<2, method = 'linear'; else method = varargin{2}; end
      sizeT   = size(T{1});
      vx_inc  = round(vx_vola./(ones(size(vx_vola))*max(0.75,min(vx_vola))));
      vx_voli = vx_vola./vx_inc;
      [Rx,Ry,Rz] = meshgrid(single(0.5+0.5/vx_inc(2):1/vx_inc(2):sizeT(2)),...
                            single(0.5+0.5/vx_inc(1):1/vx_inc(1):sizeT(1)),...
                            single(0.5+0.5/vx_inc(3):1/vx_inc(3):sizeT(3)));
                              
      for i=1:numel(T)
        if max(0.75,min(vx_vola))*2 <= max(vx_vola)
          T{i}=cat_vol_interp3f(single(T{i}),Rx,Ry,Rz,method);
        else
          vx_voli=vx_vola; 
        end

        varargout{i} = T{i};
      end
      varargout{i+1}.vx_vola = vx_vola;
      varargout{i+1}.vx_voli = vx_voli;
      varargout{i+1}.sizeTa  = sizeT;
      
      
    case 'iso2aniso'
      sTa = varargin{1}.sizeTa;
      sTi = size(T{1}); 
      if numel(varargin)<2, method = 'linear'; else method = varargin{2}; end
      
      vx_red = varargin{1}.vx_vola./varargin{1}.vx_voli;
      [Rx,Ry,Rz] = meshgrid(single(1+0.5*(vx_red(2)-1):vx_red(2):sTi(2)),...
                            single(1+0.5*(vx_red(1)-1):vx_red(1):sTi(1)),...
                            single(1+0.5*(vx_red(3)-1):vx_red(3):sTi(3)));
      
      for i=1:numel(T)
        I = cat_vol_interp3f(single(T{i}),Rx,Ry,Rz,method);
        varargout{i}=zeros(sTa,'single'); sTa=min(sTa,size(I));
        varargout{i}(1:sTa(1),1:sTa(2),1:sTa(3))=I;
      end   

        
      
    % REDUCEBRAIN & DEREDUCEBRAIN 
    % __________________________________________________________________
    case 'reducebrain'
      if numel(varargin)<1, vx_vol=[1,1,1]; else vx_vol=varargin{1}; end      
      
      if numel(varargin)<2 || isempty(varargin{2}), d=1; else d=varargin{2}; end
      if numel(d)==1, d=d(1).*(ones(1,6)); 
      elseif numel(d)~=6, error('ERROR:reduceBrain: d has to have one or six elements.'); 
      elseif any(d([2,4,6])>(size(T{1})/2)), BB=d; d=[1 1 1 1 1 1]; % ????
      else
        error('ERROR:reduceBrain: unknown error using d.');
      end
      d = round(d./[vx_vol vx_vol]);
      
      if numel(varargin)>2 && ndims(varargin{3})==3
        M = varargin{3};
      elseif numel(varargin)>2 && ndims(varargin{3})==2
        if numel(varargin{3})==6, BB = varargin{3}; else error('BB has a wrong number of elements'); end
      elseif exist('BB','var') 
      else
        [X,M] = cat_vol_iscale(T{1},'findhead',vx_vol,4); 
      end
      
      if ~exist('BB','var') 
        if sum(M(:))>0
          SSUM=sum(sum(M,3),2); BB(1)=max(1,find(SSUM>0,1,'first')-d(1)); BB(2)=min(size(M,1),find(SSUM>0,1,'last')+d(2));
          SSUM=sum(sum(M,3),1); BB(3)=max(1,find(SSUM>0,1,'first')-d(3)); BB(4)=min(size(M,2),find(SSUM>0,1,'last')+d(4));
          SSUM=sum(sum(M,2),1); BB(5)=max(1,find(SSUM>0,1,'first')-d(5)); BB(6)=min(size(M,3),find(SSUM>0,1,'last')+d(6));
        else
          BB(1)=1; BB(2)=size(T,1);
          BB(3)=1; BB(4)=size(T,2);
          BB(5)=1; BB(6)=size(T,3);
        end
      end
      
      for i=1:numel(T)
        varargout{i} = T{i}(BB(1):BB(2),BB(3):BB(4),BB(5):BB(6)); 
      end
      varargout{i+1}.BB     = BB;
      varargout{i+1}.sizeT  = size(T{1});
      varargout{i+1}.sizeTr = size(varargout{1});
      
      
    case 'dereducebrain'
      BB = varargin{1}.BB;
      for i=1:numel(T) 
        if ndims(T{i})==3
          if islogical(T{i})
            TO = false(varargin{1}.sizeT);
          else
            TO = zeros(varargin{1}.sizeT,class(T{i}));
          end
          varargout{i}=TO; varargout{i}(BB(1):BB(2),BB(3):BB(4),BB(5):BB(6)) = T{i}(:,:,:); 
        elseif ndims(T{i})==4
          TO = zeros([varargin{1}.sizeT,size(T{i},4)],class(T{i})); 
          varargout{i}=TO; varargout{i}(BB(1):BB(2),BB(3):BB(4),BB(5):BB(6),:) = T{i}(:,:,:,:); 
        end
      end
    
    case 'interp'
      if numel(varargin)>0, hdr = varargin{1}; end
      if numel(varargin)>1, res = varargin{2}; end
      if numel(varargin)>2, method = varargin{3}; end
       
      if ~exist('hdr','var') || isempty(hdr),
        hdr.mat=[1 0 0 1;0 1 0 1; 0 0 1 1; 0 0 0 1]; 
        hdr.dim=size(T); 
      end
      if numel(res)==1, res(2:3)=res; end
      if numel(res)==2, res=[res(1),res]; end
      if ~exist('method','var'), method='linear'; end
      
      T = single(T{1});

      res    = round(res*100)/100;
      resV   = round(sqrt(sum(hdr.mat(1:3,1:3).^2))*100)/100;
      hdrO   = hdr; 
      sizeO  = size(T);

      if all(res>0)
        if 1% any(resV>res)
          % final size of the interpolated image
          [Dx,Dy,Dz]=meshgrid(single(res(2) / resV(2) : res(2)/resV(2) : size(T,2)),...
                              single(res(1) / resV(1) : res(1)/resV(1) : size(T,1)),...
                              single(res(3) / resV(3) : res(3)/resV(3) : size(T,3))); 
          %T = spm_sample_vol(T,Dx,Dy,Dz,method);
          T = cat_vol_interp3f(T,Dx,Dy,Dz,method);
          clear Dx Dy Dz;

          hdr.dim=size(T);
          %hdr.mat([1,2,3,5,6,7,9,10,11]) = hdr.mat([1,2,3,5,6,7,9,10,11]) .* res(1);
        else
          d = single(res./resV);
          [Rx,Ry,Rz]=meshgrid(d(1):d(1):size(T,2),d(2):d(2):size(T,1),d(3):d(3):size(T,3));
          T = cat_vol_interp3f(T,Rx,Ry,Rz,method);
          clear Rx Ry Rz;

          hdr.dim=size(T);
         % hdr.mat([1,2,3,5,6,7,9,10,11]) = hdr.mat([1,2,3,5,6,7,9,10,11]) .* res(1);
      %    hdr.mat([1,6,11])  = sign(hdr.mat([1,6,11]))  .* res;
        end
      end
      vmat = spm_imatrix(hdr.mat);
      vmat(7:9) = sign(vmat(7:9)).*res(1:3);
      hdr.mat = spm_matrix(vmat);
      
      varargout{1}       = T;
      varargout{2}.hdrO  = hdrO;
      varargout{2}.hdrN  = hdr;
      varargout{2}.sizeO = sizeO;
      varargout{2}.sizeN = size(T);
      varargout{2}.resO  = resV;
      varargout{2}.resN  = res;
      
      
    case 'deinterp'
      res   = varargin{1}.resO;
      resV  = varargin{1}.resN;
      sizeO = varargin{1}.resO;
      if numel(varargin)>1, method = varargin{2}; else method='linear'; end
      
      
      T = single(T{1});
      if strcmp(method,'masked')
        % for interpolation of partial defined maps like the cortical
        % thickness... finally 'nearest' interpolation is often good 
        % enought and much faster 
        if all(res>0) %&& any(round(abs(res)*100)/100>resV)
          d = single(res./resV);
          %[Rx,Ry,Rz]=meshgrid(0.5:d(1):size(D,2),0.5:d(2):size(D,1),0.5:d(3):size(D,3));
          [Rx,Ry,Rz]=meshgrid(d(2):d(2):size(T,2),d(1):d(1):size(T,1),d(3):d(3):size(T,3));
          M  = T>0.5; MM = cat_vol_morph(M,'d',2);
          [D,I] = cat_vbdist(T,MM); T=T(I); clear D I; 
          %Ts = smooth3(T); MM=cat_vol_morph(M,'e'); T(~MM)=Ts(~MM); clear MM; 
          M = cat_vol_interp3f(single(M),Rx,Ry,Rz,'linear')>0.5;
          T = cat_vol_interp3f(T,Rx,Ry,Rz,'linear');
          T = T .* M; 
          clear Rx Ry Rz;
        end
      else
        if all(res>0) %&& any(round(abs(res)*100)/100>resV)
          d = single(res./resV);
          %[Rx,Ry,Rz]=meshgrid(0.5:d(1):size(D,2),0.5:d(2):size(D,1),0.5:d(3):size(D,3));
          [Rx,Ry,Rz]=meshgrid(d(2):d(2):size(T,2),d(1):d(1):size(T,1),d(3):d(3):size(T,3));
          T = cat_vol_interp3f(T,Rx,Ry,Rz,method);
          clear Rx Ry Rz;
        end
      end      
      
      varargout{1} = zeros(varargin{1}.sizeO,'single');
      varargout{1}(1:min(size(T,1),varargin{1}.sizeO(1)), ...
                   1:min(size(T,2),varargin{1}.sizeO(2)), ...
                   1:min(size(T,3),varargin{1}.sizeO(3))) = ...
                 T(1:min(size(T,1),varargin{1}.sizeO(1)), ...
                   1:min(size(T,2),varargin{1}.sizeO(2)), ...
                   1:min(size(T,3),varargin{1}.sizeO(3))); 
      
        
    % OTHERWISE
    % __________________________________________________________________
    otherwise
      error('ERROR: cat_vol_resolution: unknown operation "%s"!\n',operation);
  end
end  
function D2=imageExpand(D,d)
  if nargin<2, d=1; end
  if d>1, D=ImageExpand(D,d-1); end
  
  D2=zeros(size(D)+1,class(D));
  D2(1:end-1,1:end-1,1:end-1) = D; clear D; 
  for i=1:2
    D2(1:end,1:end,end) = D2(1:end,1:end,end-1);
    D2(1:end,end,1:end) = D2(1:end,end-1,1:end);
    D2(end,1:end,1:end) = D2(end-1,1:end,1:end);
  end  
end

