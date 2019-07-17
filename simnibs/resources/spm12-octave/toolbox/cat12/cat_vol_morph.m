function vol = cat_vol_morph(vol,action,n,vx_vol)
% ______________________________________________________________________
% Morphological operations for a volume vol based on a 26-neighborhood 
% (erode, dilate, open, close) or a distance transformation (disterode,
% distdilate, distopen, distclose). Furthermore, 3 labeling operation
% (lab, labopen, labclose) that mask the largest cluster (after an 
% distopen/disterode) are available. The voxel dimensions vx_vol are 
% only available for distancebased transformations, where n depend on 
% the distance. 
%
% out = cat_vol_morph(in,action[,n,vx_vol])
%
% in     = input volume that will be thresholded at 0.5
% action = {'d'|'e'|'c'|'o'|'dd'|'de'|'dc'|'do'|'l'|'lo'|'lc'}
% n      = 1x1 double (default=1), will be rounded for standard 
%          morphological operations, but not for distancebased operations.
% vx_vol = 1x1 or 1x3 double (default=1)
% out    = volume with the same class like the input volume
%
% Actions:
%   Morphological operations with 26-neighborhood (cube):
%    - d??| dilate 
%    - e  | erode  
%    - c  | close  
%    - o  | open   
%
%   Morphological operations with distance opereration (sphere):
%    - dd | distdilate
%    - de | disterode
%    - dc | distclose
%    - do | distopen
%    - l  | lab          n(1) largest object/cluster with at least 
%                        n(2) absolute voxels for negative n(2)
%                             or relative voxels for positive n(2)
%    - lo | labopen      (disterode  + distdilate + lab)
%    - lc | labclose     (distdilate + disterode  + lab)
%
%   Special operation:
%    - st | selftest     [in development]
% ______________________________________________________________________
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group
% University Jena 
% $Id: cat_vol_morph.m 1265 2018-02-01 10:10:55Z gaser $

% ______________________________________________________________________
%
% ToDo:
% Large n can increase computation times strongly. Oftenly, the where 
% used only for a low quality correction i.e. to create a smoothe
% complex hull of an object. 
% For a future release the cat_vol_resize function and further actions
% will allow a faster processing for images, where not the highest
% quality is necessary. 
% For fast estimation the prefix 'f' should be added to the action.
% ______________________________________________________________________


  if nargin < 4, vx_vol = 1; end
  if nargin < 3, n      = 1; end
  if nargin < 2, action = ''; end
  if nargin < 1, error('MATLAB:cat_vol_morph:NoAction','No volume given.\n'); end

  classVol = class(vol); 
  
  if iscell(vol) || ndims(vol)~=3 || isempty(vol)
    error('MATLAB:cat_vol_morph:Empty','Only nonempty 3D volumes!\n'); 
  end
  
  vol = vol>0.5; vol(isnan(vol)) = 0;
  
  if numel(vx_vol)==1, vx_vol=repmat(vx_vol,1,3); end
  if any(size(vx_vol)~=[1,3]), 
    error('MATLAB:cat_vol_morph:vx_vol', ...
      'Wrong vx_vol size. It has to be a 1x3 matrix.\n'); 
  end
  
  no=n; n=double(n); n(1)=round(n(1)); 
  switch lower(action)
    case {'l' 'lc' 'lo' 'labclose' 'labopen'}
      % not return in this case
    otherwise
      if n==0, return; end 
  end
  
  switch lower(action)
    case {'dilate' 'd'}
      % use of single input for convn is faster and less memory demanding
      vol = convn(single(vol),ones(2*n+1,2*n+1,2*n+1),'same') > 0; 

    case {'erode' 'e'}
      vol = ~cat_vol_morph(~vol,'dilate',n,vx_vol); 

    case {'close' 'c'}
      sz = size(vol);
      vol2 = zeros(sz(1)+(2*n),sz(2)+(2*n),sz(3)+(2*n),'uint8');
      vol2(n+1:sz(1)+n,n+1:sz(2)+n,n+1:sz(3)+n) = uint8(vol);
      vol2=cat_vol_morph(vol2,'dilate',n,vx_vol); 
      vol2=cat_vol_morph(vol2,'erode' ,n,vx_vol); 
      vol = vol2(n+1:sz(1)+n,n+1:sz(2)+n,n+1:sz(3)+n)>0;
      
    case {'open' 'o'}
      vol=~cat_vol_morph(~vol,'close' ,n,vx_vol); 

    case {'labclose' 'lc'}
      vol = cat_vol_morph(vol,'close',n,vx_vol); 
      vol = ~cat_vol_morph(~vol,'lab',n,vx_vol); % removing of background within the object

    case {'labopen' 'lo'}
      vol = cat_vol_morph(vol,'open',n,vx_vol); 
      vol = cat_vol_morph(vol,'lab',n,vx_vol); % removing of other objects

    %===================================================================
    case {'lab' 'l'}
    % try to catch errors, if there is no object
      try  
        [ROI,num]  = spm_bwlabel(uint8(vol),6);
        
        if numel(num)>0
          num        = hist( ROI( ROI(:)>0 ) , 1:num);
          [num,numi] = sort(num,'descend');
          vol        = ROI==numi(1);	
          
          if exist('n','var') && n(1)>1
            vol = single(vol); classVol = 'single'; 
            
            if numel(n)==1, n(2)=0; end
            for ni=2:min(numel(num),n(1))
              if (n(2)<0 && num(ni)>(-n(2))) || ... % absolute vs. 
                 (n(2)>0 && num(ni)/num(1)>(n(2)))  % relative
                vol(ROI==numi(ni)) = numi(ni);	
              end
            end
            
          end
        end
      catch %#ok<CTCH>
        %vol = [];
        %warning('MATLAB:cat_vol_morph:NoObject','WARNING: cat_vol_morph - lab - no object!');
      end 

      
    
    %===================================================================
    % You have to use the original resolution, because fine structure 
    % are bad represented for lower resolutions and lead to unaccurate 
    % results.
    case {'distdilate' 'dd'}
      vol = cat_vbdist(single(vol),true(size(vol)),vx_vol)<=no;
  
    case {'disterode' 'de'}
      vol = ~cat_vol_morph(~vol,'distdilate',n,vx_vol); 

    case {'distclose' 'dc'}
      sz   = size(vol);
      vol2 = zeros(sz(1)+(2*n),sz(2)+(2*n),sz(3)+(2*n),'single');
      vol2(n+1:sz(1)+n,n+1:sz(2)+n,n+1:sz(3)+n) = single(vol);
      vol2 = cat_vbdist(vol2,true(size(vol2)),vx_vol)<no;
      vol2 = cat_vbdist(single(~vol2),true(size(vol2)),vx_vol)<no;
      vol  = vol | ~vol2(n+1:sz(1)+n,n+1:sz(2)+n,n+1:sz(3)+n);

    case {'distopen' 'do'}
      vol = ~cat_vol_morph(~vol,'distclose',n,vx_vol); 
      
    case {'labdistclose' 'ldc'}
      vol = cat_vol_morph(vol,'distclose',n,vx_vol); 
      vol = ~cat_vol_morph(~vol,'lab',n,vx_vol); % removing of background within the object

    case {'labdistopen' 'ldo'}
      vol = cat_vol_morph(vol,'distopen',n,vx_vol); 
      vol = cat_vol_morph(vol,'lab',n,vx_vol); % removing of other objects


    %===================================================================
    case {'selftest' 'st'}
      % a=zeros(7,11,3); a(4,4,2)=1; a(4,8,2)=1; % two dots
      
      voltypes  = {'1','2','2c','2ce'};
      volclass  = {'cube','sphere'};
      method{1} = {'erode'      'e'
                   'dilate'     'd'
                   'open'       'o'
                   'close'      'c'};
      method{2} = {'disterode'  'de'
                   'distdilate' 'dd'
                   'distopen'   'do'
                   'distclose'  'dc'};
      method{3} = {'lab'        'l'
                   'labopen'    'lo'
                   'labclose'   'lc'};

      dist = 8; %[0:0.5:3 10 20];

      vol = cell(1,numel(voltypes)); 
      for vc=1:numel(volclass)
        for vt=1:numel(voltypes)
          vol{vt}.O = cat_tst_phantoms(volclass{vc},voltypes{vt});

          for cl=1:numel(method)
            for mt=1:size(method{cl},1)
              for dt=1:numel(dist)
                vol{vt}.(method{cl}{mt,2}){dt} = cat_vol_morph(vol{vt}.O,method{cl}{mt,1},dist(dt));
              end
            end
          end
        end
      end
      
%ds('d2','',[1 1 1],vol{1}.O,vol{1}.O + vol{1}.d{1} + vol{1}.e{1},vol{1}.O,vol{1}.O + vol{1}.dd{1} + vol{1}.de{1},50)

      
    %===================================================================
    % case do nothing
    case ''
      
    otherwise
      error('MATLAB:cat_vol_morph:UnknownAction','Unknown action ''%s ''',action);
  end
  
  eval(sprintf('vol = %s(vol);',classVol));
  if isa(classVol,'uint8'); vol = 255*vol; end
end