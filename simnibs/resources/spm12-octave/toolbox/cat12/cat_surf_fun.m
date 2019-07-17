function varargout = cat_surf_fun(action,S,varargin)
% Function to collect surface functions.
% 
% varargout = cat_surf_fun(action,S)
%
%   D   = cat_surf_fun('dist',S);     % Estimate the distance between the 
%                                       vertices of the faces of S.
%   A   = cat_surf_fun('area',S);     % estimate surface area (faces)
%   V   = cat_surf_fun(S,F);          % map facedata to vertices
%   HS  = cat_surf_fun('hull',S);     % estimate (optimized) hull surface
%   V   = cat_surf_fun('surf2vol',S,varargin{1}); % render surface in a volume
%   E   = cat_surf_fun('graph2edge',T); % get edges of a triangulation T
%
%  * cdata mapping:
%    Mapping of textures from S2 to S1.
%    C   = cat_surf_fun('cdatamapping',S1,S2,cdata[,opt]); 
%
% ______________________________________________________________________
% Robert Dahnke 
% Structural Brain Mapping Group
% University Jena
% ______________________________________________________________________
% $Id: cat_surf_fun.m 1219 2017-11-19 23:28:59Z gaser $ 

  switch action
    case {'dist','distance'}
      varargout{1} = cat_surf_dist(S);
    case 'area'
      varargout{1} = cat_surf_area(S);
    case 'hull'
      if nargout==1, varargout{1} = cat_surf_hull(S); end
      if nargout==2, [varargout{1},varargout{2}] = cat_surf_hull(S); end
    case {'inner','outer'}
      if numel(varargin)==1
        switch nargout % surface & texture input
          case 0, cat_surf_GMboundarySurface(action,S,varargin{1});
          case 1, varargout{1} = cat_surf_GMboundarySurface(action,S,varargin{1}); 
          case 2, [varargout{1},varargout{2}] = cat_surf_GMboundarySurface(action,S,varargin{1}); 
        end
      else % file input
        switch nargout
          case 0, cat_surf_GMboundarySurface(action,S);
          case 1, varargout{1} = cat_surf_GMboundarySurface(action,S); 
          case 2, [varargout{1},varargout{2}] = cat_surf_GMboundarySurface(action,S); 
        end
      end
    case 'surf2vol'
      if nargin>2
        [varargout{1},varargout{2},varargout{3}] = cat_surf_surf2vol(S,varargin);
      else
        [varargout{1},varargout{2},varargout{3}] = cat_surf_surf2vol(S);
      end
    case 'graph2edge'
      varargout{1} = cat_surf_edges(S); 
    case 'cdatamappingtst'
      cat_surf_cdatamappingtst;
    case 'cdatamapping' 
      if nargin<3, varargin{3} = ''; end
      if nargin<4, varargin{4} = struct(); end
      if nargout>1
        [varargout{1},varargout{2}] = cat_surf_cdatamapping(S,varargin{1},varargin{2},varargin{3});
      else
        varargout{1} = cat_surf_cdatamapping(S,varargin{1},varargin{2},varargin{3});
      end  
  end
    
end

function varargout = cat_surf_GMboundarySurface(type,varargin)
  switch type
    case 'inner', direction = -0.5;
    case 'outer', direction =  0.5;
  end
  
  if nargin==2
    %% use filenames
    [pp,ff,ee] = spm_fileparts(varargin{1});
    
    if strcmp(ee,'')
      Praw = cat_io_FreeSurfer('fs2gii',varargin{1}); 
      Praw = Praw{1};
    else
      Praw   = varargin{1};
    end
    Pthick = cat_io_strrep(Praw,{'central','.gii'},{'thickness',''});
    Ptype  = cat_io_strrep(Praw,'central',type);
    
    cmd = sprintf('CAT_Central2Pial "%s" "%s" "%s" %0.2f',Praw,Pthick,Ptype,direction); 
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,1);

    if strcmp(ee,'')
      Ptype = cat_io_FreeSurfer('gii2fs',Ptype); 
    end
    
    % filename
    varargout{1} = Ptype; 
  else
    % write temp files ...
    Praw   = 'central.';
    Pthick = strrep(Praw,'central','thickness');
    Ptype  = strrep(Praw,'central',type);
   
    cmd = sprintf('CAT_Central2Pial "%s" "%s" %0.2f',Praw,Pthick,Ptype,direction); 
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,1);
    
    % load surface 
    varargout{1} = gifti(Ptype); 
    
    % delete temp files
    delete(Praw,Pthick,Ptype);
  end
end

function cat_surf_cdatamappingtst

%% Testdata
   Psubcentral  = ['/Volumes/vbmDB/MRData/vbm12tst/results/deffiles/cg_vbm_defaults_template/template_NKI/'...
     'surf/lh.central.NKI_HC_NKI_1013090_T1_SD000000-RS00.gii'];
   PsubsphereA  = strrep(Psubcentral,'central','sphere.reg');              
   %Psubthick    = strrep(strrep(Psubcentral,'central','thickness'),'.gii','');               
   Psubthickres = strrep(strrep(Psubcentral,'central','thickness.resampled'),'lh.','s15mm.lh.'); 
   Psubtmp      = strrep(Psubcentral,'central','tmp'); 
   Pavgtmp      = strrep(strrep(Psubcentral,'central','tmp.resampled'),'lh.','s15mm.lh.'); 
 
   %Pavgcentral  = '/Users/dahnke/Neuroimaging/spm12/toolbox/cat12/templates_surfaces/lh.central.freesurfer.gii'; 
   PavgsphereA  = '/Users/dahnke/Neuroimaging/spm12/toolbox/cat12/templates_surfaces/lh.sphere.freesurfer.gii'; 
   PavgDKT40    = '/Users/dahnke/Neuroimaging/spm12/toolbox/cat12/atlases_surfaces/lh.aparc_DKT40JT.freesurfer.annot';
   
   
%% Test 1 - avg2sub - ok
   Ssub = gifti(PsubsphereA);
   Savg = gifti(PavgsphereA); 
   [vertices, label, colortable]  = cat_io_FreeSurfer('read_annotation',PavgDKT40); 
   Savg.cdata = label; 
   
   S3 = gifti(Psubcentral); 
   S3.cdata = cat_surf_fun('cdatamapping',Ssub,Savg,'nearest');
   save(gifti(S3),Psubtmp);
   
%% Test2 - sub2avg - ok
   Savg = gifti(PavgsphereA); 
   Ssub = gifti(PsubsphereA);
   %Ssub.cdata = cat_io_FreeSurfer('read_surf_data',Psubthick); 
   Ssub.cdata = cat_surf_fun('area',gifti(Psubcentral));
   
   S3 = gifti(Psubthickres); 
   mapping = {'directed'}; %,'undirected'}; %'nearest',
   for mi = 1:numel(mapping)
     S3.cdata  = cat_surf_fun('cdatamapping',Savg,Ssub,mapping{mi},1);
     S3.cdata  = spm_mesh_smooth(struct('vertices',S3.vertices,'faces',S3.faces),double(S3.cdata'),5);
     fprintf('mapping = %10s: A(sub) = %0.2f vs. A(avg) = %0.2f\n',mapping{mi},sum(Ssub.cdata(:)),sum(S3.cdata(:))); 
     save(gifti(S3),Pavgtmp); cat_surf_display(Pavgtmp)
   end
   
end
% nearest connection between to surfaces
function varargout = cat_surf_cdatamapping(S1,S2,cdata,opt) 
  if ischar(S1), S1 = gifti(S1); end
  if ischar(S2), S2 = gifti(S2); end
  if ischar(cdata),
    Pcdata = cdata;
    [pp,ff,ee] = spm_fileparts(cdata); 
    switch ee
      case '.annot'
        [vertices, cdata]  = cat_io_FreeSurfer('read_annotation',Pcdata); 
        clear vertices
      case '.gii'
        Scdata = gifti(S2); 
        if isfield(Scdata,'cdata')
          cdata = SX.cdata;
        else
          error('cat_surf_fun:cdatamapping:noTexture','No texture found in "%s"!\n',Pcdata);
        end
      otherwise
        cdata =  cat_io_FreeSurfer('read_surf_data',Pcdata);   
    end
  end
  
  if ~exist('cdata','var') || isempty(cdata)
    if isfield(S2,'cdata'), cdata = S2.cdata; end
  end
  
  if ~exist('opt','var'), opt = struct(); end
  def.method = 'nearest';
  def.verb   = 0; 
  def.smooth = 0; 
  opt        = cat_io_checkinopt(opt,def);
  
  if opt.verb, stime1 = cat_io_cmd(sprintf('Data-mapping (%s)',method)); fprintf('\n'); end
  
  % prepare vertices
  S1.vertices = S1.vertices ./ repmat(max(S1.vertices),size(S1.vertices,1),1)*1.1; % *100 
  S2.vertices = S2.vertices ./ repmat(max(S2.vertices),size(S2.vertices,1),1); 
  verticesS1  = double(S1.vertices - repmat(mean(S1.vertices),size(S1.vertices,1),1)); 
  verticesS2  = double(S2.vertices - repmat(mean(S2.vertices),size(S2.vertices,1),1)); 
  
  
  % estimate mapping
  switch opt.method
    case {'nearest'}
      [varargout{2},varargout{3}] = dsearchn([verticesS2;inf(1,3)],double([S2.faces ones(size(S2.faces,1),1)*(size(S2.vertices,1)+1)]),verticesS1);
      varargout{1} = cdata(varargout{2}); 
    case {'undirected','directed'}
      %% use the surface as delauny graph
      switch opt.method 
        case 'directed'
          if opt.verb,  stime = cat_io_cmd('  Edge-Estimation (Nearest)','g5',''); end
          nextS2fromS1 = dsearchn([verticesS2;inf(1,3)],double([S2.faces ones(size(S2.faces,1),1)*(size(S2.vertices,1)+1)]),verticesS1);
          nextS1fromS2 = dsearchn([verticesS1;inf(1,3)],double([S1.faces ones(size(S1.faces,1),1)*(size(S1.vertices,1)+1)]),verticesS2);
          tmp = nextS1fromS2; nextS1fromS2 = nextS2fromS1; nextS2fromS1 = tmp;
          nearestedges = [ (1:numel(nextS2fromS1))', nextS2fromS1; nextS1fromS2 , (1:numel(nextS1fromS2))' ]; 
          nearestedges = unique(nearestedges,'rows');
        case 'undirected'
          if opt.verb,  stime = cat_io_cmd('  Edge-Estimation (Delaunay','g5',''); end
          % nearest is required too
          nextS2fromS1 = dsearchn([verticesS2;inf(1,3)],double([S2.faces ones(size(S2.faces,1),1)*(size(S2.vertices,1)+1)]),verticesS1);
          nextS1fromS2 = dsearchn([verticesS1;inf(1,3)],double([S1.faces ones(size(S1.faces,1),1)*(size(S1.vertices,1)+1)]),verticesS2);
          tmp = nextS1fromS2; nextS1fromS2 = nextS2fromS1; nextS2fromS1 = tmp;
          nearestedges  = [ (1:numel(nextS2fromS1))', nextS2fromS1; nextS1fromS2 , (1:numel(nextS1fromS2))' ]; 
          nearestedges1 = unique(nearestedges,'rows');
          % delauany
          triangulation = delaunayn([verticesS2;verticesS1]);              % delaunay triangulation
          nearestedges  = cat_surf_fun('graph2edge',triangulation);        % get edges 
          nearestedges(sum(nearestedges<=size(verticesS2,1),2)~=1,:)=[];   % only edges between S1 and S2
          nearestedges(:,2) = nearestedges(:,2) - size(verticesS2,1); 
          nearestedges = unique([nearestedges;nearestedges1],'rows');
      end
      if opt.verb, stime = cat_io_cmd('  Weighting','g5','',1,stime); end
      
      if 0
        %% my little testset
        nextS1fromS2 = [1; 1; 3; 4; 4; 4; 5; 5]; 
        nextS2fromS1 = [1; 3; 3; 5; 8; 8];
        cdata        = [1 1 1 1 1 1]';
        nearestedges = [ (1:numel(nextS2fromS1))', nextS2fromS1; nextS1fromS2 , (1:numel(nextS1fromS2))' ]; 
        nearestedges = unique(nearestedges,'rows');
      end
      
      
      %% simplify edges 1
      if 0
        % simpler, but much slower 
        nearestedges = [nearestedges, ones(size(nearestedges,1),1)]; % default weight
        [NeighborsS1,NidS1]  = hist(nearestedges(:,1),1:1:max(nearestedges(:,1)));
        for ni=NidS1(NeighborsS1>1)
          NumNi = nearestedges(:,1)==ni; 
          nearestedges(NumNi,3) =  nearestedges(NumNi,3) ./ sum(NumNi);
        end
      else
        % faster 
        %nearestedges = [nearestedges, ones(size(nearestedges,1),1)]; % default weight
        dist = sum( (S2.vertices(nearestedges(:,1),:) - S1.vertices(nearestedges(:,2),:)).^2 , 2) .^ 0.5; 
        nearestedges = [nearestedges, dist]; % default weight
        list = [1; find(nearestedges(1:end-1,1)~=nearestedges(2:end,1))+1; size(nearestedges,1)]; 
        for ni=1:numel(list)-1
          %nearestedges(list(ni):list(ni+1)-1,3) = nearestedges(list(ni):list(ni+1)-1,3) ./ (list(ni+1) - list(ni)); 
          nearestedges(list(ni):list(ni+1)-1,3) = nearestedges(list(ni):list(ni+1)-1,3) ./ sum(nearestedges(list(ni):list(ni+1)-1,3)); 
        end
      end
      if opt.verb, stime = cat_io_cmd('  Mapping','g5','',1,stime); end

      %%
      if 0
        % correct & simple, but very slow
        varargout{1} = zeros(1,max(nearestedges(:,2)));
        for ni=1:size(nearestedges,1)
          varargout{1}(nearestedges(ni,2)) = varargout{1}(nearestedges(ni,2)) + ...
            cdata(nearestedges(ni,1))' .* nearestedges(ni,3)';
        end
      else
        varargout{1} = zeros(1,max(nearestedges(:,2)));
        if 0
          list = [1; find(nearestedges(1:end-1,2)~=nearestedges(2:end,2))+1; size(nearestedges,1)+1]; 
          for ni=1:numel(list)-1
            varargout{1}(nearestedges(list(ni),2)) = varargout{1}(nearestedges(list(ni),2)) + ...
              sum(cdata(nearestedges(list(ni):list(ni+1)-1,1)) .*  nearestedges(list(ni):list(ni+1)-1,3));
          end
        else
          nearestedges2 = sortrows([nearestedges(:,2) nearestedges(:,1) nearestedges(:,3)]);  
          list = [1; find(nearestedges2(1:end-1,1)~=nearestedges2(2:end,1))+1; size(nearestedges2,1)+1]; 
          for ni=1:numel(list)-1
            varargout{1}(nearestedges2(list(ni),1)) = varargout{1}(nearestedges2(list(ni),1)) + ...
              sum(cdata(nearestedges2(list(ni):list(ni+1)-1,2)) .*  nearestedges2(list(ni):list(ni+1)-1,3));
          end
        end
      end
      if numel(varargout{1})<20, disp(varargout{1}); end
      if opt.verb, cat_io_cmd(' ','g5','',1,stime); end
  end
  
  % default smoothing???
  if opt.smooth
    varargout{1}  = spm_mesh_smooth(struct('vertices',S3.vertices,'faces',S3.faces),double(varargout{1}'),opt.smooth);
  end
  
  if isfield(opt,'fname')
    save(gifti(struct('vertices',S1.vertices,'faces',S1.faces,'cdata',varargout{1})),opt.fname); 
  end
  
  if opt.verb, cat_io_cmd('','','',1,stime1); end
end

function E = cat_surf_edges(T)
  if isstruct(T) && isfield(T,'faces')
    T = T.faces;
  end

  T = sort(T,2); E = []; 
  for i=1:size(T,2)-1
    E = [E; T(:,[i i+1])]; %#ok<AGROW>
  end
  E = unique(E,'rows');
end

function D = cat_surf_dist(S)
% Estimate the distance between the vertices of the faces of S.
% D = [c,a,b] = [d(AB),d(BC),d(CA)]

  D = [sum( (S.vertices(S.faces(:,1),:) - S.vertices(S.faces(:,2),:)).^2 , 2) .^ 0.5, ...
       sum( (S.vertices(S.faces(:,2),:) - S.vertices(S.faces(:,3),:)).^2 , 2) .^ 0.5, ...
       sum( (S.vertices(S.faces(:,3),:) - S.vertices(S.faces(:,1),:)).^2 , 2) .^ 0.5]; 
     
end

function [AV,AF] = cat_surf_area(S)
% Calculate surface area of the faces AF (Horonsche Form) and map it to the
% vertices AV.

  % facearea (Horonsche Form)
  D = cat_surf_dist(S);
  facesp = sum(D,2) / 2;  % s = (a + b + c) / 2;
  AF = (facesp .* (facesp - D(:,1)) .* (facesp - D(:,2)) .* (facesp - D(:,3))).^0.5; % area=sqrt(s*(s-a)*(s-b)*(s-c));
  
  % numerical (to small point diffences) and mapping problems (crossing of streamlines)
  % -> correction because this is theoretical not possible (laplace field theory)
  AF(AF==0) = eps; % to small values
  AF = abs(AF);    % streamline-crossing
    
  AV = cat_surf_F2V(S,AF);
end

function data = cat_surf_F2V(S,odata)
%% mapping of facedata to vertices

  data   = zeros(size(S.vertices,1),1);
  [v,f]  = sort(S.faces(:)); 
  [f,fj] = ind2sub(size(S.faces),f);  %#ok<ASGLU>
  far = odata(f);
  for i=1:numel(S.faces), data(v(i)) = data(v(i)) + far(i); end

  data = data / size(S.vertices,2); % Schwerpunkt... besser Voronoi, aber wie bei ner Oberfl?che im Raum???
end

function [SH,V] = cat_surf_hull(S)
%% hull creation

  % render surface points
  V = false( round(max(S.vertices,[],1) - min(S.vertices))+10 );     
  I = sub2ind(size(V),round(S.vertices(:,1) - min(S.vertices(:,1)) + 5),...
                      round(S.vertices(:,2) - min(S.vertices(:,2)) + 5),...
                      round(S.vertices(:,3) - min(S.vertices(:,3)) + 5));
  V(I) = 1; clear I; 
  
  % 
  V  = cat_vol_morph(V,'lc',mean(size(V))/6); % closing 
  V  = cat_vol_smooth3X(V,2);    % smoothing
  SH = isosurface(V,0.4);        % create hull 
  V  = V>0.4;
  
  SH.vertices = [SH.vertices(:,2) SH.vertices(:,1) SH.vertices(:,3)]; % matlab flip
  SH.vertices = SH.vertices + repmat(min(S.vertices),size(SH.vertices,1),1) - 5;
end

function [V,vmat,vmati] = cat_surf_surf2vol(S,opt)
%% render inner surface area 
%  Render the volume V with V==1 within the surface. 
%  Use type=1 to render also the surface area with 0.5.
%  The transformation imat to create 
%  SH.vertices = [SH.vertices(:,2) SH.vertices(:,1) SH.vertices(:,3)];     % matlab flip
%  SH.vertices = SH.vertices + imat;

  if ~exist('opt','var'), opt = struct(); end
  def.debug  = 1;
  def.type   = 0; 
  def.refine = 0.8;
  def.bdist  = 5; 
  def.res    = 1; % not yet ...
  
  opt = cat_io_checkinopt(opt,def);
  
  % save a temporary version of S and refine it
  Praw = [tempname '.gii'];
  save(gifti(S),Praw);
  
  cmd = sprintf('CAT_RefineMesh "%s" "%s" %0.2f',Praw,Praw,opt.refine); 
  [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.debug);

  S = gifti(Praw);
  delete(Praw);
  
  %% render surface points
  V    = false( round(max(S.vertices,[],1) - min(S.vertices))+10 );     
  vmat = -[min(S.vertices(:,1)) min(S.vertices(:,2)) min(S.vertices(:,3))] + opt.bdist; 
  I    = sub2ind(size(V),...
        max(1,min(size(V,1),round(S.vertices(:,1) + vmat(1)))),...
        max(1,min(size(V,2),round(S.vertices(:,2) + vmat(2)))),...
        max(1,min(size(V,3),round(S.vertices(:,3) + vmat(3)))));
  V(I) = 1; 
  
  V    = cat_vol_morph(V,'lc',1);  % closeing 
  V(I) = 0;
  V    = cat_vol_morph(V,'lab');
  if opt.type==1
    Vd = cat_vol_morph(V,'d',1); 
    V  = single(V);
    V(intersect(I,find(Vd>0))) = 0.5;
    V  = cat_vol_smooth3X(V,0.6);    % smoothing
  end
  
  vmati = repmat(min(S.vertices),size(S.vertices,1),1) - 5; 
  %%
  %SH = isosurface(V,0.6);        % create hull 
  %SH.vertices = [SH.vertices(:,2) SH.vertices(:,1) SH.vertices(:,3)]; % matlab flip
  %SH.vertices = SH.vertices + repmat(min(S.vertices),size(SH.vertices,1),1) - 5;
end
