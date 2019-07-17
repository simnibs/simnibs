function [varargout] = cat_surf_info(P,readsurf,gui,verb)
% ______________________________________________________________________
% Extact surface information from filename.
%
% sinfo = cat_surf_info(P,readsurf,gui,verb)
%
%   P         .. surface filename
%   readsurf  .. read gifti or Freesurfer file to get more information
%   gui       .. interactive hemisphere selection
%   verb      .. verbose
%
% sinfo(i). 
%   fname     .. full filename
%   pp        .. filepath
%   ff        .. filename
%   ee        .. filetype
%   exist     .. exist file?
%   fdata     .. structure from dir command
%   ftype     .. filetype [0=no surface,1=gifti,2=freesurfer]
%   statready .. ready for statistic (^s#.*.gii) [0|1]
%   side      .. hemisphere [lh|rh|lc|rc|mesh] 
%   name      .. subject/template name
%   datatype  .. [-1=unknown|0=nosurf|1=mesh|2=data|3=surf]
%                only defined for readsurf==1 and surf=mesh+sidata
%   dataname  .. datafieldname [central|thickness|intensity...]
%   texture   .. textureclass [central|sphere|thickness|...]
%   label     .. labelmap
%   resampled .. resampled data [0|1] 
%   template  .. template or individual mesh [0|1] 
%   name      .. name of the dataset
%   roi       .. roi data
%   nvertices .. number vertices
%   nfaces    .. number faces
%   Pmesh     .. underlying meshfile
%   Psphere   .. sphere mesh
%   Pspherereg.. registered sphere mesh
%   Pdefects  .. topology defects mesh
%   Pdata     .. datafile
%   preside   .. prefix before hemi info (i.e. after smoothing)
%   posside   .. string after hemi info
%   smoothed  .. smoothing size
%   Phull     .. convex hull mesh
% ______________________________________________________________________
% Robert Dahnke
% $Id: cat_surf_info.m 1270 2018-02-07 13:32:14Z gaser $

%#ok<*RGXP1>

  if ~exist('P','var'), P=''; end
  if strcmp(P,'selftest')
    pps = {
      fullfile(spm('dir'),'toolbox','cat12','templates_surfaces')
      fullfile('User','08.15','T1 T2','subs','mri')
      };
    ffs = {
      'lh.central.freesurfer'
      'lh.mymask'
      'output'
      'lh.texture.sub1.sub2'
      'lh.tex1.tex2.resampled.sub1.sub2'
      '08.15'
      's15mm.lh.tex03.33.resampled.S01.mri'
      's5mm.lh.t1.t2-3_3.resampled.S01_-.kdk.mri'
      'rh.s33mmtexture.S01.native.mri'
      'rh'
      'rh.'
      'rh.sphere.reg.sub1'
      'rc.defects.038.37.477'
      'lc.s33mmtexture.S01.native.mri'
      'rc.texture.sub1.sub2'
      };
    ees = {
      ''
      ... '.gii'
      ... '.annot'
      };
    varargout = cell(numel(pps),numel(ffs),numel(ees)); 
    for ppsi = 1:numel(pps)
      for ffsi = 1:numel(ffs)
        for eesi = 1:numel(ees)
          varargout{1}(ppsi,ffsi,eesi) = cat_surf_info(fullfile(pps{ppsi},[ffs{ffsi} ees{eesi}]),0,0,1);
        end
      end
    end
    return; 
  end
  
  
  
  if nargin<2, readsurf = 0; end
  if nargin<3, gui  = 0; end
  if nargin<4, verb = 0; end

  P = cellstr(P);
  
  sinfo = struct(...
    'fname','',...      % full filename
    'pp','',...         % filepath
    'ff','',...         % filename
    'ee','',...         % filetype
    'exist','',...      % exist
    'fdata','',...      % datainfo (filesize)
    'ftype','',...      % filetype [0=no surface,1=gifti,2=freesurfer]
    'statready',0,...   % ready for statistic (^s#.*.gii)
    'side','',...       % hemishphere
    'name','',...       % subject/template name
    'datatype','',...   % datatype [0=nosurf/file|1=mesh|2=data|3=surf] with surf=mesh+data
    'dataname','',...   % datafieldname [central|thickness|s3thickness...]
    'texture','',...    % textureclass [central|sphere|thickness|...]
    'label','',...      % labelmap
    'resampled','',...  % dataspace
    'template','',...   % individual surface or tempalte
    'roi','',...        % roi data
    'nvertices',[],...  % number vertices
    'nfaces',[],...     % number faces
    'Pmesh','',...      % meshfile
    'Psphere','',...    % meshfile
    'Pspherereg','',... % meshfile
    'Pdefects','',...   % meshfile
    'Pdata','',...      % datafile
    'preside','', ...
    'posside','' ...
  );

  if isempty(P), varargout{1}=sinfo; return; end
  
  for i=1:numel(P)
    [pp,ff,ee] = spm_fileparts(P{i});
    sinfo(i).fdata = dir(P{i});
    
    sinfo(i).fname = P{i};
    sinfo(i).exist = exist(P{i},'file') > 0; 
    sinfo(i).pp = pp;
    switch ee
      case {'.xml','.txt','.html','.csv'}
        sinfo(i).ff = ff;
        sinfo(i).ee = ee;
        sinfo(i).ftype = 0;
        continue
      case '.gii'
        sinfo(i).ff = ff;
        sinfo(i).ee = ee;
        sinfo(i).ftype = 1;
        if sinfo(i).exist && readsurf
          S = gifti(P{i});
        end
      case '.annot'
        sinfo(i).ff = ff;
        sinfo(i).ee = ee;
        sinfo(i).ftype = 1;
        sinfo(i).label = 1; 
        if sinfo(i).exist && readsurf
          clear S; 
          try
            S = cat_io_FreeSurfer('read_annotation',P{1}); 
          end
        end
        if exist('S','var')
          sinfo(i).ftype = 2;
        end
      otherwise
        sinfo(i).ff = [ff ee];
        sinfo(i).ee = '';
        sinfo(i).ftype = 0;
        if sinfo(i).exist && readsurf
          clear S; 
          try
            S = cat_io_FreeSurfer('read_surf',P{1}); 
            if size(S.faces,2)~=3 || size(S.faces,1)<10000
              clear S; 
            end
          end
          try
            S.cdata = cat_io_FreeSurfer('read_surf_data',P{1}); 
            if size(S.face,2)==3 || size(S.face,1)<10000
              S = rmfield(S,'cdata'); 
            end
          end
        end
        if exist('S','var')
          sinfo(i).ftype = 2;
        end
    end
    
    
    noname = sinfo(i).ff; 
    
    % smoothed data
    sinfo(i).statready = ~isempty(regexp(noname,'^s(?<smooth>\d+)\..*')); 
    
    % side
    if     strfind(noname,'lh.'),   sinfo(i).side='lh';   sidei = strfind(noname,'lh.');
    elseif strfind(noname,'rh.'),   sinfo(i).side='rh';   sidei = strfind(noname,'rh.');
    elseif strfind(noname,'mesh.'), sinfo(i).side='mesh'; sidei = strfind(noname,'mesh.');
    elseif strfind(noname,'lc.'),   sinfo(i).side='lc';   sidei = strfind(noname,'lc.');
    elseif strfind(noname,'rc.'),   sinfo(i).side='rc';   sidei = strfind(noname,'rc.');
    else
      % if SPM.mat exist use that for side information
      if exist(fullfile(pp,'SPM.mat'),'file')
        load(fullfile(pp,'SPM.mat'));
        [pp2,ff2]   = spm_fileparts(SPM.xY.VY(1).fname);        
      
        % find mesh string
        hemi_ind = strfind(ff2,'mesh.');
        if ~isempty(hemi_ind)
          sinfo(i).side = ff2(hemi_ind(1):hemi_ind(1)+3);
        else
          % find lh|rh string
          hemi_ind = [strfind(ff2,'lh.') strfind(ff2,'rh.') strfind(ff2,'lc.') strfind(ff2,'rc.')];
          sinfo(i).side = ff2(hemi_ind(1):hemi_ind(1)+1);
        end

        sidei=[];
      else
        if gui
          if cat_get_defaults('extopts.expertgui')
            sinfo(i).side = spm_input('Hemisphere',1,'lh|rh|lc|rc|mesh');
          else
            sinfo(i).side = spm_input('Hemisphere',1,'lh|rh|mesh');
          end
        else
          sinfo(i).side = ''; 
        end
        sidei = strfind(noname,[sinfo(i).side '.']);
      end
    end
    if isempty(sidei), sidei = strfind(noname,sinfo(i).side); end
    if sidei>0
      sinfo(i).preside = noname(1:sidei-1);
      sinfo(i).posside = noname(sidei+numel(sinfo(i).side)+1:end);
    else
      sinfo(i).posside = noname;
    end
    
    % smoothed
    if isempty(sinfo(i).preside)
      sinfo(i).smoothed = 0; 
    else
      sinfo(i).smoothed = max([0,double(cell2mat(textscan(sinfo(i).preside,'s%dmm.')))]);
    end

    % datatype
    if sinfo(i).exist && readsurf
      switch num2str([isfield(S,'vertices'),isfield(S,'cdata')],'%d%d')
        case '00',  sinfo(i).datatype  = 0;
        case '01',  sinfo(i).datatype  = 1;
        case '10',  sinfo(i).datatype  = 2;
        case '11',  sinfo(i).datatype  = 3;
      end
    else
      sinfo(i).datatype = -1;
    end
   
    
    % resampled
    sinfo(i).resampled = ~isempty(strfind(sinfo(i).posside,'.resampled'));
    % template
    sinfo(i).template  = ~isempty(strfind(lower(sinfo(i).ff),'.template')); 
    if sinfo(i).template,  sinfo(i).resampled = 1; end
    

    % name / texture
    %  -----------------------------------------------------------------
    % ... name extraction is a problem, because the name can include points
    % and also the dataname / texture can include points ...
    resi = [strfind(sinfo(i).posside,'template.'),... 
            strfind(sinfo(i).posside,'resampled.'),...
            strfind(sinfo(i).posside,'sphere.reg.')]; 
    if ~isempty(resi)
      sinfo(i).name = cat_io_strrep(sinfo(i).posside(max(resi):end),...
        {'template.','resampled.','sphere.reg'},''); %sinfo(i).posside,
      if ~isempty(sinfo(i).name) && sinfo(i).name(1)=='.', sinfo(i).name(1)=[]; end
      sinfo(i).texture = sinfo(i).posside(1:min(resi)-2);
    else
      % without no template/resampled string
      doti = strfind(sinfo(i).posside,'.');
      if numel(doti)==0 
      % if not points exist that the string is the name
        sinfo(i).name    = '';
        sinfo(i).texture = sinfo(i).posside;
      elseif numel(doti)==1 
      % if one point exist that the first string is the dataname and the second the subject name 
        sinfo(i).name    = sinfo(i).posside(doti+1:end);
        sinfo(i).texture = sinfo(i).posside(1:doti-1);
      else
      % this is bad
        sinfo(i).name    = sinfo(i).posside(min(doti)+1:end);
        sinfo(i).texture = sinfo(i).posside(1:min(doti)-1);
      end
    end
    if verb
      fprintf('%50s: s%04.1f %2s ',sinfo(i).ff,sinfo(i).smoothed,sinfo(i).side);
      cat_io_cprintf([0.2 0.2 0.8],'%15s',sinfo(i).texture);
      cat_io_cprintf([0.0 0.5 0.2],'%15s',sinfo(i).name);
      fprintf('%4s\n',sinfo(i).ee);
    end
    % dataname
    sinfo(i).dataname  = cat_io_strrep(sinfo(i).posside,{sinfo(i).name,'template.','resampled.'},''); 
    if ~isempty(sinfo(i).dataname) && sinfo(i).dataname(end)=='.', sinfo(i).dataname(end)=[]; end
    
    % ROI
    sinfo(i).roi = ~isempty(strfind(sinfo(i).posside,'.ROI'));
    
    
    
    % find Mesh and Data Files
    %  -----------------------------------------------------------------
    sinfo(i).Pmesh = '';
    sinfo(i).Pdata = '';
    % here we know that the gifti is a surf
    if sinfo(i).statready 
      sinfo(i).Pmesh = sinfo(i).fname;
      sinfo(i).Pdata = sinfo(i).fname;
    end
    % if we have read the gifti than we can check for the fields
    if isempty(sinfo(i).Pmesh) && sinfo(i).exist && readsurf && isfield(S,'vertices')
      sinfo(i).Pmesh = sinfo(i).fname; 
    end
    if isempty(sinfo(i).Pdata) && sinfo(i).exist && readsurf && isfield(S,'cdata')
      sinfo(i).Pdata = sinfo(i).fname;
    end
    % if the dataname is central we got a mesh or surf datafile
    if isempty(sinfo(i).Pdata) || isempty(sinfo(i).Pmesh) 
      switch sinfo(i).texture
        case {'defects'} % surf
          sinfo(i).Pmesh = sinfo(i).fname;
          sinfo(i).Pdata = sinfo(i).fname;
        case {'central','inner','outer','sphere','hull'} % only mesh
          sinfo(i).Pmesh = sinfo(i).fname;
          sinfo(i).Pdata = '';
        case {'thickness','gyrification','frac','logsulc','GWMdepth','WMdepth','CSFdepth',...
             'depthWM','depthGWM','depthCSF','depthWMg',...
             'gyruswidth','gyruswidthWM','sulcuswidth'} % only thickness
          sinfo(i).Pdata = sinfo(i).fname;
      end
    end
    % if we still dont know what kind of datafile, we can try to find a
    % mesh surface
    if isempty(sinfo(i).Pmesh) 
      if strcmp(ee,'.gii') && isempty(sinfo(i).side)
        sinfo(i).Pmesh = sinfo(i).fname;
        sinfo(i).Pdata = sinfo(i).fname;
      else
        % template mesh handling !!!
        Pmesh = char(cat_surf_rename(sinfo(i),'dataname','central','ee','.gii'));
        if exist(Pmesh,'file')
          sinfo(i).Pmesh = Pmesh;
          sinfo(i).Pdata = sinfo(i).fname;
        end
      end
    end
    % if we got still no mesh than we can use SPM.mat information or average mesh
    % ...
    if isempty(sinfo(i).Pmesh) %&& sinfo(i).ftype==1
      try 
        if ischar(SPM.xVol.G)
          % data or analysis moved or data are on a different computer?
          if ~exist(SPM.xVol.G,'file')
            [pp2,ff2,xx2] = spm_fileparts(SPM.xVol.G);
            if strfind(ff2,'.central.freesurfer')
            disp('FS')
              if strfind(pp2,'templates_surfaces_32k')
                SPM.xVol.G = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces_32k',[ff2 xx2]);
              else
                SPM.xVol.G = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces',[ff2 xx2]);
              end
            end
          end

          sinfo(i).Pmesh = SPM.xVol.G;
        else
          % 32k mesh?
          if SPM.xY.VY(1).dim(1) == 32492 || SPM.xY.VY(1).dim(1) == 64984
            sinfo(i).Pmesh = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces_32k',...
              [sinfo(i).side '.central.freesurfer.gii']);
          else
            sinfo(i).Pmesh = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces',...
              [sinfo(i).side '.central.freesurfer.gii']);
          end
        end
      catch
        % 32k mesh? 
        switch sinfo(i).ee
          case '.gii'
            if sinfo(i).exist && ~readsurf
              S = gifti(P{i});
            end
          case '.annot'
            if sinfo(i).exist && ~readsurf
              clear S; 
              try
                S = cat_io_FreeSurfer('read_annotation',P{1});
              end
            end
        end
        
        if isfield(S,'cdata') && (length(S.cdata) == 32492 || length(S.cdata) == 64984)
          sinfo(i).Pmesh = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces_32k',...
            [sinfo(i).side '.central.freesurfer.gii']);
        elseif isfloat(S) && (length(S) == 32492 || length(S) == 64984)
          sinfo(i).Pmesh = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces_32k',...
            [sinfo(i).side '.central.freesurfer.gii']);
        else
          sinfo(i).Pmesh = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces',...
            [sinfo(i).side '.central.freesurfer.gii']);
        end
      end
      sinfo(i).Pdata = sinfo(i).fname;
    end
    
    [ppm,ffm,eem]        = fileparts(sinfo(i).Pmesh);
    sinfo(i).Phull       = fullfile(ppm,strrep(strrep([ffm eem],'.central.','.hull.'),'.gii',''));
    sinfo(i).Psphere     = fullfile(ppm,strrep([ffm eem],'.central.','.sphere.'));
    sinfo(i).Pspherereg  = fullfile(ppm,strrep([ffm eem],'.central.','.sphere.reg.'));
    sinfo(i).Pdefects    = fullfile(ppm,strrep([ffm eem],'.central.','.defects.'));
    if ~exist(sinfo(i).Psphere ,'file'), sinfo(i).Psphere  = ''; end
    if ~exist(sinfo(i).Pdefects,'file'), sinfo(i).Pdefects = ''; end

    
    if sinfo(i).exist && readsurf
      if isfield(S,'vertices'), 
        sinfo(i).nvertices = size(S.vertices,1);
      else
        if ~isempty(sinfo(i).Pmesh) && exist(sinfo(i).Pmesh,'file')
          S2 = gifti(sinfo(i).Pmesh);
          if ~isstruct(S), clear S; end
          if isfield(S2,'vertices'), S.vertices = S2.vertices; else S.vertices = []; end
          if isfield(S2,'faces'),    S.faces    = S2.faces;    else S.faces = []; end
        end
        if isfield(S,'vertices'),
          sinfo(i).nvertices = size(S.vertices,1);
        elseif isfield(S,'cdata'),
          sinfo(i).nvertices = size(S.cdata,1);
        else 
          sinfo(i).nvertices = nan;
        end
      end
      if isfield(S,'faces'),    sinfo(i).nfaces    = size(S.faces,1); end
      if isfield(S,'cdata'),    sinfo(i).ncdata    = size(S.cdata,1); end
    end
    
    sinfo(i).catxml = fullfile(pp,['cat_' sinfo(i).name '*.xml']);
    if ~exist(sinfo(i).catxml,'file'), sinfo(i).catxml = ''; end 
    
    if nargout>1
      varargout{2}{i} = S; 
    else
      clear S
    end
  end
  varargout{1} = sinfo; 
end