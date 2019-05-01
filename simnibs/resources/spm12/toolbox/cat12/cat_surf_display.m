function varargout = cat_surf_display(varargin)
% ______________________________________________________________________
% Function to display surfaces. Wrapper to cat_surf_render.
%
% [Psdata] = cat_surf_display(job)
% 
% job.data      .. (lh|rh|mesh).* surfaces 
% job.colormap  .. colormap
% job.caxis     .. range of the colormap
% job.multisurf .. load both sides, if possible (default = 0)
%                   1 - load other side of same structure
%                   2 - load all structures of the same side 
%                   3 - load all structures of both sides  
% job.usefsaverage .. use average surface (for resampled data only)
%                  (default = 0)
% job.view      .. view 
%                   l=left, r=right
%                   a=anterior, p=posterior
%                   s=superior, i=inferior
% job.verb      .. SPM command window report (default = 1)
% job.readsurf  .. get full surface informtion by loading the image
%                  (default = 1; see cat_surf_info)
% job.parent    .. axis handle to print in other (sub)figures
%
% job.imgprint.do   .. print image (default = 0)
% job.imgprint.type .. render image type (default = png)
% job.dpi           .. print resolution of the image (default = 600 dpi)
%
% Examples: 
%  - Open both hemispheres of one subject S01:
%   cat_surf_display(struct('data','lh.thickness.S01.gii','multisurf',1))
%  - Use another scaling of the intensities
%   cat_surf_display(struct('caxis',[0 10]))
% ______________________________________________________________________
% Robert Dahnke
% $Id: cat_surf_display.m 1271 2018-02-08 14:28:51Z gaser $

  SVNid = '$Rev: 1271 $';
  if nargout>0, varargout{1}{1} = []; end   

  if nargin>0
    if isstruct(varargin{1})
      job = varargin{1};
      if ~isfield(job,'data') || isempty(job.data)
        if cat_get_defaults('extopts.expertgui')
          job.data = spm_select([1 24],'any','Select surfaces or textures','','','(lh|rh|lc|rc|mesh).*');
        else
          job.data = spm_select([1 24],'any','Select surfaces or textures','','','.*gii');
        end
        job.imgprint.do    = 0;
        job.imgprint.close = 0;  
      end
    else
      job.data = varargin{1};
    end
  else
    if cat_get_defaults('extopts.expertgui')
        job.data = spm_select([1 24],'any','Select surfaces or textures','','','(lh|rh|lc|rc|mesh).*');
    else
        job.data = spm_select([1 24],'any','Select surfaces or textures','','','.*gii');
    end
    job.imgprint.do    = 0;
    job.imgprint.close = 0;  
  end
  if isempty(job.data), return; end
  job.data = cellstr(job.data);
  
  % scaling options for textures
  def.colormap = '';
  def.usefsaverage = 0; 
  def.caxis    = []; % default/auto, range
  def.expert   = cat_get_defaults('extopts.expertgui'); 
  
  % print options ... just a quick output > cat_surf_print as final function 
  def.imgprint.type  = 'png';
  def.imgprint.dpi   = 600;
  def.imgprint.fdpi  = @(x) ['-r' num2str(x)];
  def.imgprint.ftype = @(x) ['-d' num2str(x)];
  def.imgprint.do    = 0;
  def.imgprint.close = 0;
  def.imgprint.dir   = '';
  
  % multi-surface output for one subject 
  def.multisurf = 0; % 0 - no; 1 - both hemispheres;
  def.verb      = 1;
  def.readsurf  = 0;  % readsurf=1 for individual average surface (e.g. apes); readsurf=0 for group average surface 
  
  job = cat_io_checkinopt(job,def);
  
  %% ... need further development 
  sinfo = cat_surf_info(job.data,job.readsurf,job.usefsaverage); 

  if job.verb
    spm('FnBanner',mfilename,SVNid); 
  end

  for i=1:numel(job.data)
    if job.usefsaverage
    
      if ~isempty(strfind(fileparts(sinfo(i).Pmesh),'_32k'))
        templates_surfaces = 'templates_surfaces_32k';
      else
        templates_surfaces = 'templates_surfaces';
      end
    
      job.fsaverage    = {
        fullfile(spm('dir'),'toolbox','cat12',templates_surfaces,'lh.central.freesurfer.gii');  
        fullfile(spm('dir'),'toolbox','cat12',templates_surfaces,'lh.inflated.freesurfer.gii');  
        fullfile(spm('dir'),'toolbox','cat12',templates_surfaces,'lh.central.Template_T1_IXI555_MNI152_GS.gii');  
        };
      sinfo(i).Pmesh = cat_surf_rename(job.fsaverage{job.usefsaverage},'side',sinfo(i).side); 
    end
        
    % load multiple surfaces
    % 3 - load all structures of both sides  
    % 2 - load all structures of the same side 
    % 1 - load other side of same structure
    if job.multisurf
      if strcmp('r',sinfo(i).side(1)), oside = ['l' sinfo(i).side(2)]; else oside = ['r' sinfo(i).side(2)]; end
      if job.multisurf==3
        Pmesh = [ ...
          cat_surf_rename(sinfo(i).Pmesh,'side','lh') cat_surf_rename(sinfo(i).Pmesh,'side','rh') ...
          cat_surf_rename(sinfo(i).Pmesh,'side','lc') cat_surf_rename(sinfo(i).Pmesh,'side','rc')]; 
        Pdata = [ ...
          cat_surf_rename(sinfo(i).Pdata,'side','lh') cat_surf_rename(sinfo(i).Pdata,'side','rh') ...
          cat_surf_rename(sinfo(i).Pdata,'side','lc') cat_surf_rename(sinfo(i).Pdata,'side','rc')];
      elseif job.multisurf==2
        if strcmp('h',sinfo(i).side(2)), oside = [sinfo(i).side(1) 'c']; else oside = [sinfo(i).side(1) 'h']; end
        Pmesh = [sinfo(i).Pmesh cat_surf_rename(sinfo(i).Pmesh,'side',oside)]; 
        Pdata = [sinfo(i).Pdata cat_surf_rename(sinfo(i).Pdata,'side',oside)]; 
      else
        Pmesh = [sinfo(i).Pmesh cat_surf_rename(sinfo(i).Pmesh,'side',oside)]; 
        Pdata = [sinfo(i).Pdata cat_surf_rename(sinfo(i).Pdata,'side',oside)]; 
      end
      for im=numel(Pmesh):-1:1
        if ~exist(Pmesh{im},'file'), Pmesh(im) = []; end
        if ~isempty(Pdata) && ~exist(Pdata{im},'file'), Pdata(im) = []; end
      end
      if numel(Pmesh)==1; Pmesh=char(Pmesh); end
      if numel(Pdata)==1; Pdata=char(Pdata); end
    else
      Pmesh = sinfo(i).Pmesh;
      Pdata = sinfo(i).Pdata; 
    end
    
    if job.verb
      fprintf('Display %s\n',spm_file(job.data{i},'link','cat_surf_display(''%s'')'));
    end

    if job.expert>1
      fprintf('Developer display mode!\n');
    end

    if ~isempty(Pdata) && ~all(strcmp(Pmesh,Pdata)) 
      % only gifti surface without texture
      if isfield(job,'parent')
        if job.expert<2
          h = cat_surf_render('disp',Pmesh,'Pcdata',Pdata,'parent',job.parent);
        else
          h = cat_surf_render2('disp',Pmesh,'Pcdata',Pdata,'parent',job.parent);
        end
      else
        if job.expert<2
          h = cat_surf_render('disp',Pmesh,'Pcdata',Pdata);
        else
          h = cat_surf_render2('disp',Pmesh,'Pcdata',Pdata);
        end
      end  
    else
      % only gifti surface without texture
      if isfield(job,'parent')
        if job.expert<2
          h = cat_surf_render(Pmesh,'parent',job.parent);
        else
          h = cat_surf_render2(Pmesh,'parent',job.parent);
        end
      else
        if job.expert<2
          h = cat_surf_render(Pmesh);
        else
          h = cat_surf_render2(Pmesh);
        end
      end
    end

    set(h.figure,'MenuBar','none','Toolbar','none','Name',spm_file(job.data{i},'short60'),'NumberTitle','off');

    % shift each figure slightly
    if i==1
        pos = get(h.figure,'Position');
    else
        pos = pos - [20 20 0 0];
        set(h.figure,'Position',pos);
    end

    if sinfo(i).label, continue; end      
      
   % try   
    %% textur handling
      if job.expert<2 
        cat_surf_render('ColourBar',h.axis,'on');
      else
        cat_surf_render2('ColourBar',h.axis,'on');
      end
      
      if ~job.multisurf && strcmp(sinfo(i).side,'rh'), view(h.axis,[90 0]); end
      
      
      % temporary colormap
      if any(strcmpi({'neuromorphometrics','lpba40','ibsr','hammers','mori','aal'},sinfo(i).dataname))
        %%
        switch lower(sinfo(i).dataname)
          case 'neuromorphometrics', rngid=3; 
          case 'lpba40',             rngid=12; 
          case 'ibsr',               rngid=1; 
          case 'hammers',            rngid=5;  
          case 'mori',               rngid=3; 
          case 'aal',                rngid=11; 
          otherwise,                 rngid=1; 
        end
        
        sideids = ceil(max(h.cdata(:))/2)*2;  
        if exist('rng','builtin') == 5
          rng('default')
          rng(rngid)
        else
          rand('state',rngid);
        end

        cmap = colorcube(ceil((sideids/2) * 8/7)); % greater to avoid grays
        cmap(ceil(sideids/2):end,:)=[]; % remove grays
        cmap(sum(cmap,2)<0.3,:) = min(1,max(0.1,cmap(sum(cmap,2)<0.3,:)+0.2)); % not to dark
        cmap = cmap(randperm(size(cmap,1)),:); % random 
        cmap = reshape(repmat(cmap',2,1),3,size(cmap,1)*2)'; 
       
        if job.expert<2
          cat_surf_render('ColourMap',h.axis,cmap);
        else
          cat_surf_render2('ColourMap',h.axis,cmap);
        end
        %%
        continue
      else
        if isempty(job.colormap)
          if job.expert<2
            h = cat_surf_render('ColourMap',h.axis,jet(256)); 
          else
            h = cat_surf_render2('ColourMap',h.axis,jet(256)); 
          end
        else
          if job.expert<2
            h = cat_surf_render('ColourMap',h.axis,eval(job.colormap));
          else
            h = cat_surf_render2('ColourMap',h.axis,eval(job.colormap));
          end
        end
      end
      
      % scaling
      if isempty(job.caxis)
        switch sinfo(i).texture
          case {'ROI'}
            %%
            if     strfind(sinfo(i).posside,'-Igm.ROI'),  clim = [2/3 2/3] .* [0.9 1.1];        % balanced PVE
            elseif strfind(sinfo(i).posside,'-Iwm.ROI'),  clim = [0.85 1.05];                   % below 1 because of a lot of GM/WM PVE
            elseif strfind(sinfo(i).posside,'-Icsf.ROI'), clim = [1.33/3 1.33/3] .* [0.8 1.2];  % higher 1/3 because of a lot of GM/CSF PVE
            else                                          clim = cat_vol_iscaling(h.cdata);
            end
            if job.expert<2
              cat_surf_render('clim',h.axis,clim);
            else
              cat_surf_render2('clim',h.axis,clim);
            end
          case {'defects','sphere'}
            % no texture
          case {'central'}
            % default curvature
            set(h.patch,'AmbientStrength',0.2,'DiffuseStrength',0.8,'SpecularStrength',0.1)
          case ''
            % no texture name
            if ~isempty(h.cdata)
              clim = cat_vol_iscaling(h.cdata);
              if clim(1)<0
                clim = [-max(abs(clim)) max(abs(clim))];
                if job.expert<2
                  cat_surf_render('clim',h.axis,clim);
                else
                  cat_surf_render2('ColourMap',h.axis,cat_io_colormaps('BWR',128));
                  cat_surf_render2('clim',h.axis,clim);
                end
              else
                if job.expert<2
                  cat_surf_render('clim',h.axis,clim);
                else
                  cat_surf_render2('ColourMap',h.axis,cat_io_colormaps('hotinv',128)); 
                  cat_surf_render2('clim',h.axis,clim);
                end
              end
              
            end
          otherwise
            %%
            ranges = {
              ... name single group
              'thickness'         [0.5  5.0]  [0.5  5.0]
              'gyruswidthWM'      [0.5  8.0]  [1.0  7.0]
              'gyruswidth'        [1.0 12.0]  [1.5 11.0]
              'fractaldimension'  [0.0  4.0]  [1.0  4.0]
              'sulcuswidth'       [0.0  3.0]  [0.0  3.0]
              'gyrification'      [ 15   35]  [ 15   35]
              'sqrtsulc'          [0.0  5.0]  [0.0  5.0]
              'WMdepth'           [1.0  6.0]  [1.0  5.0]
              'GWMdepth'          [1.5 10.0]  [1.5  9.0]
              'CSFdepth'          [0.5  2.0]  [0.5  2.0]
              'depthWM'           [0.0  4.0]  [0.0  3.0]
              'depthWMg'          [0.0  1.0]  [0.0  0.5]
              'depthGWM'          [0.5  5.0]  [2.5  6.0]
              'depthCSF'          [0.5  2.0]  [0.5  2.0]  
            };

            texturei = find(cellfun('isempty',strfind(ranges(:,1),sinfo(i).texture))==0,1,'first');

            if ~isempty(texturei)
              if job.expert<2
                cat_surf_render('clim',h.axis,ranges{texturei,3});
              else
                cat_surf_render2('clim',h.axis,ranges{texturei,3});
              end
            else
              if ~isempty(h.cdata)
                clim = cat_vol_iscaling(h.cdata);  
                if job.expert<2
                  cat_surf_render('clim',h.axis,clim);
                else
                  cat_surf_render2('clim',h.axis,clim);
                end
              end
            end
        end    
      else
        if job.expert<2
          cat_surf_render('clim',h.axis,job.caxis);
        else
          cat_surf_render2('clim',h.axis,job.caxis);
        end
      end
%     catch %#ok<CTCH>
%       if ~exist('h','var')
%         try
%           cat_io_cprintf('err',sprintf('Texture error. Display surface only.\n'));
%           h = cat_surf_render(job.data{i});
%         catch %#ok<CTCH>
%           cat_io_cprintf('err',sprintf('ERROR: Can''t display surface %s\n',job.data{i})); 
%         end
%       end
%       continue
%     end
    
    
    
    
    %% view
    
    if ~isfield(job,'view')
      if strcmp(sinfo(i).side,'lh') && ~job.multisurf
        job.view = 'left'; 
      elseif strcmp(sinfo(i).side,'rh') && ~job.multisurf
        job.view = 'right';
      else
        job.view = 'top';
      end
    end
    
    switch lower(job.view)
      case {'r','right'},                 cat_surf_render('view',h,[  90   0]); viewname = '.r';
      case {'l','left'},                  cat_surf_render('view',h,[ -90   0]); viewname = '.l'; 
      case {'t','s','top','superior'},    cat_surf_render('view',h,[   0  90]); viewname = '.s';
      case {'b','i','bottom','inferior'}, cat_surf_render('view',h,[-180 -90]); viewname = '.i'; 
      case {'f','a','front','anterior'},  cat_surf_render('view',h,[-180   0]); viewname = '.a';
      case {'p','back','posterior'},      cat_surf_render('view',h,[   0   0]); viewname = '.p';
      otherwise
        if isnumeric(job.view) && size(job.view)==2
          view(job.view); viewname = sprintf('.%04dx%04d',mod(job.view,360));
        else
          error('Unknown view.\n')
        end
    end
    
    
    
    
    %% print
    if job.imgprint.do 
      %%
      if isempty(job.imgprint.dir), ppp = sinfo(i).pp; else  ppp=job.imgprint.dir;  end
      if ~exist(ppp,'dir'), mkdir(ppp); end
      pfname = fullfile(ppp,sprintf('%s%s.%s',sinfo(i).ff,viewname,job.imgprint.type));
      print(h.figure , def.imgprint.ftype(job.imgprint.type) , job.imgprint.fdpi(job.imgprint.dpi) , pfname ); 
      
      if job.imgprint.close
        close(h.figure);
      end
    end
    
    if nargout>0
      varargout{1}{i} = h;
    end   
  end
end



